# coding=utf-8
import dash
from dash import html, dcc
from dash_extensions.enrich import  DashProxy, MultiplexerTransform
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc
import requests
from flask import Flask, send_file
from flask_caching import Cache
import layout
import back
import predictor
import pandas as pd
from datetime import datetime
import os

from flask import Flask



FASTA = []
cargaFasta = True
MAX_NUM_SEQUENCES = 500



# Initialize the Dash app
external_stylesheets = [dbc.themes.SANDSTONE, "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.15.1/css/all.min.css"]

server = Flask(__name__)
app = DashProxy(
    __name__
    ,server=server
    , external_stylesheets=external_stylesheets
    , prevent_initial_callbacks=True
    , transforms=[MultiplexerTransform()]
)
#server = app.server




#Para manejo dinamico de callbacks
app.config.suppress_callback_exceptions = True

cache = Cache(app.server, config={
    #'CACHE_TYPE': 'redis',
    # Note that filesystem cache doesn't work on systems with ephemeral
    # filesystems like Heroku.
    'CACHE_TYPE': 'filesystem',
    'CACHE_DIR': 'cache-directory',

    # should be equal to maximum number of users on the app at a single time
    # higher numbers will store more data in the filesystem / redis cache
    'CACHE_THRESHOLD': 1
})


#Descargar filtro FASTA
@server.route("/download/secuenciasFasta/<path>")
def downloadFiltroFasta(path):
    return send_file("cache-directory/secuencias"+str(path)+".FASTA", as_attachment=True)

#Descargar Clasificacion
@server.route("/download/clasificacionExcel/<path>")
def downloadClasificacionExcel(path):
    return send_file("cache-directory/clasificacion"+str(path)+".xlsx", as_attachment=True)


app.layout = layout.principal()

# Navbar. Cambiar entre "Resultados" y " Analisis Estadistico"
@app.callback([Output('right-column-content-loading', 'children'),
               Output('resultados-button', 'active'),
               Output('analisis-button', 'active')],
              [Input('resultados-button', 'n_clicks'),
               Input('analisis-button', 'n_clicks')],
               State('session-id','children'),
              )
def display_page(resultados_clicks, analisis_clicks, session_id):
    if not os.path.exists("cache-directory/clasificacion"+session_id+".xlsx"):
        return html.Div('Aun no se ha generado ningun archivo', style={'color': '#FF0000'}), False, False

    ctx = dash.callback_context
    button_id = ctx.triggered[0]['prop_id'].split('.')[0]

    if button_id == "resultados-button":
            clasificacion = pd.read_excel("cache-directory/clasificacion"+session_id+".xlsx")
            tabClasificacion = layout.tabsResultadoClasificacion(clasificacion, functionalities)
            return tabClasificacion, True, False
    elif button_id == "analisis-button":
        tabAnalisisEstadistico = layout.analisisEstadistico
        return tabAnalisisEstadistico, False, True
    else:
        return html.Div(), False, False

#GRAFICA ANALISIS CLASIFICACION
@app.callback(
    Output('output-grafica-analisis-clasificacion', 'children'),
    Input('propiedad-grafica-analisis-clasificacion', 'value'),
    State('session-id', 'children')
)
def actualizar_grafica_analisis(propiedad_seleccionada, session_id):
    """
    Actualiza la gráfica de análisis de clasificación basada en la propiedad seleccionada.
    
    Args:
        propiedad_seleccionada: Propiedad seleccionada para analizar
        session_id: Identificador único de la sesión
    
    Returns:
        Componente de gráfica KDE
    """
    # Cargar datos de clasificación
    ruta_archivo = f"cache-directory/clasificacion{session_id}.xlsx"
    clasificacion = pd.read_excel(ruta_archivo)
    
    # Ordenar datos por Bosque Aleatorio
    clasificacion = clasificacion.sort_values('Bosque Aleatorio')
    
    # Generar y retornar la gráfica
    return layout.tablaKDEDoble(clasificacion, propiedad_seleccionada)


#GENERAR
@app.callback([Output('output-descarga','children'),
               Output('right-column-content-loading', 'children')],

              [Input('generar', 'n_clicks')],

              [State('numero','value'),
               State('metodo-generacion-input','value'),
               State('funcionalidades-input', 'value'),
               State('aminoacidos_input_text', 'value'),
               State('tamano_secuencias_slider','value'),
               State('probabilidad_funcional_input', 'value'),
               State('session-id','children')
               ]
              )
def generar(clicks, numero, metodo_generacion, checked_values, aminoacidos_excluir_list, tamano_secuencias_slider, probabilidad_funcional_input, session_id):
    if clicks != None:
        # Validacion numero de secuencias
        if(numero == None):
            return html.H4(""),html.Div('Debe definir el numero de peptidos a generar', style={'color': '#FF0000'})

        if(int(numero) <= 0):
            return html.H4(""),html.Div('El número de secuencias a generar debe ser un número positivo', style={'color': '#FF0000'})

        if(int(numero) > MAX_NUM_SEQUENCES):
            return html.H4(""),html.Div(f'El número de secuencias a generar debe ser menor que {MAX_NUM_SEQUENCES}', style={'color': '#FF0000'})

        #Validación metodo de generacion
        if(len(metodo_generacion) <= 0):
            return html.H4(""),html.Div('Debe seleccionar al menos un metodo de generación', style={'color': '#FF0000'})

        #Validación de funcionalidades
        if(len(checked_values) <= 0):
            return html.H4(""),html.Div('Debe seleccionar al menos una funcionalidad', style={'color': '#FF0000'})


        # Validacion umbral de generacion
        if(float(probabilidad_funcional_input) < 0 or float(probabilidad_funcional_input) > 1):
            return html.H4(""),html.Div('El umbral de generación debe ser un valor decimal entre 0 y 1', style={'color': '#FF0000'})



        # API query
        api_endpoint = 'https://multipepgen-api.medellin.unal.edu.co/generate/'
        global functionalities
        functionalities = [value for value in checked_values]
        payload = {
            'n_gens': numero,
            'functionalities': functionalities
        }
        response = requests.post(api_endpoint, json=payload)
        data = response.json()   # json file
        print("data------------------------")
        print(data)


        # Tomar peptidos de cada metodo de generacion
        peptides = []
        if len(metodo_generacion) == 1:
            if metodo_generacion[0] == "generated_better_prediction":
                peptides = data['generated_better_prediction']
            else:
                peptides = data['generated_more_stable']
        else:
            peptides_1 = data['generated_better_prediction']
            peptides_2 = data['generated_more_stable']
            peptides = peptides_1 + peptides_2


        # Validar formato aminoacidos a excluir
        if aminoacidos_excluir_list != None and aminoacidos_excluir_list != '':
            for aminoacido in aminoacidos_excluir_list.split(","):
                if len(aminoacido) > 1 or aminoacido.isupper() == False:
                    return html.H4(""),html.Div('El formato de ingreso para excluir aminoácidos no cumple con las especificaciones. Por ejemplo, un formato válido para excluir los aminoácidos A, C y K sería: A, C, K',         style={'color': '#FF0000'})


        # Filtrar secuencias por aminoacidos y tamaño
        filter_peptides = []
        print(' Tamaño de secuencias: ', tamano_secuencias_slider[0], tamano_secuencias_slider[1])
        for p in peptides:
            if back.sizeFilter(p, tamano_secuencias_slider[0], tamano_secuencias_slider[1]) == True:
                if aminoacidos_excluir_list != None:
                    if len(aminoacidos_excluir_list) > 0 and aminoacidos_excluir_list[0] != "":
                        if back.aminoacidosFilter(p, aminoacidos_excluir_list.split(",")) == True:
                            filter_peptides.append(p)
                    else:
                        filter_peptides.append(p)
                else:
                    filter_peptides.append(p)

        # Validar peptidos filtrarPeptidos
        if len(filter_peptides) == 0:
            return html.H4(""),html.Div('No se han encontrado secuencias con los filtros seleccionados', style={'color': '#FF0000'})


        # Generar identificadores para las secuencias
        ids = []
        current_datetime = datetime.now()
        formatted_datetime = current_datetime.strftime("%Y-%m-%d_%H:%M:%S")
        id_counter = 1
        for p in filter_peptides:
            if p in data['generated_better_prediction']:
                ids += ["MULTIPEPGEN_SEQUENCE_BETTER_PREDICTION_{0}_{1:04d}".format(formatted_datetime, id_counter)]
                id_counter += 1
            elif p in data['generated_more_stable']:
                ids += ["MULTIPEPGEN_SEQUENCE_MORE_STABLE_{0}_{1:04d}".format(formatted_datetime, id_counter)]
                id_counter += 1


        back.writeFasta(ids,filter_peptides,"cache-directory/secuencias" + session_id + ".FASTA")
        clasificacion = predictor.predecir(filter_peptides)


        # Filtrar por probabilidad funcional
        threshold = int(100*(float(probabilidad_funcional_input)))
        clasificacion = clasificacion[(clasificacion['Regresión Lógistica'] >= threshold) &
                                      (clasificacion['Red Neuronal'] >= threshold) &
                                      (clasificacion['Arbol de Decisión'] >= threshold) &
                                      (clasificacion['Bosque Aleatorio'] >= threshold) &
                                      (clasificacion['XGboost'] >= threshold)]

        if(len(clasificacion) == 0):
            return html.H4(""),html.H4(""),html.Div('No se han encontrado secuencias con los filtros seleccionados', style={'color': '#FF0000'})


        metodo_generacion_list = []
        for p in clasificacion['Peptido'].tolist():
            if p in data['generated_better_prediction']: metodo_generacion_list.append("better")
            elif p in data['generated_more_stable']: metodo_generacion_list.append("more stable")

        clasificacion['metodo generacion'] = metodo_generacion_list


        clasificacion.to_excel("cache-directory/clasificacion"+session_id+".xlsx")

        descarga_resultados_card = dbc.Card(
            children=[
                dbc.CardBody(
                    [
                        dbc.Button("Descargar EXCEL",
                                   outline=True,
                                   color="success",
                                   class_name="mt-3",
                                   href='/download/clasificacionExcel/' + str(session_id),
                                   download='generacion.xlsx',
                                   external_link=True),
                        dbc.Button("Descargar FASTA",
                                   outline=True,
                                   color="success",
                                   class_name="mt-3",
                                   href='/download/secuenciasFasta/' + str(session_id),
                                   download='secuencias.FASTA',
                                   external_link=True)
                    ]
                )
            ],
            class_name="ml-3"
        )

        #outputPlotClasificacion = layout.tabsResultadoClasificacion(clasificacion)

        clasificacion = pd.read_excel("cache-directory/clasificacion"+session_id+".xlsx")
        tabClasificacion = layout.tabsResultadoClasificacion(clasificacion, functionalities)

        return dcc.Loading(descarga_resultados_card), dcc.Loading(tabClasificacion)
        #return str(data), None
        #return html.H4("cache-directory/clasificacion" + session_id + ".xlsx"), html.H4("c")
    else:
        return html.H4(""), html.H4("")


# Run the app
if __name__ == "__main__":
    app.run(
        debug=True
        ,port=8100
    )
