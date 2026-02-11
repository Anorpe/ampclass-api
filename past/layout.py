
from dash import dcc
from dash import html
from dash import dash_table
import dash_bootstrap_components as dbc
import plotly.graph_objs as go
import numpy as np
import uuid
import scipy.stats as stats


minWidthGrafico = '300px'
widthGrafico = '31%'
heightGrafico = '400px'
#heightGrafico = '20wv'
margen = '1%'

TITULOS = {'longitud':'Longitud',
           'carga':'Carga',
           'momento hidrofobico':'Momento Hidrofóbico',
           'porcentaje hidrofobico':'% Aminoácidos Hidrofóbicos',
           'punto isoelectrico':'Punto Isoeléctrico',
           'indice de boman':'Índice de Boman',
           'polar angles':'Ángulo Polar',
           'wimley':'Wimley White'}

# Description
description_section = dbc.Row(
    children=[
        dbc.Col(
            dbc.Card(
                id="description-card",
                outline=False,
                children=[
                    dbc.CardBody(
                        [
                            dbc.Row(
                                [
                                    dbc.Col(
                                        [
                                            html.Img(
                                                src=("/assets/MultiPepGen_logo.png"),
                                                width="200px",
                                            )
                                        ],
                                        md="auto",
                                    ),
                                    dbc.Col(
                                        children=[
                                            html.P(
                                                "MultiPepGen es una red neuronal Conditional GAN con celdas LSTM para la generación de secuencias sintéticas de péptidos antimicrobianos con la capacidad de generar secuencias con las siguientes funcionalidades específicas: antimicrobiana, antibacteriana, anti gramnegativa, anti grampositiva, antifúngica, antiviral y anticancerígena. Este método se evaluó a partir de una estrategia de validación enfocada en medir la calidad y la diversidad de las secuencias sintéticas generadas. A partir de una serie de entrenamientos realizados al modelo se seleccionaron dos modelos sobresalientes:",
                                            ),
                                            html.P(
                                                children=[
                                                    html.Strong("More Stable: "),
                                                    html.Span("Tuvo resultados estables o superiores en la mayoría de las medidas de validación, por lo que podríamos decir que es un modelo equilibrado tanto en términos de calidad como de la diversidad de las secuencias sintéticas  generadas.")
                                                ]
                                            ),
                                            html.P(
                                                children=[
                                                    html.Strong("Better at Prediction: "),
                                                    html.Span("Obtuvo resultados superiores en cuanto a la capacidad antimicrobiana de las secuencias generadas; Sin embargo, las propiedades fisicoquímicas y de composición de las secuencias sintéticas generadas no son similares a las secuencias de entrenamiento.")
                                                ]
                                            ),
                                        ],
                                        md=True
                                    ),
                                ]
                            ),
                        ]
                    ),
                ],
            )
        )
    ],
    class_name="mt-3",
)

numero_secuencias_input = dbc.FormFloating(
    [
        dbc.Input(type="number", id="numero", placeholder="Número de secuencias", min=1, max=10000, value=10),
        dbc.Label("Cantidad de secuencias")
    ],
    class_name="mb-3"
)

metodo_generacion_input = dbc.Card(
    children=[
        dbc.CardHeader("Metodo de generación"),
        dbc.CardBody(
            dbc.Checklist(
                options=[
                    {"label": "Better prediction", "value": "generated_better_prediction"},
                    {"label": "More Stable", "value": "generated_more_stable"},
                ],
                value = ["generated_better_prediction", "generated_more_stable"],
                id="metodo-generacion-input",
                switch=True
            )
        )
    ],
    class_name="mb-3"
)

funcionalidades_input = dbc.Card(
    children=[
        dbc.CardHeader("Funcionalidades"),
        dbc.CardBody(
            dbc.Checklist(
                options=[
                    {"label": "Antimicrobiano", "value": "antimicrobiano"},
                    {"label": "Antibacteriano", "value": "antibacteriano"},
                    {"label": "Antigramnegativo", "value": "antigramnegativo"},
                    {"label": "Antigrampositivo", "value": "antigrampositivo"},
                    {"label": "Antifúngico", "value": "antifungico"},
                    {"label": "Antiviral", "value": "antiviral"},
                    {"label": "Anticáncer", "value": "anticancer"},
                ],
                value=["antimicrobiano"],
                id="funcionalidades-input",
                switch=True,

            ),
        ),

    ],
    class_name="ml-3",
)

aminoacidos_input = dbc.Card(
    children=[
        dbc.CardHeader(
            children=[
                "Excluir Aminoácidos",
                html.H4([html.I(className="fas fa-question-circle", id="question-icon")]),
                dbc.Tooltip("Inserte los aminoacidos a excluir en mayúscula, separados por coma y sin espacios", target="question-icon", placement="left"),
            ]
        ),
        dbc.CardBody(
            children=[
                dbc.FormFloating(
                    children=[
                        dbc.Input(type="text", id="aminoacidos_input_text",placeholder="Aminoácidos a excluir"),
                        dbc.Label("Aminoácidos a excluir")
                    ],
                    class_name="mb-3"
                )
            ]
        )
    ]
)

tamano_secuencias_input = dbc.Card(
    children=[
        dbc.CardHeader("Longitud de las secuencias"),
        dbc.CardBody(
            children=[
                dcc.RangeSlider(
                    id='tamano_secuencias_slider',
                    min=0,
                    max=35,
                    step=1,
                    marks={i: str(i) for i in range(10, 36, 5)},
                    value=[0, 35],
                    tooltip={"always_visible": True, "placement": "top"},
                )
            ]
        )
    ]
)

probabilidad_input = dbc.Card(
    children=[
        dbc.CardHeader(
            children=[
                "Probabilidad funcional",
                #html.H4([html.I(className="fas fa-question-circle", id="question-icon-2")]),
                #dbc.Tooltip("Ingrese un valor decimal correspondiente al umbral para la generación de la secuencias", target="question-icon-2", placement="top"),
            ]
        ),
        dbc.CardBody(
            children=[
                dcc.Slider(
                    id="probabilidad_funcional_input",
                    min=0,
                    max=1,
                    step=0.05,
                    marks={i / 5: str(i / 5) for i in range(0, 6)},
                    value=0.5,
                    tooltip={"always_visible": True, "placement": "top"},
                )

            ]
        )
    ]
)

generar_button = dbc.Button(
    "Generar",
    id="generar",
    outline=True,
    color="primary",
    size="lg",
    class_name="mt-3 mb-3"
)

right_column_navbar = dbc.Nav(
    card=True,
    navbar=True,
    pills=True,
    children=[
        dbc.NavItem(dbc.NavLink("Resultados", id="resultados-button", active=False, href="#", n_clicks=0)),
        dbc.NavItem(dbc.NavLink("Análisis Estadístico", id="analisis-button", active=False, href="#", n_clicks=0)),
    ],
)


left_column = [
    dbc.Card(
        id="parameters-card",
        #color="light",
        class_name="mt-3",
        style={"margin-left":"15px", "text-align": "center"},
        children=[
            dbc.CardHeader("Parámetros"),
            dbc.CardBody(
                [
                    numero_secuencias_input,
                    metodo_generacion_input,
                    funcionalidades_input,
                ]
            ),
        ],
    ),
    dbc.Card(
        id="filters-card",
        class_name="mt-3",
        style={"margin-left":"15px", "text-align": "center"},
        children=[
            dbc.CardHeader("Filtros"),
            dbc.CardBody(
                [
                    aminoacidos_input,
                    tamano_secuencias_input,
                    probabilidad_input,
                    generar_button
                ]
            ),
            dcc.Loading(
                id="output-descarga-loading",
                children=[html.Div(
                            id="output-descarga",
                            children=[]
                        )],
                type="circle",
                overlay_style={"visibility":"visible", "filter": "blur(2px)"},
            #fullscreen = True,
            )

        ]
    )
]

right_column = dbc.Col(
    md=9,
    children=[
        dbc.Row(description_section),
        dbc.Row(dbc.Col(


            dbc.Card(
                    id="right-card",
                    class_name="mt-3",
                    children=[
                        dbc.CardHeader(right_column_navbar),
                        dcc.Loading(
                            id="right-column-content-loading",
                            children=[
                                dbc.CardBody(
                                    id="right-column-content",
                                    children=[]
                                ),
                            ],
                            type="circle",
                            overlay_style={"visibility": "visible", "filter": "blur(2px)"},
                            # fullscreen = True,
                        )
                    ],
                )




        )),
    ]
)

vista_principal = html.Div(
    [
        #dcc.Location(id='url', refresh=False),
        #dbc.Container
        html.Div(
            [
                dbc.Row(
                    children=[dbc.Col(left_column, md=3), right_column],
                ),
                dbc.Row(
                    html.Div(
                            id="output-generar",
                            children=[]
                        )
                ),
                dbc.Row(
                    html.Div(
                        id="output-plot-generar",
                        children=[]
                    )
                )
            ],
            #fluid=True,
        ),
    ]
)

def principal():
    session_id = str(uuid.uuid4())
    return html.Div([

    dcc.Location(id='url', refresh=True),
    html.Div(session_id, id='session-id', style={'display': 'none'}),


    html.Div(id='principal',children = [

        vista_principal

    ])
    #], type="default")
])

def tablaClasificacion_1(datos, funcionalidades):

    for i in range(len(funcionalidades)):
        if(funcionalidades[i] == "antifungico"): funcionalidades[i] = "Antifúngico"
        if(funcionalidades[i] == "anticancer"): funcionalidades[i] = "Anticáncer"


    # Tomar funcionalidades seleccionadas para generar el nombre de la tabla
    titulo = "% Probabilidad de ser " + ", ".join([i.capitalize() for i in funcionalidades])

    return dash_table.DataTable(
                            columns=[
                                {"name": [titulo, 'Secuencia'], "id": 'Peptido'},
                                {"name": [titulo, 'Regresión Logística'],"id": 'Regresión Lógistica'},
                                {"name": [titulo, 'Red Neuronal'],"id": 'Red Neuronal'},
                                {"name": [titulo, 'Arbol de Decisión'],"id": 'Arbol de Decisión'},
                                {"name": [titulo, 'Bosque Aleatorio'],"id": 'Bosque Aleatorio'},
                                {"name": [titulo, 'XGboost'],"id": 'XGboost'},
                            ],
                            merge_duplicate_headers=True,
                            data=datos.to_dict('records'),
                            style_as_list_view=True,

                            sort_action='native',

                            style_table={
                                'overflowY': 'scroll',
                                'border': 'thin lightgrey solid'
                            },
                            style_cell={
                                'padding': '5px',
                                'textAlign': 'center',
                                'overflow': 'hidden',
                                'textOverflow': 'ellipsis',
                                'maxWidth': 0,
                                'whiteSpace': 'normal'
                            },
                            style_header={
                                'whiteSpace': 'normal',
                                'height': 'auto',

                                'backgroundColor': '#D6D8DC',
                                'fontWeight': 'bold'
                            },
                            style_data_conditional=[
                                {
                                    'if': {'column_id': 'Clasificación XGboost'},
                                    'backgroundColor': '#D1D8F7',
                                    'color': 'black',
                                },
                                {
                                'if': {'column_id': 'Clasificación Bosque Aleatorio'},
                                'backgroundColor': '#D1D8F7',
                                'color': 'black',
                                },
                                {
                                    'if': {'column_id': 'Clasificación Arbol de Decisión'},
                                    'backgroundColor': '#D1D8F7',
                                    'color': 'black',
                                },
                                {
                                    'if': {'column_id': 'Clasificación Red Neuronal'},
                                    'backgroundColor': '#D1D8F7',
                                    'color': 'black',
                                },
                                {
                                    'if': {'column_id': 'Clasificación Regresión Lógistica'},
                                    'backgroundColor': '#D1D8F7',
                                    'color': 'black',
                                },
                                {
                                    'if': {
                                        'filter_query': '{metodo generacion} = "better"',
                                    },
                                    'backgroundColor': '#edeff7',
                                    'color': 'black',

                                }

                            ]

                        )

def tablaClasificacion_2(datos, funcionalidades):

    for i in range(len(funcionalidades)):
        if(funcionalidades[i] == "antifungico"): funcionalidades[i] = "Antifúngico"
        if(funcionalidades[i] == "anticancer"): funcionalidades[i] = "Anticáncer"

    # Tomar funcionalidades seleccionadas para generar el nombre de la tabla
    titulo = "Actividad " + ", ".join([i.capitalize() for i in funcionalidades])

    return dash_table.DataTable(
                            columns=[
                                {"name": [titulo, 'Secuencia'], "id": 'Peptido'},
                                {"name": [titulo, 'Regresión Logística'],"id": 'Clasificación Regresión Lógistica'},
                                {"name": [titulo, 'Red Neuronal'],"id": 'Clasificación Red Neuronal'},
                                {"name": [titulo, 'Arbol de Decisión'],"id": 'Clasificación Arbol de Decisión'},
                                {"name": [titulo, 'Bosque Aleatorio'],"id": 'Clasificación Bosque Aleatorio'},
                                {"name": [titulo, 'XGboost'],"id": 'Clasificación XGboost'},
                            ],
                            merge_duplicate_headers=True,
                            data=datos.to_dict('records'),
                            style_as_list_view=True,

                            sort_action='native',

                            style_table={
                                'overflowY': 'scroll',
                                'border': 'thin lightgrey solid'
                            },
                            style_cell={
                                'padding': '5px',
                                'textAlign': 'center',
                                'overflow': 'hidden',
                                'textOverflow': 'ellipsis',
                                'maxWidth': 0,
                                'whiteSpace': 'normal'
                            },
                            style_header={
                                'whiteSpace': 'normal',
                                'height': 'auto',

                                'backgroundColor': '#D6D8DC',
                                'fontWeight': 'bold'
                            },
                            style_data_conditional=[
                                {
                                    'if': {'column_id': 'Clasificación XGboost'},

                                    'color': 'black',
                                },
                                {
                                    'if': {'column_id': 'Clasificación Bosque Aleatorio'},

                                    'color': 'black',
                                },
                                {
                                    'if': {'column_id': 'Clasificación Arbol de Decisión'},

                                    'color': 'black',
                                },
                                {
                                    'if': {'column_id': 'Clasificación Red Neuronal'},

                                    'color': 'black',
                                },
                                {
                                    'if': {'column_id': 'Clasificación Regresión Lógistica'},

                                    'color': 'black',
                                },
                                {
                                    'if': {
                                        'filter_query': '{metodo generacion} = "better"',
                                    },
                                    'backgroundColor': '#edeff7',
                                    'color': 'black',

                                }

                            ]

                        )

def tabsResultadoClasificacion(clasificacion, funcionalidades):
    return dcc.Loading(children=[html.Div(
            [
                dbc.Row(tablaClasificacion_1(clasificacion, funcionalidades), style={"margin-down": "15px"}),
                dbc.Row(tablaClasificacion_2(clasificacion, funcionalidades), style={"margin-top": "15px"}),
            ],
            className="pretty_container")], type="default",style={"transform": "scale(6)"})


analisisEstadistico = html.Div([

        html.Div([html.H4('Estimación de Densidad del Núcleo por Propiedad',
                        style={'text-align': 'center'})]),
        html.P('En esta sección se presenta la estimación de la FDP (función de densidad de probabilidad) mediante el método KDE (Estimación de Densidad de Núcleo) y el histograma para cada propiedad de los péptidos discriminados por su clasificación dependiendo de la técnica de clasificación seleccionada.',
            style = {'margin':'15px'}),

        html.Label('Método de Clasificación:',style={'font-weight':'bold','margin':'15px'}),
        dcc.Dropdown(
            options=[
                {'label': 'Regresión Logística', 'value': 'Clasificación Regresión Lógistica'},
                {'label': 'Red Neuronal', 'value': 'Clasificación Red Neuronal'},
                {'label': 'Árbol de Decisión', 'value': 'Clasificación Arbol de Decisión'},
                {'label': 'Bosque Aleatorio', 'value': 'Clasificación Bosque Aleatorio'},
                {'label': 'XGboost', 'value': 'Clasificación XGboost'},
            ],
            id='propiedad-grafica-analisis-clasificacion', style={'width': 400,'margin':'10px'}
        ),
        dcc.Loading(html.Div(id='output-grafica-analisis-clasificacion'))
    ])


def KDEDoble(data,propiedad,metodo):
    """
    Dado un dataframe, una columna, y escoger un metodo de aprendizaje automatico genera un KDE doble de los valores arrojados por dicho
    metodo a dichos datos, separados por los datos de prediccion positiva y prediccion negativa
    :param data: Dataframe con los datos
    :param propiedad: string con columna del dataframe
    :param metodo: string con el metodo de aprendizaje automatico
    :return: grafico Doble, grafico de las predicciones positivas, y de las predicciones negativas
    """
    try:
        positivos = data[data[metodo]=='Si']
        positivos = np.array(positivos[propiedad])
        negativos = data[data[metodo]=='No']
        negativos = np.array(negativos[propiedad])

        x_kde = np.linspace(data[propiedad].min(), data[propiedad].max())


        if np.unique(positivos).size > 1:
            kdePositivos = stats.gaussian_kde(positivos.T)
            y_kde_pos = kdePositivos.evaluate(x_kde)
        else:
            kdePositivos= []
            y_kde_pos=[]


        if np.unique(negativos).size > 1:
            kdeNegativos = stats.gaussian_kde(negativos.T)
            y_kde_neg = kdeNegativos.evaluate(x_kde)
        else:
            kdeNegativos=[]
            y_kde_neg =[]

        n_bins = int(1 + 3.3 * np.log10(max(1,len(positivos))))
        hist = np.histogram(positivos,bins = n_bins)
        y_hist_pos = hist[0]/len(positivos)
        x_hist_pos = hist[1]

        n_bins = int(1 + 3.3 * np.log10(max(1,len(negativos))))
        hist = np.histogram(negativos,bins = n_bins)
        y_hist_neg = hist[0]/len(negativos)
        x_hist_neg = hist[1]

        #(x1,y1,etiqueta1, x2, y2,etiqueta2,titulo=None)
        #return graficoDoble(dominio,rangoPositivos,'Antimicrobiano',dominio,rangoNegativos,'No Antimicrobiano',propiedad)

        grafico = dcc.Graph(
            figure=go.Figure(
                data=[
                    go.Scatter(
                        x=x_kde,
                        y=y_kde_pos,
                        name='FDP Antimicrobiano'

                    ),
                    go.Bar(
                        x=x_hist_pos,
                        y=y_hist_pos,
                        marker=go.bar.Marker(
                            color='rgb(51, 113, 255)'
                        ),
                        name='Histograma Antimicrobiano'
                    ),
                    go.Scatter(
                        x=x_kde,
                        y=y_kde_neg,
                        name='FDP No antimicrobiano'

                    ),
                    go.Bar(
                        x=x_hist_neg,
                        y=y_hist_neg,
                        marker=go.bar.Marker(
                            color='rgb(35,208,18)'
                        ),
                        name='Histograma No Antimicrobiano'
                    ),

                ],
                layout=go.Layout(
                    title=TITULOS[propiedad],
                    showlegend=True,
                    legend=go.layout.Legend(
                        x=0,
                        y=1.0
                    ),
                    margin=go.layout.Margin(l=40, r=0, t=40, b=30)
                )
            ),
            style={'width': widthGrafico, 'height': heightGrafico,'margin':margen,'min-width':minWidthGrafico},
        )

        return grafico
    except:
        return html.Div([html.H4('Estimación de Densidad del Núcleo por Propiedad', style={'text-align': 'center'}),
                  html.Div(className="bootstrap-wrapper", children=[

                      html.Div(className="row", children=[

                          html.Div('ANo se puede realizar la estimación de densidad del núcleo',
                                   style={'color': '#FF0000'})

                      ])

                  ])
                  ])


def tablaKDEDoble(clasificacion,value):
    if(len(clasificacion) < 3):
        return html.Div([
            html.Div(className="bootstrap-wrapper", children=[

            html.Div(className="row",children=[

                html.Div('No se puede realizar la estimación de densidad del núcleo', style={'color': '#FF0000'})

            ])

            ])
        ])
    else:
        return html.Div([

            html.Div(className="bootstrap-wrapper", children=[

                dcc.Loading(children=[html.Div(className="row",children=[

                    KDEDoble(clasificacion, 'longitud', value),
                    KDEDoble(clasificacion, 'carga', value),
                    KDEDoble(clasificacion, 'momento hidrofobico', value),


                    KDEDoble(clasificacion, 'porcentaje hidrofobico', value),
                    KDEDoble(clasificacion, 'punto isoelectrico', value),
                    KDEDoble(clasificacion, 'indice de boman', value),

                    KDEDoble(clasificacion, 'polar angles', value),
                    KDEDoble(clasificacion, 'wimley', value),
                    # html.Td(estadistica.KDE(filtro['helices transmembrana']))
                ])
                ], type="default")


        ])

        ])
