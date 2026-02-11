# coding=utf-8
from modlamp.descriptors import GlobalDescriptor
from Bio import SeqIO
import pandas as pd
import numpy as np
import propiedades
import base64
import io
import re
from time import time

def sizeFilter(peptide, tamano_min, tamano_max):
    """
    Verifica si la longitud de una secuencia se encuentra
    dentro del intervalo Dado

    :param peptide: secuencias
    :param tamano_min: longitud minima de la secuencia
    :param tamano_max: longitud maxima de la secuencia
    :return: verdadero o falso
    """
    if (len(peptide) >= tamano_min and len(peptide) <= tamano_max): return True
    else: return False

def aminoacidosFilter(peptide, aminoacidos_excluir_list):
    """
    Verifica si la secuencia no contiene alguno aminoacido
    que se desea excluir

    :param peptide: secuencias
    :pram aminoacidos_excluir_list: lista con aminoacidos que se desean excluir
    :return: verdadero o falso
    """
    for aminoacido in aminoacidos_excluir_list:
        if aminoacido in peptide:
            print("aaaaaaa", aminoacido, peptide)
            return False

    return True


def writeFasta(ids, secuencias, ruta):
    """
    Escribe una lista de secuencias con sus identificadores en formato FASTA

    :param ids: lista de identificadores de las secuencias
    :param secuencias: lista de secuencias
    :param ruta:
    :return:
    """
    idsaux = []
    if (ids == None):
        for i in range(len(secuencias)):
            idsaux.append("")
        ids = idsaux
    escrito = ""
    for i in range(len(secuencias)):
        escrito = escrito + ">" + ids[i] + "\n" + secuencias[i] + "\n"

    with open(ruta, "w") as fasta:
        fasta.write(escrito)


def filtrarPetidos(peptidos):
    """
    Funcion que dado una lista de secuencias de peptidos, filtra aquellos caracteres que no representan
    un aminoacido definido
    :param peptidos: Lista de string de secuencias
    :return: Lista de string de secuencias filtradas
    """
    filtrados = []
    mayusculas = re.compile('\A[A-Z]*\Z')
    #Regex con aminoacidos que no se van a trabajar en la aplicacion
    excluidos = re.compile('\A[^BJOUXZ]*\Z')
    for peptido in peptidos:
        if (mayusculas.match(peptido) != None and excluidos.match(peptido) != None):
            filtrados.append(peptido)

    return filtrados


def filtrarPetido(peptido):
    """
    Funcion que dado un peptido, devuelve verdadero o falso, si los caracteres del peptido cumplen con aminoacidos definidos
    :param peptido: string de la secuencia del peptido
    :return: Booleano
    """
    mayusculas = re.compile('\A[A-Z]*\Z')
    excluidos = re.compile('\A[^BJOUXZ]*\Z')
    if (mayusculas.match(peptido) != None and excluidos.match(peptido) != None):
        return True
    else:
        return False


TERMINACIONES_FASTA = ['.fasta','.fa','.fas','.fna','.faa','.frn','.ffn','.txt']

#FUNCION CARGA FASTA
def intersecciones(TERMINACIONES_FASTA, filename):
    for fasta in  TERMINACIONES_FASTA:
        if fasta in filename:
            return True
    return False


def parse_contents(contents, filename, date):
    if contents == None:
        cargaFasta  =False
        return '', cargaFasta
    content_type, content_string = contents.split(',')

    decoded = base64.b64decode(content_string)
    try:
        if intersecciones(TERMINACIONES_FASTA,filename):
            # Assume that the user uploaded a fasta file
            FASTA = []
            longitud = []
            origen = []
            for seq_record in SeqIO.parse(io.StringIO(decoded.decode('utf-8')), 'fasta'):
                # print(seq_record.id)
                # print(repr(seq_record.seq))
                # print(seq_record.seq)
                # print()
                FASTA.append(str(seq_record.seq))
                longitud.append(len(seq_record))
                origen.append(str(seq_record.name))

                # print(len(seq_record))

            data = pd.DataFrame({
                'origen':origen,
                'longitud':longitud,
                'secuencia':FASTA})
            data.to_csv('secuencias.csv',index = False,sep=",")
            cargaFasta = True
            return filename,cargaFasta
        else:
            cargaFasta = False
            return "Debe ser un archivo tipo FASTA",cargaFasta


    except Exception as e:
        print(e)
        cargaFasta = False
        return "Hubo un error procesando el archivo"




def cargar_secuencias(secuencias):
    secuencias = secuencias.split(',')
    longitudes = []
    origenes = []
    for secuencia in secuencias:
        longitudes.append(len(secuencia))
        origenes.append('No definido')
    data = pd.DataFrame({
        'origen':origenes ,
        'longitud': longitudes,
        'secuencia': secuencias})
    return data


def cargar_fasta(contents, filename, date):
    dataVacio = pd.DataFrame({
        'origen': [],
        'longitud': [],
        'secuencia': []})


    if contents == None:
        cargaFasta=False
        print("vacio")
        return dataVacio,'Archivo no contiene información', cargaFasta



    content_type, content_string = contents.split(',')

    decoded = base64.b64decode(content_string)
    try:
        if intersecciones(TERMINACIONES_FASTA,filename):
            # Assume that the user uploaded a fasta file
            FASTA = []
            longitud = []
            origen = []
            for seq_record in SeqIO.parse(io.StringIO(decoded.decode('utf-8')), 'fasta'):
                # print(seq_record.id)
                # print(repr(seq_record.seq))
                # print(seq_record.seq)
                # print()
                if(filtrarPetido(str(seq_record.seq))):
                    FASTA.append(str(seq_record.seq))
                    longitud.append(len(seq_record))
                    origen.append(str(seq_record.name))

                # print(len(seq_record))

            data = pd.DataFrame({
                'origen':origen,
                'longitud':longitud,
                'secuencia':FASTA})
            #data.to_csv('secuencias.csv',index = False,sep=",")
            cargaFasta = True
            print("correcto")
            return data,filename,cargaFasta


        else:
            cargaFasta = False
            print("no tipo fasta")
            return dataVacio,"Debe ser un archivo tipo FASTA",cargaFasta





    except Exception as e:
        print(e)
        cargaFasta = False
        print("archivo vacio")
        return dataVacio,"Hubo un error procesando el archivo",False






def validarIntervalos(entradaMinimo,entradaMaximo,minimo,maximo):
    """
    Funcion que dado dos entradas y un intervalo, devuelve verdadero si dichas entradas se encuentra en el intervalo,
    de lo contrario devuelve falso.

    :param entradaMinimo: float o int de la entrada minima
    :param entradaMaximo: float o int de la entrada maxima
    :param minimo: float o int del valor minimo del intervalo
    :param maximo: float o int del valor maximo del intervalo
    :return: Boolean con el valor de verdad si las entrada pertenecen al intervalo
    """
    if (entradaMinimo == None and entradaMaximo == None):
        return False
    if(entradaMinimo != None and entradaMaximo != None):

        if (entradaMaximo < entradaMinimo or entradaMinimo < minimo or entradaMinimo> maximo or entradaMaximo < minimo or entradaMaximo > maximo):
            return True
    if(entradaMinimo != None):
        if (entradaMinimo < minimo or entradaMinimo> maximo ):
            return True

    if(entradaMaximo != None):
        if (entradaMaximo < minimo or entradaMaximo > maximo):
            return True
    if(entradaMinimo == None and entradaMaximo == None):
        return True


def interseccion(string1, string2):
    """
    Dado dos objetos iteradores devuelve verdadero si encuentra almenos una coincidencia
    :param string1: String a
    :param string2: String b
    :return: Boolean con el valor de verdad
    """
    for i in string1:
        if (i in string2):
            return True

    return False


def busqueda(cadena, minimo, maximo, excluir):
    """
    Busca los pedazos de cadena de longitudes de minimo a maximo, excluyendo las letras
    dadas en excluir

    :param cadena: String donde se va a buscar
    :param minimo: Int de longitud minima de busqueda
    :param maximo: Int de longitud maxima de busqueda
    :param excluir: String con caracteres que se van a excluir
    :return: Lista de string encontrados,lista de posiciones del string en la cadena original,longitud del string encontrado
    """
    posiciones = []
    combinaciones = []
    longitudes = []
    aux = []

    for i in range(minimo, maximo + 1):
        for j in range(len(cadena) - i + 1):
            aux = cadena[j:j + i]

            if (j == j + i - 1):
                posicion = str(j + 1)
            else:
                posicion = str(j + 1) + "-" + str(j + i)

            if (not interseccion(aux, excluir)):
                combinaciones.append(aux)
                posiciones.append(posicion)
                longitudes.append(len(aux))

    return combinaciones,posiciones,longitudes



CARGAMINIMA=-10
CARGAMAXIMA=10

HIDROFOBICOMINIMO=0
HIDROFOBICOMAXIMO=100

MOMENTOMINIMO=-1000
MOMENTOMAXIMO=1000

ISOELECTRICOMINIMO= 1
ISOELECTRICOMAXIMO= 14

BOMANMINIMO = -4.92
BOMANMAXIMO = 14.92

WIMLEYMINIMO = -1000
WIMLEYMAXIMO = 1000


APLICARHELICES = False
METODOHELICES = 'eisenberg'
ESCALA = 'kyte'

EXCLUIRDEFECTO = []

POLARANGLESMINIMO = 0
POLARANGLESMAXIMO = 360


def filtro(fasta,
            longitudMinima, longitudMaxima,
            cargaMinima=None, cargaMaxima=None,
            momentoMinimo=None, momentoMaximo=None,
            isoelectricoMinimo=None, isoelectricoMaximo=None,
            hidrofobicoMinimo=None, hidrofobicoMaximo=None,
            bomanMinimo=None, bomanMaximo=None,
            escala=None,
            aplicarHelices=None,
            metodoHelices=None,
            wimleyMinimo=None, wimleyMaximo=None,
            polarAnglesMinimo=None, polarAnglesMaximo=None,

            excluir=None):
    """
    Funcion que busca secuencias parciales de peptidos en secuencias de proteinas, calcula propiedades, y filtra
    por rango de propiedades

    :param fasta: CSV de una columna con las aminoacidos listados
    :param longitudMinima: Parametro Obligatorio, tamaño minimo de peptidos a buscar
    :param longitudMaxima: Parametro Obligatorio, tamaño maximo de peptidos a buscar
    los demas argumentos son los filtros a las propiedades
    :param excluir: Aminoacidos a excluir en la busqueda
    :return: Dataframe con las secuencias encontradas
    """

    # Parametros por defecto
    if (cargaMinima == None):
        cargaMinima = CARGAMINIMA
    if (cargaMaxima == None):
        cargaMaxima = CARGAMAXIMA

    if (momentoMinimo == None):
        momentoMinimo = MOMENTOMINIMO
    if (momentoMaximo == None):
        momentoMaximo = MOMENTOMAXIMO

    if (isoelectricoMinimo == None):
        isoelectricoMinimo = ISOELECTRICOMINIMO
    if (isoelectricoMaximo == None):
        isoelectricoMaximo = ISOELECTRICOMAXIMO

    if (hidrofobicoMinimo == None):
        hidrofobicoMinimo = HIDROFOBICOMINIMO
    if (hidrofobicoMaximo == None):
        hidrofobicoMaximo = HIDROFOBICOMAXIMO

    if (bomanMinimo == None):
        bomanMinimo = BOMANMINIMO
    if (bomanMaximo == None):
        bomanMaximo = BOMANMAXIMO

    if (escala == None):
        escala = ESCALA

    if (excluir == None):
        excluir = EXCLUIRDEFECTO

    if (aplicarHelices == None):
        aplicarHelices = APLICARHELICES
    if (aplicarHelices == 'si'):
        aplicarHelices = True
    if (aplicarHelices == 'no'):
        aplicarHelices = False

    if (metodoHelices == None):
        metodoHelices = METODOHELICES

    if (wimleyMinimo == None):
        wimleyMinimo = WIMLEYMINIMO

    if (wimleyMaximo == None):
        wimleyMaximo = WIMLEYMAXIMO
    if (polarAnglesMinimo == None):
        polarAnglesMinimo = POLARANGLESMINIMO
    if (polarAnglesMaximo == None):
        polarAnglesMaximo = POLARANGLESMAXIMO

    secuencias = list(fasta['secuencia'])
    fastaOrigen = list(fasta['origen'])
    fastaLongitud = list(fasta['longitud'])
    peptidos = []
    origen = []
    longitudOrigen = []
    longitudes = []
    posiciones = []
    contador = 0

    t1 = time()

    # BUSCA POR SECUENCIA LOS PEPTIDOS
    for secuencia in secuencias:
        auxbusqueda, auxposiciones, auxlongitudes = busqueda(secuencia, longitudMinima, longitudMaxima, excluir)

        # se halla el origen y se multiplica para que en todas las encontradas aparezca el origen
        origen = origen + [fastaOrigen[contador]] * len(auxbusqueda)
        # se hala la longitud del origen y se multiplica para que aparezca en todas las encontradas
        longitudOrigen = longitudOrigen + [fastaLongitud[contador]] * len(auxbusqueda)
        posiciones = posiciones + auxposiciones
        peptidos = peptidos + auxbusqueda
        longitudes = longitudes + auxlongitudes

    t2 = time()
    print("Tiempo busqueda de secuencias:", t2 - t1)

    if (len(peptidos) == 0):
        return pd.DataFrame()

    dataPeptidos = pd.DataFrame(
        {
            'origen': origen,
            'longitud Origen': longitudOrigen,
            'posicion': posiciones,
            'longitud': longitudes,
            'peptido': peptidos
        })

    t0 = time()
    # wimley_white
    wimleyWhite = propiedades.wimley_white(peptidos)
    t1 = time()
    print("Tiempo calculo de wimleywhite:", t1 - t0)
    dataPeptidos['wimley'] = wimleyWhite
    dataPeptidos = dataPeptidos[dataPeptidos['wimley'] >= wimleyMinimo]
    dataPeptidos = dataPeptidos[dataPeptidos['wimley'] <= wimleyMaximo]
    if (len(dataPeptidos['peptido']) == 0):
        return pd.DataFrame()

    # Indice de Boman
    indiceBoman = propiedades.indice_boman(list(dataPeptidos['peptido']))
    t2 = time()
    print("Tiempo calculo de boman:", t2 - t1)
    dataPeptidos['indice de boman'] = indiceBoman
    dataPeptidos = dataPeptidos[dataPeptidos['indice de boman'] >= bomanMinimo]
    dataPeptidos = dataPeptidos[dataPeptidos['indice de boman'] <= bomanMaximo]
    if (len(dataPeptidos['peptido']) == 0):
        return pd.DataFrame()

    # Porcentaje hidrofobico
    porcentajeHidrofobico = list(
        np.array(propiedades.porcentaje_hidrofobico(list(dataPeptidos['peptido']), escala)).round(2))
    t3 = time()
    print("Tiempo calculo de porcerntaje hidrofobic:", t3 - t2)
    dataPeptidos['porcentaje hidrofobico'] = porcentajeHidrofobico
    dataPeptidos = dataPeptidos[dataPeptidos['porcentaje hidrofobico'] >= hidrofobicoMinimo]
    dataPeptidos = dataPeptidos[dataPeptidos['porcentaje hidrofobico'] <= hidrofobicoMaximo]
    if (len(dataPeptidos['peptido']) == 0):
        return pd.DataFrame()

    # Momento hidrofobico
    # descriptor = PeptideDescriptor(peptidos,'kytedoolittle')
    # descriptor.calculate_moment()
    # momentoHidrofobico = descriptor.descriptor
    # momentoHidrofobico = list(np.array(momentoHidrofobico.ravel()).round(2))

    momentoHidrofobico = list(
        np.array(propiedades.momento_hidrofobico_relativo(list(dataPeptidos['peptido']))).round(2))
    t4 = time()
    print("Tiempo calculo de momento:", t4 - t3)
    dataPeptidos['momento hidrofobico'] = momentoHidrofobico
    dataPeptidos = dataPeptidos[dataPeptidos['momento hidrofobico'] >= momentoMinimo]
    dataPeptidos = dataPeptidos[dataPeptidos['momento hidrofobico'] <= momentoMaximo]
    if (len(dataPeptidos['peptido']) == 0):
        return pd.DataFrame()

    # Helices Transmembrana
    helicesTransmembrana = propiedades.helices_transmembrana(list(dataPeptidos['peptido']), metodoHelices)
    dataPeptidos['helices transmembrana'] = helicesTransmembrana
    t5 = time()
    print("Tiempo de helices trasnmembrana:", t5 - t4)
    if (aplicarHelices):
        dataPeptidos = dataPeptidos[dataPeptidos['helices transmembrana'] == 'si']
    if (len(dataPeptidos['peptido']) == 0):
        return pd.DataFrame()

    # Carga
    carga = propiedades.carga(list(dataPeptidos['peptido']))
    t6 = time()
    print("Tiempo calculo de carga:", t6 - t5)
    dataPeptidos['carga'] = carga
    dataPeptidos = dataPeptidos[dataPeptidos['carga'] >= cargaMinima]
    dataPeptidos = dataPeptidos[dataPeptidos['carga'] <= cargaMaxima]
    if (len(dataPeptidos['peptido']) == 0):
        return pd.DataFrame()

    # Polar Angles
    polarAngles = propiedades.polar_angles(list(dataPeptidos['peptido']))
    t7 = time()
    print("Tiempo calculo de polarangles:", t7 - t6)
    dataPeptidos['polar angles'] = polarAngles
    dataPeptidos = dataPeptidos[dataPeptidos['polar angles'] >= polarAnglesMinimo]
    dataPeptidos = dataPeptidos[dataPeptidos['polar angles'] <= polarAnglesMaximo]
    if (len(dataPeptidos['peptido']) == 0):
        return pd.DataFrame()

    # Punto Isoelectrico
    glob = GlobalDescriptor(list(dataPeptidos['peptido']))
    glob.isoelectric_point()
    puntoIsoelectrico = glob.descriptor
    puntoIsoelectrico = list(np.array(puntoIsoelectrico.ravel()).round(2))
    dataPeptidos['punto isoelectrico'] = puntoIsoelectrico

    # puntoIsoelectrico = propiedades.punto_isoelectrico(dataPeptidos['peptido'])

    # puntoIsoelectrico = propiedades.punto_isoelectricoManual(peptidos)
    t8 = time()
    print("Tiempo calculo de isoelectrico:", t8 - t7)

    dataPeptidos = dataPeptidos[dataPeptidos['punto isoelectrico'] >= isoelectricoMinimo]
    dataPeptidos = dataPeptidos[dataPeptidos['punto isoelectrico'] <= isoelectricoMaximo]
    if (len(dataPeptidos['peptido']) == 0):
        return pd.DataFrame()

    t9 = time()
    print("Tiempo filtrado de propiedades:", t9 - t1)
    return dataPeptidos
