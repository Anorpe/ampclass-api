# -*- coding: utf-8 -*-

import warnings
warnings.filterwarnings('ignore')
from collections import Counter



#from biovec.models import prot_vec
import collections

import pandas as pd

from modlamp.descriptors import GlobalDescriptor
from Bio.SeqUtils.ProtParam import ProteinAnalysis

#Se importa la librería pydpi de Python 2.
#Para que esto funcione, se debe crear una carpeta "tmp" en la raíz del disco donde se esté corriendo el programa.
# from past import autotranslate
# autotranslate(['pydpi'])

import time
from pydpi.protein import AAComposition
from pydpi.protein import CTD
from pydpi.protein import Autocorrelation  # 109, 1023, 1024
from pydpi.protein import QuasiSequenceOrder # 468, 470
from pydpi.protein import PseudoAAC
#from sklearn.externals.joblib import load
from joblib import load
import numpy as np
import propiedades
import re


CARAC = [
    'GTPC_positivecharge.aliphatic.aliphatic',
    'GTPC_aliphatic.aliphatic.positivecharge',
    '_SolventAccessibilityC2',
    'QSOgrant35',
    'MoreauBrotoAuto_Mutability14',
    'QSOgrant9',
    '_ChargeD3100',
    'E',
    'GearyAuto_AvFlexibility14',
    '_PolarityD1001',
    '_SecondaryStrD1050',
    'CTriad_g5.g5.g2',
    'isoelectric_point',
    '_HydrophobicityD3025',
    'MoreauBrotoAuto_Steric18',
    'MoranAuto_Hydrophobicity2',
    '_SecondaryStrD3075',
    'CTriad_g2.g2.g5',
    'taugrant25',
    'taugrant28',
    'MoranAuto_Hydrophobicity11',
    'GearyAuto_Hydrophobicity25',
    '_SolventAccessibilityD3050',
    'taugrant2',
    'GearyAuto_Hydrophobicity26',
    'GearyAuto_Hydrophobicity30',
    'MoranAuto_ResidueASA14',
    'tausw15',
    '_NormalizedVDWVD1050',
    '_SecondaryStrD1001',
    'L',
    'QSOgrant22',
    'MoreauBrotoAuto_Mutability2',
    'CTriad_g2.g5.g5',
    'GearyAuto_Hydrophobicity29',
    '_PolarityD3075',
    '_NormalizedVDWVD1025',
    '_HydrophobicityD2001',
    'MoreauBrotoAuto_Mutability18',
    'CTriad_g5.g1.g7',
    'GAAC_uncharged',
    'taugrant13',
    'QSOgrant24',
    'R',
    '_PolarizabilityC3',
    'GearyAuto_Hydrophobicity4',
    'QSOgrant31',
    'QSOgrant12',
    '_ChargeC3',
    '_SecondaryStrC1',
    '_SolventAccessibilityC3',
    'A',
    'taugrant17',
    'molecular_weight',
    '_PolarityD2050',
    'GearyAuto_Mutability1',
    'MoreauBrotoAuto_Mutability7',
    'taugrant14',
    '_PolarizabilityD3001',
    '_PolarizabilityT23',
    'taugrant30',
    '_SolventAccessibilityC1',
    'MoreauBrotoAuto_ResidueVol12',
    '_NormalizedVDWVD1001',
    'GDPC_negativecharge.aliphatic',
    'AA',
    'MoreauBrotoAuto_Mutability3',
    'GearyAuto_Mutability9',
    'MoreauBrotoAuto_Mutability1',
    'MoreauBrotoAuto_Mutability12',
    '_NormalizedVDWVC3',
    'taugrant27',
    'MoreauBrotoAuto_AvFlexibility12',
    'MoreauBrotoAuto_Polarizability12',
    'charge',
    'MoreauBrotoAuto_Steric15',
    'GearyAuto_Hydrophobicity27',
    'QSOSW22',
    '_HydrophobicityD3001',
    'CTriad_g2.g5.g2',
    '_PolarityC3',
    '_SolventAccessibilityT12',
    'taugrant29',
    'GearyAuto_Hydrophobicity28',
    '_SecondaryStrD1075',
    'QSOgrant10',
    'tausw17',
    'G',
    '_PolarizabilityD3075',
    'charge_density',
    'GearyAuto_AvFlexibility10',
    'P',
    'MoreauBrotoAuto_Mutability10',
    '_PolarityT13',
    'CTriad_g2.g5.g3',
    'MoreauBrotoAuto_Mutability13',
    '_SolventAccessibilityD3075',
    '_PolarityD2025',
    'MoreauBrotoAuto_ResidueVol16',
    '_HydrophobicityT13',
    '_SecondaryStrD1100',
    'K',
    '_ChargeD3075',
    'aromaticity',
    '_SolventAccessibilityD1001',
    'I',
    'GearyAuto_Hydrophobicity24',
    'GearyAuto_FreeEnergy10',
    'GearyAuto_Hydrophobicity15',
    'MoreauBrotoAuto_Steric11',
    'PR',
    '_PolarizabilityD3025',
    'MoranAuto_FreeEnergy1',
    'QSOgrant8',
    'QSOgrant26',
    'KL',
    '_SolventAccessibilityD3001',
    'GDPC_positivecharge.positivecharge',
    'QSOSW5',
    'MoreauBrotoAuto_Steric13',
    'GTPC_aliphatic.positivecharge.positivecharge',
    '_SecondaryStrD3001',
    'MoranAuto_Hydrophobicity4',
    'MoreauBrotoAuto_AvFlexibility13',
    'W',
    '_SecondaryStrD2001',
    '_ChargeT23',
    '_ChargeC1',
    '_SolventAccessibilityD3025',
    'GearyAuto_Mutability4',
    'C',
    'M',
    'H',
    '_ChargeD1001',
    'V',
    'GDPC_aliphatic.positivecharge',
    'tausw13',
    '_ChargeD2001',
    'boman_index',
    'TW',
    'TF',
    'MoreauBrotoAuto_ResidueVol17',
    '_SolventAccessibilityD3100',
    '_NormalizedVDWVD3001',
    'QSOSW12',
    'AR',
    'MoreauBrotoAuto_ResidueVol18',
    '_SecondaryStrD3100',
    '_ChargeD3001',
    'GearyAuto_Mutability3',
    'MoreauBrotoAuto_Steric19',
    '_SecondaryStrD3050',
    'MoreauBrotoAuto_Polarizability14',
    'T',
    'MoreauBrotoAuto_Steric16',
    'MoreauBrotoAuto_Hydrophobicity12',
    'GearyAuto_Hydrophobicity17',
    'MoreauBrotoAuto_Mutability5',
    'D',
    'AN',
    'MoreauBrotoAuto_Mutability15',
    'MoreauBrotoAuto_ResidueVol13',
    'MoreauBrotoAuto_Mutability4',
    '_NormalizedVDWVD3050',
    'GearyAuto_Mutability5',
    '_PolarityD2075',
    'MoreauBrotoAuto_Hydrophobicity18',
    'taugrant15',
    'F',
    'GAAC_aromatic',
    'QSOgrant1',
    'MoranAuto_FreeEnergy6',
    'GAAC_aliphatic',
    'GAAC_positivecharge',
    'MoreauBrotoAuto_Polarizability13',
    'MoreauBrotoAuto_Mutability11',
    'taugrant5',
    'QSOgrant14',
    'GearyAuto_ResidueVol7',
    'MoreauBrotoAuto_Hydrophobicity14',
    'hydrophobic_ratio',
    'GearyAuto_Mutability6',
    'GearyAuto_AvFlexibility15',
    'QSOgrant13',
    '_ChargeD1075',
    'MoreauBrotoAuto_Hydrophobicity13',
    '_PolarizabilityD1001',
    'QSOSW26',
    'MoranAuto_Hydrophobicity7',
    'QSOSW25',
    '_SecondaryStrD3025',
    'Y',
    'tausw2',
    'KK',
    'Q',
    'GearyAuto_Mutability8',
    '_NormalizedVDWVD3075',
    'GearyAuto_Mutability7',
    'MoreauBrotoAuto_Mutability8',
    'QSOSW31',
    'AD',
    '_ChargeD1100',
    'GDPC_positivecharge.aliphatic',
    'MoreauBrotoAuto_Steric17',
    'QSOgrant25',
    'GearyAuto_AvFlexibility18',
    'N',
    'QSOSW13',
    'tausw14',
    'GearyAuto_Mutability10',
    'GAAC_negativecharge',
    'GTPC_positivecharge.positivecharge.aliphatic',
    'LK',
    '_ChargeD3025',
    'taugrant24',
    'instability_index',
    'MoreauBrotoAuto_ResidueVol14',
    'MoreauBrotoAuto_Mutability16',
    'QSOSW35',
    '_PolarityD2100',
    '_PolarityD3100',
    'aliphatic_index',
    'MoreauBrotoAuto_Mutability6',
    'QSOgrant2',
    'QSOgrant6',
    '_PolarityD2001',
    'MoreauBrotoAuto_Polarizability17',
    'GearyAuto_Hydrophobicity7',
    'GearyAuto_Mutability2',
    '_PolarizabilityD1025',
    'gravy',
    'TP',
    '_HydrophobicityC3',
    'MoreauBrotoAuto_Steric12',
    'S',
    'MoreauBrotoAuto_AvFlexibility14',
    'taugrant26',
    '_ChargeD3050',
    '_ChargeT12',
    'MoreauBrotoAuto_Steric14',
    '_NormalizedVDWVD3025',
    '_HydrophobicityC1',
    'MoreauBrotoAuto_FreeEnergy13',
    'MoreauBrotoAuto_ResidueASA14'
]

CARAC_MERA = [
    'charge',
 'isoelectric_point',
 'R',
 'D',
 'C',
 'E',
 'L',
 'K',
 'M',
 'GAAC_positivecharge',
 'GAAC_negativecharge',
 'GAAC_uncharged',
 'GDPC_aliphatic.positivecharge',
 'GDPC_aliphatic.uncharged',
 'GDPC_positivecharge.aliphatic',
 'GDPC_uncharged.uncharged',
 'GTPC_aliphatic.aliphatic.positivecharge',
 'GTPC_positivecharge.aliphatic.aliphatic',
 '_PolarizabilityC1',
 '_PolarizabilityC3',
 '_SolventAccessibilityC3',
 '_SecondaryStrC3',
 '_ChargeC2',
 '_PolarityC2',
 '_NormalizedVDWVC1',
 '_HydrophobicityC2',
 '_PolarizabilityT12',
 '_SolventAccessibilityT12',
 '_SolventAccessibilityT13',
 '_SecondaryStrT23',
 '_PolarityT13',
 '_NormalizedVDWVT12',
 '_HydrophobicityT13',
 '_PolarizabilityD3001',
 '_PolarizabilityD3050',
 '_SolventAccessibilityD1001',
 '_SolventAccessibilityD1050',
 '_SolventAccessibilityD3001',
 '_SolventAccessibilityD3100',
 '_SecondaryStrD1001',
 '_SecondaryStrD2001',
 '_SecondaryStrD2100',
 '_SecondaryStrD3001',
 '_SecondaryStrD3075',
 '_SecondaryStrD3100',
 '_ChargeD1001',
 '_ChargeD1075',
 '_ChargeD1100',
 '_ChargeD2001',
 '_ChargeD3100',
 '_PolarityD1001',
 '_PolarityD2100',
 '_NormalizedVDWVD2050',
 '_HydrophobicityD2050',
 'MoreauBrotoAuto_Steric21',
 'MoreauBrotoAuto_Steric22',
 'MoreauBrotoAuto_Mutability1',
 'MoranAuto_Hydrophobicity2',
 'MoranAuto_Hydrophobicity4',
 'MoranAuto_Hydrophobicity7',
 'MoranAuto_Hydrophobicity11',
 'GearyAuto_Hydrophobicity7',
 'GearyAuto_Hydrophobicity11',
 'GearyAuto_Steric8',
 'GearyAuto_Steric9',
 'GearyAuto_Mutability3',
 'GearyAuto_Mutability5',
 'GearyAuto_Mutability11',
 'tausw7',
 'tausw10',
 'QSOgrant2',
 'QSOgrant6',
 'QSOgrant13',
 'QSOgrant31']

scaler = load('models/scaler.pkl')
robust = load('models/robust.pkl')
## Regresión logistica deprecada

regresionLogistica = load('models/regresion_logistica.pkl')
arbolDecision = load('models/arbol_decision.pkl')
bosqueAleatorio = load('models/bosque_aleatorio.pkl')
xgboost = load('models/xgbc.pkl')
#Organizar Red Neuronal
redNeuronal = load('models/xgbc.pkl')





UMBRAL = 50


def filtrarPeptidos(peptidos):
    """
    Funcion que dado una lista de secuencias de peptidos, filtra aquellos caracteres que no representan
    un aminoacido definido
    :param peptidos: Lista de string de secuencias
    :return: Lista de string de secuencias filtradas
    """

    filtrados = []
    mayusculas = re.compile(r"\A[A-Z]*\Z")
    excluidos = re.compile(r"\A[^JBOXZU]*\Z")

    for peptido in peptidos:
        if (mayusculas.match(peptido) != None and excluidos.match(peptido) != None and len(peptido)>=7 and len(peptido)<=35):
            filtrados.append(peptido)

    return filtrados


def predecir(secuencias,umbral=None,origenes=None,longitudes=None,posiciones=None):
    """
    Dado una lista de secuencias de peptido genera las predicciones y probabilidades de cada una de ser
    un peptido antimicrobiano
    :param secuencias: Lista de secuencias
    :param umbral: Umbral de clasificacion
    :param origenes: LIsta de origenes de las secuencias
    :param longitudes: Lista de longitudes de las secuencias
    :param posiciones: Lista de posiciones en el origen de las secuencias
    :return: Dataframe con predicciones
    """
    #if origenes is None:
    secuencias = filtrarPeptidos(secuencias)

    if umbral == None:
        umbral = UMBRAL

    probabilidadesLogistica = []
    probabilidadesRed = []
    probabilidadesArbol = []
    probabilidadesBosque = []
    probabilidadesXgboost = []


    peptidos = []
    for secuencia in secuencias:
        print("Calculando propiedades de la secuencia:",secuencia)
        if (len(secuencia) >= 7):


            descriptores = get_features(secuencia)

            descriptores = pd.DataFrame({clave: [valor] for clave, valor in descriptores.items()})

            descriptores1 = descriptores[CARAC_MERA]
            descriptores2 = descriptores[CARAC]

            descriptores1 = robust.transform(descriptores1)
            descriptores2 = scaler.transform(descriptores2)



            prediccionLogistica = int(round(regresionLogistica.predict_proba(descriptores1)[0][1] * 100, 2))

            prediccionRed = int(round(redNeuronal.predict_proba(descriptores2)[0][1] * 100, 2))
            prediccionArbol = int(round(arbolDecision.predict_proba(descriptores2)[0][1] * 100, 2))
            prediccionBosque = int(round(bosqueAleatorio.predict_proba(descriptores2)[0][1] * 100, 2))
            prediccionXgboost = int(round(xgboost.predict_proba(descriptores2)[0][1] * 100, 2))




            probabilidadesLogistica.append(prediccionLogistica)
            probabilidadesRed.append(prediccionRed)
            probabilidadesArbol.append(prediccionArbol)
            probabilidadesBosque.append(prediccionBosque)
            probabilidadesXgboost.append(prediccionXgboost)
            peptidos.append(secuencia)
            # except:
            #     probabilidadesRed.append(100)
            #     probabilidadesLogistica.append(100)
            #     probabilidadesXgboost.append(100)
            #     probabilidadesBosque.append(100)
            #     probabilidadesArbol.append(100)
            #     peptidos.append(secuencia)
            #     print("Hubo un error con el peptido: ", secuencia)
            #     continue

    if len(peptidos) == 0:
        return []

    #CLASIFICACIONES SEGUN UMBRAL
    clasificacionLogistica = umbral_clasificacion(probabilidadesLogistica,umbral)
    clasificacionRed = umbral_clasificacion(probabilidadesRed,umbral)

    clasificacionArbol = umbral_clasificacion(probabilidadesArbol,umbral)
    clasificacionBosque = umbral_clasificacion(probabilidadesBosque,umbral)
    clasificacionXgboost = umbral_clasificacion(probabilidadesXgboost,umbral)




    #longitud
    longitud = propiedades.longitud(peptidos)

    # Momento hidrofobico
    momentoHidrofobico = list(np.array(propiedades.momento_hidrofobico_relativo(peptidos)).round(2))

    # Porcentaje hidrofobico
    porcentajeHidrofobico = list(np.array(propiedades.porcentaje_hidrofobico(peptidos, 'kyte')).round(2))

    # Carga
    carga = propiedades.carga(peptidos)

    # Punto Isoelectrico
    puntoIsoelectrico = propiedades.punto_isoelectricoManual(peptidos)

    # Indice de Boman
    indiceBoman = propiedades.indice_boman(peptidos)

    # Polar Angles
    polarAngles = propiedades.polar_angles(peptidos)

    # wimley_white
    wimleyWhite = propiedades.wimley_white(peptidos)

    # Helices Transmembrana
    helicesTransmembrana = propiedades.helices_transmembrana(peptidos, 'eisenberg')
    dataClasificacion = pd.DataFrame({
        'Peptido': peptidos,
        'Regresión Lógistica': probabilidadesLogistica,
        'Red Neuronal': probabilidadesRed,
        'Arbol de Decisión': probabilidadesArbol,
        'Bosque Aleatorio': probabilidadesBosque,
        'XGboost': probabilidadesXgboost,

        'Clasificación Regresión Lógistica': clasificacionLogistica,
        'Clasificación Red Neuronal': clasificacionRed,
        'Clasificación Arbol de Decisión': clasificacionArbol,
        'Clasificación Bosque Aleatorio': clasificacionBosque,
        'Clasificación XGboost':clasificacionXgboost,

        'longitud': longitud,
        'carga': carga,
        'momento hidrofobico': momentoHidrofobico,
        'porcentaje hidrofobico': porcentajeHidrofobico,
        'punto isoelectrico': puntoIsoelectrico,
        'indice de boman': indiceBoman,
        'polar angles': polarAngles,
        'wimley': wimleyWhite,
        'helices transmembrana': helicesTransmembrana
    })
    print("Clasificación:",dataClasificacion.head())
    return dataClasificacion



def umbral_clasificacion(probabilidades, umbral):
    """
    Dado la prediccion de la probabilidad y el umbral de clasificacion, genera la prediccion
    :param probabilidades: Lista de la prediccion de la probabilidades
    :param umbral: Int o float de umbral de clasificacion
    :return: lista con las clasificaciones ('Si' o 'No')
    """
    clasificacion = []
    for probabilidad in probabilidades:
        if probabilidad >= umbral:
            clasificacion.append('Si')
        else:
            clasificacion.append('No')
    return clasificacion


def get_GAAC(seq):
    group = {
        'GAAC_aliphatic': 'GAVLMI',
        'GAAC_aromatic': 'FYW',
        'GAAC_positivecharge': 'KRH',
        'GAAC_negativecharge': 'DE',
        'GAAC_uncharged': 'STCPNQ'
    }

    fts = {}

    groupKey = group.keys()
    count = Counter(seq)
    myDict = {}

    for key in groupKey:
        for aa in group[key]:
            myDict[key] = myDict.get(key, 0) + count[aa]
    for key in groupKey:
        fts[key] = (myDict[key] / len(seq))
    return fts


def get_GDPC(seq):
    group = {
        'aliphatic': 'GAVLMI',
        'aromatic': 'FYW',
        'positivecharge': 'KRH',
        'negativecharge': 'DE',
        'uncharged': 'STCPNQ'
    }

    groupKey = group.keys()
    baseNum = len(groupKey)
    dipeptide = [g1 + '.' + g2 for g1 in groupKey for g2 in groupKey]

    index = {}
    for key in groupKey:
        for aa in group[key]:
            index[aa] = key

    fts = {}

    myDict = {}
    for t in dipeptide:
        myDict[t] = 0

    sum = 0
    for j in range(len(seq) - 2 + 1):
        myDict[index[seq[j]] + '.' + index[seq[j + 1]]] = myDict[index[seq[j]] + '.' + index[
            seq[j + 1]]] + 1
        sum = sum + 1

    if sum == 0:
        for t in dipeptide:
            fts['GDPC_' + t] = 0
            # code.append(0)
    else:
        for t in dipeptide:
            fts['GDPC_' + t] = (myDict[t] / sum)

    return fts


def get_GTPC(seq):
    group = {
        'aliphatic': 'GAVLMI',
        'aromatic': 'FYW',
        'positivecharge': 'KRH',
        'negativecharge': 'DE',
        'uncharged': 'STCPNQ'
    }

    groupKey = group.keys()
    baseNum = len(groupKey)
    triple = [g1 + '.' + g2 + '.' + g3 for g1 in groupKey for g2 in groupKey for g3 in groupKey]

    index = {}
    for key in groupKey:
        for aa in group[key]:
            index[aa] = key

    fts = {}
    myDict = {}
    for t in triple:
        myDict[t] = 0
    sum = 0
    for j in range(len(seq) - 3 + 1):
        myDict[index[seq[j]] + '.' + index[seq[j + 1]] + '.' + index[seq[j + 2]]] = myDict[index[seq[j]] + '.' + index[
            seq[j + 1]] + '.' + index[seq[j + 2]]] + 1
        sum = sum + 1
    if sum == 0:
        for t in triple:
            fts['GTPC_' + t] = 0
    else:
        for t in triple:
            fts['GTPC_' + t] = (myDict[t] / sum)

    return fts


def get_grouped_aa_features(features):
    gaac = get_GAAC(features['sequence'])
    gdpc = get_GDPC(features['sequence'])
    gdtp = get_GTPC(features['sequence'])
    features.update(gaac)
    features.update(gdpc)
    features.update(gdtp)
    return features


def get_global_features(features):
    seq = features['sequence']
    desc = GlobalDescriptor(seq)
    desc.calculate_MW(amide=False)
    features['molecular_weight'] = desc.descriptor[0][0]
    desc.calculate_charge(ph=7.0)
    features['charge'] = desc.descriptor[0][0]
    desc.charge_density(ph=7.0)
    features['charge_density'] = desc.descriptor[0][0]
    # desc.isoelectric_point(amide=False)
    # features['isoelectric_point'] = desc.descriptor[0][0]
    biop_analysis = ProteinAnalysis(seq)
    features['isoelectric_point'] = biop_analysis.isoelectric_point()
    # features['flexibility'] = biop_analysis.flexibility()
    features['gravy'] = biop_analysis.gravy()
    desc.instability_index()
    features['instability_index'] = desc.descriptor[0][0]
    desc.aromaticity()
    features['aromaticity'] = desc.descriptor[0][0]
    desc.aliphatic_index()
    features['aliphatic_index'] = desc.descriptor[0][0]
    desc.boman_index()
    features['boman_index'] = desc.descriptor[0][0]
    desc.hydrophobic_ratio()
    features['hydrophobic_ratio'] = desc.descriptor[0][0]
    return features


def get_composition_features(features):
    seq = features['sequence']
    aa_composition = AAComposition.CalculateAAComposition(seq)
    dipep_composition = AAComposition.CalculateDipeptideComposition(seq)
    features.update(aa_composition)
    features.update(dipep_composition)
    return features


# Composición-Transición-Distribución (CTD)
# Cálcula los descriptores correspondientes Composition-Transition-Distribution#
def get_ctd_features(features):
    # hardProteinSequence=string.replace(hardProteinSequence,index,k)
    # cambiar por
    # hardProteinSequence=hardProteinSequence.replace(index,k)
    # module 'string' has no attribute 'find'
    seq = features['sequence']
    ctd = CTD.CalculateCTD(seq)
    features.update(ctd)
    return features


# Autocorrelación
# Cálcula los descriptores de autocorrelación basados en Moran, Geary y Broto  #
def get_autocorr_features(features):
    seq = features['sequence']
    auto = Autocorrelation.CalculateAutoTotal(seq)
    features.update(auto)
    return features


def get_pseudoaa_features(features):
    seq = features['sequence']
    Hydrophobicity = PseudoAAC._Hydrophobicity
    hydrophilicity = PseudoAAC._hydrophilicity
    residuemass = PseudoAAC._residuemass
    pK1 = PseudoAAC._pK1
    pK2 = PseudoAAC._pK2
    pI = PseudoAAC._pI
    PseAACI1 = PseudoAAC.GetPseudoAAC1(seq, lamda=3, weight=0.05,
                                       AAP=[Hydrophobicity,
                                            hydrophilicity,
                                            residuemass, pK1, pK2,
                                            pI])
    PseAACI2 = PseudoAAC.GetPseudoAAC2(seq, lamda=3, weight=0.05,
                                       AAP=[Hydrophobicity,
                                            hydrophilicity,
                                            residuemass, pK1, pK2,
                                            pI])

    PseAACII1 = PseudoAAC.GetAPseudoAAC1(seq, lamda=3, weight=0.5)
    pre_PseAACII2 = PseudoAAC.GetAPseudoAAC2(seq, lamda=3, weight=0.5)
    PseAACII2 = {}
    for key in pre_PseAACII2.keys():
        PseAACII2["A" + key] = pre_PseAACII2[key]
    features.update(PseAACI1)
    features.update(PseAACI2)
    features.update(PseAACII1)
    features.update(PseAACII2)
    return features


def get_quasi_sec_features(features):
    seq = features['sequence']
    ocnt = QuasiSequenceOrder.GetSequenceOrderCouplingNumberTotal(seq)
    q_grant1 = QuasiSequenceOrder.GetQuasiSequenceOrder1Grant(seq)
    q_grant2 = QuasiSequenceOrder.GetQuasiSequenceOrder2Grant(seq)
    q_sw1 = QuasiSequenceOrder.GetQuasiSequenceOrder1SW(seq)
    q_sw2 = QuasiSequenceOrder.GetQuasiSequenceOrder2SW(seq)
    features.update(ocnt)
    features.update(q_grant1)
    features.update(q_grant2)
    features.update(q_sw1)
    features.update(q_sw2)
    return features


def get_ctriad_features(features):
    # https://github.com/Superzchen/iLearn/blob/master/descproteins/CTriad.py
    def CalculateKSCTriad(sequence, gap, features, AADict):
        res = []
        for g in range(gap + 1):
            myDict = {}
            for f in features:
                myDict[f] = 0

            for i in range(len(sequence)):
                if i + gap + 1 < len(sequence) and i + 2 * gap + 2 < len(sequence):
                    fea = AADict[sequence[i]] + '.' + AADict[sequence[i + gap + 1]] + '.' + AADict[
                        sequence[i + 2 * gap + 2]]
                    myDict[fea] = myDict[fea] + 1

            maxValue, minValue = max(myDict.values()), min(myDict.values())
            for f in features:
                res.append((myDict[f] - minValue) / maxValue)

        return res

    AAGroup = {
        'g1': 'AGV',
        'g2': 'ILFP',
        'g3': 'YMTS',
        'g4': 'HNQW',
        'g5': 'RK',
        'g6': 'DE',
        'g7': 'C'
    }
    myGroups = sorted(AAGroup.keys())
    AADict = {}
    for g in myGroups:
        for aa in AAGroup[g]:
            AADict[aa] = g
    feats = [f1 + '.' + f2 + '.' + f3 for f1 in myGroups for f2 in myGroups for f3 in myGroups]
    seq = features['sequence']
    desc = CalculateKSCTriad(seq, 0, feats, AADict)
    feats = ['CTriad_' + feat for feat in feats]
    ctriad_dic = {feats[i]: desc[i] for i in range(len(desc))}
    features.update(ctriad_dic)
    return features


# def get_embed_features(features):
#     seq = features['sequence']
#     embedding = pv2.to_vecs(seq)
#     for i in range(3):
#         for j in range(100):
#             features["embed_{}_{}".format(i, j)] = embedding[i][j]
#     return features


def get_features(seq):
    """This function receives a seqIO sequence container as inpunt and returns a feature
    dictionary."""
    features = collections.OrderedDict()
    features['sequence'] = seq
    features['length'] = len(seq)
    features = get_global_features(features)
    features = get_composition_features(features)
    features = get_grouped_aa_features(features)
    features = get_ctd_features(features)
    features = get_autocorr_features(features)
    features = get_quasi_sec_features(features)
    # features = get_pseudoaa_features(features)
    features = get_ctriad_features(features)
    #features = get_embed_features(features)
    return features


def get_features_timed(seq):
    """This function receives a seqIO sequence container as inpunt and returns a feature
    dictionary."""
    features = collections.OrderedDict()
    features['sequence'] = seq
    features['length'] = len(seq)
    t1 = time.time()
    features = get_global_features(features)
    t2 = time.time()
    print(len(features))
    print(t2 - t1)
    features = get_composition_features(features)
    t3 = time.time()
    print(len(features))
    print(t3 - t2)
    features = get_ctd_features(features)
    t4 = time.time()
    print(len(features))
    print(t4 - t3)
    features = get_autocorr_features(features)
    t5 = time.time()
    print(len(features))
    print(t5 - t4)
    features = get_quasi_sec_features(features)
    t6 = time.time()
    print(len(features))
    print(t6 - t5)
    features = get_pseudoaa_features(features)
    t7 = time.time()
    print(len(features))
    print("aaa", t7 - t6)
    features = get_ctriad_features(features)
    t8 = time.time()
    print(len(features))
    print(t8 - t7)
    features = get_grouped_aa_features(features)
    t9 = time.time()
    print(len(features))
    print(t9 - t8)
    #features = get_embed_features(features)
    t10 = time.time()
    print(len(features))
    print(t10 - t9)
    print(t10 - t1)
    return features
