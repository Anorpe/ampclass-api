# -*- coding: utf-8 -*-
import math
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from modlamp.descriptors import GlobalDescriptor




#AAs Alifaticos
def cargas(peptidos):
    """

    :param peptidos:
    :return:
    """
    centro = []
    extremo = []
    aleatorio = []
    for peptido in peptidos:

        ventana = int(len(peptido) / 3)
        z = []
        t = [0, 0, 0]
        z.append(peptido[:ventana])
        z.append(peptido[ventana:ventana + ventana])
        z.append(peptido[ventana * 2:])

        for i in range(3):
            for aminoacido in z[i]:
                if aminoacido in ['K', 'R', 'D', 'E']:
                    t[i] += 1

        if (t[1] > ventana / 2):
            if (t[0] < ventana / 2):
                if (t[2] < ventana / 2):
                    centro.append(peptido)
        elif (t[0] > ventana / 2):
            if (t[1] < ventana / 2):
                extremo.append(peptido)

        elif (t[2] > ventana / 2):
            if (t[1] < ventana / 2):
                extremo.append(peptido)
        else:
            aleatorio.append(peptido)

    return centro, extremo, aleatorio



def alifaticos(peptidos):
    retorno = []

    for peptido in peptidos:
        alifaticos = ['G','A','V','L','I']
        for alifatico in alifaticos:
            if( alifatico in peptido):
                retorno.append(peptido)

    return retorno

def longitud(peptidos):
    longitudes = []
    for peptido in peptidos:
        longitudes.append(len(peptido))
    return longitudes

def indice_boman(peptidos):
    glob = GlobalDescriptor(peptidos)
    glob.boman_index()
    indiceBoman = glob.descriptor
    indiceBoman = list(np.array(indiceBoman.ravel()).round(2))

    return indiceBoman


def  carga(peptidos):
    glob = GlobalDescriptor(peptidos)
    glob.calculate_charge(ph=7)
    carga = glob.descriptor
    carga = list(np.array(carga.ravel()).round(0))

    return carga


def aromaticos(peptidos):
    retorno = []
    for peptido in peptidos:
        aromaticos = ['F','W','Y']
        for aromatico in aromaticos:
            if(aromatico in peptido):
                retorno.append(peptido)


    return retorno


def no_alifaticos(peptidos):
    retorno = []

    for peptido in peptidos:
        alifaticos = ['G','A','V','L','I']
        for alifatico in alifaticos:
            if(not alifatico in peptido):
                retorno.append(peptido)

    return retorno




def no_aromaticos(peptidos):
    retorno = []
    for peptido in peptidos:
        aromaticos = ['F','W','Y']
        for aromatico in aromaticos:
            if(not aromatico in peptido):
                retorno.append(peptido)


    return retorno





# Punto Isoelectrico
def punto_isoelectrico(peptidos):
    puntoIsoelectrico = []
    for peptido in peptidos:
        biop_analysis = ProteinAnalysis(peptido)
        isoelectric_point = biop_analysis.isoelectric_point()
        puntoIsoelectrico.append(isoelectric_point)

    puntoIsoelectrico  = list(np.array(puntoIsoelectrico).round(2))
    return puntoIsoelectrico

def punto_isoelectricoManual(peptidos):
    puntoIsoelectrico = []
    for peptido in peptidos:
        CNi = 0
        CNj = 0
        PI = 0
        CN = 0
        i = 0

        while (i <= 14):

            for aminoacido in peptido:
                if (aminoacido == 'R'):
                    CNi = CNi + ((10 ** 12.48) / ((10 ** i) + (10 ** 12.48)))

                elif (aminoacido == 'K'):
                    CNi = CNi + ((10 ** 10.79) / ((10 ** i) + (10 ** 10.79)))

                elif (aminoacido == 'H'):
                    CNi = CNi + ((10 ** 6.04) / ((10 ** i) + (10 ** 6.04)))

                elif (aminoacido == 'D'):
                    CNj = CNj + ((10 ** i) / ((10 ** i) + (10 ** 3.86)))

                elif (aminoacido == 'E'):
                    CNj = CNj + ((10 ** i) / ((10 ** i) + (10 ** 4.25)))

                elif (aminoacido == 'C'):
                    CNj = CNj + ((10 ** i) / ((10 ** i) + (10 ** 8.33)))

                elif (aminoacido == 'Y'):
                    CNj = CNj + ((10 ** i) / ((10 ** i) + (10 ** 10.07)))

                else:
                    CNi = CNi
                    CNj = CNj

            CNi = CNi + ((10 ** 9.69) / ((10 ** i) + (10 ** 9.69)))

            CNj = CNj + ((10 ** i) / ((10 ** i) + (10 ** 2.34)))

            CN = CNi - CNj
            CN = round(CN, 1)

            if (CN >= -0.1):

                if (CN <= 0.1):
                    PI = round(i, 1)
            else:
                i = i

            CNi = 0
            CNj = 0

            i += 0.01

        puntoIsoelectrico.append(PI)

    return puntoIsoelectrico


# Momento hidrofobico relativo
# revisar
def momento_hidrofobico_relativo(peptidos):
    momentos = []
    for peptido in peptidos:
        L = len(peptido)
        delta = (5 / 9) * math.pi
        # Numero de aminoacidos segun escala Kyte & Doolitle
        hidrofobicidad = {
            'A': 1.8,
            'R': -4.5,
            'N': 3.5,
            'D': -3.5,
            'C': 2.5,
            'Q': -3.5,
            'E': -3.5,
            'G': -0.4,
            'H': -3.2,
            'I': 4.5,
            'L': 3.8,
            'K': -3.9,
            'M': 1.9,
            'F': 2.8,
            'P': -1.6,
            'S': -0.8,
            'T': -0.7,
            'W': -0.9,
            'Y': -1.3,
            'V': 4.2
        }
        uH = 0
        i = 1

        sum1 = 0
        sum2 = 0
        for aminoacido in peptido:
            sum1 = sum1 + hidrofobicidad[aminoacido] * math.cos(i * delta)
            sum2 = sum2 + hidrofobicidad[aminoacido] * math.sin(i * delta)

            i = i + 1

        uH = sum1 ** 2 + sum2 ** 2
        uH = math.sqrt(uH)
        uH = uH / L

        uH = (uH / 2.88) * 100


        momentos.append(uH)
    return momentos



def wimley_white(peptidos):
    wimleyWhites = []
    for peptido in peptidos:
        L = len(peptido)

        wimley = {'A':0.33, 'I':-0.81, 'L':-0.69, 'W':-0.24, 'F':-0.58, 'V':-0.53, 'M':-0.44, 'Y':0.23, 'P':-0.31, 'T':0.11,
        'S':0.33, 'C':0.22, 'G':1.14, 'N':0.43, 'D':2.41, 'Q':0.19, 'E':1.61, 'H':-0.06, 'K':1.81, 'R':1.00}



        wimleyWhite = 0
        for aminoacido in peptido:
            wimleyWhite = wimleyWhite + wimley[aminoacido]


        wimleyWhites.append(wimleyWhite)

    return wimleyWhites


def porcentaje_hidrofobico(peptidos, tipo):
    porcentajes = []
    for peptido in peptidos:
        L = len(peptido)
        # Numero de aminoacidos segun escala Kyte & Doolitle

        kyte = {'A': 1.8, 'R': -4.5, 'N': 3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2,
                'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9,
                'Y': -1.3, 'V': 4.2}
        eisenberg = {'A': 0.25, 'R': -1.76, 'N': -0.64, 'D': -0.72, 'C': 0.04, 'Q': -0.69, 'E': -0.62, 'G': 0.16,
                     'H': -0.4, 'I': 0.73, 'L': 0.53, 'K': -1.1, 'M': 0.26, 'F': 0.61, 'P': -0.07, 'S': -0.26,
                     'T': -0.18, 'W': 0.37, 'Y': 0.02, 'V': 0.54}
        engleman = {'A': 1.6, 'R': -12.3, 'N': -4.8, 'D': -9.2, 'C': 2, 'Q': -4.1, 'E': -8.2, 'G': 1, 'H': -3, 'I': 3.1,
                    'L': 2.8, 'K': -8.8, 'M': 3.4, 'F': 3.7, 'P': -0.2, 'S': 0.6, 'T': 1.2, 'W': 1.9, 'Y': -0.7,
                    'V': 2.6}
        hoops = {'A': -0.5, 'R': 3, 'N': 0.2, 'D': 3, 'C': - 1, 'Q': 0.2, 'E': 3, 'G': 0, 'H': - 0.5, 'I': - 1.8,
                 'L': - 1.8, 'K': 3, 'M': - 1.3, 'F': - 2.5, 'P': 0, 'S': 0.3, 'T': - 0.4, 'W': - 3.4, 'Y': - 2.3,
                 'V': - 1.5}
        janin = {'A': 0.3, 'R': - 1.4, 'N': - 0.5, 'D': - 0.6, 'C': 0.9, 'Q': - 0.7, 'E': - 0.7, 'G': 0.3, 'H': - 0.1,
                 'I': 0.7, 'L': 0.5, 'K': - 1.8, 'M': 0.4, 'F': 0.5, 'P': - 0.3, 'S': - 0.1, 'T': - 0.2, 'W': 0.3,
                 'Y': - 0.4, 'V': 0.6}
        wimley = {'A':-0.17, 'I':0.31, 'L':0.56, 'W':1.85, 'F':1.13, 'V':-0.07, 'M':0.23, 'Y':0.94, 'P':-0.45, 'T':-0.14,
        'S':-0.13, 'C':0.24, 'G':-0.01, 'N':-0.42, 'D':-1.23, 'Q':-0.58, 'E':-2.02, 'H':-0.96, 'K':-0.99, 'R':-0.81}


        escalas = {'kyte': kyte, 'einsenberg': eisenberg, 'engleman': engleman, 'hoops': hoops, 'janin': janin,'wimley':wimley}
        hidrofobicidad = escalas[tipo]
        H = 0
        for aminoacido in peptido:
            H = H + hidrofobicidad[aminoacido]

        H = H / L
        H = H * 100
        if (H < 0):
            H = 0
        if (H > 100):
            H = 100
        porcentajes.append(H)

    return porcentajes



# Helices transmembranales

def helices_transmembrana(peptidos, metodo):
    if (metodo == 'eisenberg'):

        # Metodo Eisenberg

        pi = math.pi
        helices = []

        for peptido in peptidos:

            i = 0
            Ky = 0
            Kp = 0
            sen = 0
            cos = 0
            MH = 0
            uH = 0

            for aminoacido in peptido:

                if (aminoacido == 'A'):
                    Ky = Ky + 0.310
                    sen = sen + 0.310 * math.sin(((i + 1) * 100 * pi) / 180)
                    cos = cos + 0.310 * math.cos(((i + 1) * 100 * pi) / 180)

                elif (aminoacido == 'R'):
                    Ky = Ky - 1.010
                    sen = sen - 1.010 * math.sin(((i + 1) * 100 * pi) / 180)
                    cos = cos - 1.010 * math.cos(((i + 1) * 100 * pi) / 180)


                elif (aminoacido == 'N'):
                    Ky = Ky - 0.600
                    sen = sen - 0.600 * math.sin(((i + 1) * 100 * pi) / 180)
                    cos = cos - 0.600 * math.cos(((i + 1) * 100 * pi) / 180)

                elif (aminoacido == 'D'):
                    Ky = Ky - 0.770
                    sen = sen - 0.770 * math.sin(((i + 1) * 100 * pi) / 180)
                    cos = cos - 0.770 * math.cos(((i + 1) * 100 * pi) / 180)


                elif (aminoacido == 'C'):
                    Ky = Ky + 1.540
                    sen = sen + 1.540 * math.sin(((i + 1) * 100 * pi) / 180)
                    cos = cos + 1.540 * math.cos(((i + 1) * 100 * pi) / 180)

                elif (aminoacido == 'Q'):
                    Ky = Ky - 0.220
                    sen = sen - 0.220 * math.sin(((i + 1) * 100 * pi) / 180)
                    cos = cos - 0.220 * math.cos(((i + 1) * 100 * pi) / 180)


                elif (aminoacido == 'E'):
                    Ky = Ky - 0.640
                    sen = sen - 0.640 * math.sin(((i + 1) * 100 * pi) / 180)
                    cos = cos - 0.640 * math.cos(((i + 1) * 100 * pi) / 180)

                elif (aminoacido == 'G'):
                    Ky = Ky + 0.000
                    sen = sen + 0.000 * math.sin(((i + 1) * 100 * pi) / 180)
                    cos = cos + 0.000 * math.cos(((i + 1) * 100 * pi) / 180)

                elif (aminoacido == 'H'):
                    Ky = Ky + 0.130
                    sen = sen + 0.130 * math.sin(((i + 1) * 100 * pi) / 180)
                    cos = cos + 0.130 * math.cos(((i + 1) * 100 * pi) / 180)

                elif (aminoacido == 'I'):
                    Ky = Ky + 1.800
                    sen = sen + 1.800 * math.sin(((i + 1) * 100 * pi) / 180)
                    cos = cos + 1.800 * math.cos(((i + 1) * 100 * pi) / 180)

                elif (aminoacido == 'L'):
                    Ky = Ky + 1.700
                    sen = sen + 1.700 * math.sin(((i + 1) * 100 * pi) / 180)
                    cos = cos + 1.700 * math.cos(((i + 1) * 100 * pi) / 180)


                elif (aminoacido == 'K'):
                    Ky = Ky - 0.990
                    sen = sen - 0.990 * math.sin(((i + 1) * 100 * pi) / 180)
                    cos = cos - 0.990 * math.cos(((i + 1) * 100 * pi) / 180)


                elif (aminoacido == 'M'):
                    Ky = Ky + 1.230
                    sen = sen + 1.230 * math.sin(((i + 1) * 100 * pi) / 180)
                    cos = cos + 1.230 * math.cos(((i + 1) * 100 * pi) / 180)

                elif (aminoacido == 'F'):
                    Ky = Ky + 1.790
                    sen = sen + 1.790 * math.sin(((i + 1) * 100 * pi) / 180)
                    cos = cos + 1.790 * math.cos(((i + 1) * 100 * pi) / 180)

                elif (aminoacido == 'P'):
                    Ky = Ky + 0.720
                    sen = sen + 0.720 * math.sin(((i + 1) * 100 * pi) / 180)
                    cos = cos + 0.720 * math.cos(((i + 1) * 100 * pi) / 180)

                elif (aminoacido == 'S'):
                    Ky = Ky - 0.040
                    sen = sen - 0.040 * math.sin(((i + 1) * 100 * pi) / 180)
                    cos = cos - 0.040 * math.cos(((i + 1) * 100 * pi) / 180)

                elif (aminoacido == 'T'):
                    Ky = Ky + 0.260
                    sen = sen + 0.260 * math.sin(((i + 1) * 100 * pi) / 180)
                    cos = cos + 0.260 * math.cos(((i + 1) * 100 * pi) / 180)

                elif (aminoacido == 'W'):
                    Ky = Ky + 2.250
                    sen = sen + 2.250 * math.sin(((i + 1) * 100 * pi) / 180)
                    cos = cos + 2.250 * math.cos(((i + 1) * 100 * pi) / 180)

                elif (aminoacido == 'Y'):
                    Ky = Ky + 0.960
                    sen = sen + 0.960 * math.sin(((i + 1) * 100 * pi) / 180)
                    cos = cos + 0.960 * math.cos(((i + 1) * 100 * pi) / 180)
                elif (aminoacido == 'V'):
                    Ky = Ky + 1.220
                    sen = sen + 1.220 * math.sin(((i + 1) * 100 * pi) / 180)
                    cos = cos + 1.220 * math.cos(((i + 1) * 100 * pi) / 180)

            Kp = round(Ky / len(peptido), 3)
            MH = round((1 / len(peptido)) * (((sen ** 2) + (cos ** 2)) ** 0.5), 3)
            uH = round(0.654 - (0.324 * Kp), 3)

            if (Kp > 0.75 and MH < uH):
                helices.append('si')
            else:
                helices.append('no')

        return helices




    elif (metodo == 'white'):
        # Metodo White and Wimley
        helices = []
        for peptido in peptidos:

            WL = 0
            for aminoacido in peptido:
                valores = {'A': 0.33,'R': 1,'N': 0.43,'D': 2.41,'C': 0.22,'Q': 1.61,'E': 0.19,'G': 1.14,'H': -0.06,
                    'I': -0.81,'L': -0.69,'K': 1.81,'M': -0.44,'F': -0.58,'P': -0.31,'S': 0.33,'T': 0.11,'W': -0.24,'Y': 0.23,'V': -0.53
                }

                WL += valores[aminoacido]

            WL = round(WL, 3)
            if (WL < 0):
                helices.append('si')
            else:
                helices.append('no')

        return helices


# Polar Angles


def get_polar_dic():
    """
    Retorna un diccionario de los valores de hidrofobicidad según la escala Fauchere-Pliska usada
    en Polar Angle as a Determinant of Amphipathic a-Helix-Lipid Interactions: A Model Peptide Study.
    https://web.expasy.org/protscale/pscale/Hphob.Fauchere.html"""
    dic_idx = {'A': 0.310, 'R': -1.010, 'N': -0.600, 'D': -0.770, 'C': 1.540, 'Q': -0.220, 'E': -0.640, 'G': 0.000,
               'H': 0.130,
               'I': 1.800, 'L': 1.700, 'K': -0.990, 'M': 1.230, 'F': 1.790, 'P': 0.720, 'S': -0.040, 'T': 0.260,
               'W': 2.250, 'Y': 0.960, 'V': 1.220}
    return dic_idx


def get_positions():
    theta = 0
    lista = []
    for i in range(0, 18):
        lista.append(theta)
        theta = theta + 100
        if (theta >= 360):
            theta = theta - 360
    seq = sorted(lista)
    indexes = [seq.index(v) for v in lista]
    return indexes


"""


"""


def to_helix(pep):
    """Retorna la representación como helice ideal de los valores de hidrofobicidad de un peptido."""
    positions = get_positions()
    polar_dic = get_polar_dic()
    helix = []
    idx = 0
    floor = [0] * 18

    for aa in pep:
        floor[positions[idx]] = polar_dic[aa]
        idx = idx + 1
        if (idx == 18):
            helix.append(floor)
            idx = 0
            floor = [0] * 18
    if idx != 0:
        helix.append(floor)
    return helix


def helix_to_polar_struct(helix):
    n_floors = len(helix)
    flattened_helix = [0] * 18
    for i in range(0, 18):
        for j in range(n_floors):
            flattened_helix[i] += helix[j][i]
    return flattened_helix


def get_longest_hydrophilic_sub_chain(flattened_helix):
    longest = []
    for i in range(0, 18):
        cont = 0
        if (flattened_helix[i] >= 0):
            continue
        sub_chain = []
        idx = i
        while (cont <= 2):
            aa = flattened_helix[idx]
            sub_chain.append(aa)
            if (aa >= 0):
                cont += 1
            else:
                cont = 0
            idx += 1
            if (idx == 18):
                idx = 0
            if (idx == i):
                break
        if (len(sub_chain) < 18):
            sub_chain = sub_chain[:-cont]
        if (len(sub_chain) > len(longest)):
            longest = sub_chain
    return longest


def polar_angles(peptidos):
    polarAngles = []

    for peptido in peptidos:
        helix = to_helix(peptido)

        flattened_helix = helix_to_polar_struct(helix)

        longest = get_longest_hydrophilic_sub_chain(flattened_helix)

        polar_angle = len(longest) / 18 * 360

        polarAngles.append(polar_angle)

    return polarAngles
