from calculation import Calculation, CalculationDataMC, WeightCalculation
import numpy as np

def weight(vertex, isData):
    if not isData:
        return vertex["mcWeight"]
    return np.ones(len(vertex))
calc_weight = WeightCalculation(weight, ["mcWeight"])

def vertex_mass(vertex):
    return vertex["vertex_mass"]
calc_vertex_mass = Calculation(vertex_mass, ["vertex_mass"])

def vertex_count(vertex):
    return np.ones(len(vertex)) * 0.5
calc_vertex_count = Calculation(vertex_count, [])

def vertex_chiSquared(vertex):
    return vertex["vertex_chiSquared"]
calc_vertex_chiSquared = Calculation(vertex_chiSquared, ["vertex_chiSquared"])

def vertex_Rxy(vertex):
    return vertex["vertex_Rxy"]
calc_vertex_Rxy = Calculation(vertex_Rxy, ["vertex_Rxy"])

def higher_trk_p(vertex):
    trk_p = np.zeros(len(vertex))
    trk_p[vertex["trk1_p"] > vertex["trk2_p"]] = trk1_p[vertex["trk1_p"] > vertex["trk2_p"]]
    trk_p[vertex["trk2_p"] >= vertex["trk1_p"]] = trk2_p[vertex["trk2_p"] >= vertex["trk1_p"]]
    return trk_p

