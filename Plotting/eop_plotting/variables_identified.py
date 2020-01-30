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
    return np.ones(len(vertex)) * 0.0
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

def vertex_pt(vertex):
    return vertex["vertex_pt"]

calc_vertex_pt = Calculation(vertex_pt, ["vertex_pt"])

def cos_theta(vertex):
    pv_x = vertex["primary_vertex_x"]
    pv_y = vertex["primary_vertex_y"]
    pv_z = vertex["primary_vertex_z"]

    x = vertex["vertex_x"]
    y = vertex["vertex_y"]
    z = vertex["vertex_z"]

    #this is the vector pointing from the primary vertex to the secondary vertex
    x = x - pv_x
    y = y - pv_y
    z = z - pv_z

    # get the momentum vector of the primary vertex
    px = vertex["trk1_px"] + vertex["trk2_px"] 
    py = vertex["trk1_py"] + vertex["trk2_py"] 
    pz = vertex["trk1_pz"] + vertex["trk2_pz"] 

    #calculate the cos theta star
    numerator = px*x + py*y + pz+z
    mag_pv = np.sqrt(px**2 + py**2 + pz**2)
    mag_v = np.sqrt(x**2 + y**2 + z**2)
    non_zero = (mag_pv > 0.0) & (mag_v > 0.0)
    cos_theta = np.ones(len(vertex)) * 1000.0
    cos_theta[non_zero] = numerator[non_zero] / (mag_pv[non_zero] * mag_v[non_zero])
    return cos_theta

calc_cos_theta = Calculation(cos_theta, ["primary_vertex_{}".format(x) for x in ["x", "y","z"]]\
        + ["vertex_{}".format(x) for x in ["x", "y","z"]] \
        + ["trk1_{}".format(x) for x in ["px", "py","pz"]]\
        + ["trk2_{}".format(x) for x in ["px", "py", "pz"]])

def pos_track_momentum(vertex):
    trk_momentum = np.ones(len(vertex))
    trk1_pos = (vertex["trk1_charge"] > 0.5) & (vertex["trk2_charge"] < -0.5)
    trk1_neg = (vertex["trk1_charge"] < -0.5) & (vertex["trk2_charge"] > 0.5)
    two_same_charge = (vertex["trk1_charge"]  == vertex["trk2_charge"])
    assert not np.any(two_same_charge) #check that the tracks really did both have opposite charge
    trk_momentum[trk1_pos] = vertex["trk1_p"][trk1_pos]
    trk_momentum[trk1_neg] = vertex["trk2_p"][trk1_neg]
    return trk_momentum
calc_pos_track_momentum = Calculation(pos_track_momentum, ["trk1_momentum", "trk2_momentum", "trk1_charge", "trk2_charge"])

def neg_track_momentum(vertex):
    trk_momentum = np.ones(len(vertex))
    trk1_pos = (vertex["trk1_charge"] > 0.5) & (vertex["trk2_charge"] < -0.5)
    trk1_neg = (vertex["trk1_charge"] < -0.5) & (vertex["trk2_charge"] > 0.5)
    two_same_charge = (vertex["trk1_charge"]  == vertex["trk2_charge"])
    assert not np.any(two_same_charge) #check that the tracks really did both have opposite charge
    trk_momentum[trk1_pos] = vertex["trk2_p"][trk1_pos] #track 1 positive -> trk2 is negative
    trk_momentum[trk1_neg] = vertex["trk1_p"][trk1_neg] #track 1 negative -> trk1 is negative
    return trk_momentum
calc_neg_track_momentum = Calculation(neg_track_momentum, ["trk1_momentum", "trk2_momentum", "trk1_charge", "trk2_charge"])

def pos_track_energy(vertex):
    trk_energy = np.ones(len(vertex))
    trk1_pos = (vertex["trk1_charge"] > 0.5) & (vertex["trk2_charge"] < -0.5)
    trk1_neg = (vertex["trk1_charge"] < -0.5) & (vertex["trk2_charge"] > 0.5)
    two_same_charge = (vertex["trk1_charge"]  == vertex["trk2_charge"])
    assert not np.any(two_same_charge) #check that the tracks really did both have opposite charge
    trk_energy[trk1_pos] = (vertex["trk1_ClusterEnergy_EM_200"] + vertex["trk1_ClusterEnergy_HAD_200"])[trk1_pos] #track 1 positive -> trk2 is negative
    trk_energy[trk1_neg] = (vertex["trk2_ClusterEnergy_EM_200"] + vertex["trk2_ClusterEnergy_HAD_200"])[trk1_neg] #track 1 negative -> trk1 is negative
    return trk_energy
branches = ["trk1_ClusterEnergy_EM_200", "trk1_ClusterEnergy_HAD_200"]
branches += ["trk2_ClusterEnergy_EM_200", "trk2_ClusterEnergy_HAD_200"]
calc_pos_track_energy = Calculation(pos_track_energy, ["trk1_charge", "trk2_charge"] + branches)

def neg_track_energy(vertex):
    trk_energy = np.ones(len(vertex))
    trk1_pos = (vertex["trk1_charge"] > 0.5) & (vertex["trk2_charge"] < -0.5)
    trk1_neg = (vertex["trk1_charge"] < -0.5) & (vertex["trk2_charge"] > 0.5)
    two_same_charge = (vertex["trk1_charge"]  == vertex["trk2_charge"])
    assert not np.any(two_same_charge) #check that the tracks really did both have opposite charge
    trk_energy[trk1_pos] = (vertex["trk2_ClusterEnergy_EM_200"] + vertex["trk2_ClusterEnergy_HAD_200"])[trk1_pos] #track 1 positive -> trk2 is negative
    trk_energy[trk1_neg] = (vertex["trk1_ClusterEnergy_EM_200"] + vertex["trk1_ClusterEnergy_HAD_200"])[trk1_neg] #track 1 negative -> trk1 is negative
    return trk_energy
branches = ["trk1_ClusterEnergy_EM_200", "trk1_ClusterEnergy_HAD_200"]
branches += ["trk2_ClusterEnergy_EM_200", "trk2_ClusterEnergy_HAD_200"]
calc_neg_track_energy = Calculation(neg_track_energy, ["trk1_charge", "trk2_charge"] + branches)

def pos_track_nearest_dR_EM(vertex):
    trk_energy = np.ones(len(vertex))
    trk1_pos = (vertex["trk1_charge"] > 0.5) & (vertex["trk2_charge"] < -0.5)
    trk1_neg = (vertex["trk1_charge"] < -0.5) & (vertex["trk2_charge"] > 0.5)
    two_same_charge = (vertex["trk1_charge"]  == vertex["trk2_charge"])
    assert not np.any(two_same_charge) #check that the tracks really did both have opposite charge
    trk_nearest[trk1_pos] = vertex["trk1_nearest_dR_EM"]#track 1 positive -> trk2 is negative
    trk_nearest[trk1_neg] = vertex["trk2_nearest_dR_EM"]#track 1 negative -> trk1 is negative
    return trk_nearest
calc_pos_track_nearest_dR_EM = Calculation(pos_track_nearest_dR_EM, ["trk1_nearest_dR_EM", "trk2_nearest_dR_EM", "trk1_charge", "trk2_charge"])

def neg_track_nearest_dR_EM(vertex):
    trk_energy = np.ones(len(vertex))
    trk1_pos = (vertex["trk1_charge"] > 0.5) & (vertex["trk2_charge"] < -0.5)
    trk1_neg = (vertex["trk1_charge"] < -0.5) & (vertex["trk2_charge"] > 0.5)
    two_same_charge = (vertex["trk1_charge"]  == vertex["trk2_charge"])
    assert not np.any(two_same_charge) #check that the tracks really did both have opposite charge
    trk_nearest[trk1_pos] = vertex["trk2_nearest_dR_EM"]#track 1 positive -> trk2 is negative
    trk_nearest[trk1_neg] = vertex["trk1_nearest_dR_EM"]#track 1 negative -> trk1 is negative
    return trk_nearest
calc_neg_track_nearest_dR_EM = Calculation(neg_track_nearest_dR_EM, ["trk1_nearest_dR_EM", "trk2_nearest_dR_EM", "trk1_charge", "trk2_charge"])

def pos_track_pt(vertex):
    trk_pt = np.ones(len(vertex))
    trk1_pos = (vertex["trk1_charge"] > 0.5) & (vertex["trk2_charge"] < -0.5)
    trk1_neg = (vertex["trk1_charge"] < -0.5) & (vertex["trk2_charge"] > 0.5)
    two_same_charge = (vertex["trk1_charge"]  == vertex["trk2_charge"])
    assert not np.any(two_same_charge) #check that the tracks really did both have opposite charge
    trk_pt[trk1_pos] = vertex["trk1_pt"][trk1_pos]
    trk_pt[trk1_neg] = vertex["trk2_pt"][trk1_neg]
    return trk_pt
calc_pos_track_pt = Calculation(pos_track_pt, ["trk1_pt", "trk2_pt", "trk1_charge", "trk2_charge"])

def neg_track_pt(vertex):
    trk_pt = np.ones(len(vertex))
    trk1_pos = (vertex["trk1_charge"] > 0.5) & (vertex["trk2_charge"] < -0.5)
    trk1_neg = (vertex["trk1_charge"] < -0.5) & (vertex["trk2_charge"] > 0.5)
    two_same_charge = (vertex["trk1_charge"]  == vertex["trk2_charge"])
    assert not np.any(two_same_charge) #check that the tracks really did both have opposite charge
    trk_pt[trk1_pos] = vertex["trk2_pt"][trk1_pos] #track 1 positive -> trk2 is negative
    trk_pt[trk1_neg] = vertex["trk1_pt"][trk1_neg] #track 1 negative -> trk1 is negative
    return trk_pt
calc_neg_track_pt = Calculation(neg_track_pt, ["trk1_pt", "trk2_pt", "trk1_charge", "trk2_charge"])

def pos_track_eop(vertex):
    return calc_pos_track_energy.eval(vertex)/calc_pos_track_momentum.eval(vertex)
calc_pos_track_eop = Calculation(pos_track_eop, calc_pos_track_momentum.branches + calc_pos_track_energy.branches)

def neg_track_eop(vertex):
    return calc_neg_track_energy.eval(vertex)/calc_neg_track_momentum.eval(vertex)
calc_neg_track_eop = Calculation(neg_track_eop, calc_neg_track_momentum.branches + calc_neg_track_energy.branches)

