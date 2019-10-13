from calculation import Calculation

## Take a look here for the run-1 eop selections: https://journals.aps.org/prd/pdf/10.1103/PhysRevD.85.012001

def chi_square_fifteen(vertex):
    return vertex["vertex_chiSquared"] < 15
sel_chi_square_fifteen = Calculation(chi_square_fifteen, ["vertex_chiSquared"])

from variables_identified import calc_cos_theta
def tight_cos_theta_ks(vertex):
    return calc_cos_theta.eval(vertex) > 0.999
sel_tight_cos_theta_ks = Calculation(tight_cos_theta_ks, calc_cos_theta.branches)

def rxy_ks(vertex):
    return (vertex["vertex_Rxy"] < 450) & (vertex["vertex_Rxy"] > 4.0)
sel_rxy_ks = Calculation(rxy_ks, ["vertex_Rxy"])

def pt_ks(vertex):
    return (vertex["vertex_pt"] > 0.1)
sel_pt_ks = Calculation(pt_ks, ["vertex_pt"])

def tight_cos_theta_lambda(vertex):
    return calc_cos_theta.eval(vertex) > 0.9998
sel_tight_cos_theta_lambda = Calculation(tight_cos_theta_lambda, calc_cos_theta.branches)

def rxy_lambda(vertex):
    return (vertex["vertex_Rxy"] < 450) & (vertex["vertex_Rxy"] > 17.0)
sel_rxy_lambda = Calculation(rxy_lambda, ["vertex_Rxy"])

def pt_lambda(vertex):
    return (vertex["vertex_pt"] > 0.5)
sel_pt_lambda = Calculation(pt_lambda, ["vertex_pt"])

def pos_track_momentum(vertex):
    trk_momentum = np.ones(len(vertex))
    trk1_pos = (vertex["trk1_charge"] > 0.5) & (vertex["trk2_charge"] < -0.5)
    trk1_neg = (vertex["trk1_charge"] < -0.5) & (vertex["trk2_charge"] > 0.5)
    two_same_charge = (vertex["trk1_charge"]  == vertex["trk2_charge"])
    assert not np.any(two_same_charge) #check that the tracks really did both have opposite charge
    trk_momentum[trk1_pos] = vertex["trk1_p"][trk1_pos]
    trk_momentum[trk1_neg] = vertex["trk2_p"][trk1_neg]
    return trk_momentum

def neg_track_momentum(vertex):
    trk_momentum = np.ones(len(vertex))
    trk1_pos = (vertex["trk1_charge"] > 0.5) & (vertex["trk2_charge"] < -0.5)
    trk1_neg = (vertex["trk1_charge"] < -0.5) & (vertex["trk2_charge"] > 0.5)
    two_same_charge = (vertex["trk1_charge"]  == vertex["trk2_charge"])
    assert not np.any(two_same_charge) #check that the tracks really did both have opposite charge
    trk_momentum[trk1_pos] = vertex["trk2_p"][trk1_pos] #track 1 positive -> trk2 is negative
    trk_momentum[trk1_neg] = vertex["trk1_p"][trk1_neg] #track 1 negative -> trk1 is negative
    return trk_momentum

#use the mass hypothesis of the higher-momentum track and it's momentum to select lambda and lambda_bars
def sel_lambda(vertex):
    #get the charge of the higher-pt track
    return (vertex["trk1_charge"] > 0.5) 

from variables_identified import calc_cos_theta
def tight_cos_theta_phi(vertex):
    return calc_cos_theta.eval(vertex) > 0.999 #This probably doesn't work for phis...
sel_tight_cos_theta_phi = Calculation(tight_cos_theta_phi, calc_cos_theta.branches)

def rxy_phi(vertex):
    return (vertex["vertex_Rxy"] < 450) & (vertex["vertex_Rxy"] > 4.0)
sel_rxy_phi = Calculation(rxy_phi, ["vertex_Rxy"])

def pt_phi(vertex):
    return (vertex["vertex_pt"] > 0.1)
sel_pt_phi = Calculation(pt_phi, ["vertex_pt"])

