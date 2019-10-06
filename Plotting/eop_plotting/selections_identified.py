from calculation import Calculation

## Take a look here for the run-1 eop selections: https://journals.aps.org/prd/pdf/10.1103/PhysRevD.85.012001

def chi_square_fifteen(vertex):
    return vertex["vertex_chiSquared"] < 15
sel_chi_square_fifteen = Calculation(chi_square_fifteen, ["vertex_chiSquared"])

from variables_identified import calc_cos_theta
def tight_cos_theta(vertex):
    return calc_cos_theta.eval(vertex) > 0.999
sel_tight_cos_theta = Calculation(tight_cos_theta, calc_cos_theta.branches)
