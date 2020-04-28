import ROOT, array, sys
import copy, math

#TODO find sample weights?
#TODO break up into multiple files
#TODO selection with no hadronic clusters
#TODO are muons a candidate for low cluster efficiency?


#Main function to create all histograms and return them as a list of HistContainer objects
def load_histograms(meta_data):
  

  eta_dimension = create_eta_dimension()


  if meta_data.isData:
    particleID_types = {'AllParticle':   lambda vd: True }
  else:
    particleID_types = {'AllParticle':   lambda vd: True,
                        'PionAll':       lambda vd: abs(vd['trk_truthPdgId'][0])==211, 
                        'PionPos':       lambda vd: vd['trk_truthPdgId'][0]==211,
                        'PionNeg':       lambda vd: vd['trk_truthPdgId'][0]==-211,
                        'KaonPos':       lambda vd: vd['trk_truthPdgId'][0]==321,
                        'KaonNeg':       lambda vd: vd['trk_truthPdgId'][0]==-321,
                        'ProtonPos':     lambda vd: vd['trk_truthPdgId'][0]==2212,
                        'ProtonNeg':     lambda vd: vd['trk_truthPdgId'][0]==-2212, 
                        'OtherParticle': lambda vd: not (abs(vd['trk_truthPdgId'][0]) in [211,321,2212])
                       }


  hist_cont_list = []
  for particle_name, particle_selection in particleID_types.items():
    ClusR = '200'
    hist_cont = HistContainer(particle_name+"_Rad200_nTRT20", clus_radius=ClusR)
  
    global_selection = lambda vd: vd['trk_nTRT'][0] >= 20
    hist_cont.add_branches([('trk_nTRT', 'i')])

    #Couldn't use dictionary values here (everyone was overwritten and only last selection in the loop was used)
    if particle_name == 'PionAll':
      hist_cont.selection_func = lambda vd: ( global_selection(vd) and abs(vd['trk_truthPdgId'][0])==211 )
    elif particle_name == 'PionPos':
      hist_cont.selection_func = lambda vd: ( global_selection(vd) and vd['trk_truthPdgId'][0]==211 )
    elif particle_name == 'PionNeg':
      hist_cont.selection_func = lambda vd: ( global_selection(vd) and vd['trk_truthPdgId'][0]==-211 )
    elif particle_name == 'PionNeg':
      hist_cont.selection_func = lambda vd: ( global_selection(vd) and vd['trk_truthPdgId'][0]==-211 )
    elif particle_name == 'KaonPos':
      hist_cont.selection_func = lambda vd: ( global_selection(vd) and vd['trk_truthPdgId'][0]==321 )
    elif particle_name == 'KaonNeg':
      hist_cont.selection_func = lambda vd: ( global_selection(vd) and vd['trk_truthPdgId'][0]==-321 )
    elif particle_name == 'ProtonPos':
      hist_cont.selection_func = lambda vd: ( global_selection(vd) and vd['trk_truthPdgId'][0]==2212 )
    elif particle_name == 'ProtonNeg':
      hist_cont.selection_func = lambda vd: ( global_selection(vd) and vd['trk_truthPdgId'][0]==-2212 )
    elif particle_name == 'OtherParticle':
      hist_cont.selection_func = lambda vd: ( global_selection(vd) and not (abs(vd['trk_truthPdgId'][0]) in [211,321,2212]) )
    else:
      hist_cont.selection_func = lambda vd: global_selection(vd)


    if not meta_data.isData:
      this_ph = PileupHist("ParticleType")
      this_ph.set_func_x( lambda vd: sign(vd['trk_truthPdgId'][0]) * ([0,211,321,2212].index(abs(vd['trk_truthPdgId'][0])) if (abs(vd['trk_truthPdgId'][0]) in [0,211,321,2212]) else 4) )
      this_ph.set_axis_x( "PDGID", [9,-4,4] )
      this_ph.add_branches( [('trk_truthPdgId', 'i')] )
      hist_cont.add(this_ph)
  
    ###--- 2D E/p vs P  ---###
    this_ph = PileupHist("EOPvP")
    this_ph.set_func_y( lambda vd: (vd['trk_ClusterEnergy_EM_'+ClusR][0] + vd['trk_ClusterEnergy_HAD_'+ClusR][0])/vd['trk_p'][0] )
    this_ph.set_func_x( lambda vd: vd['trk_p'][0] )
    this_ph.set_axis_y( "E/p", [150, -0.5, 2.5] )
    this_ph.set_axis_x( "track p [GeV]", [180, 0., 36.] )
    this_ph.add_branches( [('trk_ClusterEnergy_EM_'+ClusR,'f'), ('trk_ClusterEnergy_HAD_'+ClusR,'f'), ('trk_p','f') ] )
    this_ph.add_dimension(eta_dimension)
    hist_cont.add(this_ph)
  
    ###--- 2D Annulus Energy vs P  ---###
    this_ph = PileupHist("AnnEMvP")
    this_ph.set_func_y( lambda vd: (vd['trk_ClusterEnergy_EM_200'][0] - vd['trk_ClusterEnergy_EM_100'][0])  )
    this_ph.set_func_x( lambda vd: vd['trk_p'][0] )
    this_ph.set_axis_y( "E_{EM Annulus}", [84, -1, 20.] )
    this_ph.set_axis_x( "track p [GeV]", [350, 0., 35.] )
    this_ph.add_branches( [('trk_ClusterEnergy_EM_200','f'), ('trk_p','f'), 'trk_ClusterEnergy_EM_100'] )
    this_ph.add_dimension(eta_dimension)
    hist_cont.add(this_ph)
  
    ###--- Track kinematics for all tracks ---###
    this_ph = PileupHist("trk_p")
    this_ph.set_func_x( lambda vd: vd['trk_p'][0] )
    this_ph.set_axis_x( "track p [GeV]", [350, 0., 35.] )
    this_ph.add_dimension( eta_dimension )
    this_ph.add_branches( [('trk_p','f')] )
    hist_cont.add(this_ph)
    ###--- Track p when matched to a cluster ---###
    this_ph = copy.deepcopy(this_ph) #Make a copy
    this_ph.name = 'trkP_clusmatch'
    #At least one matched cluster to the track
    this_ph.selection_func = lambda vd: (vd['trk_nclusters_EM_'+ClusR][0] + vd['trk_nclusters_HAD_'+ClusR][0]) > 0.5
    hist_cont.add(this_ph)
    
    this_ph = PileupHist("trk_pt")
    this_ph.set_func_x( lambda vd: vd['trk_pt'][0] )
    this_ph.set_axis_x( "track p_{t} [GeV]", [35, 0., 35.] )
    this_ph.add_branches( [('trk_pt','f')] )
    hist_cont.add(this_ph)
    
    this_ph = PileupHist("trk_etaID")
    this_ph.set_func_x( lambda vd: vd['trk_etaID'][0] )
    this_ph.set_axis_x( "track #eta_{ID}", [25, -2.5, 2.5] )
    this_ph.add_branches( [('trk_etaID','f')] )
    hist_cont.add(this_ph)
    
    this_ph = PileupHist("trk_eta")
    this_ph.set_func_x( lambda vd: get_eta(vd) )
    this_ph.set_axis_x( "track #eta_{EM2}", [25, -2.5, 2.5] )
    this_ph.add_branches( ['trk_etaEMB2', 'trk_etaEME2'] )
    hist_cont.add(this_ph)
    
    this_ph = PileupHist("trk_phi")
    this_ph.set_func_x( lambda vd: get_phi(vd) )
    this_ph.set_axis_x( "track #phi_{EM2}", [50, -5., 5.] )
    this_ph.add_branches( ['trk_phiEMB2', 'trk_phiEME2'] )
    hist_cont.add(this_ph)
    
    this_ph = PileupHist("trkDef_dR")
    this_ph.set_func_x( lambda vd: get_trkDef_dR(vd) )
    this_ph.set_axis_x( "track dR b/w ID and EM2", [40, 0., 2.] )
    this_ph.add_branches( ['trk_phiEMB2', 'trk_phiEME2', 'trk_etaEMB2', 'trk_etaEME2', 'trk_etaID', 'trk_phiID'] )
    hist_cont.add(this_ph)
    
    this_ph = PileupHist("actualmu")
    this_ph.set_func_x( lambda vd: vd['trk_actualmu'][0] )
    this_ph.set_axis_x( "Actual #mu", [30, 0., 15.] )
    this_ph.add_branches( ['trk_actualmu'] )
    hist_cont.add(this_ph)
    
    this_ph = PileupHist("NPV2")
    this_ph.set_func_x( lambda vd: vd['trk_NPV_2'][0] )
    this_ph.set_axis_x( "N_{PV}", [15, 0., 15.] )
    this_ph.add_branches( ['trk_NPV_2'] )
    hist_cont.add(this_ph)
    
    this_ph = PileupHist("trk_d0")
    this_ph.set_func_x( lambda vd: vd['trk_d0'][0] )
    this_ph.set_axis_x( "track d_{0}", [120, -6., 6.] )
    this_ph.add_branches( ['trk_d0'] )
    hist_cont.add(this_ph)
    
    this_ph = PileupHist("trk_z0sintheta")
    this_ph.set_func_x( lambda vd: vd['trk_z0sintheta'][0] )
    this_ph.set_axis_x( "track z_{0}sin#theta", [100, -200., 200.] )
    this_ph.add_branches( ['trk_z0sintheta'] )
    hist_cont.add(this_ph)
    

  
    ### --- Ncluster vs Sum Cluster Energy for each good track ---###
    this_ph = PileupHist("Cluster_NvE")
    this_ph.set_func_x( lambda vd: (vd['trk_nclusters_EM_'+ClusR][0] + vd['trk_nclusters_HAD_'+ClusR][0]) )
    this_ph.set_func_y( lambda vd: (vd['trk_ClusterEnergy_EM_'+ClusR][0] + vd['trk_ClusterEnergy_HAD_'+ClusR][0]) )
    this_ph.set_axis_x( "N_{clusters}", [15, 0, 15] )
    this_ph.set_axis_y( "Sum E_{clusters} [GeV]", [110, -10, 100] )
    this_ph.add_branches( [('trk_nclusters_EM_'+ClusR, 'i'), ('trk_nclusters_HAD_'+ClusR, 'i'), 'trk_ClusterEnergy_EM_'+ClusR, 'trk_ClusterEnergy_HAD_'+ClusR] )
    this_ph.add_dimension( eta_dimension )
    hist_cont.add(this_ph)
  
    hist_cont.create_hists()
    hist_cont_list.append(hist_cont)


  return hist_cont_list

###--- Histogramming Classes ---###

#This class holds histograms of various types.
#The fill call will fill them when given a dictionary of values
#Histograms of the same distribution can be made for multiple regions of a PlotDimension (i.e. eta bins)
class PileupHist:
  def __init__(self, name):
    self.name = name
    self.dir_name = ""

    self.axis_x = []
    self.axis_name_x = ""
    self.axis_y = []
    self.axis_name_y = ""
    
    self.fill_func_x = -999
    self.fill_func_y = -999

    self.is2D = False
    self.selection_func = lambda vd: True
    self.needed_branches = [ 'trkWeight' ]
    self.region_edges = []
    self.region_funcs = []
    self.plot_dimension = None

    self.pass_one_region = True

    self.hist_list = []
    self.event_weight = lambda vd: vd['trkWeight'][0]

  def add_branches(self, branches):
    self.needed_branches += branches
  
  def get_event_weight(self, vd):
    return self.event_weight(vd)

  #This is a dimension to make histograms from.  By definition, an event can only fill one bin of the dimension.
  #e.g. add a track eta binning, creating one histogram for each track eta region
  def add_dimension(self, dimension):
    self.plot_dimension = dimension
    self.add_branches( self.plot_dimension.needed_branches )

  def set_axis_x(self, axis_name, x_vals, is_variable=False, var_type='f'):
    if is_variable:
      x_array = array.array(var_type, x_vals)
      self.axis_x = [len(x_array)-1, x_array]
    else: 
      self.axis_x = x_vals
    self.axis_name_x = axis_name

  def set_axis_y(self, axis_name, y_vals, is_variable=False, var_type='f'):
    if is_variable:
      y_array = array.array(var_type, y_vals)
      self.axis_y = [len(x_array)-1, y_array]
    else: 
      self.axis_y = y_vals
    self.axis_name_y = axis_name

  def set_func_x(self, func):
    self.fill_func_x = func
  
  def set_func_y(self, func):
    self.fill_func_y = func
    self.is2D = True

  def fill_func(self, vd):
    if self.is2D:
      return [self.fill_func_x(vd), self.fill_func_y(vd)]
    else:
      return [self.fill_func_x(vd)]


  def branch_list(self):
    return self.needed_branches
 
  def clear_hists(self):
    self.hist_list.clear()

  def create_hists(self, prepend):
    self.dir_name = prepend
    this_name = prepend+'_'+self.name
    if( self.plot_dimension == None ):
      self.create_one_hist(this_name)
    else:
      for iR, region_name in enumerate(self.plot_dimension.region_names):
        new_name = this_name+'_'+region_name
        self.create_one_hist(new_name)

  def create_one_hist(self, hist_name):
    if( self.is2D ):
      this_hist = ROOT.TH2F(hist_name, hist_name, *self.axis_x, *self.axis_y) 
    else:
      this_hist = ROOT.TH1F(hist_name, hist_name, *self.axis_x) # "splat" (unpack) the x-axis arguments

    self.hist_list.append(this_hist)


  def fill_hists(self, vd):
    if( self.selection_func(vd) == False ):
      return

    for iH, hist in enumerate(self.hist_list):
      if( self.plot_dimension == None or self.plot_dimension.region_funcs[iH](vd) == True ):
        hist.Fill( *self.fill_func(vd), self.event_weight(vd) )
        if( self.pass_one_region ):
          break #Only allow 1 region to pass!
    return
  
  def Write(self):
    for hist in self.hist_list: 
      #hist.SetDirectory(self.dir_name)
      hist.Write()

#This class holds multiple PileupHist objects for a given selection
#For example. there could be one HistContainer for each truth particle type
class HistContainer:
  def __init__(self, name, clus_radius='200'):
    self.name = name
    self.selection_func = lambda vd: True
    self.clus_radius = clus_radius
    self.PH_list = []
    self.needed_branches = []
    self.var_dict = None

  def add(self, this_PileupHist):
    self.PH_list.append(this_PileupHist)
  
  def add_branches(self, branches):
    self.needed_branches += branches

  def create_hists(self):
    for iP, pileup_hist in enumerate(self.PH_list):
      pileup_hist.create_hists(self.name)
      self.add_branches( pileup_hist.needed_branches )

    self.needed_branches = list(set(self.needed_branches)) #unique list only

  def connect_vars(self, var_dict):
    self.var_dict = var_dict

  #Call fill function on all PileupHist objects
  def fill_hists(self):
    if( self.selection_func(self.var_dict) == False ):
      return
    for iP, pileup_hist in enumerate(self.PH_list):
      pileup_hist.fill_hists(self.var_dict)
    return
  
  def Write(self, out_file):
    out_file.mkdir(self.name);
    out_file.cd(self.name)
    for iP, pileup_hist in enumerate(self.PH_list):
      pileup_hist.Write()
    out_file.cd('..')

###--- Plotting Dimensions ---###

#This class provides a dimension for which to provide a set of histograms
#i.e. we may want a histogram for various eta regions, so this will define the eta regions and the selection for each one
#Can be given to the PileupHist object to create histograms for this dimension
class PlotDimension:
  def __init__(self, name):
    self.name = name
    self.region_names = []
    self.region_funcs = []
    self.needed_branches = []

  def add_region(self, name, function):
    self.region_names.append(name)
    self.region_funcs.append(function)

def create_eta_dimension():
  eta_edges = [0, 0.2, 0.7, 1.3, 1.8, 2.5]
  eta_dimension = PlotDimension("eta")
  eta_dimension.needed_branches = ['trk_etaEMB2', 'trk_etaEME2']
  eta_dimension.add_region('etaAll', lambda vd: True )
  for iEta in range(len(eta_edges)-1):
    low_str, high_str = str(int(eta_edges[iEta]*10)), str(int(eta_edges[iEta+1]*10))
    eta_dimension.add_region('eta{0}to{1}'.format(low_str,high_str), lambda vd: abs(get_eta(vd)) >= eta_edges[iEta] and abs(get_eta(vd)) < eta_edges[iEta+1] )
  return eta_dimension


###--- Helper Functions ---###

def get_eta(vd):
  if vd['calc_trk_eta'] == -999:
    if( abs(vd["trk_etaEMB2"][0]) < 1000.0 ):
      vd['calc_trk_eta'] = vd["trk_etaEMB2"][0]
    elif( abs(vd["trk_etaEME2"][0]) < 1000.0 ):
      vd['calc_trk_eta'] = vd["trk_etaEME2"][0]
    else:
      print( "Couldn't find eta?" )
      vd['calc_trk_eta'] = -999
  
  return vd['calc_trk_eta']

def get_phi(vd):
  if vd['calc_trk_phi'] == -999:
    if( abs(vd["trk_phiEMB2"][0]) < 1000.0 ):
      vd['calc_trk_phi'] = vd["trk_phiEMB2"][0]
    elif( abs(vd["trk_phiEME2"][0]) < 1000.0 ):
      vd['calc_trk_phi'] = vd["trk_phiEME2"][0]
    else:
      print( "Couldn't find phi?" )
      vd['calc_trk_phi'] = -999

  return vd['calc_trk_phi']

def get_trkDef_dR(vd):
  dPhi = ROOT.TVector2.Phi_0_2pi( get_phi(vd) - vd['trk_phiID'][0] )
  dPhi = min(dPhi, 2*math.pi-dPhi)
  dEta = get_eta(vd) - vd['trk_etaID'][0]
  return math.sqrt(dEta**2 + dPhi**2)

class metadata():
  def __init__(self, isData=True):
    self.isData = isData

def define_and_connect_vars(in_tree, branch_list):
  in_tree.SetBranchStatus('*', 0)
  var_dict = {}
  for br_name in branch_list:
    if type(br_name) == tuple:
      br_name, br_type = br_name #If a tuple, unwrap it
      var_dict[br_name] = array.array(br_type, [0])
    else:
      var_dict[br_name] = array.array('f', [0])
    in_tree.SetBranchStatus(br_name, 1)
    in_tree.SetBranchAddress(br_name, var_dict[br_name])
  return var_dict

def sign(value):
  if( value == 0 ):
    return 1.
  else:
    return value/abs(value)
