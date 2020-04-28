import ROOT
from PileupHistogrammer import *

def process_EOP_tree(file_list, out_name, tree_name):

  if 'mc16_13TeV' in out_name:
    md = metadata(isData=False)
  else:
    md = metadata(isData=True)

  hist_cont_list = load_histograms(md)
  needed_branches = []
  for hist_container in hist_cont_list:
    needed_branches += hist_container.needed_branches
  needed_branches = list(set(needed_branches))

  in_tree = ROOT.TChain(tree_name)
  for in_file_name in file_list:
    in_tree.Add(in_file_name)

  var_dict = define_and_connect_vars(in_tree, needed_branches)

  for hist_container in hist_cont_list:
    hist_container.connect_vars(var_dict)
 
  n_events = in_tree.GetEntries()
  print( "Will being processing of", n_events, "events.")
  for iEntry in range(n_events):
    var_dict['calc_trk_eta'] = -999
    var_dict['calc_trk_phi'] = -999
    in_tree.GetEntry(iEntry, 1)
    for hist_container in hist_cont_list:
      hist_container.fill_hists()
    if(iEntry%10000==0):
      print ("Event",iEntry)
    if(iEntry > 10000):
      break

  out_file = ROOT.TFile(out_name+".root", "RECREATE")
  for hist_container in hist_cont_list:
    hist_container.Write(out_file)
  out_file.Close()



if __name__ == "__main__":

  print("Starting script")
  import argparse
  parser = argparse.ArgumentParser(description='Simplifed E/p plotting script for pileup-related studies.')
  parser.add_argument('--tree_name', default='LA_EoverP_InDetTrackParticlesSortedLooseIsolatedVertexAssociated_tree', help='Name of the input TTree')
  parser.add_argument('--job_index', type=int, default='-1', help='Index of the job when using condor, used for input file list')
  
  args = parser.parse_args()

  #Determine file-list for condor jobs
  import pickle
  file_list_pickle_name = "job_file_list.pickle"
  file_list = pickle.load( open( file_list_pickle_name, "rb" ) )
  out_name, input_file_list = file_list[args.job_index]
  print("Imported pickle")

  #Run script
  process_EOP_tree( file_list=input_file_list, out_name=out_name, tree_name = args.tree_name)
  print ("Finished process_EOP_tree")
