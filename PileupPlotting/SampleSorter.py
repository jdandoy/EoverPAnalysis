import ROOT, glob, os

dir_inclusive="/eos/atlas/atlascerngroupdisk/perf-jets/EoverP/v03_tuples/"

all_files = {}
all_files['data'] = ["user.luadamek.data17_13TeV.00341294.physics_MinBias.EoverP_Jan7_hist/",\
"user.luadamek.data17_13TeV.00341312.physics_MinBias.EoverP_Jan7_hist/",\
"user.luadamek.data17_13TeV.00341419.physics_MinBias.EoverP_Jan7_hist/",\
"user.luadamek.data17_13TeV.00341534.physics_MinBias.EoverP_Jan7_hist/",\
"user.luadamek.data17_13TeV.00341615.physics_MinBias.EoverP_Jan7_hist/",\
"user.luadamek.data17_13TeV.00341649.physics_MinBias.EoverP_Jan7_hist/"]
all_files['singlePion'] = [\
"user.luadamek.mc16_13TeV.428001.ParticleGun_single_piplus_logE0p2to2000.EoverP_Jan7_hist/",\
"user.luadamek.mc16_13TeV.428002.ParticleGun_single_piminus_logE0p2to2000.EoverP_Jan7_hist/"]
all_files['pythia'] = [\
"user.luadamek.mc16_13TeV.361020.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ0W.EoverP_Jan7_hist/",\
"user.luadamek.mc16_13TeV.361021.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1W.EoverP_Jan7_hist/",\
"user.luadamek.mc16_13TeV.361022.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ2W.EoverP_Jan7_hist/"]


def make_sample_lists(njobs, tree_name, on_eos=True):
  given_files = []
  for filetype, files in all_files.items():
    for this_file in files:
      given_files.append( dir_inclusive+this_file )
  
  
  #Get list of full files, particularly if we were given input directories
  container_files = {}
  #Dfor filetype, files in given_files.items():
  for this_file in given_files:
    container_files[this_file] = []
    #for this_file in files:
    #  assert this_file not in container_files[filetype]
    #  container_files[filetype][this_file] = [] #ROOT.TChain(tree_name)
  
    #check if this file was a directory or a file
    if os.path.isfile(this_file):
      #print("Found file with full path {0}".format(this_file))
      if on_eos:
          container_files[this_file].append('root://eosatlas.cern.ch/' + this_file)
      else:
          container_files[this_file].append(this_file)
    
    else: #this was a directory
      #go and get all of the files in the directory
      if not on_eos:
          wildcards = ["*.root", "*.root*"]
          files = []
          for wild_card in wildcards:
              files += glob.glob(os.path.join(this_file, wild_card))
          files = list(set(files))
      else:
          from XRootD import client
          from XRootD.client.flags import DirListFlags
          xrootd_client = client.FileSystem('root://eosatlas.cern.ch')
          files = [el.name for el in  xrootd_client.dirlist(this_file, DirListFlags.STAT)[1] if ".root" in os.path.split(el.name)[-1]]
          files = [os.path.join(this_file, el) if this_file not in el else el for el in files]
    
      unique_files = []
      for file_with_path in files:
          assert "//" not in file_with_path
          #print("Found file {}".format(file_with_path))
          if on_eos:
             container_files[this_file].append('root://eosatlas.cern.ch/' + file_with_path)
          else:
             container_files[this_file].append(file_with_path)


  #Find number of events of each file and container
  file_count = {}
  container_counts = {}
  for container_name, files in container_files.items():
    #Dfor file_container, files in file_list.items():
    file_count[container_name] = []
    container_counts[container_name] = 0
    for this_file_name in files:
      this_file = ROOT.TFile.Open(this_file_name, 'READ')
      this_tree = this_file.Get(tree_name)

      n_events = this_tree.GetEntries()
      file_count[container_name].append( n_events )
      container_counts[container_name] += n_events

      this_file.Close()

  if njobs < len(container_counts):
    print("Warning, number of requested jobs {0} is less than number of input containers {1}.  Set number of jobs equal to number of containers".format(njobs, len(container_count)) )
    njobs = len(container_counts)

  #Get the number of jobs per container for an even split
  nJobs_per_container = get_container_event_split(njobs, container_counts)

  #Get list of files for each job according to an even event split
  files_per_job = get_file_event_split(nJobs_per_container, file_count, container_files)

  return files_per_job

def make_file_pickles(file_name, files_per_job):
  import pickle
  if not os.path.exists(os.path.dirname(file_name)):
    os.makedirs(os.path.dirname(file_name))
  pickle.dump( files_per_job, open(file_name, "wb" ) )

#This function takes the number of jobs and a dictionary of number events in each container,
#and returns the optimal number of jobs per container
def get_container_event_split(njobs, container_count):
  n_containers = len(container_count)
  total_count = 0
  for container_key, nEvents in container_count.items():
    total_count += nEvents
  
  nJobs_per_container = {}
  nEv_per_job = {}

  #Get containers already at minimum
  nJobs_remaining = njobs 
  target_nEv_per_job = total_count / nJobs_remaining
  for container_key, nEvents in container_count.items():
    if container_count[container_key] <= target_nEv_per_job:
      nJobs_per_container[container_key] = 1
      nEv_per_job[container_key] = float(container_count[container_key])
      total_count -= container_count[container_key]
      nJobs_remaining -= 1
      n_containers -= 1

  #Divy up remaining jobs between containers, rounding down
  target_nEv_per_job = total_count / nJobs_remaining
  for container_key, nEvents in container_count.items():
    if container_key in nJobs_per_container:
      continue #Skip containers with 1 minimum job already 
    print(container_key, nEvents)
    this_nJobs = max( int( nEvents / target_nEv_per_job ), 1)  #Max protects against edge cases
    nJobs_per_container[container_key] = this_nJobs
    nJobs_remaining -= this_nJobs
    nEv_per_job[container_key] = float(nEvents) / this_nJobs
 

  #Give remaining jobs (from rounding down) to containers with the most events per job
  while( nJobs_remaining > 0 ):
    #Get sorted tuples of containers and their nEv_per_job
    sortednEvList = sorted(nEv_per_job.items() ,  key=lambda x: x[1], reverse=True)
    this_key, this_nEv = sortednEvList[0]
    nJobs_per_container[this_key] += 1
    nJobs_remaining -= 1
    nEv_per_job[this_key] = float(container_count[this_key]) / nJobs_per_container[this_key] 

  print(nJobs_per_container)
  return nJobs_per_container

def get_file_event_split(nJobs_per_container, file_count, container_files):
  files_per_job = []

  job_index = 0
  for container_name, files in container_files.items():
      
    out_name = container_name.rstrip('/').split('/')[-1]
    out_name = 'hists.'+'.'.join(out_name.split('.')[2:])
    #Get files and event counts listed in descending event count
    this_nJobs = nJobs_per_container[container_name]
    nEv_files = file_count[container_name]
    nEv_files, files = zip(*sorted(zip(nEv_files, files), reverse=True))
    nEv_files, files = list(nEv_files), list(files)
    
    this_files_per_job = []
    this_nEv_per_job = []
    if len(files) <= this_nJobs:  #If more jobs than files, just give one file per job
      this_files_per_job = files
    else:
      #First give one file per job, starting with largest files
      for iJob in range(this_nJobs):
        this_files_per_job.append( [files.pop(0)] )
        this_nEv_per_job.append( nEv_files.pop(0) )

      #Then, give largest remaning file to the job with smallest event count, constantly updating job ev count list
      for iF, this_file in enumerate(files):
        this_file_nEv = nEv_files[iF]
        this_files_per_job[-1].append( this_file )
        this_nEv_per_job[-1] += this_file_nEv
        this_nEv_per_job, this_files_per_job = zip(*sorted(zip(this_nEv_per_job, this_files_per_job), reverse=True))
        this_nEv_per_job, this_files_per_job = list(this_nEv_per_job), list(this_files_per_job)
  
    #Change to a tuple with the output name for each job at the beginning
    for iJob, file_list in enumerate(this_files_per_job):
      this_files_per_job[iJob] = (out_name+'.Job'+str(job_index), file_list)
      job_index += 1

    files_per_job += this_files_per_job #Add the file list for all jobs for this container

  return files_per_job

if __name__ == "__main__":
  files_per_job = make_sample_lists(njobs=20, tree_name='LA_EoverP_InDetTrackParticlesSortedLooseIsolatedVertexAssociated_tree')
  make_file_pickles(file_name='./job_file_list.pickle', files_per_job=files_per_job)
