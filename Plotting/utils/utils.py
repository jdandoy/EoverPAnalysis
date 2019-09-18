import os
import ROOT
import glob

#the location of the files on eos
dir_inclusive="/eos/atlas/atlascerngroupdisk/perf-jets/EoverP/v01_tuples/"
dir_identified="/eos/atlas/atlascerngroupdisk/perf-jets/EoverP/v01_secondaries/"

# a dictionary to store the directories
directories = {}
directories["inclusive"] = dir_inclusive
directories["test"] = dir_inclusive
directories["identified"] = dir_identified
directories["inclusive_hadiso"] = dir_inclusive

inclusive_hadiso_files = {}
inclusive_hadiso_files["LowMuData"] = ["user.luadamek.data17_13TeV.00341294.physics_MinBias.calibhits_v01_hist",\
"user.luadamek.data17_13TeV.00341312.physics_MinBias.calibhits_v01_hist",\
"user.luadamek.data17_13TeV.00341419.physics_MinBias.calibhits_v01_hist",\
"user.luadamek.data17_13TeV.00341534.physics_MinBias.calibhits_v01_hist",\
"user.luadamek.data17_13TeV.00341615.physics_MinBias.calibhits_v01_hist",\
"user.luadamek.data17_13TeV.00341649.physics_MinBias.calibhits_v01_hist"]
inclusive_hadiso_files["PythiaJetJet"] = [\
"user.luadamek.mc16_13TeV.361020.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ0W.calibhits_v01_hist",\
"user.luadamek.mc16_13TeV.361021.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1W.calibhits_v01_hist",\
"user.luadamek.mc16_13TeV.361022.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ2W.calibhits_v01_hist"]

# a dictionary to store all of the filenames for each channel
inclusive_files = {}
inclusive_files["LowMuData"] = ["user.luadamek.data17_13TeV.00341294.physics_MinBias.calibhits_v01_hist",\
"user.luadamek.data17_13TeV.00341312.physics_MinBias.calibhits_v01_hist",\
"user.luadamek.data17_13TeV.00341419.physics_MinBias.calibhits_v01_hist",\
"user.luadamek.data17_13TeV.00341534.physics_MinBias.calibhits_v01_hist",\
"user.luadamek.data17_13TeV.00341615.physics_MinBias.calibhits_v01_hist",\
"user.luadamek.data17_13TeV.00341649.physics_MinBias.calibhits_v01_hist"]
inclusive_files["LowMuDataTightIso"] = ["user.luadamek.data17_13TeV.00341294.physics_MinBias.calibhits_v01_hist",\
"user.luadamek.data17_13TeV.00341312.physics_MinBias.calibhits_v01_hist",\
"user.luadamek.data17_13TeV.00341419.physics_MinBias.calibhits_v01_hist",\
"user.luadamek.data17_13TeV.00341534.physics_MinBias.calibhits_v01_hist",\
"user.luadamek.data17_13TeV.00341615.physics_MinBias.calibhits_v01_hist",\
"user.luadamek.data17_13TeV.00341649.physics_MinBias.calibhits_v01_hist"]
inclusive_files["SinglePion"] = [\
"user.luadamek.mc16_13TeV.428001.ParticleGun_single_piplus_logE0p2to2000.singlepartapr30_calibhits_hist",\
"user.luadamek.mc16_13TeV.428002.ParticleGun_single_piminus_logE0p2to2000.singlepartapr30_calibhits_hist"]
inclusive_files["SinglePionTightIso"] = [\
"user.luadamek.mc16_13TeV.428001.ParticleGun_single_piplus_logE0p2to2000.singlepartapr30_calibhits_hist",\
"user.luadamek.mc16_13TeV.428002.ParticleGun_single_piminus_logE0p2to2000.singlepartapr30_calibhits_hist"]
inclusive_files["PythiaJetJet"] = [\
"user.luadamek.mc16_13TeV.361020.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ0W.calibhits_v01_hist",\
"user.luadamek.mc16_13TeV.361021.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1W.calibhits_v01_hist",\
"user.luadamek.mc16_13TeV.361022.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ2W.calibhits_v01_hist"]
inclusive_files["PythiaJetJetTightIso"] = [\
"user.luadamek.mc16_13TeV.361020.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ0W.calibhits_v01_hist",\
"user.luadamek.mc16_13TeV.361021.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1W.calibhits_v01_hist",\
"user.luadamek.mc16_13TeV.361022.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ2W.calibhits_v01_hist"]
inclusive_files["PythiaJetJetTightIso"] = [\
"user.luadamek.mc16_13TeV.361020.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ0W.calibhits_v01_hist",\
"user.luadamek.mc16_13TeV.361021.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1W.calibhits_v01_hist",\
"user.luadamek.mc16_13TeV.361022.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ2W.calibhits_v01_hist"]
inclusive_files["PythiaJetJetHardScatter"] = [\
"user.luadamek.mc16_13TeV.361020.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ0W.calibhits_v01_hist",\
"user.luadamek.mc16_13TeV.361021.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1W.calibhits_v01_hist",\
"user.luadamek.mc16_13TeV.361022.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ2W.calibhits_v01_hist"]
inclusive_files["PythiaJetJetHardScatterTightIso"] = [\
"user.luadamek.mc16_13TeV.361020.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ0W.calibhits_v01_hist",\
"user.luadamek.mc16_13TeV.361021.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1W.calibhits_v01_hist",\
"user.luadamek.mc16_13TeV.361022.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ2W.calibhits_v01_hist"]
#inclusive_files["PythiaJetJetPionsReweighted"] = [\
#"user.luadamek.mc16_13TeV.361020.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ0W.calibhits_v01_hist",\
#"user.luadamek.mc16_13TeV.361021.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1W.calibhits_v01_hist",\
#"user.luadamek.mc16_13TeV.361022.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ2W.calibhits_v01_hist"]
#inclusive_files["PythiaJetJetPionsPosReweighted"] = [\
#"user.luadamek.mc16_13TeV.361020.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ0W.calibhits_v01_hist",\
#"user.luadamek.mc16_13TeV.361021.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1W.calibhits_v01_hist",\
#"user.luadamek.mc16_13TeV.361022.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ2W.calibhits_v01_hist"]
#inclusive_files["PythiaJetJetPionsNegReweighted"] = [\
#"user.luadamek.mc16_13TeV.361020.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ0W.calibhits_v01_hist",\
#"user.luadamek.mc16_13TeV.361021.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1W.calibhits_v01_hist",\
#"user.luadamek.mc16_13TeV.361022.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ2W.calibhits_v01_hist"]
#inclusive_files["PythiaJetJetKaonsPosReweighted"] = [\
#"user.luadamek.mc16_13TeV.361020.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ0W.calibhits_v01_hist",\
#"user.luadamek.mc16_13TeV.361021.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1W.calibhits_v01_hist",\
#"user.luadamek.mc16_13TeV.361022.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ2W.calibhits_v01_hist"]
#inclusive_files["PythiaJetJetKaonsNegReweighted"] = [\
#"user.luadamek.mc16_13TeV.361020.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ0W.calibhits_v01_hist",\
#"user.luadamek.mc16_13TeV.361021.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1W.calibhits_v01_hist",\
#"user.luadamek.mc16_13TeV.361022.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ2W.calibhits_v01_hist"]
#inclusive_files["PythiaJetJetProtonsPosReweighted"] = [\
#"user.luadamek.mc16_13TeV.361020.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ0W.calibhits_v01_hist",\
#"user.luadamek.mc16_13TeV.361021.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1W.calibhits_v01_hist",\
#"user.luadamek.mc16_13TeV.361022.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ2W.calibhits_v01_hist"]
#inclusive_files["PythiaJetJetProtonsNegReweighted"] = [\
#"user.luadamek.mc16_13TeV.361020.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ0W.calibhits_v01_hist",\
#"user.luadamek.mc16_13TeV.361021.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1W.calibhits_v01_hist",\
#"user.luadamek.mc16_13TeV.361022.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ2W.calibhits_v01_hist"]
#inclusive_files["PythiaJetJetOtherReweighted"] = [\
#"user.luadamek.mc16_13TeV.361020.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ0W.calibhits_v01_hist",\
#"user.luadamek.mc16_13TeV.361021.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1W.calibhits_v01_hist",\
#"user.luadamek.mc16_13TeV.361022.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ2W.calibhits_v01_hist"]
#inclusive_files["PythiaJetJetHardScatter"] = [\
#"user.luadamek.mc16_13TeV.361020.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ0W.calibhits_v01_hist",\
#"user.luadamek.mc16_13TeV.361021.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1W.calibhits_v01_hist",\
#"user.luadamek.mc16_13TeV.361022.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ2W.calibhits_v01_hist"]
#inclusive_files["PythiaJetJetHardScatterTightIso"] = [\
#"user.luadamek.mc16_13TeV.361020.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ0W.calibhits_v01_hist",\
#"user.luadamek.mc16_13TeV.361021.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1W.calibhits_v01_hist",\
#"user.luadamek.mc16_13TeV.361022.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ2W.calibhits_v01_hist"]
identified_files = {}
identified_files["LowMuData"] = [\
"user.luadamek.data17_13TeV.00341294.physics_MinBias.calibhits_v01_ANALYSIS.root",\
"user.luadamek.data17_13TeV.00341312.physics_MinBias.calibhits_v01_ANALYSIS.root",\
"user.luadamek.data17_13TeV.00341419.physics_MinBias.calibhits_v01_ANALYSIS.root",\
"user.luadamek.data17_13TeV.00341534.physics_MinBias.calibhits_v01_ANALYSIS.root",\
"user.luadamek.data17_13TeV.00341615.physics_MinBias.calibhits_v01_ANALYSIS.root",\
"user.luadamek.data17_13TeV.00341649.physics_MinBias.calibhits_v01_ANALYSIS.root"]
identified_files["PythiaJetJet"] = [\
"user.luadamek.mc16_13TeV.361020.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ0W.calibhits_v01_ANALYSIS.root",\
"user.luadamek.mc16_13TeV.361021.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1W.calibhits_v01_ANALYSIS.root",\
"user.luadamek.mc16_13TeV.361022.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ2W.calibhits_v01_ANALYSIS.root"]

files = {}
files["inclusive"] = inclusive_files
files["test"] = {"PythiaJetJet":["user.luadamek.mc16_13TeV.361020.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ0W.calibhits_v01_hist"],\
                 "SinglePion":["user.luadamek.mc16_13TeV.428001.ParticleGun_single_piplus_logE0p2to2000.singlepartapr30_calibhits_hist"],\
                 "LowMuData":["user.luadamek.data17_13TeV.00341294.physics_MinBias.calibhits_v01_hist"],\
                }
files["identified"] = identified_files 
files["inclusive_hadiso"] = inclusive_hadiso_files 

def get_files(flavour):
    '''
    Get a list of files associated with a given type of file (e.g. the inclusive files on eos)
    '''
    assert flavour in directories
    assert flavour in files

    directory = directories[flavour]
    all_files = files[flavour]

    to_return = {}
    for key in all_files:
        to_return[key] = [os.path.join(directory, f) for f in all_files[key]]

    return to_return

def tchain_files_together(tree_name, channel_to_filelist, on_eos = True):
    '''
    Given a tree_name, and a dictionary of channel to file list, return a dictionary of channel to filename to tchain.
    '''
    trees = {}
    print("\n"*10)
    print("Chaining files together for {}".format(list(trees.keys())))
    for channel in channel_to_filelist:
        trees[channel] = {}
        files = channel_to_filelist[channel]
        print("For channel {}, found files {}".format(channel, files))
        for f in files:
            #create the tchain for these files
            assert f not in trees[channel]
            trees[channel][f] = ROOT.TChain(tree_name)

            #check if this file was a directory or a file
            if os.path.isfile(f):
                print("For channel {}, and file {}, found files {}".format(channel, f, f))
                if on_eos:
                    trees[channel][f].Add("root://eosatlas/" + f)
                else:
                    trees[channel][f].Add(f)

            else: #this was a directory
                #go and get all of the files in the directory
                wildcards = ["*.root", "*.root*"]
                files = []
                for wild_card in wildcards:
                    files += glob.glob(os.path.join(f, wild_card))


                unique_files = []
                for raw_file in files:
                    file_with_path = os.path.join(f, raw_file)
                    if file_with_path not in unique_files and os.path.isfile(file_with_path):
                        assert "//" not in file_with_path
                        print("Found file {}".format(file_with_path))
                        if on_eos:
                            trees[channel][f].Add("root://eosatlas/" + file_with_path)
                        else:
                            trees[channel][f].Add(file_with_path)
    return trees

def generate_partitions(trees, NPartitions):
    '''
    generate a dictionary of channel to file to list of tuples with information about what events to read for each partition
    '''
    partitions = {}
    for channel in trees:
        assert channel not in partitions
        partitions[channel] = {}
        for f in trees[channel]:
            assert f not in partitions[channel]
            tree = trees[channel][f] 
            entries = tree.GetEntries()
            #OK we need to create n event splits from 0 to entries
            step = int(float(entries)/float(NPartitions)) - 1 
            cuts = []
            #are there enough entires to warrant a split?
            if step > 50:
                cuts.append((0, step))
                while cuts[-1][-1] < entries:
                    last_value = cuts[-1][-1]
                    cuts.append( (last_value, last_value + step))
                cuts = cuts[:-2]
                cuts.append( (cuts[-1][-1], entries))
                assert len(cuts) == NPartitions
            else:
                cuts.append((0.0, entries))
                for i in range(1, NPartitions):
                    cuts.append( (entries, entries) )
            print("Found partitions for channel {}, and file {}, and they were {}".format(channel, f, cuts))
            partitions[channel][f] = cuts

    return partitions

