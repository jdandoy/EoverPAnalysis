# E/p secondary vertex analysis for Run 2
# Matt LeBlanc (matt.leblanc@cern.ch)

from xAH_config import xAH_config
c = xAH_config()

c.output("ANALYSIS")

# input containers
trks = "InDetTrackParticles"
trks_loose = trks + "Loose"
trks_tight = trks + "Tight"
trks_loose_isolated = trks_loose + "Isolated"
trks_tight_isolated = trks_tight + "Isolated"
trks_loose_isolated_vertex = trks_loose_isolated + "VertexAssociated"
trks_tight_isolated_vertex = trks_tight_isolated + "VertexAssociated"

# Set up all the algorithms
c.algorithm("BasicEventSelection", {"m_name": "BasicEventSelection",				    
				    "m_applyGRLCut": False,
				    "m_doPUreweighting": False,
				    "m_applyPrimaryVertexCut": True,
				    "m_applyEventCleaningCut": True,
				    "m_applyCoreFlagsCut": False,
				    "m_applyTriggerCut": False,
				    #"m_useCutflow": True,
				    "m_GRLxml": "EoverPAnalysis/data17_13TeV.periodN_DetStatus-v98-pro21-16_Unknown_PHYS_StandardGRL_All_Good_25ns_ignore_GLOBAL_LOWMU_for_specific_run341294.xml",
				    "m_lumiCalcFileNames": "EoverPAnalysis/ilumicalc_histograms_HLT_mb_sptrk_341294_OflLumi-13TeV-010.root",
				    "m_PRWFileNames": "EoverPAnalysis/ntup_prw_36102_JZ012.root",
				    "m_triggerSelection": "HLT_mb_sptrk",
				    "m_PVNTrack": 2,
				    "m_useMetaData": False,
				    "m_checkDuplicatesData": False,
				    "m_checkDuplicatesMC": False})

# Fill histograms with tracking details, passing only basic event selection
c.algorithm("TrackHistsAlgo", {"m_name": "Tracks_BasicEvtSel",
			       "m_inContainerName": trks,
			       "m_detailStr": "2D IPDetails HitCounts Chi2Details",
			       "m_msgLevel": "info"})


# track selection algorithm
c.algorithm("InDetTrackSelectionToolAlgo", {"m_name": "Sel_" + trks_loose,
					    "m_inputTrackContainer": trks,
					    "m_minPt": 0.5,
					    "m_CutLevel": "Loose",
					    "m_outputTrackContainer": trks_loose,
					    "m_msgLevel": "info"})

c.algorithm("InDetTrackSelectionToolAlgo", {"m_name": "Sel_" + trks_tight,
                                            "m_inputTrackContainer": trks,
                                            "m_minPt": 0.5,
                                            "m_CutLevel": "Tight",
                                            "m_outputTrackContainer": trks_tight,
                                            "m_msgLevel": "info"})

# Fill histograms with tracking details, after LoosePrimary selection
c.algorithm("TrackHistsAlgo", {"m_name": "TrackHist_" + trks_loose,
			       "m_inContainerName": trks_loose,
			       "m_detailStr": "2D IPDetails HitCounts Chi2Details",
			       "m_msgLevel": "info"})

c.algorithm("TrackExtrapolationIsolationTool", {"m_name": "TrackIso_" + trks_loose,
						"m_inputTrackContainer": trks_loose,
						"m_outputTrackContainer": trks_loose_isolated,
						"m_trkIsoDRmax": 0.4,
						"m_msgLevel": "info"})

# Fill histograms with tracking details, after LoosePrimary selection
c.algorithm("TrackHistsAlgo", {"m_name": "TrackHist_" + trks_loose_isolated,
			       "m_inContainerName": trks_loose_isolated,
			       "m_detailStr": "2D IPDetails HitCounts Chi2Details",
			       "m_msgLevel": "info"})

c.algorithm("InDetTrackSelectionToolAlgo", {"m_name": "Sel_" + trks_tight_isolated ,
					    "m_inputTrackContainer": trks_loose_isolated,
					    "m_outputTrackContainer": trks_tight_isolated,
					    "m_CutLevel": "TightPrimary",
					    "m_msgLevel": "info"})

# Fill histograms with tracking details, after LoosePrimary selection
c.algorithm("TrackHistsAlgo", {"m_name": "TrackHist_" + trks_tight_isolated,
			       "m_inContainerName": trks_tight_isolated,
			       "m_detailStr": "2D IPDetails HitCounts Chi2Details",
			       "m_msgLevel": "info"})


#c.setalg("TightTrackVertexAssociationToolAlgo", {"m_name":"TrackVertexAssociationTool",\
#                                            "m_inputTrackContainer": trks_loose_isolated,\
#                                            "m_outputTrackContainer": trks_loose_isolated_vertex,\
#                                            "m_dzSinTheta_cut": 0.5,
#                                            "m_d0_cut": 0.5,
#                                            "m_msgLevel": "debug",
#                                            })

#c.setalg("TightTrackVertexAssociationToolAlgo", {"m_name":"TrackVertexAssociatedTool",\
#                                            "m_inputTrackContainer": trks_tight_isolated,\
#                                            "m_outputTrackContainer": trks_tight_isolated_vertex,\
#                                            "m_dzSinTheta_cut":0.5,
#                                            "m_d0_cut":0.5,
#                                            "m_msgLevel": "debug",
#                                            })

# Take only tracks associated with secondary vertices

#for sv in ['Lambda','Ks','Phi']:

for sv in ['Lambda','Ks','Phi']:
	#trks_loose_sv = trks + "Loose" + sv
	#trks_tight_sv = trks + "Tight" + sv
	#trks_loose_isolated_sv = trks_loose_isolated + sv
	#trks_tight_isolated_sv = trks_tight_isolated + sv

	for t in [trks_loose, trks_tight, trks_loose_isolated, trks_tight_isolated]:
		c.algorithm("TrackSecondaryVertexAssociationToolAlgo", {"m_name":"TrackVertexAssociationTool"+t+sv,
									"m_inputTrackContainer": t,
									"m_outputTrackContainer": t+sv,
									"m_inputVertexContainer" : sv+"Candidates",
									"m_dzSinTheta_cut": 999.,
									"m_d0_cut": 999.,
									"m_doPV": False,
									"m_msgLevel": "info", 
									})

		c.algorithm("TrackHistsAlgo", {"m_name": "TrackHist_" + t + sv,
					       "m_inContainerName": t+sv,
					       "m_detailStr": "2D IPDetails HitCounts Chi2Details",
					       "m_msgLevel": "info"})
		
		# E/p histograms with LoosePrimary track selection
		#c.algorithm("EoverPTreeAlgo", {"m_name": "EoverP_" + t+sv,
		#"m_inTrackContainerName": t+sv,
		#"m_energyCalibList": "ClusterEnergy,CellEnergy,ClusterEnergyLCW",
		#"m_useCutFlow": True,
		#"m_msgLevel": "info"})
		
		#E/p histograms with LoosePrimary track selection
		c.algorithm("SecondariesTrees", {"m_name": "EoverP_" + t+sv,
						 "label": sv,
						 "isData": True,				
						 "MessageFrequency": 1,
						 "VertexContainer" : sv+"Candidates"
						 })

