# E/p analysis for run 2
# Joakim Olsson (joakim.olsson@cern.ch)

from xAH_config import xAH_config
c = xAH_config()

# input containers
trks = "InDetTrackParticles"

# selected version
trks_loose = trks+"Loose"

''' Set up all the algorithms '''
c.setalg("BasicEventSelection", {"m_name": "BasicEventSelection",
                                 "m_applyGRLCut": False,
                                 "m_doPUreweighting": False,
                                 "m_applyPrimaryVertexCut": True,
                                 "m_applyEventCleaningCut": False,
                                 "m_applyCoreFlagsCut": False,
                                 "m_applyTriggerCut":False,
                                 "m_triggerSelection": "HLT_mb_sptrk",
                                 "m_PVNTrack": 2,
                                 "m_useMetaData": False,
                                 "m_checkDuplicatesData": False,
                                 "m_checkDuplicatesMC": False})

''' Fill histograms with tracking details, passing only basic event selection'''
c.setalg("TrackHistsAlgo", {"m_name": "Tracks_BasicEvtSel",
                            "m_inContainerName": trks,
                            "m_detailStr": "2D IPDetails HitCounts Chi2Details",
                            "m_msgLevel": "info"})

''' Select tracks passing the "LoosePrimary" Tracking CP Recommendations (Moriond 2017)'''
c.setalg("TrackVertexSelection", {"m_name": "TrackSel_LoosePrimary",
                                  "m_inContainerName": trks,
                                  "m_decorateSelectedObjects": False,
                                  "m_createSelectedContainer": True,
                                  "m_pass_min": 1,
                                  "m_minPt": 0.5,
                                  "m_cutLevel": "Loose",
                                  "m_outContainerName": trks_loose,
                                  "m_useCutFlow": True,
                                  "m_msgLevel": "info"})

''' Fill histograms with tracking details, after LoosePrimary selection '''
c.setalg("TrackHistsAlgo", {"m_name": "Tracks_LoosePrimary",
                            "m_inContainerName": trks_loose,
                            "m_detailStr": "2D IPDetails HitCounts Chi2Details",
                            "m_msgLevel": "info"})

#''' Select tracks passing the "LoosePrimary" Tracking CP Recommendations (Moriond 2016), with nTRT > 20'''
# https://twiki.cern.ch/twiki/bin/view/AtlasProtected/TrackingCPMoriond2016
#c.setalg("TrackVertexSelection", {"m_name": "TrackSel_LoosePrimary_nTRTG20",
#                                   "m_inContainerName": trks,
#                                   "m_decorateSelectedObjects": False,
#                                   "m_createSelectedContainer": True,
#                                   "m_pass_min": 1,
#                                   "m_cutLevel": "LoosePrimary",
#                                   "m_minPt": 0.4,
#                                   "m_maxAbsEta": 2.5,
#                                   "m_maxD0": 2.0,
#                                   "m_maxZ0SinTheta": 3.0,
#                                   "m_minNTrtHits": 20,
#                                   "m_outContainerName": trks_loose_ntrtG20,
#                                   "m_useCutFlow": False,
#                                   "m_debug": False})

#''' Fill histograms with tracking details, after LoosePrimary and with TRT selection '''
#c.setalg("TrackHistsAlgo", {"m_name": "Tracks_LoosePrimaryTRT",
#                             "m_inContainerName": trks_loose_ntrtG20,
#                             "m_detailStr": "2D IPDetails HitCounts Chi2Details",
#                             "m_debug": False})

# ''' Select tracks passing the "TightPrimary" Tracking CP Recommendations (Moriond 2016)'''
# # https://twiki.cern.ch/twiki/bin/view/AtlasProtected/TrackingCPMoriond2016
# c.setalg("TrackVertexSelection", {"m_name": "TrackSel_TightPrimary",
#                                   "m_inContainerName": trks,
#                                   "m_decorateSelectedObjects": False,
#                                   "m_createSelectedContainer": True,
#                                   "m_pass_min": 1,
#                                   "m_cutLevel": "TightPrimary",
#                                   "m_minPt": 0.4,
#                                   "m_maxAbsEta": 2.5,
#                                   "m_maxD0": 2.0,
#                                   "m_maxZ0SinTheta": 3.0,
#                                   "m_minNTrtHits": -1,
#                                   "m_outContainerName": trks_tight,
#                                   "m_useCutFlow": False,
#                                   "m_debug": False})
#
# ''' Fill histograms with tracking details, after TightPrimary selection '''
# c.setalg("TrackHistsAlgo", {"m_name": "Tracks_TightPrimary",
#                             "m_inContainerName": trks_tight,
#                             "m_detailStr": "2D IPDetails HitCounts Chi2Details",
#                             "m_debug": False})
#

#### Make E/p plots

for track_container in [trks_loose]:
    for energy_calib in ["ClusterEnergy", "ClusterEnergyLCW", "CellEnergy"]:
        ''' E/p histograms with LoosePrimary track selection'''
        c.setalg("EoverPAnalysis", {"m_name": "EoverP_"+energy_calib + track_container,
                                    "m_fillOutputTree": True,
                                    "m_inTrackContainerName": track_container,
                                    "m_energyCalib": energy_calib, # ClusterEnergy, ClusterEnergyLCW, or CellEnergy
                                    "m_LarEmax": 1.0, #This selection is applied if m_applyTileCuts is true
                                    "m_applyTileCuts": False, #whether or not to cut on TileEfracMin
                                    "m_TileEfracmin": 0.7, #the fraction of energy that should be deposited in the tile calorimeter
                                    "m_doGlobalTileEfracRanges": False, #plots of different TileEfrac selections that could be applied
                                    "m_doGlobalEtaRanges": False, #plots with eta < 0.5, 0.5 < eta < 0.7, and eta < 0.7
                                    "m_doExtraEtaEnergyBinHists": False, 
                                    "m_doGlobalExtraRanges": False,
                                    "m_doGlobalEnergyRanges": False,
                                    "m_doCaloTotal": False,
                                    "m_doCaloEM": False,
                                    "m_doCaloHAD": False,
                                    "m_doBgSubtr" : False,
                                    "m_doTileLayer": False,
                                    "m_trkIsoDRmax": .4,
                                    "m_trkIsoPfrac": .0,
                                    "m_doTrkPcut": False,
                                    "m_trkPmin": 1.0,
                                    "m_trkPmax": 1e8,
                                    "m_doTrkEtacut": False,
                                    "m_trkEtamin": 0.,
                                    "m_trkEtamax": 1.7,
                                    "m_doTrkPtReweighting": False, # do_trkPtRewighting,
                                    "m_trkPtReweightingFile": "",
                                    "m_Pbins": "500, 0, 50",
                                    "m_doPbinsArray": False,
                                    "m_Etabins": "50, 0., 2.5",
                                    "m_doEtabinsArray": False,
                                    "m_doProfileEta": False,
                                    "m_doProfileP": False,
                                    "m_detailStr": "all",
                                    "m_useCutFlow": False,
                                    "m_msgLevel": "info"})

