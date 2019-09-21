#ifndef eoverpanalysis_SecondariesTrees_H
#define eoverpanalysis_SecondariesTrees_H

#include <iostream>
#include <iomanip>
#include <vector>

#include <AnaAlgorithm/AnaAlgorithm.h>
#include "EventLoop/OutputStream.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TLorentzVector.h"

#include "xAODRootAccess/TStore.h"
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"

#include "xAODTracking/TrackParticleContainer.h"
#include "xAODTruth/TruthParticle.h"

class SecondariesTrees : public EL::AnaAlgorithm
{
 public:
  SecondariesTrees (const std::string& name, ISvcLocator* pSvcLocator);
  ~SecondariesTrees();

  virtual StatusCode initialize () override;
  virtual StatusCode execute () override;
  virtual StatusCode finalize () override;
  const xAOD::TruthParticle* getTruthPtr(const xAOD::TrackParticle* trackParticle); 

  xAOD::TEvent *m_event; //!
  xAOD::TStore *m_store; //!


  // pileup reweighting
  bool m_doCustomPUreweighting = false;
  std::string m_pileupReweightingFile = "pileup_reweighting.root";
  // pt reweighting
  bool m_doTrkPtReweighting = false;
  std::string m_trkPtReweightingFile = "pt_reweighting.root";

 private:
  std::map<std::string, float> m_energyVariablesForTree;
  std::map<std::string, int> m_clusterVariablesForTree;
  std::map<std::string, float> m_trk1_extrapolEta;
  std::map<std::string, float> m_trk1_extrapolPhi;
  std::map<std::string, float> m_trk2_extrapolEta;
  std::map<std::string, float> m_trk2_extrapolPhi;
  std::vector<std::string> m_energyCalibList;
  std::vector<std::string> m_radiusCutList;
  
  // Set by config
  Int_t m_MessageFrequency;
  Bool_t m_isData;
  TString m_label;
  TString m_VertexContainer;
  TString m_TrackContainer;

  // energy calibration, either "ClusterEnergy", "LCWClusterEnergy", or "CellEnergy"
  std::string m_energyCalibCommaList = "ClusterEnergy,CellEnergy,LCWClusterEnergy";
  std::string m_radiusCutCommaList = "100,200";

  // Output variables for tree
  
  Int_t m_DSID;
  Float_t m_mcWeight;
  Float_t m_weight;
  Float_t m_mu;
  Int_t m_n_candidates;
  
  Int_t m_vertex_N;
  Float_t m_vertex_pt;
  Float_t m_vertex_eta;
  Float_t m_vertex_phi;
  Float_t m_vertex_mass;
  Float_t m_vertex_massErr;
  Float_t m_vertex_Rxy;
  Float_t m_vertex_dr;
  Float_t m_vertex_chiSquared;
  Float_t m_vertex_numberDoF; 
  Int_t m_vertex_isIsolatedPairEM;
  Int_t m_vertex_isIsolatedPairHAD;

  /*
  Float_t m_track1_pt;
  Float_t m_track2_pt;
  Float_t m_track1_etaID;
  Float_t m_track1_phiID;
  Float_t m_track2_etaID;
  Float_t m_track2_phiID;
  Float_t m_track1_m;
  Float_t m_track2_m;
  
  Int_t m_trackN;
  */
  /*
  std::vector<float> *m_trackPt = nullptr;
  std::vector<float> *m_trackEta = nullptr;
  std::vector<float> *m_trackPhi = nullptr;
  std::vector<float> *m_trackM = nullptr;
  */

  unsigned long long m_eventNumber;

  // TRACK ONE
  float trk1_etaID;
  float trk1_phiID;

  uint8_t m_trk1Index; 

  uint8_t trk1_nTRT;
  int trk1_truthPdgId;

  float trk1_truthEnergy;
  float trk1_truthP;
  float trk1_truthProb;
  char trk1_hasTruthParticle;
  TLorentzVector trk1_truthPartVec;

  float trk1_pt;
  float trk1_px;
  float trk1_py;
  float trk1_pz;
  float trk1_m;
  float trk1_p;
  float trk1_d0;
  float trk1_z0sintheta;
  float trk1_etaEMB2;
  float trk1_phiEMB2;
  float trk1_etaTileBar2;
  float trk1_phiTileBar2;
  float trk1_etaTileExt1;
  float trk1_phiTileExt1;
  float trk1_etaHEC1;
  float trk1_phiHEC1;
  float trk1_etaEME2;
  float trk1_phiEME2;
  float trk1_nearest_dR_EM;
  float trk1_nearest_dR_HAD;
  float trk1Weight;
  uint8_t trk1_NPV_2;
  uint8_t trk1_NPV_4;
  uint8_t trk1_nclusters;
  uint8_t trk1_nclusters_EM;
  uint8_t trk1_nclusters_hadlike;
  uint8_t trk1_nclusters_HAD;
  uint8_t trk1_nclusters_emlike;
  char trk1_charge;
  char trk1_truthFromPileup;
  char trk1_truthIsFake;
  float trk1_actualmu;
  float trk1_averagemu;
  float trk1_corrected_averagemu;
  float trk1_iso1_EM2;
  float trk1_iso1_HAD2;
  float trk1_iso2_EM2;
  float trk1_iso2_HAD2;

  float m_vertex_x;
  float m_vertex_y;
  float m_vertex_z;
  float m_primary_vertex_x;
  float m_primary_vertex_y;
  float m_primary_vertex_z;


  // TRACK TWO
  float trk2_etaID;
  float trk2_phiID;

  uint8_t m_trk2Index; 

  uint8_t trk2_nTRT;
  int trk2_truthPdgId;

  float trk2_truthEnergy;
  float trk2_truthP;
  float trk2_truthProb;
  char trk2_hasTruthParticle;
  TLorentzVector trk2_truthPartVec;

  float trk2_pt;
  float trk2_px;
  float trk2_py;
  float trk2_pz;
  float trk2_m;
  float trk2_p;
  float trk2_d0;
  float trk2_z0sintheta;
  float trk2_etaEMB2;
  float trk2_phiEMB2;
  float trk2_etaTileBar2;
  float trk2_phiTileBar2;
  float trk2_etaTileExt1;
  float trk2_phiTileExt1;
  float trk2_etaHEC1;
  float trk2_phiHEC1;
  float trk2_etaEME2;
  float trk2_phiEME2;
  float trk2_nearest_dR_EM;
  float trk2_nearest_dR_HAD;
  float trk2Weight;
  uint8_t trk2_NPV_2;
  uint8_t trk2_NPV_4;
  uint8_t trk2_nclusters;
  uint8_t trk2_nclusters_EM;
  uint8_t trk2_nclusters_hadlike;
  uint8_t trk2_nclusters_HAD;
  uint8_t trk2_nclusters_emlike;
  char trk2_charge;
  char trk2_truthFromPileup;
  char trk2_truthIsFake;
  float trk2_actualmu;
  float trk2_averagemu;
  float trk2_corrected_averagemu;
  float trk2_iso1_EM2;
  float trk2_iso1_HAD2;
  float trk2_iso2_EM2;
  float trk2_iso2_HAD2;

};

#endif
