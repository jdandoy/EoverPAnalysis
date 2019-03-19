#ifndef digihstruthgrooming_SecondariesTrees_H
#define digihstruthgrooming_SecondariesTrees_H

#include <iostream>
#include <iomanip>
#include <vector>

#include <AnaAlgorithm/AnaAlgorithm.h>
#include "EventLoop/OutputStream.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"

#include "xAODRootAccess/TStore.h"
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"

class SecondariesTrees : public EL::AnaAlgorithm
{
 public:
  SecondariesTrees (const std::string& name, ISvcLocator* pSvcLocator);
  ~SecondariesTrees();

  virtual StatusCode initialize () override;
  virtual StatusCode execute () override;
  virtual StatusCode finalize () override;

  xAOD::TEvent *m_event; //!

 private:
  // Set by config
  Int_t m_MessageFrequency;
  Bool_t m_isData;
  TString m_label;
  TString m_VertexContainer;

  // Output variables for tree

  Int_t m_DSID;
  Float_t m_mcWeight;
  Float_t m_weight;
  Float_t m_mu;
  Int_t m_npv;

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
  
  Float_t m_track1_pt;
  Float_t m_track2_pt;
  Float_t m_track1_etaID;
  Float_t m_track1_phiID;
  Float_t m_track2_etaID;
  Float_t m_track2_phiID;
  Float_t m_track1_m;
  Float_t m_track2_m;
  
  Int_t m_trackN;
  /*
  std::vector<float> *m_trackPt = nullptr;
  std::vector<float> *m_trackEta = nullptr;
  std::vector<float> *m_trackPhi = nullptr;
  std::vector<float> *m_trackM = nullptr;
  */
};

#endif
