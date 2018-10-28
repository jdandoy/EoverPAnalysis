// E/p analysis for run 2
// Joakim Olsson (joakim.olsson@cern.ch)

#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <EventLoop/OutputStream.h>
#include "AthContainers/ConstDataVector.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODTruth/TruthParticle.h"
#include "xAODAnaHelpers/HelperFunctions.h"
#include <EoverPAnalysis/TrackExtrapolationIsolationTool.h>
#include "AsgTools/MessageCheck.h"
#include "TMath.h"

#include "xAODTracking/TrackParticleContainer.h"
#include "AthLinks/ElementLink.h"
#include "AthContainers/AuxElement.h"

// this is needed to distribute the algorithm to the workers
ClassImp(TrackExtrapolationIsolationTool)

TrackExtrapolationIsolationTool :: TrackExtrapolationIsolationTool (std::string className) :
    Algorithm(className),
    // cutflows
    m_trk_cutflowHist_1(nullptr)
{
  m_inputTrackContainer    = "";
  m_outputTrackContainer   = "";
}

EL::StatusCode TrackExtrapolationIsolationTool :: setupJob (EL::Job& job)
{
  ANA_CHECK_SET_TYPE (EL::StatusCode);
  job.useXAOD();
  xAOD::Init("TrackExtrapolationIsolationTool").ignore();

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode TrackExtrapolationIsolationTool :: histInitialize ()
{

  Info("histInitialize()", "%s", m_name.c_str() );
  ANA_CHECK(xAH::Algorithm::algInitialize());
  // needed here and not in initalize since this is called first
  if( m_inputTrackContainer.empty()){
    Error("configure()", "One or more required configuration values are empty");
    return EL::StatusCode::FAILURE;
  }

  if(m_useCutFlow) {
    TFile *file = wk()->getOutputFile ("cutflow");
    m_trk_cutflowHist_1  = (TH1D*)file->Get("cutflow_trks_1");
  }

  return EL::StatusCode::SUCCESS;
}


EL::StatusCode TrackExtrapolationIsolationTool :: fileExecute () { return EL::StatusCode::SUCCESS; }
EL::StatusCode TrackExtrapolationIsolationTool :: changeInput (bool /*firstFile*/) { return EL::StatusCode::SUCCESS; }

EL::StatusCode TrackExtrapolationIsolationTool :: initialize ()
{
  Info("initialize()", "TrackExtrapolationIsolationTool");
  m_event = wk()->xaodEvent();
  m_store = wk()->xaodStore();

  if(m_useCutFlow) {
    // retrieve the object cutflow
    m_trk_cutflowHist_1  = (TH1D*)file->Get("cutflow_trks_1");
  }
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode TrackExtrapolationIsolationTool :: execute ()
{
  // retrieve event

  // retrieve the tracking cutflows for the plots
 
  const xAOD::TrackParticleContainer* inputTracks(nullptr);
  if (!m_event->retrieve(inputTracks, m_inputTrackContainer).isSuccess()){
     Error("execute()","couldn't retrieve the input track container");
  }
  static SG::AuxElement::Decorator< float > nearestDRDecorator ("dRToNearestTrack");

  ConstDataVector<xAOD::TrackParticleContainer>* selectedTracks(nullptr);
  selectedTracks = new ConstDataVector<xAOD::TrackParticleContainer>(SG::VIEW_ELEMENTS);
  // loop over all tracks only once
  int trk1_counter = 0;
  ANA_MSG_DEBUG("Beginning track loop");
  for(auto trk1 : *inputTracks) {
    trk1_counter += 1; //keep track of which track this is. 
    if(m_useCutFlow) m_trk_cutflowHist_1->Fill("Before Isolation Requirement", 1);

    // coordinates of the track extrapolated to the calorimeter
    // EMB2
    trk1_etaEMB2 = trk1->auxdata<float>("CALO_trkEta_EMB2");
    trk1_phiEMB2 = trk1->auxdata<float>("CALO_trkPhi_EMB2");

    // EME2
    trk1_etaEME2 = trk1->auxdata<float>("CALO_trkEta_EME2");
    trk1_phiEME2 = trk1->auxdata<float>("CALO_trkPhi_EME2");

    //is the track extrapolated to either the EME2 or EMB2?
    bool trk1_hasEME2 = (fabs(trk1_etaEME2) < (float)1000.0) && (fabs(trk1_phiEME2) < (float)1000.0);
    bool trk1_hasEMB2 = (fabs(trk1_etaEMB2) < (float)1000.0) && (fabs(trk1_phiEMB2) < (float)1000.0);

    //continue because the track doesn't have an extrapolation
    if ((!trk1_hasEME2) && (!trk1_hasEMB2)) continue;

    trk1_nearest_dR = 999999.0; //how close was the nearest track?

    //These values will be set to true if the track is not isolated
    bool trk1_not_isolated_EMB2 = false;
    bool trk1_not_isolated_EME2 = false;

    int trk2_counter = 0; //keep track of which track it is
    for (auto trk2 : *inputTracks ) {
      trk2_counter += 1;
      if (trk1_counter == trk2_counter)continue; //Don't compare a track to itself. continue the for loop. 

      //coordinates of the second track extrapolated to the calorimeter.
      //EMB
      float trk2_etaEMB2 = trk2->auxdata<float>("CALO_trkEta_EMB2");
      float trk2_phiEMB2 = trk2->auxdata<float>("CALO_trkPhi_EMB2");

      //EME
      float trk2_etaEME2 = trk2->auxdata<float>("CALO_trkEta_EME2");
      float trk2_phiEME2 = trk2->auxdata<float>("CALO_trkPhi_EME2");

      bool trk2_hasEME2 = (fabs(trk2_etaEME2) < (float)1000.0) && (fabs(trk2_phiEME2) < (float)1000.0);
      bool trk2_hasEMB2 = (fabs(trk2_etaEMB2) < (float)1000.0) && (fabs(trk2_phiEMB2) < (float)1000.0);

      //Do the tracks have an extrapolation to EMB?
      if (trk1_hasEMB2 && trk2_hasEMB2) {
        //the distance between track 1 and track 2 in the EMB
        float trk1_trk2_dR_EMB2 = deltaR(trk1_etaEMB2, trk1_phiEMB2, trk2_etaEMB2, trk2_phiEMB2);
        // was this the nearest track? is the track isolated?
        if (trk1_trk2_dR_EMB2 < trk1_nearest_dR) trk1_nearest_dR = trk1_trk2_dR_EMB2;
        if (trk1_trk2_dR_EMB2 <= m_trkIsoDRmax) trk1_not_isolated_EMB2 = true;
      } //tracks have extrapolation to EMB

      //Does track 1 have an EMB extrapolation, and track 2 an EME?
      //This is to make sure that the tracks are isolated in the calorimeter cracks between the EME and EMB
      if (trk1_hasEMB2 && trk2_hasEME2){
        float trk1_EMB_trk2_EME_dR = deltaR(trk1_etaEMB2, trk1_phiEMB2, trk2_etaEME2, trk2_phiEME2);
        // was this the nearest track? Is the track isolated?
        if (trk1_EMB_trk2_EME_dR < trk1_nearest_dR) trk1_nearest_dR = trk1_EMB_trk2_EME_dR;
        if (trk1_EMB_trk2_EME_dR < m_trkIsoDRmax) trk1_not_isolated_EMB2 = true; //the track is in the EMB, and it is not isolated
      }

      //does track 1 have an EME extrapolation, and track 2 an EMB?
      //This is to make sure that the tracks are isolated in the calorimeter cracks between the EME and EMB
      if (trk1_hasEME2 && trk2_hasEMB2){
        float trk1_EME_trk2_EMB_dR = deltaR(trk1_etaEME2, trk1_phiEME2, trk2_etaEMB2, trk2_phiEMB2);
        // was this the nearest track? Is the track isolated?
        if (trk1_EME_trk2_EMB_dR < trk1_nearest_dR) trk1_nearest_dR = trk1_EME_trk2_EMB_dR;
        if (trk1_EME_trk2_EMB_dR < m_trkIsoDRmax) trk1_not_isolated_EME2 = true;
      }

      //Do the tracks have an extrapoltion to the EME?
      if (trk1_hasEME2 && trk2_hasEME2) {
        //the distance between track 2 and track 2 in the EME
        float trk1_trk2_dR_EME2 = deltaR(trk1_etaEME2, trk1_phiEME2, trk2_etaEME2, trk2_phiEME2);
        // was this the nearest track? Is the track isolated?
        if (trk1_trk2_dR_EME2 < trk1_nearest_dR) trk1_nearest_dR = trk1_trk2_dR_EME2;
        if (trk1_trk2_dR_EME2 <= m_trkIsoDRmax) trk1_not_isolated_EME2 = true;
      } //tracks have extrapolation to EME
    } // END looping trk2

    if (trk1_not_isolated_EMB2) {ANA_MSG_DEBUG("Track failed isolation"); continue;}
    if (trk1_not_isolated_EME2) {ANA_MSG_DEBUG("Track failed isolation");  continue;}
    ANA_MSG_DEBUG("Track passed isolation cut, decorating with dR = " + std::to_string(trk1_nearest_dR));
    nearestDRDecorator(*trk1) = trk1_nearest_dR;//decorate the track with the distance to the nearest track
    selectedTracks->push_back(trk1);
    if(m_useCutFlow) m_trk_cutflowHist_1->Fill("Isolated in EM Layer 2", 1);
  } //END looping trk1
  m_store->record( selectedTracks, m_outputTrackContainer );
  return EL::StatusCode::SUCCESS;
}//Close execute function

EL::StatusCode TrackExtrapolationIsolationTool :: postExecute () { return EL::StatusCode::SUCCESS; }

EL::StatusCode TrackExtrapolationIsolationTool :: finalize () { 
    return EL::StatusCode::SUCCESS;
}

EL::StatusCode TrackExtrapolationIsolationTool :: histFinalize ()
{
    return EL::StatusCode::SUCCESS;
}

double TrackExtrapolationIsolationTool :: deltaR (double trk_eta, double trk_phi, double trk2_eta, double trk2_phi)
{
  double trk_trk2_dEta = fabs(trk2_eta - trk_eta);
  double trk_trk2_dPhi = fabs(trk2_phi - trk_phi);

  if (trk_trk2_dPhi > TMath::Pi())
    trk_trk2_dPhi = 2*TMath::Pi() - trk_trk2_dPhi;

  return sqrt( pow(trk_trk2_dEta, 2) + pow(trk_trk2_dPhi, 2) );
}

std::vector<double> TrackExtrapolationIsolationTool::str2vec(std::string str)
{
  std::vector<double> vec;
  if (str.size() > 0) {
    str.erase (std::remove (str.begin(), str.end(), ' '), str.end()); // remove whitespace
    std::stringstream ss(str);
    std::string token; // split string at ','
    while ( std::getline(ss, token, ','))
      vec.push_back(std::stof(token));
  }
  return vec;
}
