#include "EoverPAnalysis/EoverPHistsTrks.h"
#include "xAODAnaHelpers/HelperFunctions.h"
#include <regex>
#include <math.h>

EoverPHistsTrks :: EoverPHistsTrks (std::string name, std::string detailStr, float trkIsoDRmax, float trkIsoPfrac, bool doTrkPcut, float trkPmin, float trkPmax, bool doTrkEtacut, float trkEtamin, float trkEtamax, bool doTrkIsocut) :
  HistogramManager(name, detailStr)
{
  m_trkIsoDRmax = trkIsoDRmax; 
  m_trkIsoPfrac = trkIsoPfrac; 
  m_doTrkPcut = doTrkPcut;
  m_trkPmin = trkPmin;
  m_trkPmax = trkPmax; 
  m_doTrkEtacut = doTrkEtacut;
  m_trkEtamin = trkEtamin;
  m_trkEtamax = trkEtamax; 
  m_doTrkIsocut = doTrkIsocut;
}

EoverPHistsTrks :: ~EoverPHistsTrks () {}

StatusCode EoverPHistsTrks::initialize()
{

  // number of bins and ranges for histograms
  unsigned int nBinsMu = 50;           float minMu = 0.0;           float maxMu = 50.0;
  unsigned int nBinsMu_many = 500;
  float minMu_shift = -0.5; float maxMu_shift = 49.5;
  unsigned int nBinsNPV = 50;          float minNPV = -0.5;         float maxNPV = 49.5;
  unsigned int nBinsTrkN = 200;        float minTrkN = -0.5;        float maxTrkN = 199.5;
  unsigned int nBinsP = 500;           float minP = 0;              float maxP = 50;
  unsigned int nBinsDR = 60;           float minDR = 0;             float maxDR = 3;
  unsigned int nBinsPhi = 32;          float minPhi = -TMath::Pi(); float maxPhi = TMath::Pi(); 
  unsigned int nBinsEta = 100;         float minEta = -2.5;         float maxEta = 2.5;

  //// Book histograms

  // event level plots
  m_mu = book(m_name, "mu", "#mu", nBinsMu, minMu, maxMu); 
  m_mu_avg = book(m_name, "mu_avg", "<#mu>", nBinsMu, minMu, maxMu); 
  m_mu_avg_many = book(m_name, "mu_avg_many", "<#mu>", nBinsMu_many, minMu, maxMu); 
  m_mu_avg_shift = book(m_name, "mu_avg_shift", "<#mu>", nBinsMu, minMu_shift, maxMu_shift); 
  m_mu_avg_vs_npv = book(m_name, "mu_avg_vs_npv", "NPV", nBinsNPV, minNPV, maxNPV, "<#mu>", nBinsMu, minMu, maxMu); 
  m_mu_avg_vs_trk_n_nocut = book(m_name, "mu_avg_vs_trk_n_nocut", "N_{trk}", nBinsTrkN, minTrkN, maxTrkN, "<#mu>", nBinsMu, minMu, maxMu); 
  m_npv = book(m_name, "npv", "NPV", nBinsNPV, minNPV, maxNPV); 
  m_npv_trks = book(m_name, "npv_trks", "N_{trks} associated with a PV", nBinsTrkN, minTrkN, maxTrkN); 

  // track plots
  m_trk_n_nocut       = book(m_name, "trk_n_nocut", "N_{trk}", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_p_noiso       = book(m_name, "trk_p_noiso", "p_{trk} [GeV]", nBinsP, minP, maxP); 
  m_trk_etaID_noiso   = book(m_name, "trk_etaID_noiso", "#eta_{trk}", nBinsEta, minEta, maxEta); 
  m_trk_etaEMB2_noiso = book(m_name, "trk_etaEMB2_noiso", "#eta_{trk}", nBinsEta, minEta, maxEta); 
  m_trk_etaEME2_noiso = book(m_name, "trk_etaEME2_noiso", "#eta_{trk}", nBinsEta, minEta, maxEta); 
  m_trk_phiID_noiso   = book(m_name, "trk_phiID_noiso", "#phi_{trk}", nBinsPhi, minPhi, maxPhi); 
  m_trk_phiEMB2_noiso = book(m_name, "trk_phiEMB2_noiso", "#phi_{trk}", nBinsPhi, minPhi, maxPhi); 
  m_trk_phiEME2_noiso = book(m_name, "trk_phiEME2_noiso", "#phi_{trk}", nBinsPhi, minPhi, maxPhi); 

  m_trk_etaEMB2_vs_etaID_noiso = book(m_name, "trk_etaEMB2_vs_etaID_noiso", "#eta_{trk,ID}", nBinsEta, minEta, maxEta, "#eta_{trk,EMB2}", nBinsEta, minEta, maxEta); 
  m_trk_etaEME2_vs_etaID_noiso = book(m_name, "trk_etaEME2_vs_etaID_noiso", "#eta_{trk,ID}", nBinsEta, minEta, maxEta, "#eta_{trk,EME2}", nBinsEta, minEta, maxEta); 
  m_trk_etaEME2_vs_etaEMB2_noiso = book(m_name, "trk_etaEME2_vs_etaEMB2_noiso", "#eta_{trk,EMB2}", nBinsEta, minEta, maxEta, "#eta_{trk,EME2}", nBinsEta, minEta, maxEta); 

  m_trk_trk2_dR_ID = book(m_name, "trk_trk2_dR_ID", "#DeltaR(trk,trk2)", nBinsDR, minDR, maxDR);
  m_trk_trk2_dR_ID_vs_trk_p = book(m_name, "trk_trk2_dR_ID_vs_trk_p", "p_{trk} [GeV]", nBinsP, minP, maxP, "#DeltaR(trk,trk2)", nBinsDR, minDR, maxDR);
  m_trk_trk2_dR_EMB2 = book(m_name, "trk_trk2_dR_EMB2", "#DeltaR(trk,trk2)", nBinsDR, minDR, maxDR);
  m_trk_trk2_dR_EMB2_vs_trk_p = book(m_name, "trk_trk2_dR_EMB2_vs_trk_p", "p_{trk} [GeV]", nBinsP, minP, maxP, "#DeltaR(trk,trk2)", nBinsDR, minDR, maxDR);
  m_trk_trk2_dR_EME2 = book(m_name, "trk_trk2_dR_EME2", "#DeltaR(trk,trk2)", nBinsDR, minDR, maxDR);
  m_trk_trk2_dR_EME2_vs_trk_p = book(m_name, "trk_trk2_dR_EME2_vs_trk_p", "p_{trk} [GeV]", nBinsP, minP, maxP, "#DeltaR(trk,trk2)", nBinsDR, minDR, maxDR);
  m_trk_trk2_dR_EMB2_EME2 = book(m_name, "trk_trk2_dR_EMB2_EME2", "#DeltaR(trk_EMB2,trk2_EME2)", nBinsDR, minDR, maxDR);
  m_trk_trk2_dR_EMB2_EME2_vs_trk_p = book(m_name, "trk_trk2_dR_EMB2_EME2_vs_trk_p", "p_{trk} [GeV]", nBinsP, minP, maxP, "#DeltaR(trk_EMB2,trk2_EME2)", nBinsDR, minDR, maxDR);
  m_trk_trk2_dR_EME2_EMB2 = book(m_name, "trk_trk2_dR_EME2_EMB2", "#DeltaR(trk_EME2,trk2_EMB2)", nBinsDR, minDR, maxDR);
  m_trk_trk2_dR_EME2_EMB2_vs_trk_p = book(m_name, "trk_trk2_dR_EME2_EMB2_vs_trk_p", "p_{trk} [GeV]", nBinsP, minP, maxP, "#DeltaR(trk_EME2,trk2_EMB2)", nBinsDR, minDR, maxDR);

  m_trk_ntrks_maxDR01 = book(m_name, "trk_ntrks_maxDR01", "N_{trks} within #Delta R<0.1 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR02 = book(m_name, "trk_ntrks_maxDR02", "N_{trks} within #Delta R<0.2 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR03 = book(m_name, "trk_ntrks_maxDR03", "N_{trks} within #Delta R<0.3 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR04 = book(m_name, "trk_ntrks_maxDR04", "N_{trks} within #Delta R<0.4 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR05 = book(m_name, "trk_ntrks_maxDR05", "N_{trks} within #Delta R<0.5 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR06 = book(m_name, "trk_ntrks_maxDR06", "N_{trks} within #Delta R<0.6 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR07 = book(m_name, "trk_ntrks_maxDR07", "N_{trks} within #Delta R<0.7 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR08 = book(m_name, "trk_ntrks_maxDR08", "N_{trks} within #Delta R<0.8 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR09 = book(m_name, "trk_ntrks_maxDR09", "N_{trks} within #Delta R<0.9 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR10 = book(m_name, "trk_ntrks_maxDR10", "N_{trks} within #Delta R<1.0 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 

  m_trk_ntrks_maxDR01_vs_trk_p = book(m_name, "trk_ntrks_maxDR01_vs_trk_p", "p_{trk} [GeV]", nBinsP, minP, maxP, "N_{trks} within #Delta R<0.1 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR02_vs_trk_p = book(m_name, "trk_ntrks_maxDR02_vs_trk_p", "p_{trk} [GeV]", nBinsP, minP, maxP, "N_{trks} within #Delta R<0.2 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR03_vs_trk_p = book(m_name, "trk_ntrks_maxDR03_vs_trk_p", "p_{trk} [GeV]", nBinsP, minP, maxP, "N_{trks} within #Delta R<0.3 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR04_vs_trk_p = book(m_name, "trk_ntrks_maxDR04_vs_trk_p", "p_{trk} [GeV]", nBinsP, minP, maxP, "N_{trks} within #Delta R<0.4 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR05_vs_trk_p = book(m_name, "trk_ntrks_maxDR05_vs_trk_p", "p_{trk} [GeV]", nBinsP, minP, maxP, "N_{trks} within #Delta R<0.5 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR06_vs_trk_p = book(m_name, "trk_ntrks_maxDR06_vs_trk_p", "p_{trk} [GeV]", nBinsP, minP, maxP, "N_{trks} within #Delta R<0.6 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR07_vs_trk_p = book(m_name, "trk_ntrks_maxDR07_vs_trk_p", "p_{trk} [GeV]", nBinsP, minP, maxP, "N_{trks} within #Delta R<0.7 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR08_vs_trk_p = book(m_name, "trk_ntrks_maxDR08_vs_trk_p", "p_{trk} [GeV]", nBinsP, minP, maxP, "N_{trks} within #Delta R<0.8 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR09_vs_trk_p = book(m_name, "trk_ntrks_maxDR09_vs_trk_p", "p_{trk} [GeV]", nBinsP, minP, maxP, "N_{trks} within #Delta R<0.9 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_maxDR10_vs_trk_p = book(m_name, "trk_ntrks_maxDR10_vs_trk_p", "p_{trk} [GeV]", nBinsP, minP, maxP, "N_{trks} within #Delta R<1.0 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 

  m_trk_ntrks_trkIsoDRmax_vs_sum_surr_trk_p = book(m_name, "trk_ntrks_trkIsoDRmax_vs_sum_surr_trk_p", "#Sigma_{i#neq0} p_{trk,#DeltaR<0.4}_{i} [GeV]", nBinsP, minP, maxP, "N_{trks} within #Delta R<0.4 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_trkIsoDRmax_vs_avg_surr_trk_p = book(m_name, "trk_ntrks_trkIsoDRmax_vs_avg_surr_trk_p", "p_{trk,#DeltaR<0.4}^{avg} [GeV]", nBinsP, minP, maxP, "N_{trks} within #Delta R<0.4 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 
  m_trk_ntrks_trkIsoDRmax_vs_leading_surr_trk_p = book(m_name, "trk_ntrks_trkIsoDRmax_vs_leading_surr_trk_p", "p_{trk,#DeltaR<0.4}^{leading} [GeV]", nBinsP, minP, maxP, "N_{trks} within #Delta R<0.4 of selected trk", nBinsTrkN, minTrkN, maxTrkN); 

  m_trk_ntrks_trkIsoDRmax_sum_surr_trk_p_vs_trk_p = book(m_name, "trk_ntrks_trkIsoDRmax_sum_surr_trk_p_vs_trk_p", "p_{trk} [GeV]", nBinsP, minP, maxP, "#Sigma_{i#neq0} p_{trk,#DeltaR<0.4}_{i} [GeV]", nBinsP, minP, maxP); 
  m_trk_ntrks_trkIsoDRmax_avg_surr_trk_p_vs_trk_p = book(m_name, "trk_ntrks_trkIsoDRmax_avg_surr_trk_p_vs_trk_p", "p_{trk} [GeV]", nBinsP, minP, maxP, "p_{trk,#DeltaR<0.4}^{avg} [GeV]", nBinsP, minP, maxP); 
  m_trk_ntrks_trkIsoDRmax_leading_surr_trk_p_vs_trk_p = book(m_name, "trk_ntrks_trkIsoDRmax_leading_surr_trk_p_vs_trk_p", "p_{trk} [GeV]", nBinsP, minP, maxP, "p_{trk,#DeltaR<0.4}^{leading} [GeV]", nBinsP, minP, maxP); 

  m_trk_ntrks_trkIsoDRmax_sum_surr_trk_p_over_trk_p = book(m_name, "trk_ntrks_trkIsoDRmax_sum_surr_trk_p_over_trk_p", "#Sigma_{i#neq0} p_{trk,#DeltaR<0.4}_{i} / p_{trk}", 100, 0, 5.0); 
  m_trk_ntrks_trkIsoDRmax_avg_surr_trk_p_over_trk_p = book(m_name, "trk_ntrks_trkIsoDRmax_avg_surr_trk_p_over_trk_p", "p_{trk,#DeltaR<0.4}^{avg} / p_{trk}", 100, 0, 5.0); 
  m_trk_ntrks_trkIsoDRmax_leading_surr_trk_p_over_trk_p = book(m_name, "trk_ntrks_trkIsoDRmax_leading_surr_trk_p_over_trk_p", "p_{trk,#DeltaR<0.4}^{leading} / p_{trk}", 100, 0, 5.0); 


  // if worker is passed to the class add histograms to the output
  return StatusCode::SUCCESS;
}

StatusCode EoverPHistsTrks::execute( const xAOD::TrackParticleContainer* trks, const xAOD::VertexContainer* vtxs, const xAOD::EventInfo* eventInfo, float eventWeight )
{

  // get pileup
  float mu(-1.);
  if( eventInfo->isAvailable< float >( "actualInteractionsPerCrossing" ) ) {
    mu = eventInfo->actualInteractionsPerCrossing();
  }
  float mu_avg(-1.);
  if( eventInfo->isAvailable< float >( "corrected_averageInteractionsPerCrossing" ) &&
      !eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION )) 
    mu_avg = eventInfo->auxdata< float >( "corrected_averageInteractionsPerCrossing" );
  else if( eventInfo->isAvailable< float >( "averageInteractionsPerCrossing" ) )
    mu_avg = eventInfo->averageInteractionsPerCrossing();

  m_mu -> Fill(mu, eventWeight);
  m_mu_avg -> Fill(mu_avg, eventWeight);
  m_mu_avg_many -> Fill(mu_avg, eventWeight);
  m_mu_avg_shift -> Fill(mu_avg, eventWeight);

  // get number of primary vtxs 
  float npv = HelperFunctions::countPrimaryVertices(vtxs, 2);
  m_npv -> Fill(npv, eventWeight);
  for( auto vtx_itr : *vtxs)
    m_npv_trks->Fill(vtx_itr->nTrackParticles());
  m_mu_avg_vs_npv -> Fill(npv, mu_avg, eventWeight);

  // number of tracks
  // m_trk_n_nocut-> Fill(trks->size());

  // track-cluster plots

  // loop over tracks and clusters around a selected track, plot how many tracks (clusters) fall within a certain dR
  float trk_n = 0.;
  float maxDR_ranges[] = {.1, .2, .3, .4, .5, .6, .7, .8, .9, 1.0};

  xAOD::TrackParticleContainer::const_iterator trk_itr = trks->begin();
  xAOD::TrackParticleContainer::const_iterator trk_end = trks->end();
  xAOD::TrackParticleContainer::const_iterator trk2_itr = trks->begin();
  xAOD::TrackParticleContainer::const_iterator trk2_end = trks->end();

  // track-track and track-cluster (outer loop over all tracks)
  for( ; trk_itr != trk_end; ++trk_itr ) {

    // count the number of tracks passing all selections
    trk_n += 1;
    const xAOD::TrackParticle* trk = (*trk_itr);
    float trk_p = 0;
    if (TMath::Abs(trk->qOverP())>0.) trk_p = (1./TMath::Abs(trk->qOverP()))/1e3; 
    m_trk_p_noiso->Fill(trk_p, eventWeight);

    float trk_etaID = trk->eta();
    float trk_phiID = trk->phi();
    // coordinates of the track extrapolated to the calorimeter
    // EMB2
    float trk_etaEMB2 = trk->auxdata<float>("CALO_trkEta_EMB2");
    float trk_phiEMB2 = trk->auxdata<float>("CALO_trkPhi_EMB2");
    // EME2
    float trk_etaEME2 = trk->auxdata<float>("CALO_trkEta_EME2");
    float trk_phiEME2 = trk->auxdata<float>("CALO_trkPhi_EME2");

    // check that the track is extrapolated to either EMB2 or EME2
    // (if not then trk_eta = trk_phi = -999999999)
    if (TMath::Abs(trk_etaEMB2) > 100.0 && TMath::Abs(trk_etaEME2) > 100.0) { 
      continue;
    }

    m_trk_etaID_noiso->Fill(trk_etaID, eventWeight);
    m_trk_etaEMB2_noiso->Fill(trk_etaEMB2, eventWeight);
    m_trk_etaEME2_noiso->Fill(trk_etaEME2, eventWeight);
    m_trk_phiID_noiso->Fill(trk_phiID, eventWeight);
    m_trk_phiEMB2_noiso->Fill(trk_phiEMB2, eventWeight);
    m_trk_phiEME2_noiso->Fill(trk_phiEME2, eventWeight);

    m_trk_etaEMB2_vs_etaID_noiso->Fill(trk_etaID, trk_etaEMB2, eventWeight);
    m_trk_etaEME2_vs_etaID_noiso->Fill(trk_etaID, trk_etaEME2, eventWeight);
    m_trk_etaEME2_vs_etaEMB2_noiso->Fill(trk_etaEMB2, trk_etaEME2, eventWeight);

    // check track p requirement
    if (m_doTrkPcut) {
      if (trk_p < m_trkPmin) continue;
      if (trk_p > m_trkPmax) continue;
    }

    // check track eta requirement
    if (m_doTrkEtacut) {
      if (TMath::Abs(trk_etaID) < m_trkEtamin) continue;
      if (TMath::Abs(trk_etaID) > m_trkEtamax) continue;
    }

    // keep track of the number of tracks within a certain dR from the selected track
    int trk_ntrks_maxDR[10] = {0};
    float trk2_p = 0.;
    float surr_trk_sum_p = 0.;
    float surr_trk_leading_p = 0.;
    trk2_itr = trks->begin();
    for( ; trk2_itr != trk2_end; ++trk2_itr ) {
      if (trk_itr != trk2_itr) { // do not double count the selected track 
        const xAOD::TrackParticle* trk2 = (*trk2_itr);

        float trk2_etaID  = trk2->eta();
        float trk2_phiID  = trk2->phi();
        float trk_trk2_dR_ID = deltaR(trk_etaID, trk_phiID, trk2_etaID, trk2_phiID);
        m_trk_trk2_dR_ID -> Fill(trk_trk2_dR_ID, eventWeight);
        m_trk_trk2_dR_ID_vs_trk_p -> Fill(trk_p, trk_trk2_dR_ID, eventWeight);

        //EMB2
        float trk2_etaEMB2 = trk2->auxdata<float>("CALO_trkEta_EMB2");
        float trk2_phiEMB2 = trk2->auxdata<float>("CALO_trkPhi_EMB2");
        //EME2
        float trk2_etaEME2 = trk2->auxdata<float>("CALO_trkEta_EME2");
        float trk2_phiEME2 = trk2->auxdata<float>("CALO_trkPhi_EME2");

        // Calculate all possible permutations of dR between tracks extrapolated to EMB2 or EME2
        // Modification: skip the case when one track is on EMB2 and the other in EME2
        float trk_trk2_dR[4] = {1e8, 1e8, 1e8, 1e8};
        trk_trk2_dR[0] = deltaR(trk_etaEMB2, trk_phiEMB2, trk2_etaEMB2, trk2_phiEMB2);
        // trk_trk2_dR[1] = deltaR(trk_etaEMB2, trk_phiEMB2, trk2_etaEME2, trk2_phiEME2);
        // trk_trk2_dR[2] = deltaR(trk_etaEME2, trk_phiEME2, trk2_etaEMB2, trk2_phiEMB2);
        trk_trk2_dR[3] = deltaR(trk_etaEME2, trk_phiEME2, trk2_etaEME2, trk2_phiEME2);

        float trk_trk2_dR_min = trk_trk2_dR[0];
        for (int i = 1; i < 4; ++i) {
          if (trk_trk2_dR[i] < trk_trk2_dR_min)
            trk_trk2_dR_min = trk_trk2_dR[i];
        }

        m_trk_trk2_dR_EMB2 -> Fill(trk_trk2_dR[0], eventWeight);
        m_trk_trk2_dR_EMB2_vs_trk_p -> Fill(trk_p, trk_trk2_dR[0], eventWeight);

        m_trk_trk2_dR_EME2 -> Fill(trk_trk2_dR[3], eventWeight);
        m_trk_trk2_dR_EME2_vs_trk_p -> Fill(trk_p, trk_trk2_dR[3], eventWeight);

        m_trk_trk2_dR_EMB2_EME2 -> Fill(trk_trk2_dR[1], eventWeight);
        m_trk_trk2_dR_EMB2_EME2_vs_trk_p -> Fill(trk_p, trk_trk2_dR[1], eventWeight);

        m_trk_trk2_dR_EME2_EMB2 -> Fill(trk_trk2_dR[2], eventWeight);
        m_trk_trk2_dR_EME2_EMB2_vs_trk_p -> Fill(trk_p, trk_trk2_dR[2], eventWeight);

        // count the number of track wihin a given radius
        for (int i = 0; i < 10; ++i) {
          if (trk_trk2_dR_ID < maxDR_ranges[i]) {
            trk_ntrks_maxDR[i]++;
          }
        }

        // track isolation
        if (trk_trk2_dR_min < m_trkIsoDRmax) { // check if trk2 falls within DRmax of trk 

          // calculate the leading and avg p of the surrounding tracks
          if (TMath::Abs(trk2->qOverP())>0.) trk2_p = (1./TMath::Abs(trk2->qOverP()))/1e3; 
          surr_trk_sum_p += trk2_p;
          if (trk2_p > surr_trk_leading_p) {
            surr_trk_leading_p = trk2_p;
          }
        }
      }
    } // END looping trk2

    // plots of the number of close-by tracks within a certain radius
    m_trk_ntrks_maxDR01 -> Fill(trk_ntrks_maxDR[0], eventWeight);
    m_trk_ntrks_maxDR02 -> Fill(trk_ntrks_maxDR[1], eventWeight);
    m_trk_ntrks_maxDR03 -> Fill(trk_ntrks_maxDR[2], eventWeight);
    m_trk_ntrks_maxDR04 -> Fill(trk_ntrks_maxDR[3], eventWeight);
    m_trk_ntrks_maxDR05 -> Fill(trk_ntrks_maxDR[4], eventWeight);
    m_trk_ntrks_maxDR06 -> Fill(trk_ntrks_maxDR[5], eventWeight);
    m_trk_ntrks_maxDR07 -> Fill(trk_ntrks_maxDR[6], eventWeight);
    m_trk_ntrks_maxDR08 -> Fill(trk_ntrks_maxDR[7], eventWeight);
    m_trk_ntrks_maxDR09 -> Fill(trk_ntrks_maxDR[8], eventWeight);
    m_trk_ntrks_maxDR10 -> Fill(trk_ntrks_maxDR[9], eventWeight);

    // plots of the number of close-by tracks vs p of the selected track
    m_trk_ntrks_maxDR01_vs_trk_p -> Fill(trk_p, trk_ntrks_maxDR[0], eventWeight);
    m_trk_ntrks_maxDR02_vs_trk_p -> Fill(trk_p, trk_ntrks_maxDR[1], eventWeight);
    m_trk_ntrks_maxDR03_vs_trk_p -> Fill(trk_p, trk_ntrks_maxDR[2], eventWeight);
    m_trk_ntrks_maxDR04_vs_trk_p -> Fill(trk_p, trk_ntrks_maxDR[3], eventWeight);
    m_trk_ntrks_maxDR05_vs_trk_p -> Fill(trk_p, trk_ntrks_maxDR[4], eventWeight);
    m_trk_ntrks_maxDR06_vs_trk_p -> Fill(trk_p, trk_ntrks_maxDR[5], eventWeight);
    m_trk_ntrks_maxDR07_vs_trk_p -> Fill(trk_p, trk_ntrks_maxDR[6], eventWeight);
    m_trk_ntrks_maxDR08_vs_trk_p -> Fill(trk_p, trk_ntrks_maxDR[7], eventWeight);
    m_trk_ntrks_maxDR09_vs_trk_p -> Fill(trk_p, trk_ntrks_maxDR[8], eventWeight);
    m_trk_ntrks_maxDR10_vs_trk_p -> Fill(trk_p, trk_ntrks_maxDR[9], eventWeight);

    // plots of leading p and avg p of the tracks surrounding the selected track 
    float surr_trk_avg_p = surr_trk_sum_p/(trks->size()-1); // calculate avg p of the surrounding tracks 
    m_trk_ntrks_trkIsoDRmax_vs_sum_surr_trk_p -> Fill(surr_trk_sum_p, trk_ntrks_maxDR[3], eventWeight);
    m_trk_ntrks_trkIsoDRmax_vs_avg_surr_trk_p -> Fill(surr_trk_avg_p, trk_ntrks_maxDR[3], eventWeight);
    m_trk_ntrks_trkIsoDRmax_vs_leading_surr_trk_p -> Fill(surr_trk_leading_p, trk_ntrks_maxDR[3], eventWeight);
    m_trk_ntrks_trkIsoDRmax_sum_surr_trk_p_vs_trk_p -> Fill(trk_p, surr_trk_sum_p, eventWeight);
    m_trk_ntrks_trkIsoDRmax_avg_surr_trk_p_vs_trk_p -> Fill(trk_p, surr_trk_avg_p, eventWeight);
    m_trk_ntrks_trkIsoDRmax_leading_surr_trk_p_vs_trk_p -> Fill(trk_p, surr_trk_leading_p, eventWeight);
    m_trk_ntrks_trkIsoDRmax_sum_surr_trk_p_over_trk_p -> Fill(surr_trk_sum_p/trk_p, eventWeight);
    m_trk_ntrks_trkIsoDRmax_avg_surr_trk_p_over_trk_p -> Fill(surr_trk_avg_p/trk_p, eventWeight);
    m_trk_ntrks_trkIsoDRmax_leading_surr_trk_p_over_trk_p -> Fill(surr_trk_leading_p/trk_p, eventWeight);

  } // END loop tracks 

  m_trk_n_nocut-> Fill(trk_n);
  m_mu_avg_vs_trk_n_nocut -> Fill(trk_n, mu_avg, eventWeight);

  return StatusCode::SUCCESS;
}

float EoverPHistsTrks::deltaR (float trk_eta, float trk_phi, float trk2_eta, float trk2_phi)
{
  float trk_trk2_dEta = TMath::Abs(trk2_eta - trk_eta);
  float trk_trk2_dPhi = TMath::Abs(trk2_phi - trk_phi);
  if (trk_trk2_dPhi > TMath::Pi())
    trk_trk2_dPhi = 2*TMath::Pi() - trk_trk2_dPhi;
  return sqrt( pow(trk_trk2_dEta, 2) + pow(trk_trk2_dPhi, 2) );
}
