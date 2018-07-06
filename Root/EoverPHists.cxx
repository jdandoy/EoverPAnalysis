// E/p analysis for run 2
// Joakim Olsson (joakim.olsson@cern.ch)

#include "EoverPAnalysis/EoverPHists.h"
#include "xAODAnaHelpers/HelperFunctions.h"
#include <regex>
#include <sstream>
#include <cmath>

EoverPHists :: EoverPHists (std::string name, std::string detailStr, std::string energyCalib, bool doCaloTotal, bool doCaloEM, bool doCaloHAD, bool doBgSubtr, bool doTileLayer, std::string Pbins, bool doPbinsArray, std::string PbinsArray, std::string Etabins, bool doEtabinsArray, std::string EtabinsArray, bool doProfileEta, bool doProfileP, bool doExtraEtaEnergyBinHists) : HistogramManagerD(name, detailStr)
{
  m_energyCalib = energyCalib;
  m_doCaloTotal = doCaloTotal;
  m_doCaloEM = doCaloEM;
  m_doCaloHAD = doCaloHAD;
  m_doBgSubtr = doBgSubtr;
  m_doTileLayer = doTileLayer;
  m_Pbins = Pbins;
  m_doPbinsArray = doPbinsArray;
  m_PbinsArray = PbinsArray; 
  m_Etabins = Etabins; 
  m_doEtabinsArray = doEtabinsArray;
  m_EtabinsArray = EtabinsArray; 
  m_doProfileEta= doProfileEta;
  m_doProfileP= doProfileP;
  m_doExtraEtaEnergyBinHists = doExtraEtaEnergyBinHists;
}

EoverPHists :: ~EoverPHists () {}

StatusCode EoverPHists::initialize()
{
  // number of bins and ranges for histograms
  unsigned int nBinsMu = 50;          double minMu = 0;                  double maxMu = 50;
  unsigned int nBinsNPV = 50;         double minNPV = -0.5;              double maxNPV = 49.5;
  unsigned int nBinsDR = 60;          double minDR = 0;                  double maxDR = 3;
  unsigned int nBinsPhi = 32;         double minPhi = -TMath::Pi();      double maxPhi = TMath::Pi(); 
  unsigned int nBinsPhiExtra = 1024;  double minPhiExtra = -TMath::Pi(); double maxPhiExtra = TMath::Pi(); 
  unsigned int nBinsPhiExtra2 = 800;  double minPhiExtra2 = -4.0;        double maxPhiExtra2 = 4.0; 
  unsigned int nBinsEop = 300;        double minEop = -4;                double maxEop = 20;
  unsigned int nBinsEop_l = 255;      double minEop_l = -100;            double maxEop_l = 5000;
  unsigned int nBinsE = 700;          double minE = -20;                 double maxE = 50;

  unsigned int nBinsP = 500;          double minP = 0;                   double maxP = 50;
  std::vector<double> Pbins = str2vec(m_Pbins);
  if (Pbins.size() > 2) {
    nBinsP = (int) Pbins[0];
    minP = (double) Pbins[1];
    maxP = (double) Pbins[2];
  }

  nPbinsArray = 0;
  if (m_doPbinsArray){
    PbinsArray = str2vec(m_PbinsArray);
    nPbinsArray = PbinsArray.size()-1;
  }
  if (nPbinsArray <= 1)
    m_doPbinsArray = false;

  unsigned int nBinsEta = 100;      double minEta = -2.5;           double maxEta = 2.5;
  std::vector<double> Etabins = str2vec(m_Etabins);
  if (Etabins.size() > 2) {
    nBinsEta = (int) Etabins[0];
    minEta = (double) Etabins[1];
    maxEta = (double) Etabins[2];
  }
  if (minEta < 0) m_doEtaAbs = false;
  else m_doEtaAbs = true;

  nEtabinsArray = 0;
  if (m_doEtabinsArray){
    EtabinsArray = str2vec(m_EtabinsArray);
    nEtabinsArray = EtabinsArray.size()-1;
  }
  if (nEtabinsArray <= 1)
    m_doEtabinsArray = false;

  // track kinematics
  m_trk_p = book(m_name, "trk_p", "p_{trk} [GeV]", nBinsP, minP, maxP); 
  if (m_doPbinsArray) m_trk_p_array = book(m_name, "trk_p_PbinsArray", "p_{trk} [GeV]", nPbinsArray, &PbinsArray[0]); 
  m_trk_pt = book(m_name, "trk_pt", "p_{T,trk} [GeV]", nBinsP, minP, maxP); 
  if (m_doPbinsArray) m_trk_pt_array = book(m_name, "trk_pt_PbinsArray", "p_{T,trk} [GeV]", nPbinsArray, &PbinsArray[0]); 
  m_trk_etaID = book(m_name, "trk_etaID", "#eta_{trk}", nBinsEta, minEta, maxEta); 
  m_trk_etaEMB2 = book(m_name, "trk_etaEMB2", "#eta_{trk}", nBinsEta, minEta, maxEta); 
  m_trk_etaEME2 = book(m_name, "trk_etaEME2", "#eta_{trk}", nBinsEta, minEta, maxEta); 
  if (m_doEtabinsArray) m_trk_eta_array = book(m_name, "trk_eta_EtabinsArray", "#eta_{trk}", nEtabinsArray, &EtabinsArray[0]); 
  m_trk_phiID = book(m_name, "trk_phiID", "#phi_{trk}", nBinsPhi, minPhi, maxPhi); 
  m_trk_phiEMB2 = book(m_name, "trk_phiEMB2", "#phi_{trk}", nBinsPhi, minPhi, maxPhi); 
  m_trk_phiEME2 = book(m_name, "trk_phiEME2", "#phi_{trk}", nBinsPhi, minPhi, maxPhi); 
  m_trk_phi_extra = book(m_name, "trk_phi_extra", "#phi_{trk}", nBinsPhiExtra, minPhiExtra, maxPhiExtra); 
  m_trk_phi_extra2 = book(m_name, "trk_phi_extra2", "#phi_{trk}", nBinsPhiExtra2, minPhiExtra2, maxPhiExtra2); 

  m_trk_p_vs_eta = book(m_name, "trk_p_vs_eta", "#eta_{trk}", nBinsEta, minEta, maxEta, "p_{trk} [GeV]", nBinsP, minP, maxP); 
  m_trk_p_vs_etaEMB2 = book(m_name, "trk_p_vs_etaEMB2", "#eta_{trk}", nBinsEta, minEta, maxEta, "p_{trk} [GeV]", nBinsP, minP, maxP); 
  m_trk_p_vs_etaEME2 = book(m_name, "trk_p_vs_etaEME2", "#eta_{trk}", nBinsEta, minEta, maxEta, "p_{trk} [GeV]", nBinsP, minP, maxP); 
  if (m_doPbinsArray && m_doEtabinsArray) m_trk_p_vs_eta_array = book(m_name, "trk_p_vs_eta_pEtabinsArray", "#eta_{trk}", nEtabinsArray, &EtabinsArray[0], "p_{T,trk} [GeV]", nBinsP, minP, maxP); 

  m_trk_DR_EMB2_ID = book(m_name, std::string("trk_DR_EMB2_ID"), std::string("#Delta R_{trk}(EMB2, ID)"), nBinsDR, minDR, maxDR); 
  m_trk_DEta_EMB2_ID = book(m_name, std::string("trk_DEta_EMB2_ID"), std::string("#Delta #eta_{trk}(EMB2, ID)"), nBinsEta, minEta, maxEta); 
  m_trk_DPhi_EMB2_ID = book(m_name, std::string("trk_DPhi_EMB2_ID"), std::string("#Delta #phi_{trk}(EMB2, ID)"), nBinsPhi, minPhi, maxPhi); 
  m_trk_DR_EMB2_ID_vs_trk_p = book(m_name, std::string("trk_DR_EMB2_ID_vs_trk_p"), "p_{trk} [GeV]", nBinsP, minP, maxP, std::string("#Delta R_{trk}(EMB2, ID)"), minDR, maxDR); 
  m_trk_DEta_EMB2_ID_vs_trk_p = book(m_name, std::string("trk_DEta_EMB2_ID_vs_trk_p"), "p_{trk} [GeV]", nBinsP, minP, maxP, std::string("#Delta #eta_{trk}(EMB2, ID)"), minEta, maxEta); 
  m_trk_DPhi_EMB2_ID_vs_trk_p = book(m_name, std::string("trk_DPhi_EMB2_ID_vs_trk_p"), "p_{trk} [GeV]", nBinsP, minP, maxP, std::string("#Delta #phi_{trk}(EMB2, ID)"), minPhi, maxPhi); 

  m_trk_d0 = book(m_name, "trk_d0", "d0[mm]", 100,-5.0, 5.0 );
  m_trk_d0_s = book(m_name, "trk_d0_s" , "d0[mm]", 100,  -1.0, 1.0 );

  m_trk_z0 = book(m_name, "trk_z0","z0[mm]", 100,-5.0, 5.0 );
  m_trk_z0_s = book(m_name, "trk_z0_s", "z0[mm]", 100,-1.0, 1.0 );
  m_trk_z0sinT = book(m_name, "trk_z0sinT", "z0xsin(#theta)[mm]", 100, -5.0, 5.0 );

  m_trk_chi2Prob = book(m_name, "trk_chi2Prob", "chi2Prob", 100,   -0.01,     1.0);
  m_trk_charge = book(m_name, "trk_charge", "charge",   3,  -1.5,  1.5   );

  m_trk_nSi         = book(m_name, "trk_nSi",         "nSi",         30,   -0.5, 29.5 );
  m_trk_nSiAndDead  = book(m_name, "trk_nSiAndDead",  "nSi(+Dead)",  30,   -0.5, 29.5 );
  m_trk_nSiDead     = book(m_name, "trk_nSiDead",     "nSiDead",     10,   -0.5,  9.5 );
  m_trk_nSCT        = book(m_name, "trk_nSCT",        "nSCTHits",    20,   -0.5, 19.5 );
  m_trk_nPix        = book(m_name, "trk_nPix",        "nPix",        10,   -0.5,  9.5 );
  m_trk_nPixHoles   = book(m_name, "trk_nPixHoles",   "nPixHoles",   10,   -0.5,  9.5 );
  m_trk_nBL         = book(m_name, "trk_nBL",         "nBL",          3,   -0.5,  2.5 );
  m_trk_nTRT        = book(m_name, "trk_nTRT",        "nTRT",        50,   -0.5, 49.5 );
  m_trk_nTRT_vs_p   = book(m_name, "trk_nTRT_vs_p",   "p_{trk} [GeV]", nBinsP, minP, maxP, "nTRT", 50, -0.5, 49.5 );
  m_trk_nTRT_vs_eta = book(m_name, "trk_nTRT_vs_eta", "#eta_{trk}", nBinsEta, minEta, maxEta, "nTRT", 50, -0.5, 49.5 );

  // Tile energy fractions
  m_trk_TileEfrac_100 = book (m_name, std::string("trk_"+m_energyCalib+"_TileEfrac_100"), "E(tile)/E(total)", 100, 0, 5.0); 
  m_trk_TileEfrac_100_vs_trk_p = book (m_name, std::string("trk_"+m_energyCalib+"_TileEfrac_100_vs_trk_p"), "p_{trk} [GeV]", nBinsP, minP, maxP, "E(tile)/E(total)", 100, 0, 5.0); 
  m_trk_TileEfrac_200 = book (m_name, std::string("trk_"+m_energyCalib+"_TileEfrac_200"), "E(tile)/E(total)", 100, 0, 5.0); 
  m_trk_TileEfrac_200_vs_trk_p = book (m_name, std::string("trk_"+m_energyCalib+"_TileEfrac_200_vs_trk_p"), "p_{trk} [GeV]", nBinsP, minP, maxP, "E(tile)/E(total)", 100, 0, 5.0); 

  // basic "sanity" checks of the E/p xAOD derivation
  m_trk_sumE_Tile_100 = book (m_name, std::string("trk_"+m_energyCalib+"_sumE_Tile_100"), "E(Tile)", nBinsE, minE, maxE); 
  m_trk_sumE_Tile_200 = book (m_name, std::string("trk_"+m_energyCalib+"_sumE_Tile_200"), "E(Tile)", nBinsE, minE, maxE); 
  m_trk_sumE_Lar_100 = book (m_name, std::string("trk_"+m_energyCalib+"_sumE_Lar_100"), "E(Lar)", nBinsE, minE, maxE); 
  m_trk_sumE_Lar_200 = book (m_name, std::string("trk_"+m_energyCalib+"_sumE_Lar_200"), "E(Lar)", nBinsE, minE, maxE); 
  m_trk_SumTileLayers_over_HAD_100 = book (m_name, std::string("trk_"+m_energyCalib+"_SumTileLayers_over_HAD_100"), "#Sigma Tile layers / HAD", 100, 0, 5.0); 
  m_trk_SumLarLayers_over_EM_100 = book (m_name, std::string("trk_"+m_energyCalib+"_SumLarLayers_over_EM_100"), "#Sigma Lar layers / EM", 100, 0, 5.0); 
  m_trk_SumHADLayers_over_HAD_100 = book (m_name, std::string("trk_"+m_energyCalib+"_SumHADLayers_over_HAD_100"), "#Sigma Tile+EM layers / EM", 100, 0, 5.0); 
  m_trk_SumEMLayers_over_EM_100 = book (m_name, std::string("trk_"+m_energyCalib+"_SumEMLayers_over_EM_100"), "#Sigma EM layers / EM", 100, 0, 5.0); 
  m_trk_EMandHAD_over_Total_100 = book (m_name, std::string("trk_"+m_energyCalib+"_EMandHAD_over_Total_100"), "(EM+HAD) / Total", 100, 0, 5.0); 
  m_trk_SumAllLayers_over_Total_100 = book (m_name, std::string("trk_"+m_energyCalib+"_SumAllLayers_over_Total_100"), "#Sigma all layers / Total", 100, 0, 5.0); 
  m_trk_SumTileLayers_over_HAD_200 = book (m_name, std::string("trk_"+m_energyCalib+"_SumTileLayers_over_HAD_200"), "#Sigma Tile layers / HAD", 100, 0, 5.0); 
  m_trk_SumLarLayers_over_EM_200 = book (m_name, std::string("trk_"+m_energyCalib+"_SumLarLayers_over_EM_200"), "#Sigma Lar layers / EM", 100, 0, 5.0); 
  m_trk_SumHADLayers_over_HAD_200 = book (m_name, std::string("trk_"+m_energyCalib+"_SumHADLayers_over_HAD_200"), "#Sigma Tile+EM layers / EM", 100, 0, 5.0); 
  m_trk_SumEMLayers_over_EM_200 = book (m_name, std::string("trk_"+m_energyCalib+"_SumEMLayers_over_EM_200"), "#Sigma EM layers / EM", 100, 0, 5.0); 
  m_trk_EMandHAD_over_Total_200 = book (m_name, std::string("trk_"+m_energyCalib+"_EMandHAD_over_Total_200"), "(EM+HAD) / Total", 100, 0, 5.0); 
  m_trk_SumAllLayers_over_Total_200 = book (m_name, std::string("trk_"+m_energyCalib+"_SumAllLayers_over_Total_200"), "#Sigma all layers / Total", 100, 0, 5.0); 

  // calo layer with the highest cluster energy
  m_trk_E_100_highElayer = book (m_name, std::string("trk_"+m_energyCalib+"_E_100_highElayer"), "Calorimeter layer with highest E", 21, 0, 21); 
  m_trk_E_200_highElayer = book (m_name, std::string("trk_"+m_energyCalib+"_E_200_highElayer"), "Calorimeter layer with highest E", 21, 0, 21); 
  m_trk_E_100_highElayer_vs_E = book (m_name, std::string("trk_"+m_energyCalib+"_E_100_highElayer_vs_E"), "E [GeV]", nBinsE, minE, maxE, "Calorimeter layer", 21, 0, 21); 
  m_trk_E_200_highElayer_vs_E = book (m_name, std::string("trk_"+m_energyCalib+"_E_200_highElayer_vs_E"), "E [GeV]", nBinsE, minE, maxE, "Calorimeter layer", 21, 0, 21); 

  // sum of the energy of clusters matched to tracks vs. calo layer
  m_trk_E_100_vs_layer = book(m_name, std::string("trk_"+m_energyCalib+"_E_100_vs_layer"), "Calorimeter layer", 21, 0, 21, "E (#DeltaR<0.4)", nBinsE, minE, maxE); 
  m_trk_E_200_vs_layer = book(m_name, std::string("trk_"+m_energyCalib+"_E_200_vs_layer"), "Calorimeter layer", 21, 0, 21, "E (#DeltaR<0.4)", nBinsE, minE, maxE); 

  // Zero fraction histograms
  if (m_doPbinsArray && m_doEtabinsArray) {
    m_trk_n_E_200 = book(m_name, "trk_n_E_200", "p_{trk} [GeV]", nPbinsArray, &PbinsArray[0], "|#eta_{trk}|", nEtabinsArray, &EtabinsArray[0]); 
    m_trk_n_E_200_l0 = book(m_name, "trk_n_E_200_l0", "p_{trk} [GeV]", nPbinsArray, &PbinsArray[0], "|#eta_{trk}|", nEtabinsArray, &EtabinsArray[0]); 
    m_trk_n_E_200_eq0 = book(m_name, "trk_n_E_200_eq0", "p_{trk} [GeV]", nPbinsArray, &PbinsArray[0], "|#eta_{trk}|", nEtabinsArray, &EtabinsArray[0]); 
    m_trk_n_E_200_leq0 = book(m_name, "trk_n_E_200_leq0", "p_{trk} [GeV]", nPbinsArray, &PbinsArray[0], "|#eta_{trk}|", nEtabinsArray, &EtabinsArray[0]); 
  }

  // create E/p histograms

  // total calorimeter energy
  if (m_doCaloTotal) {
    // dR(trk,cluster) < 0.1
    m_trk_E_Total_100 = book(m_name, std::string("trk_E_Total_"+m_energyCalib+"_0_100"), "E", nBinsE, minE, maxE);
    m_trk_E_Total_100_vs_mu_avg = book(m_name, std::string("trk_E_Total_"+m_energyCalib+"_0_100_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "E", minE, maxE);

    m_eop_Total_100 = book(m_name, std::string("eop_Total_"+m_energyCalib+"_0_100"), "E/p", nBinsEop, minEop, maxEop);
    if (m_doPbinsArray) m_eop_Total_100_vs_trkP = book(m_name, std::string("eop_Total_"+m_energyCalib+"_0_100_vs_trkP"), "p_{trk}", nPbinsArray, &PbinsArray[0], "<E/p>", minEop_l, maxEop_l);
    else m_eop_Total_100_vs_trkP = book(m_name, std::string("eop_Total_"+m_energyCalib+"_0_100_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "<E/p>", minEop_l, maxEop_l);
    if(m_doEtabinsArray) m_eop_Total_100_vs_trkEta = book(m_name, std::string("eop_Total_"+m_energyCalib+"_0_100_vs_trkEta"), "|#eta_{trk}|", nEtabinsArray, &EtabinsArray[0], "<E/p>", minEop_l, maxEop_l);
    else m_eop_Total_100_vs_trkEta = book(m_name, std::string("eop_Total_"+m_energyCalib+"_0_100_vs_trkEta"), "|#eta_{trk}|", nBinsEta, minEta, maxEta, "<E/p>", minEop_l, maxEop_l);
    m_eop_Total_100_vs_trkPhi = book(m_name, std::string("eop_Total_"+m_energyCalib+"_0_100_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "<E/p>", minEop_l, maxEop_l);
    m_eop_Total_100_vs_trkPhiID = book(m_name, std::string("eop_ID_Total_"+m_energyCalib+"_0_100_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "<E/p>", minEop_l, maxEop_l);
    m_eop_Total_100_vs_trkPhi_extra = book(m_name, std::string("eop_extra_Total_"+m_energyCalib+"_0_100_vs_trkPhi"), "#phi_{trk}", nBinsPhiExtra, minPhiExtra, maxPhiExtra, "<E/p>", minEop_l, maxEop_l);
    m_eop_Total_100_vs_trkPhi_extra2 = book(m_name, std::string("eop_extra2_Total_"+m_energyCalib+"_0_100_vs_trkPhi"), "#phi_{trk}", nBinsPhiExtra2, minPhiExtra2, maxPhiExtra2, "<E/p>", minEop_l, maxEop_l);
    m_eop_Total_100_vs_mu = book(m_name, std::string("eop_Total_"+m_energyCalib+"_0_100_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "<E/p>", minEop_l, maxEop_l);
    m_eop_Total_100_vs_mu_avg = book(m_name, std::string("eop_Total_"+m_energyCalib+"_0_100_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "<E/p>", minEop_l, maxEop_l);
    m_eop_Total_100_vs_npv = book(m_name, std::string("eop_Total_"+m_energyCalib+"_0_100_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "<E/p>", minEop_l, maxEop_l);
    m_eop_Total_100_vs_highE_layer = book(m_name, std::string("eop_Total_"+m_energyCalib+"_0_100_vs_highE_layer"), "Calorimeter layer with the highest energy", 21, 0, 21, "<E/p>", minEop_l, maxEop_l);
    if (m_doEtabinsArray && m_doPbinsArray) {
      if (m_doProfileEta) {
        // TProfile binned in p (one profile for each eta bin)
        m_eop_Total_100_profileEta = std::vector<TProfile* >();
        char buffer [200];
        for (unsigned int i=0; i<nEtabinsArray; i++) {
          std::snprintf(buffer, 200, "etaG%dL%d", (int)std::round((EtabinsArray[i]*10)), (int)std::round((EtabinsArray[i+1]*10)));
          TProfile* tmp_profile = book(m_name, std::string(std::string(buffer)+"_eop_Total_"+m_energyCalib+"_0_100"), "p_{trk}", nPbinsArray, &PbinsArray[0], "<E/p>", minEop_l, maxEop_l);
          m_eop_Total_100_profileEta.push_back(tmp_profile);
        }
      }
      if (m_doProfileP) {
        // TProfile binned in eta (one profile for each p bin)
        m_eop_Total_100_profileP = std::vector<TProfile* >();
        char buffer [200];
        for (unsigned int i=0; i<nPbinsArray; i++) {
          std::snprintf(buffer, 200, "pG%dL%d", (int)std::round((PbinsArray[i]*1000)), (int)std::round((PbinsArray[i+1]*1000)));
          TProfile* tmp_profile = book(m_name, std::string(std::string(buffer)+"_eop_Total_"+m_energyCalib+"_0_100"), "eta_{trk}", nEtabinsArray, &EtabinsArray[0], "<E/p>", minEop_l, maxEop_l);
          m_eop_Total_100_profileP.push_back(tmp_profile);
        }
      }
      if (m_doExtraEtaEnergyBinHists) {
        // 1D histograms for E/p distributions
        m_eop_Total_100_EtaEnergyRanges = std::vector<std::vector<TH1D*> >(nPbinsArray); 
        char buffer [200];
        for (unsigned int i=0; i<nPbinsArray; i++) {
          for (unsigned int j=0; j<nEtabinsArray; j++) {
            std::snprintf(buffer, 200, "pG%dL%d_etaG%dL%d", (int)std::round((PbinsArray[i]*1000)), (int)std::round((PbinsArray[i+1]*1000)), (int)std::round((EtabinsArray[j]*10)), (int)std::round((EtabinsArray[j+1]*10)));
            TH1D* tmp_hist = book(m_name, std::string(std::string(buffer)+"_eop_Total_"+m_energyCalib+"_0_100"), "E/p", nBinsEop, minEop, maxEop);
            m_eop_Total_100_EtaEnergyRanges[i].push_back(tmp_hist);
          }
        }
      }
    }
    
    // dR(trk,cluster) < 0.2
    m_trk_E_Total_200 = book(m_name, std::string("trk_E_Total_"+m_energyCalib+"_0_200"), "E", nBinsE, minE, maxE);
    m_trk_E_Total_200_vs_mu_avg = book(m_name, std::string("trk_E_Total_"+m_energyCalib+"_0_200_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "E", minE, maxE);
    m_eop_Total_200 = book(m_name, std::string("eop_Total_"+m_energyCalib+"_0_200"), "E/p", nBinsEop, minEop, maxEop);
    m_eop_Total_200_l = book(m_name, std::string("eop_Total_"+m_energyCalib+"_0_200_l"), "E/p", nBinsEop_l, minEop_l, maxEop_l);
    if (m_doPbinsArray) m_eop_Total_200_vs_trkP = book(m_name, std::string("eop_Total_"+m_energyCalib+"_0_200_vs_trkP"), "p_{trk}", nPbinsArray, &PbinsArray[0], "<E/p>", minEop_l, maxEop_l);
    else m_eop_Total_200_vs_trkP = book(m_name, std::string("eop_Total_"+m_energyCalib+"_0_200_vs_trkP"), "p_{trk}", nBinsP, minE, maxE, "<E/p>", minEop_l, maxEop_l);
    if(m_doEtabinsArray) m_eop_Total_200_vs_trkEta = book(m_name, std::string("eop_Total_"+m_energyCalib+"_0_200_vs_trkEta"), "|#eta_{trk}|", nEtabinsArray, &EtabinsArray[0], "<E/p>", minEop_l, maxEop_l);
    else m_eop_Total_200_vs_trkEta = book(m_name, std::string("eop_Total_"+m_energyCalib+"_0_200_vs_trkEta"), "|#eta_{trk}|", nBinsEta, minEta, maxEta, "<E/p>", minEop_l, maxEop_l);
    m_eop_Total_200_vs_trkPhi = book(m_name, std::string("eop_Total_"+m_energyCalib+"_0_200_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "<E/p>", minEop_l, maxEop_l);
    m_eop_Total_200_vs_trkPhiID = book(m_name, std::string("eop_ID_Total_"+m_energyCalib+"_0_200_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "<E/p>", minEop_l, maxEop_l);
    m_eop_Total_200_vs_trkPhi_extra = book(m_name, std::string("eop_extra_Total_"+m_energyCalib+"_0_200_vs_trkPhi"), "#phi_{trk}", nBinsPhiExtra, minPhiExtra, maxPhiExtra, "<E/p>", minEop_l, maxEop_l);
    m_eop_Total_200_vs_trkPhi_extra2 = book(m_name, std::string("eop_extra2_Total_"+m_energyCalib+"_0_200_vs_trkPhi"), "#phi_{trk}", nBinsPhiExtra2, minPhiExtra2, maxPhiExtra2, "<E/p>", minEop_l, maxEop_l);
    m_eop_Total_200_vs_mu = book(m_name, std::string("eop_Total_"+m_energyCalib+"_0_200_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "<E/p>", minEop_l, maxEop_l);
    m_eop_Total_200_vs_mu_avg = book(m_name, std::string("eop_Total_"+m_energyCalib+"_0_200_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "<E/p>", minEop_l, maxEop_l);
    m_eop_Total_200_vs_npv = book(m_name, std::string("eop_Total_"+m_energyCalib+"_0_200_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "<E/p>", minEop_l, maxEop_l);
    m_eop_Total_200_vs_highE_layer = book(m_name, std::string("eop_Total_"+m_energyCalib+"_0_200_vs_highE_layer"), "Calorimeter layer with the highest energy", 21, 0, 21, "<E/p>", minEop_l, maxEop_l);

    if (m_doEtabinsArray && m_doPbinsArray) {
      if (m_doProfileEta) {
        // TProfile binned in p (one profile for each eta bin)
        m_eop_Total_200_profileEta = std::vector<TProfile* >();
        char buffer [200];
        for (unsigned int i=0; i<nEtabinsArray; i++) {
          std::snprintf(buffer, 200, "etaG%dL%d", (int)std::round((EtabinsArray[i]*10)), (int)std::round((EtabinsArray[i+1]*10)));
          TProfile* tmp_profile = book(m_name, std::string(std::string(buffer)+"_eop_Total_"+m_energyCalib+"_0_200"), "p_{trk}", nPbinsArray, &PbinsArray[0], "<E/p>", minEop_l, maxEop_l);
          m_eop_Total_200_profileEta.push_back(tmp_profile);
        }
      }
      if (m_doProfileP) {
        // TProfile binned in eta (one profile for each p bin)
        m_eop_Total_200_profileP = std::vector<TProfile* >();
        char buffer [200];
        for (unsigned int i=0; i<nPbinsArray; i++) {
          std::snprintf(buffer, 200, "pG%dL%d", (int)std::round((PbinsArray[i]*1000)), (int)std::round((PbinsArray[i+1]*1000)));
          TProfile* tmp_profile = book(m_name, std::string(std::string(buffer)+"_eop_Total_"+m_energyCalib+"_0_200"), "eta_{trk}", nEtabinsArray, &EtabinsArray[0], "<E/p>", minEop_l, maxEop_l);
          m_eop_Total_200_profileP.push_back(tmp_profile);
        }
      }
      if (m_doExtraEtaEnergyBinHists) {
        // 1D histograms for E/p distributions
        m_eop_Total_200_EtaEnergyRanges = std::vector<std::vector<TH1D*> >(nPbinsArray); 
        char buffer [200];
        for (unsigned int i=0; i<nPbinsArray; i++) {
          for (unsigned int j=0; j<nEtabinsArray; j++) {
            std::snprintf(buffer, 200, "pG%dL%d_etaG%dL%d", (int)std::round((PbinsArray[i]*1000)), (int)std::round((PbinsArray[i+1]*1000)), (int)std::round((EtabinsArray[j]*10)), (int)std::round((EtabinsArray[j+1]*10)));
            TH1D* tmp_hist = book(m_name, std::string(std::string(buffer)+"_eop_Total_"+m_energyCalib+"_0_200"), "E/p", nBinsEop, minEop, maxEop);
            m_eop_Total_200_EtaEnergyRanges[i].push_back(tmp_hist);
          }
        }
      }
    }

  } // doCaloTotal

  // EM calorimeter (EMB+EMEC)
  if (m_doCaloEM) {
    // dR(trk,cluster) < 0.1
    m_trk_E_EM_100 = book(m_name, std::string("trk_E_EM_"+m_energyCalib+"_0_100"), "E", nBinsE, minE, maxE);
    m_eop_EM_100 = book(m_name, std::string("eop_EM_"+m_energyCalib+"_0_100"), "E/p", nBinsEop, minEop, maxEop);
    if (m_doPbinsArray) m_eop_EM_100_vs_trkP = book(m_name, std::string("eop_EM_"+m_energyCalib+"_0_100_vs_trkP"), "p_{trk}", nPbinsArray, &PbinsArray[0], "<E/p>", minEop_l, maxEop_l);
    else m_eop_EM_100_vs_trkP = book(m_name, std::string("eop_EM_"+m_energyCalib+"_0_100_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "<E/p>", minEop_l, maxEop_l);
    if(m_doEtabinsArray) m_eop_EM_100_vs_trkEta = book(m_name, std::string("eop_EM_"+m_energyCalib+"_0_100_vs_trkEta"), "|#eta_{trk}|", nEtabinsArray, &EtabinsArray[0], "<E/p>", minEop_l, maxEop_l);
    else m_eop_EM_100_vs_trkEta = book(m_name, std::string("eop_EM_"+m_energyCalib+"_0_100_vs_trkEta"), "|#eta_{trk}|", nBinsEta, minEta, maxEta, "<E/p>", minEop_l, maxEop_l);
    m_eop_EM_100_vs_trkPhi = book(m_name, std::string("eop_EM_"+m_energyCalib+"_0_100_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "<E/p>", minEop_l, maxEop_l);
    m_eop_EM_100_vs_mu = book(m_name, std::string("eop_EM_"+m_energyCalib+"_0_100_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "<E/p>", minEop_l, maxEop_l);
    m_eop_EM_100_vs_mu_avg = book(m_name, std::string("eop_EM_"+m_energyCalib+"_0_100_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "<E/p>", minEop_l, maxEop_l);
    m_eop_EM_100_vs_npv = book(m_name, std::string("eop_EM_"+m_energyCalib+"_0_100_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "<E/p>", minEop_l, maxEop_l);
    m_eop_EM_100_vs_highE_layer = book(m_name, std::string("eop_EM_"+m_energyCalib+"_0_100_vs_highE_layer"), "Calorimeter layer with the highest energy", 21, 0, 21, "<E/p>", minEop_l, maxEop_l);
    if (m_doEtabinsArray && m_doPbinsArray) { 
      if (m_doProfileEta) {
        // TProfile binned in p (one profile for each eta bin)
        m_eop_EM_100_profileEta = std::vector<TProfile* >();
        char buffer [200];
        for (unsigned int i=0; i<nEtabinsArray; i++) {
          std::snprintf(buffer, 200, "etaG%dL%d", (int)std::round((EtabinsArray[i]*10)), (int)std::round((EtabinsArray[i+1]*10)));
          TProfile* tmp_profile = book(m_name, std::string(std::string(buffer)+"_eop_EM_"+m_energyCalib+"_0_100"), "p_{trk}", nPbinsArray, &PbinsArray[0], "<E/p>", minEop_l, maxEop_l);
          m_eop_EM_100_profileEta.push_back(tmp_profile);
        }
      }
      if (m_doProfileP) {
        // TProfile binned in eta (one profile for each p bin)
        m_eop_EM_100_profileP = std::vector<TProfile* >();
        char buffer [200];
        for (unsigned int i=0; i<nPbinsArray; i++) {
          std::snprintf(buffer, 200, "pG%dL%d", (int)std::round((PbinsArray[i]*1000)), (int)std::round((PbinsArray[i+1]*1000)));
          TProfile* tmp_profile = book(m_name, std::string(std::string(buffer)+"_eop_EM_"+m_energyCalib+"_0_100"), "eta_{trk}", nEtabinsArray, &EtabinsArray[0], "<E/p>", minEop_l, maxEop_l);
          m_eop_EM_100_profileP.push_back(tmp_profile);
        }
      }
      if (m_doExtraEtaEnergyBinHists) {
        // 1D histograms for E/p distributions
        m_eop_EM_100_EtaEnergyRanges = std::vector<std::vector<TH1D*> >(nPbinsArray); 
        char buffer [200];
        for (unsigned int i=0; i<nPbinsArray; i++) {
          for (unsigned int j=0; j<nEtabinsArray; j++) {
            std::snprintf(buffer, 200, "pG%dL%d_etaG%dL%d", (int)std::round((PbinsArray[i]*1000)), (int)std::round((PbinsArray[i+1]*1000)), (int)std::round((EtabinsArray[j]*10)), (int)std::round((EtabinsArray[j+1]*10)));
            TH1D* tmp_hist = book(m_name, std::string(std::string(buffer)+"_eop_EM_"+m_energyCalib+"_0_100"), "E/p", nBinsEop, minEop, maxEop);
            m_eop_EM_100_EtaEnergyRanges[i].push_back(tmp_hist);
          }
        }
      }
    }
    // dR(trk,cluster) < 0.2
    m_trk_E_EM_200 = book(m_name, std::string("trk_E_EM_"+m_energyCalib+"_0_200"), "E", nBinsE, minE, maxE);
    m_eop_EM_200 = book(m_name, std::string("eop_EM_"+m_energyCalib+"_0_200"), "E/p", nBinsEop, minEop, maxEop);
    if (m_doPbinsArray) m_eop_EM_200_vs_trkP = book(m_name, std::string("eop_EM_"+m_energyCalib+"_0_200_vs_trkP"), "p_{trk}", nPbinsArray, &PbinsArray[0], "<E/p>", minEop_l, maxEop_l);
    else m_eop_EM_200_vs_trkP = book(m_name, std::string("eop_EM_"+m_energyCalib+"_0_200_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "<E/p>", minEop_l, maxEop_l);
    if(m_doEtabinsArray) m_eop_EM_200_vs_trkEta = book(m_name, std::string("eop_EM_"+m_energyCalib+"_0_200_vs_trkEta"), "|#eta_{trk}|", nEtabinsArray, &EtabinsArray[0], "<E/p>", minEop_l, maxEop_l);
    else m_eop_EM_200_vs_trkEta = book(m_name, std::string("eop_EM_"+m_energyCalib+"_0_200_vs_trkEta"), "|#eta_{trk}|", nBinsEta, minEta, maxEta, "<E/p>", minEop_l, maxEop_l);
    m_eop_EM_200_vs_trkPhi = book(m_name, std::string("eop_EM_"+m_energyCalib+"_0_200_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "<E/p>", minEop_l, maxEop_l);
    m_eop_EM_200_vs_mu = book(m_name, std::string("eop_EM_"+m_energyCalib+"_0_200_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "<E/p>", minEop_l, maxEop_l);
    m_eop_EM_200_vs_mu_avg = book(m_name, std::string("eop_EM_"+m_energyCalib+"_0_200_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "<E/p>", minEop_l, maxEop_l);
    m_eop_EM_200_vs_npv = book(m_name, std::string("eop_EM_"+m_energyCalib+"_0_200_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "<E/p>", minEop_l, maxEop_l);
    m_eop_EM_200_vs_highE_layer = book(m_name, std::string("eop_EM_"+m_energyCalib+"_0_200_vs_highE_layer"), "Calorimeter layer with the highest energy", 21, 0, 21, "<E/p>", minEop_l, maxEop_l);
    if (m_doEtabinsArray && m_doPbinsArray) { 
      if (m_doProfileEta) {
        // TProfile binned in p (one profile for each eta bin)
        m_eop_EM_200_profileEta = std::vector<TProfile* >();
        char buffer [200];
        for (unsigned int i=0; i<nEtabinsArray; i++) {
          std::snprintf(buffer, 200, "etaG%dL%d", (int)std::round((EtabinsArray[i]*10)), (int)std::round((EtabinsArray[i+1]*10)));
          TProfile* tmp_profile = book(m_name, std::string(std::string(buffer)+"_eop_EM_"+m_energyCalib+"_0_200"), "p_{trk}", nPbinsArray, &PbinsArray[0], "<E/p>", minEop_l, maxEop_l);
          m_eop_EM_200_profileEta.push_back(tmp_profile);
        }
      }
      if (m_doProfileP) {
        // TProfile binned in eta (one profile for each p bin)
        m_eop_EM_200_profileP = std::vector<TProfile* >();
        char buffer [200];
        for (unsigned int i=0; i<nPbinsArray; i++) {
          std::snprintf(buffer, 200, "pG%dL%d", (int)std::round((PbinsArray[i]*1000)), (int)std::round((PbinsArray[i+1]*1000)));
          TProfile* tmp_profile = book(m_name, std::string(std::string(buffer)+"_eop_EM_"+m_energyCalib+"_0_200"), "eta_{trk}", nEtabinsArray, &EtabinsArray[0], "<E/p>", minEop_l, maxEop_l);
          m_eop_EM_200_profileP.push_back(tmp_profile);
        }
      }
      if (m_doExtraEtaEnergyBinHists) {
        // 1D histograms for E/p distributions
        m_eop_EM_200_EtaEnergyRanges = std::vector<std::vector<TH1D*> >(nPbinsArray); 
        char buffer [200];
        for (unsigned int i=0; i<nPbinsArray; i++) {
          for (unsigned int j=0; j<nEtabinsArray; j++) {
            std::snprintf(buffer, 200, "pG%dL%d_etaG%dL%d", (int)std::round((PbinsArray[i]*1000)), (int)std::round((PbinsArray[i+1]*1000)), (int)std::round((EtabinsArray[j]*10)), (int)std::round((EtabinsArray[j+1]*10)));
            TH1D* tmp_hist = book(m_name, std::string(std::string(buffer)+"_eop_EM_"+m_energyCalib+"_0_200"), "E/p", nBinsEop, minEop, maxEop);
            m_eop_EM_200_EtaEnergyRanges[i].push_back(tmp_hist);
          }
        }
      }
    }

    // No clusters in HAD matched to trk 
    m_eop_Total_200_noHAD = book(m_name, std::string("eop_Total_noHAD_"+m_energyCalib+"_0_200"), "E/p", nBinsEop, minEop, maxEop);
    m_eop_EM_200_noHAD = book(m_name, std::string("eop_EM_noHAD_"+m_energyCalib+"_0_200"), "E/p", nBinsEop, minEop, maxEop);
    m_eop_HAD_200_noHAD = book(m_name, std::string("eop_HAD_noHAD_"+m_energyCalib+"_0_200"), "E/p", nBinsEop, minEop, maxEop);

  } // doCaloEM

  // HAD calorimeter (HEC+TileBarrel+TileGap+TileExtBarrel)
  if (m_doCaloHAD) {
    // dR(trk,cluster) < 0.1
    m_trk_E_HAD_100 = book(m_name, std::string("trk_E_HAD_"+m_energyCalib+"_0_100"), "E", nBinsE, minE, maxE);
    m_eop_HAD_100 = book(m_name, std::string("eop_HAD_"+m_energyCalib+"_0_100"), "E/p", nBinsEop, minEop, maxEop);
    if (m_doPbinsArray) m_eop_HAD_100_vs_trkP = book(m_name, std::string("eop_HAD_"+m_energyCalib+"_0_100_vs_trkP"), "p_{trk}", nPbinsArray, &PbinsArray[0], "<E/p>", minEop_l, maxEop_l);
    else m_eop_HAD_100_vs_trkP = book(m_name, std::string("eop_HAD_"+m_energyCalib+"_0_100_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "<E/p>", minEop_l, maxEop_l);
    if(m_doEtabinsArray) m_eop_HAD_100_vs_trkEta = book(m_name, std::string("eop_HAD_"+m_energyCalib+"_0_100_vs_trkEta"), "|#eta_{trk}|", nEtabinsArray, &EtabinsArray[0], "<E/p>", minEop_l, maxEop_l);
    else m_eop_HAD_100_vs_trkEta = book(m_name, std::string("eop_HAD_"+m_energyCalib+"_0_100_vs_trkEta"), "|#eta_{trk}|", nBinsEta, minEta, maxEta, "<E/p>", minEop_l, maxEop_l);
    m_eop_HAD_100_vs_trkPhi = book(m_name, std::string("eop_HAD_"+m_energyCalib+"_0_100_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "<E/p>", minEop_l, maxEop_l);
    m_eop_HAD_100_vs_mu = book(m_name, std::string("eop_HAD_"+m_energyCalib+"_0_100_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "<E/p>", minEop_l, maxEop_l);
    m_eop_HAD_100_vs_mu_avg = book(m_name, std::string("eop_HAD_"+m_energyCalib+"_0_100_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "<E/p>", minEop_l, maxEop_l);
    m_eop_HAD_100_vs_npv = book(m_name, std::string("eop_HAD_"+m_energyCalib+"_0_100_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "<E/p>", minEop_l, maxEop_l);
    m_eop_HAD_100_vs_highE_layer = book(m_name, std::string("eop_HAD_"+m_energyCalib+"_0_100_vs_highE_layer"), "Calorimeter layer with the highest energy", 21, 0, 21, "<E/p>", minEop_l, maxEop_l);
    if (m_doEtabinsArray && m_doPbinsArray) { 
      if (m_doProfileEta) {
        // TProfile binned in p (one profile for each eta bin)
        m_eop_HAD_100_profileEta = std::vector<TProfile* >();
        char buffer [200];
        for (unsigned int i=0; i<nEtabinsArray; i++) {
          std::snprintf(buffer, 200, "etaG%dL%d", (int)std::round((EtabinsArray[i]*10)), (int)std::round((EtabinsArray[i+1]*10)));
          TProfile* tmp_profile = book(m_name, std::string(std::string(buffer)+"_eop_HAD_"+m_energyCalib+"_0_100"), "p_{trk}", nPbinsArray, &PbinsArray[0], "<E/p>", minEop_l, maxEop_l);
          m_eop_HAD_100_profileEta.push_back(tmp_profile);
        }
      }
      if (m_doProfileP) {
        // TProfile binned in eta (one profile for each p bin)
        m_eop_HAD_100_profileP = std::vector<TProfile* >();
        char buffer [200];
        for (unsigned int i=0; i<nPbinsArray; i++) {
          std::snprintf(buffer, 200, "pG%dL%d", (int)std::round((PbinsArray[i]*1000)), (int)std::round((PbinsArray[i+1]*1000)));
          TProfile* tmp_profile = book(m_name, std::string(std::string(buffer)+"_eop_HAD_"+m_energyCalib+"_0_100"), "eta_{trk}", nEtabinsArray, &EtabinsArray[0], "<E/p>", minEop_l, maxEop_l);
          m_eop_HAD_100_profileP.push_back(tmp_profile);
        }
      }
      if (m_doExtraEtaEnergyBinHists) {
        // 1D histograms for E/p distributions
        m_eop_HAD_100_EtaEnergyRanges = std::vector<std::vector<TH1D*> >(nPbinsArray); 
        char buffer [200];
        for (unsigned int i=0; i<nPbinsArray; i++) {
          for (unsigned int j=0; j<nEtabinsArray; j++) {
            std::snprintf(buffer, 200, "pG%dL%d_etaG%dL%d", (int)std::round((PbinsArray[i]*1000)), (int)std::round((PbinsArray[i+1]*1000)), (int)std::round((EtabinsArray[j]*10)), (int)std::round((EtabinsArray[j+1]*10)));
            TH1D* tmp_hist = book(m_name, std::string(std::string(buffer)+"_eop_HAD_"+m_energyCalib+"_0_100"), "E/p", nBinsEop, minEop, maxEop);
            m_eop_HAD_100_EtaEnergyRanges[i].push_back(tmp_hist);
          }
        }
      }
    }
    // dR(trk,cluster) < 0.2
    m_trk_E_HAD_200 = book(m_name, std::string("trk_E_HAD_"+m_energyCalib+"_0_200"), "E", nBinsE, minE, maxE);
    m_eop_HAD_200 = book(m_name, std::string("eop_HAD_"+m_energyCalib+"_0_200"), "E/p", nBinsEop, minEop, maxEop);
    if (m_doPbinsArray) m_eop_HAD_200_vs_trkP = book(m_name, std::string("eop_HAD_"+m_energyCalib+"_0_200_vs_trkP"), "p_{trk}", nPbinsArray, &PbinsArray[0], "<E/p>", minEop_l, maxEop_l);
    else m_eop_HAD_200_vs_trkP = book(m_name, std::string("eop_HAD_"+m_energyCalib+"_0_200_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "<E/p>", minEop_l, maxEop_l);
    if(m_doEtabinsArray) m_eop_HAD_200_vs_trkEta = book(m_name, std::string("eop_HAD_"+m_energyCalib+"_0_200_vs_trkEta"), "|#eta_{trk}|", nEtabinsArray, &EtabinsArray[0], "<E/p>", minEop_l, maxEop_l);
    else m_eop_HAD_200_vs_trkEta = book(m_name, std::string("eop_HAD_"+m_energyCalib+"_0_200_vs_trkEta"), "|#eta_{trk}|", nBinsEta, minEta, maxEta, "<E/p>", minEop_l, maxEop_l);
    m_eop_HAD_200_vs_trkPhi = book(m_name, std::string("eop_HAD_"+m_energyCalib+"_0_200_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "<E/p>", minEop_l, maxEop_l);
    m_eop_HAD_200_vs_mu = book(m_name, std::string("eop_HAD_"+m_energyCalib+"_0_200_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "<E/p>", minEop_l, maxEop_l);
    m_eop_HAD_200_vs_mu_avg = book(m_name, std::string("eop_HAD_"+m_energyCalib+"_0_200_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "<E/p>", minEop_l, maxEop_l);
    m_eop_HAD_200_vs_npv = book(m_name, std::string("eop_HAD_"+m_energyCalib+"_0_200_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "<E/p>", minEop_l, maxEop_l);
    m_eop_HAD_200_vs_highE_layer = book(m_name, std::string("eop_HAD_"+m_energyCalib+"_0_200_vs_highE_layer"), "Calorimeter layer with the highest energy", 21, 0, 21, "<E/p>", minEop_l, maxEop_l);
    if (m_doEtabinsArray && m_doPbinsArray) { 
      if (m_doProfileEta) {
        // TProfile binned in p (one profile for each eta bin)
        m_eop_HAD_200_profileEta = std::vector<TProfile* >();
        char buffer [200];
        for (unsigned int i=0; i<nEtabinsArray; i++) {
          std::snprintf(buffer, 200, "etaG%dL%d", (int)std::round((EtabinsArray[i]*10)), (int)std::round((EtabinsArray[i+1]*10)));
          TProfile* tmp_profile = book(m_name, std::string(std::string(buffer)+"_eop_HAD_"+m_energyCalib+"_0_200"), "p_{trk}", nPbinsArray, &PbinsArray[0], "<E/p>", minEop_l, maxEop_l);
          m_eop_HAD_200_profileEta.push_back(tmp_profile);
        }
      }
      if (m_doProfileP) {
        // TProfile binned in eta (one profile for each p bin)
        m_eop_HAD_200_profileP = std::vector<TProfile* >();
        char buffer [200];
        for (unsigned int i=0; i<nPbinsArray; i++) {
          std::snprintf(buffer, 200, "pG%dL%d", (int)std::round((PbinsArray[i]*1000)), (int)std::round((PbinsArray[i+1]*1000)));
          TProfile* tmp_profile = book(m_name, std::string(std::string(buffer)+"_eop_HAD_"+m_energyCalib+"_0_200"), "eta_{trk}", nEtabinsArray, &EtabinsArray[0], "<E/p>", minEop_l, maxEop_l);
          m_eop_HAD_200_profileP.push_back(tmp_profile);
        }
      }
      if (m_doExtraEtaEnergyBinHists) {
        // 1D histograms for E/p distributions
        m_eop_HAD_200_EtaEnergyRanges = std::vector<std::vector<TH1D*> >(nPbinsArray); 
        char buffer [200];
        for (unsigned int i=0; i<nPbinsArray; i++) {
          for (unsigned int j=0; j<nEtabinsArray; j++) {
            std::snprintf(buffer, 200, "pG%dL%d_etaG%dL%d", (int)std::round((PbinsArray[i]*1000)), (int)std::round((PbinsArray[i+1]*1000)), (int)std::round((EtabinsArray[j]*10)), (int)std::round((EtabinsArray[j+1]*10)));
            TH1D* tmp_hist = book(m_name, std::string(std::string(buffer)+"_eop_HAD_"+m_energyCalib+"_0_200"), "E/p", nBinsEop, minEop, maxEop);
            m_eop_HAD_200_EtaEnergyRanges[i].push_back(tmp_hist);
          }
        }
      }
    }

    // MIP requirement for HAD calorimeter
    m_eop_Total_200_MIP = book(m_name, std::string("eop_Total_MIP_"+m_energyCalib+"_0_200"), "E/p", nBinsEop, minEop, maxEop);
    m_eop_EM_200_MIP = book(m_name, std::string("eop_EM_MIP_"+m_energyCalib+"_0_200"), "E/p", nBinsEop, minEop, maxEop);
    m_eop_HAD_200_MIP = book(m_name, std::string("eop_HAD_MIP_"+m_energyCalib+"_0_200"), "E/p", nBinsEop, minEop, maxEop);

  } // doCaloHAD

  // background subtraction, the way it was done in run 1
  if (m_doBgSubtr) {
    m_eop_EM_BG = book(m_name, std::string("eop_EM_BG_"+m_energyCalib), "E/p Background", nBinsEop, minEop, maxEop);
    if (m_doPbinsArray) m_eop_EM_BG_vs_trkP = book(m_name, std::string("eop_EM_BG_"+m_energyCalib+"_vs_trkP"), "p_{trk}", nPbinsArray, &PbinsArray[0], "<E/p>", minEop_l, maxEop_l);
    else m_eop_EM_BG_vs_trkP = book(m_name, std::string("eop_EM_BG_"+m_energyCalib+"_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "<E/p>", minEop_l, maxEop_l);
    if(m_doEtabinsArray) m_eop_EM_BG_vs_trkEta = book(m_name, std::string("eop_EM_BG_"+m_energyCalib+"_vs_trkEta"), "|#eta_{trk}|", nEtabinsArray, &EtabinsArray[0], "<E/p>", minEop_l, maxEop_l);
    else m_eop_EM_BG_vs_trkEta = book(m_name, std::string("eop_EM_BG_"+m_energyCalib+"_vs_trkEta"), "|#eta_{trk}|", nBinsEta, minEta, maxEta, "<E/p>", minEop_l, maxEop_l);
    m_eop_EM_BG_vs_trkPhi = book(m_name, std::string("eop_EM_BG_"+m_energyCalib+"_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "<E/p>", minEop_l, maxEop_l);
    m_eop_EM_BG_vs_mu = book(m_name, std::string("eop_EM_BG_"+m_energyCalib+"_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "<E/p>", minEop_l, maxEop_l);
    m_eop_EM_BG_vs_mu_avg = book(m_name, std::string("eop_EM_BG_"+m_energyCalib+"_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "<E/p>", minEop_l, maxEop_l);
    m_eop_EM_BG_vs_npv = book(m_name, std::string("eop_EM_BG_"+m_energyCalib+"_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "<E/p>", minEop_l, maxEop_l);
    if (m_doEtabinsArray && m_doPbinsArray) {
      if (m_doProfileEta) {
        // TProfile binned in p (one profile for each eta bin)
        m_eop_EM_BG_profileEta = std::vector<TProfile* >();
        char buffer [200];
        for (unsigned int i=0; i<nEtabinsArray; i++) {
          std::snprintf(buffer, 200, "etaG%dL%d", (int)std::round((EtabinsArray[i]*10)), (int)std::round((EtabinsArray[i+1]*10)));
          TProfile* tmp_profile = book(m_name, std::string(std::string(buffer)+"_eop_EM_BG_"+m_energyCalib), "p_{trk}", nPbinsArray, &PbinsArray[0], "<E/p>", minEop_l, maxEop_l);
          m_eop_EM_BG_profileEta.push_back(tmp_profile);
        }
      }
      if (m_doProfileP) {
        // TProfile binned in eta (one profile for each p bin)
        m_eop_EM_BG_profileP = std::vector<TProfile* >();
        char buffer [200];
        for (unsigned int i=0; i<nPbinsArray; i++) {
          std::snprintf(buffer, 200, "pG%dL%d", (int)std::round((PbinsArray[i]*1000)), (int)std::round((PbinsArray[i+1]*1000)));
          TProfile* tmp_profile = book(m_name, std::string(std::string(buffer)+"_eop_EM_BG_"+m_energyCalib), "eta_{trk}", nEtabinsArray, &EtabinsArray[0], "<E/p>", minEop_l, maxEop_l);
          m_eop_EM_BG_profileP.push_back(tmp_profile);
        }
      }
      if (m_doExtraEtaEnergyBinHists) {
        // 1D histograms for E/p distributions
        m_eop_EM_BG_EtaEnergyRanges = std::vector<std::vector<TH1D*> >(nPbinsArray); 
        char buffer [200];
        for (unsigned int i=0; i<nPbinsArray; i++) {
          for (unsigned int j=0; j<nEtabinsArray; j++) {
            std::snprintf(buffer, 200, "pG%dL%d_etaG%dL%d", (int)std::round((PbinsArray[i]*1000)), (int)std::round((PbinsArray[i+1]*1000)), (int)std::round((EtabinsArray[j]*10)), (int)std::round((EtabinsArray[j+1]*10)));
            TH1D* tmp_hist = book(m_name, std::string(std::string(buffer)+"_eop_EM_BG_"+m_energyCalib), "E/p", nBinsEop, minEop, maxEop);
            m_eop_EM_BG_EtaEnergyRanges[i].push_back(tmp_hist);
          }
        }
      }
    }
  } // END doBgSubr

  // E/p where E is maximum in either Tile layer A, BC, or D
  if (m_doTileLayer) {
    // dR(trk,cluster) < 0.1
    // max in Tile Layer A
    m_trk_E_highTileA_100 = book(m_name, std::string("trk_E_highTileA_"+m_energyCalib+"_0_100"), "E", nBinsE, minE, maxE);
    m_eop_highTileA_100 = book(m_name, std::string("eop_highTileA_"+m_energyCalib+"_0_100"), "E/p", nBinsEop, minEop, maxEop);
    if (m_doPbinsArray) m_eop_highTileA_100_vs_trkP = book(m_name, std::string("eop_highTileA_"+m_energyCalib+"_0_100_vs_trkP"), "p_{trk}", nPbinsArray, &PbinsArray[0], "<E/p>", minEop_l, maxEop_l);
    else m_eop_highTileA_100_vs_trkP = book(m_name, std::string("eop_highTileA_"+m_energyCalib+"_0_100_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "<E/p>", minEop_l, maxEop_l);
    if(m_doEtabinsArray) m_eop_highTileA_100_vs_trkEta = book(m_name, std::string("eop_highTileA_"+m_energyCalib+"_0_100_vs_trkEta"), "|#eta_{trk}|", nEtabinsArray, &EtabinsArray[0], "<E/p>", minEop_l, maxEop_l);
    else m_eop_highTileA_100_vs_trkEta = book(m_name, std::string("eop_highTileA_"+m_energyCalib+"_0_100_vs_trkEta"), "|#eta_{trk}|", nBinsEta, minEta, maxEta, "<E/p>", minEop_l, maxEop_l);
    m_eop_highTileA_100_vs_trkPhi = book(m_name, std::string("eop_highTileA_"+m_energyCalib+"_0_100_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "<E/p>", minEop_l, maxEop_l);
    m_eop_highTileA_100_vs_mu = book(m_name, std::string("eop_highTileA_"+m_energyCalib+"_0_100_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "<E/p>", minEop_l, maxEop_l);
    m_eop_highTileA_100_vs_mu_avg = book(m_name, std::string("eop_highTileA_"+m_energyCalib+"_0_100_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "<E/p>", minEop_l, maxEop_l);
    m_eop_highTileA_100_vs_npv = book(m_name, std::string("eop_highTileA_"+m_energyCalib+"_0_100_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "<E/p>", minEop_l, maxEop_l);
    m_eop_highTileA_100_vs_TileA_E = book(m_name, std::string("eop_highTileA_"+m_energyCalib+"_0_100_vs_TileA_E"), "E(Tile, layer A)", nBinsE, minE, maxE, "<E/p>", minEop_l, maxEop_l);
    if (m_doEtabinsArray && m_doPbinsArray) { 
      if (m_doProfileEta) {
        // TProfile binned in p (one profile for each eta bin)
        m_eop_highTileA_100_profileEta = std::vector<TProfile* >();
        char buffer [200];
        for (unsigned int i=0; i<nEtabinsArray; i++) {
          std::snprintf(buffer, 200, "etaG%dL%d", (int)std::round((EtabinsArray[i]*10)), (int)std::round((EtabinsArray[i+1]*10)));
          TProfile* tmp_profile = book(m_name, std::string(std::string(buffer)+"_eop_highTileA_"+m_energyCalib+"_0_100"), "p_{trk}", nPbinsArray, &PbinsArray[0], "<E/p>", minEop_l, maxEop_l);
          m_eop_highTileA_100_profileEta.push_back(tmp_profile);
        }
      }
      if (m_doProfileP) {
        // TProfile binned in eta (one profile for each p bin)
        m_eop_highTileA_100_profileP = std::vector<TProfile* >();
        char buffer [200];
        for (unsigned int i=0; i<nPbinsArray; i++) {
          std::snprintf(buffer, 200, "pG%dL%d", (int)std::round((PbinsArray[i]*1000)), (int)std::round((PbinsArray[i+1]*1000)));
          TProfile* tmp_profile = book(m_name, std::string(std::string(buffer)+"_eop_highTileA_"+m_energyCalib+"_0_100"), "eta_{trk}", nEtabinsArray, &EtabinsArray[0], "<E/p>", minEop_l, maxEop_l);
          m_eop_highTileA_100_profileP.push_back(tmp_profile);
        }
      }
      if (m_doExtraEtaEnergyBinHists) {
        // 1D histograms for E/p distributions
        m_eop_highTileA_100_EtaEnergyRanges = std::vector<std::vector<TH1D*> >(nPbinsArray); 
        char buffer [200];
        for (unsigned int i=0; i<nPbinsArray; i++) {
          for (unsigned int j=0; j<nEtabinsArray; j++) {
            std::snprintf(buffer, 200, "pG%dL%d_etaG%dL%d", (int)std::round((PbinsArray[i]*1000)), (int)std::round((PbinsArray[i+1]*1000)), (int)std::round((EtabinsArray[j]*10)), (int)std::round((EtabinsArray[j+1]*10)));
            TH1D* tmp_hist = book(m_name, std::string(std::string(buffer)+"_eop_highTileA_"+m_energyCalib+"_0_100"), "E/p", nBinsEop, minEop, maxEop);
            m_eop_highTileA_100_EtaEnergyRanges[i].push_back(tmp_hist);
          }
        }
      }
    }
    // max in Tile Layer BC
    m_trk_E_highTileB_100 = book(m_name, std::string("trk_E_highTileB_"+m_energyCalib+"_0_100"), "E", nBinsE, minE, maxE);
    m_eop_highTileB_100 = book(m_name, std::string("eop_highTileB_"+m_energyCalib+"_0_100"), "E/p", nBinsEop, minEop, maxEop);
    if (m_doPbinsArray) m_eop_highTileB_100_vs_trkP = book(m_name, std::string("eop_highTileB_"+m_energyCalib+"_0_100_vs_trkP"), "p_{trk}", nPbinsArray, &PbinsArray[0], "<E/p>", minEop_l, maxEop_l);
    else m_eop_highTileB_100_vs_trkP = book(m_name, std::string("eop_highTileB_"+m_energyCalib+"_0_100_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "<E/p>", minEop_l, maxEop_l);
    if(m_doEtabinsArray) m_eop_highTileB_100_vs_trkEta = book(m_name, std::string("eop_highTileB_"+m_energyCalib+"_0_100_vs_trkEta"), "|#eta_{trk}|", nEtabinsArray, &EtabinsArray[0], "<E/p>", minEop_l, maxEop_l);
    else m_eop_highTileB_100_vs_trkEta = book(m_name, std::string("eop_highTileB_"+m_energyCalib+"_0_100_vs_trkEta"), "|#eta_{trk}|", nBinsEta, minEta, maxEta, "<E/p>", minEop_l, maxEop_l);
    m_eop_highTileB_100_vs_trkPhi = book(m_name, std::string("eop_highTileB_"+m_energyCalib+"_0_100_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "<E/p>", minEop_l, maxEop_l);
    m_eop_highTileB_100_vs_mu = book(m_name, std::string("eop_highTileB_"+m_energyCalib+"_0_100_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "<E/p>", minEop_l, maxEop_l);
    m_eop_highTileB_100_vs_mu_avg = book(m_name, std::string("eop_highTileB_"+m_energyCalib+"_0_100_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "<E/p>", minEop_l, maxEop_l);
    m_eop_highTileB_100_vs_npv = book(m_name, std::string("eop_highTileB_"+m_energyCalib+"_0_100_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "<E/p>", minEop_l, maxEop_l);
    m_eop_highTileB_100_vs_TileB_E = book(m_name, std::string("eop_highTileB_"+m_energyCalib+"_0_100_vs_TileB_E"), "E(Tile, layer B)", nBinsE, minE, maxE, "<E/p>", minEop_l, maxEop_l);
    if (m_doEtabinsArray && m_doPbinsArray) { 
      if (m_doProfileEta) {
        // TProfile binned in p (one profile for each eta bin)
        m_eop_highTileB_100_profileEta = std::vector<TProfile* >();
        char buffer [200];
        for (unsigned int i=0; i<nEtabinsArray; i++) {
          std::snprintf(buffer, 200, "etaG%dL%d", (int)std::round((EtabinsArray[i]*10)), (int)std::round((EtabinsArray[i+1]*10)));
          TProfile* tmp_profile = book(m_name, std::string(std::string(buffer)+"_eop_highTileB_"+m_energyCalib+"_0_100"), "p_{trk}", nPbinsArray, &PbinsArray[0], "<E/p>", minEop_l, maxEop_l);
          m_eop_highTileB_100_profileEta.push_back(tmp_profile);
        }
      }
      if (m_doProfileP) {
        // TProfile binned in eta (one profile for each p bin)
        m_eop_highTileB_100_profileP = std::vector<TProfile* >();
        char buffer [200];
        for (unsigned int i=0; i<nPbinsArray; i++) {
          std::snprintf(buffer, 200, "pG%dL%d", (int)std::round((PbinsArray[i]*1000)), (int)std::round((PbinsArray[i+1]*1000)));
          TProfile* tmp_profile = book(m_name, std::string(std::string(buffer)+"_eop_highTileB_"+m_energyCalib+"_0_100"), "eta_{trk}", nEtabinsArray, &EtabinsArray[0], "<E/p>", minEop_l, maxEop_l);
          m_eop_highTileB_100_profileP.push_back(tmp_profile);
        }
      }
      if (m_doExtraEtaEnergyBinHists) {
        // 1D histograms for E/p distributions
        m_eop_highTileB_100_EtaEnergyRanges = std::vector<std::vector<TH1D*> >(nPbinsArray); 
        char buffer [200];
        for (unsigned int i=0; i<nPbinsArray; i++) {
          for (unsigned int j=0; j<nEtabinsArray; j++) {
            std::snprintf(buffer, 200, "pG%dL%d_etaG%dL%d", (int)std::round((PbinsArray[i]*1000)), (int)std::round((PbinsArray[i+1]*1000)), (int)std::round((EtabinsArray[j]*10)), (int)std::round((EtabinsArray[j+1]*10)));
            TH1D* tmp_hist = book(m_name, std::string(std::string(buffer)+"_eop_highTileB_"+m_energyCalib+"_0_100"), "E/p", nBinsEop, minEop, maxEop);
            m_eop_highTileB_100_EtaEnergyRanges[i].push_back(tmp_hist);
          }
        }
      }
    }
    // max in Tile Layer D
    m_trk_E_highTileD_100 = book(m_name, std::string("trk_E_highTileD_"+m_energyCalib+"_0_100"), "E", nBinsE, minE, maxE);
    m_eop_highTileD_100 = book(m_name, std::string("eop_highTileD_"+m_energyCalib+"_0_100"), "E/p", nBinsEop, minEop, maxEop);
    if (m_doPbinsArray) m_eop_highTileD_100_vs_trkP = book(m_name, std::string("eop_highTileD_"+m_energyCalib+"_0_100_vs_trkP"), "p_{trk}", nPbinsArray, &PbinsArray[0], "<E/p>", minEop_l, maxEop_l);
    else m_eop_highTileD_100_vs_trkP = book(m_name, std::string("eop_highTileD_"+m_energyCalib+"_0_100_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "<E/p>", minEop_l, maxEop_l);
    if(m_doEtabinsArray) m_eop_highTileD_100_vs_trkEta = book(m_name, std::string("eop_highTileD_"+m_energyCalib+"_0_100_vs_trkEta"), "|#eta_{trk}|", nEtabinsArray, &EtabinsArray[0], "<E/p>", minEop_l, maxEop_l);
    else m_eop_highTileD_100_vs_trkEta = book(m_name, std::string("eop_highTileD_"+m_energyCalib+"_0_100_vs_trkEta"), "|#eta_{trk}|", nBinsEta, minEta, maxEta, "<E/p>", minEop_l, maxEop_l);
    m_eop_highTileD_100_vs_trkPhi = book(m_name, std::string("eop_highTileD_"+m_energyCalib+"_0_100_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "<E/p>", minEop_l, maxEop_l);
    m_eop_highTileD_100_vs_mu = book(m_name, std::string("eop_highTileD_"+m_energyCalib+"_0_100_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "<E/p>", minEop_l, maxEop_l);
    m_eop_highTileD_100_vs_mu_avg = book(m_name, std::string("eop_highTileD_"+m_energyCalib+"_0_100_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "<E/p>", minEop_l, maxEop_l);
    m_eop_highTileD_100_vs_npv = book(m_name, std::string("eop_highTileD_"+m_energyCalib+"_0_100_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "<E/p>", minEop_l, maxEop_l);
    m_eop_highTileD_100_vs_TileD_E = book(m_name, std::string("eop_highTileD_"+m_energyCalib+"_0_100_vs_TileD_E"), "E(Tile, layer D)", nBinsE, minE, maxE, "<E/p>", minEop_l, maxEop_l);
    if (m_doEtabinsArray && m_doPbinsArray) { 
      if (m_doProfileEta) {
        // TProfile binned in p (one profile for each eta bin)
        m_eop_highTileD_100_profileEta = std::vector<TProfile* >();
        char buffer [200];
        for (unsigned int i=0; i<nEtabinsArray; i++) {
          std::snprintf(buffer, 200, "etaG%dL%d", (int)std::round((EtabinsArray[i]*10)), (int)std::round((EtabinsArray[i+1]*10)));
          TProfile* tmp_profile = book(m_name, std::string(std::string(buffer)+"_eop_highTileD_"+m_energyCalib+"_0_100"), "p_{trk}", nPbinsArray, &PbinsArray[0], "<E/p>", minEop_l, maxEop_l);
          m_eop_highTileD_100_profileEta.push_back(tmp_profile);
        }
      }
      if (m_doProfileP) {
        // TProfile binned in eta (one profile for each p bin)
        m_eop_highTileD_100_profileP = std::vector<TProfile* >();
        char buffer [200];
        for (unsigned int i=0; i<nPbinsArray; i++) {
          std::snprintf(buffer, 200, "pG%dL%d", (int)std::round((PbinsArray[i]*1000)), (int)std::round((PbinsArray[i+1]*1000)));
          TProfile* tmp_profile = book(m_name, std::string(std::string(buffer)+"_eop_highTileD_"+m_energyCalib+"_0_100"), "eta_{trk}", nEtabinsArray, &EtabinsArray[0], "<E/p>", minEop_l, maxEop_l);
          m_eop_highTileD_100_profileP.push_back(tmp_profile);
        }
      }
      if (m_doExtraEtaEnergyBinHists) {
        // 1D histograms for E/p distributions
        m_eop_highTileD_100_EtaEnergyRanges = std::vector<std::vector<TH1D*> >(nPbinsArray); 
        char buffer [200];
        for (unsigned int i=0; i<nPbinsArray; i++) {
          for (unsigned int j=0; j<nEtabinsArray; j++) {
            std::snprintf(buffer, 200, "pG%dL%d_etaG%dL%d", (int)std::round((PbinsArray[i]*1000)), (int)std::round((PbinsArray[i+1]*1000)), (int)std::round((EtabinsArray[j]*10)), (int)std::round((EtabinsArray[j+1]*10)));
            TH1D* tmp_hist = book(m_name, std::string(std::string(buffer)+"_eop_highTileD_"+m_energyCalib+"_0_100"), "E/p", nBinsEop, minEop, maxEop);
            m_eop_highTileD_100_EtaEnergyRanges[i].push_back(tmp_hist);
          }
        }
      }
    }
    // dR(trk,cluster) < 0.2
    // max in Tile Layer A
    m_trk_E_highTileA_200 = book(m_name, std::string("trk_E_highTileA_"+m_energyCalib+"_0_200"), "E", nBinsE, minE, maxE);
    m_eop_highTileA_200 = book(m_name, std::string("eop_highTileA_"+m_energyCalib+"_0_200"), "E/p", nBinsEop, minEop, maxEop);
    if (m_doPbinsArray) m_eop_highTileA_200_vs_trkP = book(m_name, std::string("eop_highTileA_"+m_energyCalib+"_0_200_vs_trkP"), "p_{trk}", nPbinsArray, &PbinsArray[0], "<E/p>", minEop_l, maxEop_l);
    else m_eop_highTileA_200_vs_trkP = book(m_name, std::string("eop_highTileA_"+m_energyCalib+"_0_200_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "<E/p>", minEop_l, maxEop_l);
    if(m_doEtabinsArray) m_eop_highTileA_200_vs_trkEta = book(m_name, std::string("eop_highTileA_"+m_energyCalib+"_0_200_vs_trkEta"), "|#eta_{trk}|", nEtabinsArray, &EtabinsArray[0], "<E/p>", minEop_l, maxEop_l);
    else m_eop_highTileA_200_vs_trkEta = book(m_name, std::string("eop_highTileA_"+m_energyCalib+"_0_200_vs_trkEta"), "|#eta_{trk}|", nBinsEta, minEta, maxEta, "<E/p>", minEop_l, maxEop_l);
    m_eop_highTileA_200_vs_trkPhi = book(m_name, std::string("eop_highTileA_"+m_energyCalib+"_0_200_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "<E/p>", minEop_l, maxEop_l);
    m_eop_highTileA_200_vs_mu = book(m_name, std::string("eop_highTileA_"+m_energyCalib+"_0_200_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "<E/p>", minEop_l, maxEop_l);
    m_eop_highTileA_200_vs_mu_avg = book(m_name, std::string("eop_highTileA_"+m_energyCalib+"_0_200_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "<E/p>", minEop_l, maxEop_l);
    m_eop_highTileA_200_vs_npv = book(m_name, std::string("eop_highTileA_"+m_energyCalib+"_0_200_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "<E/p>", minEop_l, maxEop_l);
    m_eop_highTileA_200_vs_TileA_E = book(m_name, std::string("eop_highTileA_"+m_energyCalib+"_0_200_vs_TileA_E"), "E(Tile, layer A)", nBinsE, minE, maxE, "<E/p>", minEop_l, maxEop_l);
    if (m_doEtabinsArray && m_doPbinsArray) { 
      if (m_doProfileEta) {
        // TProfile binned in p (one profile for each eta bin)
        m_eop_highTileA_200_profileEta = std::vector<TProfile* >();
        char buffer [200];
        for (unsigned int i=0; i<nEtabinsArray; i++) {
          std::snprintf(buffer, 200, "etaG%dL%d", (int)std::round((EtabinsArray[i]*10)), (int)std::round((EtabinsArray[i+1]*10)));
          TProfile* tmp_profile = book(m_name, std::string(std::string(buffer)+"_eop_highTileA_"+m_energyCalib+"_0_200"), "p_{trk}", nPbinsArray, &PbinsArray[0], "<E/p>", minEop_l, maxEop_l);
          m_eop_highTileA_200_profileEta.push_back(tmp_profile);
        }
      }
      if (m_doProfileP) {
        // TProfile binned in eta (one profile for each p bin)
        m_eop_highTileA_200_profileP = std::vector<TProfile* >();
        char buffer [200];
        for (unsigned int i=0; i<nPbinsArray; i++) {
          std::snprintf(buffer, 200, "pG%dL%d", (int)std::round((PbinsArray[i]*1000)), (int)std::round((PbinsArray[i+1]*1000)));
          TProfile* tmp_profile = book(m_name, std::string(std::string(buffer)+"_eop_highTileA_"+m_energyCalib+"_0_200"), "eta_{trk}", nEtabinsArray, &EtabinsArray[0], "<E/p>", minEop_l, maxEop_l);
          m_eop_highTileA_200_profileP.push_back(tmp_profile);
        }
      }
      if (m_doExtraEtaEnergyBinHists) {
        // 1D histograms for E/p distributions
        m_eop_highTileA_200_EtaEnergyRanges = std::vector<std::vector<TH1D*> >(nPbinsArray);
        char buffer [200];
        for (unsigned int i=0; i<nPbinsArray; i++) {
          for (unsigned int j=0; j<nEtabinsArray; j++) {
            std::snprintf(buffer, 200, "pG%dL%d_etaG%dL%d", (int)std::round((PbinsArray[i]*1000)), (int)std::round((PbinsArray[i+1]*1000)), (int)std::round((EtabinsArray[j]*10)), (int)std::round((EtabinsArray[j+1]*10)));
            TH1D* tmp_hist = book(m_name, std::string(std::string(buffer)+"_eop_highTileA_"+m_energyCalib+"_0_200"), "E/p", nBinsEop, minEop, maxEop);
            m_eop_highTileA_200_EtaEnergyRanges[i].push_back(tmp_hist);
          }
        }
      }
    }
    // max in Tile Layer B
    m_trk_E_highTileB_200 = book(m_name, std::string("trk_E_highTileB_"+m_energyCalib+"_0_200"), "E", nBinsE, minE, maxE);
    m_eop_highTileB_200 = book(m_name, std::string("eop_highTileB_"+m_energyCalib+"_0_200"), "E/p", nBinsEop, minEop, maxEop);
    if (m_doPbinsArray) m_eop_highTileB_200_vs_trkP = book(m_name, std::string("eop_highTileB_"+m_energyCalib+"_0_200_vs_trkP"), "p_{trk}", nPbinsArray, &PbinsArray[0], "<E/p>", minEop_l, maxEop_l);
    else m_eop_highTileB_200_vs_trkP = book(m_name, std::string("eop_highTileB_"+m_energyCalib+"_0_200_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "<E/p>", minEop_l, maxEop_l);
    if(m_doEtabinsArray) m_eop_highTileB_200_vs_trkEta = book(m_name, std::string("eop_highTileB_"+m_energyCalib+"_0_200_vs_trkEta"), "|#eta_{trk}|", nEtabinsArray, &EtabinsArray[0], "<E/p>", minEop_l, maxEop_l);
    else m_eop_highTileB_200_vs_trkEta = book(m_name, std::string("eop_highTileB_"+m_energyCalib+"_0_200_vs_trkEta"), "|#eta_{trk}|", nBinsEta, minEta, maxEta, "<E/p>", minEop_l, maxEop_l);
    m_eop_highTileB_200_vs_trkPhi = book(m_name, std::string("eop_highTileB_"+m_energyCalib+"_0_200_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "<E/p>", minEop_l, maxEop_l);
    m_eop_highTileB_200_vs_mu = book(m_name, std::string("eop_highTileB_"+m_energyCalib+"_0_200_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "<E/p>", minEop_l, maxEop_l);
    m_eop_highTileB_200_vs_mu_avg = book(m_name, std::string("eop_highTileB_"+m_energyCalib+"_0_200_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "<E/p>", minEop_l, maxEop_l);
    m_eop_highTileB_200_vs_npv = book(m_name, std::string("eop_highTileB_"+m_energyCalib+"_0_200_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "<E/p>", minEop_l, maxEop_l);
    m_eop_highTileB_200_vs_TileB_E = book(m_name, std::string("eop_highTileB_"+m_energyCalib+"_0_200_vs_TileB_E"), "E(Tile, layer B)", nBinsE, minE, maxE, "<E/p>", minEop_l, maxEop_l);
    if (m_doEtabinsArray && m_doPbinsArray) { 
      if (m_doProfileEta) {
        // TProfile binned in p (one profile for each eta bin)
        m_eop_highTileB_200_profileEta = std::vector<TProfile* >();
        char buffer [200];
        for (unsigned int i=0; i<nEtabinsArray; i++) {
          std::snprintf(buffer, 200, "etaG%dL%d", (int)std::round((EtabinsArray[i]*10)), (int)std::round((EtabinsArray[i+1]*10)));
          TProfile* tmp_profile = book(m_name, std::string(std::string(buffer)+"_eop_highTileB_"+m_energyCalib+"_0_200"), "p_{trk}", nPbinsArray, &PbinsArray[0], "<E/p>", minEop_l, maxEop_l);
          m_eop_highTileB_200_profileEta.push_back(tmp_profile);
        }
      }
      if (m_doProfileP) {
        // TProfile binned in eta (one profile for each p bin)
        m_eop_highTileB_200_profileP = std::vector<TProfile* >();
        char buffer [200];
        for (unsigned int i=0; i<nPbinsArray; i++) {
          std::snprintf(buffer, 200, "pG%dL%d", (int)std::round((PbinsArray[i]*1000)), (int)std::round((PbinsArray[i+1]*1000)));
          TProfile* tmp_profile = book(m_name, std::string(std::string(buffer)+"_eop_highTileB_"+m_energyCalib+"_0_200"), "eta_{trk}", nEtabinsArray, &EtabinsArray[0], "<E/p>", minEop_l, maxEop_l);
          m_eop_highTileB_200_profileP.push_back(tmp_profile);
        }
      }
      if (m_doExtraEtaEnergyBinHists) {
        // 1D histograms for E/p distributions
        m_eop_highTileB_200_EtaEnergyRanges = std::vector<std::vector<TH1D*> >(nPbinsArray);
        char buffer [200];
        for (unsigned int i=0; i<nPbinsArray; i++) {
          for (unsigned int j=0; j<nEtabinsArray; j++) {
            std::snprintf(buffer, 200, "pG%dL%d_etaG%dL%d", (int)std::round((PbinsArray[i]*1000)), (int)std::round((PbinsArray[i+1]*1000)), (int)std::round((EtabinsArray[j]*10)), (int)std::round((EtabinsArray[j+1]*10)));
            TH1D* tmp_hist = book(m_name, std::string(std::string(buffer)+"_eop_highTileB_"+m_energyCalib+"_0_200"), "E/p", nBinsEop, minEop, maxEop);
            m_eop_highTileB_200_EtaEnergyRanges[i].push_back(tmp_hist);
          }
        }
      }
    }
    // max in Tile Layer D
    m_trk_E_highTileD_200 = book(m_name, std::string("trk_E_highTileD_"+m_energyCalib+"_0_200"), "E", nBinsE, minE, maxE);
    m_eop_highTileD_200 = book(m_name, std::string("eop_highTileD_"+m_energyCalib+"_0_200"), "E/p", nBinsEop, minEop, maxEop);
    if (m_doPbinsArray) m_eop_highTileD_200_vs_trkP = book(m_name, std::string("eop_highTileD_"+m_energyCalib+"_0_200_vs_trkP"), "p_{trk}", nPbinsArray, &PbinsArray[0], "<E/p>", minEop_l, maxEop_l);
    else m_eop_highTileD_200_vs_trkP = book(m_name, std::string("eop_highTileD_"+m_energyCalib+"_0_200_vs_trkP"), "p_{trk}", nBinsE, minE, maxE, "<E/p>", minEop_l, maxEop_l);
    if(m_doEtabinsArray) m_eop_highTileD_200_vs_trkEta = book(m_name, std::string("eop_highTileD_"+m_energyCalib+"_0_200_vs_trkEta"), "|#eta_{trk}|", nEtabinsArray, &EtabinsArray[0], "<E/p>", minEop_l, maxEop_l);
    else m_eop_highTileD_200_vs_trkEta = book(m_name, std::string("eop_highTileD_"+m_energyCalib+"_0_200_vs_trkEta"), "|#eta_{trk}|", nBinsEta, minEta, maxEta, "<E/p>", minEop_l, maxEop_l);
    m_eop_highTileD_200_vs_trkPhi = book(m_name, std::string("eop_highTileD_"+m_energyCalib+"_0_200_vs_trkPhi"), "#phi_{trk}", nBinsPhi, minPhi, maxPhi, "<E/p>", minEop_l, maxEop_l);
    m_eop_highTileD_200_vs_mu = book(m_name, std::string("eop_highTileD_"+m_energyCalib+"_0_200_vs_mu"), "#mu", nBinsMu, minMu, maxMu, "<E/p>", minEop_l, maxEop_l);
    m_eop_highTileD_200_vs_mu_avg = book(m_name, std::string("eop_highTileD_"+m_energyCalib+"_0_200_vs_mu_avg"), "<#mu>", nBinsMu, minMu, maxMu, "<E/p>", minEop_l, maxEop_l);
    m_eop_highTileD_200_vs_npv = book(m_name, std::string("eop_highTileD_"+m_energyCalib+"_0_200_vs_npv"), "NPV", nBinsNPV, minNPV, maxNPV, "<E/p>", minEop_l, maxEop_l);
    m_eop_highTileD_200_vs_TileD_E = book(m_name, std::string("eop_highTileD_"+m_energyCalib+"_0_200_vs_TileD_E"), "E(Tile, layer D)", nBinsE, minE, maxE, "<E/p>", minEop_l, maxEop_l);
    if (m_doEtabinsArray && m_doPbinsArray) { 
      if (m_doProfileEta) {
        // TProfile binned in p (one profile for each eta bin)
        m_eop_highTileD_200_profileEta = std::vector<TProfile* >();
        char buffer [200];
        for (unsigned int i=0; i<nEtabinsArray; i++) {
          std::snprintf(buffer, 200, "etaG%dL%d", (int)std::round((EtabinsArray[i]*10)), (int)std::round((EtabinsArray[i+1]*10)));
          TProfile* tmp_profile = book(m_name, std::string(std::string(buffer)+"_eop_highTileD_"+m_energyCalib+"_0_200"), "p_{trk}", nPbinsArray, &PbinsArray[0], "<E/p>", minEop_l, maxEop_l);
          m_eop_highTileD_200_profileEta.push_back(tmp_profile);
        }
      }
      if (m_doProfileP) {
        // TProfile binned in eta (one profile for each p bin)
        m_eop_highTileD_200_profileP = std::vector<TProfile* >();
        char buffer [200];
        for (unsigned int i=0; i<nPbinsArray; i++) {
          std::snprintf(buffer, 200, "pG%dL%d", (int)std::round((PbinsArray[i]*1000)), (int)std::round((PbinsArray[i+1]*10)));
          TProfile* tmp_profile = book(m_name, std::string(std::string(buffer)+"_eop_highTileD_"+m_energyCalib+"_0_200"), "eta_{trk}", nEtabinsArray, &EtabinsArray[0], "<E/p>", minEop_l, maxEop_l);
          m_eop_highTileD_200_profileP.push_back(tmp_profile);
        }
      }
      if (m_doExtraEtaEnergyBinHists) {
        // 1D histograms for E/p distributions
        m_eop_highTileD_200_EtaEnergyRanges = std::vector<std::vector<TH1D*> >(nPbinsArray);
        char buffer [200];
        for (unsigned int i=0; i<nPbinsArray; i++) {
          for (unsigned int j=0; j<nEtabinsArray; j++) {
            std::snprintf(buffer, 200, "pG%dL%d_etaG%dL%d", (int)std::round((PbinsArray[i]*1000)), (int)std::round((PbinsArray[i+1]*1000)), (int)std::round((EtabinsArray[j]*10)), (int)std::round((EtabinsArray[j+1]*10)));
            TH1D* tmp_hist = book(m_name, std::string(std::string(buffer)+"_eop_highTileD_"+m_energyCalib+"_0_200"), "E/p", nBinsEop, minEop, maxEop);
            m_eop_highTileD_200_EtaEnergyRanges[i].push_back(tmp_hist);
          }
        }
      }
    }
  } // END doTileLayer

  // if worker is passed to the class add histograms to the output
  return StatusCode::SUCCESS;
}

StatusCode EoverPHists::execute( const xAOD::TrackParticle* trk, const xAOD::VertexContainer *vtxs, const xAOD::EventInfo* eventInfo, double eventWeight )
{

  // get pileup
  double mu(-1.);
  if( eventInfo->isAvailable< float >( "actualInteractionsPerCrossing" ) ) {
    mu = eventInfo->actualInteractionsPerCrossing();
  }
  double mu_avg(-1.);
  if( eventInfo->isAvailable< float >( "corrected_averageInteractionsPerCrossing" ) &&
      !eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION )) 
    mu_avg = eventInfo->auxdata< float >( "corrected_averageInteractionsPerCrossing" );
  else if( eventInfo->isAvailable< float >( "averageInteractionsPerCrossing" ) )
    mu_avg = eventInfo->averageInteractionsPerCrossing();

  // get number of primary vtxs 
  double npv = HelperFunctions::countPrimaryVertices(vtxs, 2);

  double trk_p = 0;
  if (TMath::Abs(trk->qOverP())>0.) trk_p = (1./TMath::Abs(trk->qOverP()))/1e3; 
  double trk_pt = trk->pt()/1e3;
  // coordinates of the track in the ID
  double trk_etaID = trk->eta();
  double trk_phiID = trk->phi();
  // coordinates of the track extrapolated to the specified calorimeter layer
  double trk_etaEMB2 = trk->auxdata<float>("CALO_trkEta_EMB2");
  double trk_phiEMB2 = trk->auxdata<float>("CALO_trkPhi_EMB2");
  double trk_etaEME2 = trk->auxdata<float>("CALO_trkEta_EME2");
  double trk_phiEME2 = trk->auxdata<float>("CALO_trkPhi_EME2");
  double dEta_EMB2_ID = TMath::Abs(trk_etaEMB2 - trk_etaID);
  double dPhi_EMB2_ID = TMath::Abs(trk_phiEMB2 - trk_phiID);
  if (dPhi_EMB2_ID > TMath::Pi())
    dPhi_EMB2_ID = 2*TMath::Pi() - dPhi_EMB2_ID;
  double dR_EMB2_ID = sqrt( pow(dEta_EMB2_ID, 2) + pow(dPhi_EMB2_ID, 2) );

  if (m_doEtaAbs) {
    trk_etaID = TMath::Abs(trk_etaID);
    trk_etaEMB2 = TMath::Abs(trk_etaEMB2);
  }

  // fill histos with basic track properties
  m_trk_p -> Fill(trk_p, eventWeight); 
  if (m_doPbinsArray) m_trk_p_array -> Fill(trk_p, eventWeight); 
  m_trk_pt -> Fill(trk_pt, eventWeight); 
  if (m_doPbinsArray) m_trk_pt_array -> Fill(trk_pt, eventWeight); 
  m_trk_etaID -> Fill(trk_etaID, eventWeight); 
  m_trk_etaEMB2 -> Fill(trk_etaEMB2, eventWeight); 
  m_trk_etaEME2 -> Fill(trk_etaEME2, eventWeight); 
  if (m_doEtabinsArray) m_trk_eta_array -> Fill(trk_etaID, eventWeight); 
  m_trk_phiID -> Fill(trk_phiID, eventWeight); 
  m_trk_phiEMB2 -> Fill(trk_phiEMB2, eventWeight); 
  m_trk_phiEME2 -> Fill(trk_phiEME2, eventWeight); 
  m_trk_phi_extra -> Fill(trk_phiID, eventWeight); 
  m_trk_phi_extra2 -> Fill(trk_phiID, eventWeight); 

  m_trk_p_vs_eta -> Fill(trk_etaID, trk_p, eventWeight); 
  m_trk_p_vs_etaEMB2 -> Fill(trk_etaEMB2, trk_p, eventWeight); 
  m_trk_p_vs_etaEME2 -> Fill(trk_etaEME2, trk_p, eventWeight); 
  if (m_doPbinsArray && m_doEtabinsArray) m_trk_p_vs_eta_array -> Fill(trk_etaID, trk_p, eventWeight); 

  m_trk_DR_EMB2_ID -> Fill(dR_EMB2_ID, eventWeight); 
  m_trk_DEta_EMB2_ID -> Fill(dEta_EMB2_ID, eventWeight); 
  m_trk_DPhi_EMB2_ID -> Fill(dPhi_EMB2_ID, eventWeight); 
  m_trk_DR_EMB2_ID_vs_trk_p -> Fill(trk_p, dR_EMB2_ID, eventWeight);
  m_trk_DEta_EMB2_ID_vs_trk_p -> Fill(trk_p, dEta_EMB2_ID, eventWeight);
  m_trk_DPhi_EMB2_ID_vs_trk_p -> Fill(trk_p, dPhi_EMB2_ID, eventWeight);


  double chi2 = trk->chiSquared();
  double ndof = trk->numberDoF();
  double chi2Prob = TMath::Prob(chi2,ndof);
  double d0 = trk->d0();
  const xAOD::Vertex* pvx = HelperFunctions::getPrimaryVertex(vtxs);
  double pvz = HelperFunctions::getPrimaryVertexZ(pvx);
  double z0 = trk->z0() + trk->vz() - pvz;
  double sinT = sin(trk->theta());
  m_trk_d0 -> Fill( d0, eventWeight );
  m_trk_d0_s -> Fill( d0, eventWeight );
  m_trk_z0 -> Fill( z0, eventWeight );
  m_trk_z0_s -> Fill( z0, eventWeight );
  m_trk_z0sinT -> Fill(z0*sinT, eventWeight );
  m_trk_chi2Prob-> Fill( chi2Prob, eventWeight );
  m_trk_charge -> Fill( trk->charge(), eventWeight );

  uint8_t nBL       = -1;
  uint8_t nPix      = -1;
  uint8_t nPixDead  = -1;
  uint8_t nPixHoles = -1;
  uint8_t nSCT      = -1;
  uint8_t nSCTDead  = -1;
  uint8_t nTRT      = -1;

  if(!trk->summaryValue(nBL,       xAOD::numberOfBLayerHits))       Error("TrackHists::execute()", "BLayer hits not filled");
  if(!trk->summaryValue(nPix,      xAOD::numberOfPixelHits))        Error("TrackHists::execute()", "Pix hits not filled");
  if(!trk->summaryValue(nPixDead,  xAOD::numberOfPixelDeadSensors)) Error("TrackHists::execute()", "Pix Dead not filled");
  if(!trk->summaryValue(nPixHoles, xAOD::numberOfPixelHoles))       Error("TrackHists::execute()", "Pix holes not filled");
  if(!trk->summaryValue(nSCT,      xAOD::numberOfSCTHits))          Error("TrackHists::execute()", "SCT hits not filled");
  if(!trk->summaryValue(nSCTDead,  xAOD::numberOfSCTDeadSensors))   Error("TrackHists::execute()", "SCT Dead not filled");
  if(!trk->summaryValue(nTRT,      xAOD::numberOfTRTHits))          Error("TrackHists::execute()", "TRT hits not filled");

  uint8_t nSi     = nPix     + nSCT;
  uint8_t nSiDead = nPixDead + nSCTDead;
  m_trk_nBL         -> Fill( nBL          , eventWeight );
  m_trk_nSi         -> Fill( nSi          , eventWeight );
  m_trk_nSiAndDead  -> Fill( nSi+nSiDead  , eventWeight );
  m_trk_nSiDead     -> Fill( nSiDead      , eventWeight );
  m_trk_nSCT        -> Fill( nSCT         , eventWeight );
  m_trk_nPix        -> Fill( nPix         , eventWeight );
  m_trk_nPixHoles   -> Fill( nPixHoles    , eventWeight );
  m_trk_nTRT        -> Fill( nTRT         , eventWeight );
  m_trk_nTRT_vs_p   -> Fill( trk_p,     nTRT, eventWeight );
  m_trk_nTRT_vs_eta -> Fill( trk_etaID, nTRT, eventWeight );

  // cluster energy associated with the track
  double trk_E_Total_100_all = trk->auxdata<float>(std::string("CALO_Total_"+m_energyCalib+"_0_100"))/1e3; 
  double trk_E_Total_200_all = trk->auxdata<float>(std::string("CALO_Total_"+m_energyCalib+"_0_200"))/1e3; 
  double trk_E_EM_100 = trk->auxdata<float>(std::string("CALO_EM_"+m_energyCalib+"_0_100"))/1e3; 
  double trk_E_EM_200 = trk->auxdata<float>(std::string("CALO_EM_"+m_energyCalib+"_0_200"))/1e3; 
  double trk_E_HAD_100 = trk->auxdata<float>(std::string("CALO_HAD_"+m_energyCalib+"_0_100"))/1e3; 
  double trk_E_HAD_200 = trk->auxdata<float>(std::string("CALO_HAD_"+m_energyCalib+"_0_200"))/1e3; 

  // To agree with Millie
  double trk_E_Total_100 = trk_E_EM_100 + trk_E_HAD_100;
  double trk_E_Total_200 = trk_E_EM_200 + trk_E_HAD_200;

  // calculate energy fractions in different layers
  double trk_sumE_Tile_100 = 0; 
  double trk_sumE_Tile_200 = 0; 
  double trk_sumE_Lar_100 = 0; 
  double trk_sumE_Lar_200 = 0; 
  double trk_sumE_HAD_100 = 0; 
  double trk_sumE_HAD_200 = 0; 
  double trk_sumE_EM_100 = 0; 
  double trk_sumE_EM_200 = 0; 
  for (unsigned int i=0; i<m_layer_tile.size(); i++) {
    double trk_E_100_tmp = trk->auxdata<float>(std::string("CALO_"+m_energyCalib+"_"+m_layer_tile[i]+"_100"))/1e3; 
    if (trk_E_100_tmp > 0.) // only include E > 0 (i.e. not calorimeter noise)
      trk_sumE_Tile_100 += trk_E_100_tmp; 
    double trk_E_200_tmp = trk->auxdata<float>(std::string("CALO_"+m_energyCalib+"_"+m_layer_tile[i]+"_200"))/1e3; 
    if (trk_E_200_tmp > 0.)
      trk_sumE_Tile_200 += trk_E_200_tmp; 
  }
  for (unsigned int i=0; i<m_layer_lar.size(); i++) {
    trk_sumE_Lar_100 += trk->auxdata<float>(std::string("CALO_"+m_energyCalib+"_"+m_layer_lar[i]+"_100"))/1e3; 
    trk_sumE_Lar_200 += trk->auxdata<float>(std::string("CALO_"+m_energyCalib+"_"+m_layer_lar[i]+"_200"))/1e3; 
  }
  for (unsigned int i=0; i<m_layer_had.size(); i++) {
    trk_sumE_HAD_100 += trk->auxdata<float>(std::string("CALO_"+m_energyCalib+"_"+m_layer_had[i]+"_100"))/1e3; 
    trk_sumE_HAD_200 += trk->auxdata<float>(std::string("CALO_"+m_energyCalib+"_"+m_layer_had[i]+"_200"))/1e3; 
  }
  for (unsigned int i=0; i<m_layer_em.size(); i++) {
    trk_sumE_EM_100 += trk->auxdata<float>(std::string("CALO_"+m_energyCalib+"_"+m_layer_em[i]+"_100"))/1e3; 
    trk_sumE_EM_200 += trk->auxdata<float>(std::string("CALO_"+m_energyCalib+"_"+m_layer_em[i]+"_200"))/1e3; 
  }

  // basic "sanity" checks of the E/p xAOD derivation
  m_trk_sumE_Tile_100 -> Fill(trk_sumE_Tile_100, eventWeight);
  m_trk_sumE_Tile_200 -> Fill(trk_sumE_Tile_200, eventWeight);
  m_trk_sumE_Lar_100 -> Fill(trk_sumE_Lar_100, eventWeight);
  m_trk_sumE_Lar_200 -> Fill(trk_sumE_Lar_200, eventWeight);
  m_trk_SumTileLayers_over_HAD_100 -> Fill(trk_sumE_Tile_100/trk_E_HAD_100, eventWeight);
  m_trk_SumTileLayers_over_HAD_200 -> Fill(trk_sumE_Tile_200/trk_E_HAD_200, eventWeight);
  m_trk_SumLarLayers_over_EM_100 -> Fill(trk_sumE_Lar_100/trk_E_EM_100, eventWeight);
  m_trk_SumLarLayers_over_EM_200 -> Fill(trk_sumE_Lar_200/trk_E_EM_200, eventWeight);
  m_trk_SumHADLayers_over_HAD_100 -> Fill(trk_sumE_HAD_100/trk_E_HAD_100, eventWeight);
  m_trk_SumHADLayers_over_HAD_200 -> Fill(trk_sumE_HAD_200/trk_E_HAD_200, eventWeight);
  m_trk_SumEMLayers_over_EM_100 -> Fill(trk_sumE_EM_100/trk_E_EM_100, eventWeight);
  m_trk_SumEMLayers_over_EM_200 -> Fill(trk_sumE_EM_200/trk_E_EM_200, eventWeight);
  m_trk_EMandHAD_over_Total_100 -> Fill((trk_E_EM_100+trk_E_HAD_100)/trk_E_Total_100_all, eventWeight);
  m_trk_EMandHAD_over_Total_200 -> Fill((trk_E_EM_200+trk_E_HAD_200)/trk_E_Total_200_all, eventWeight);

  // determine which layer has the highest energy deposit 
  int trk_highE_100_layer = 0;
  int trk_highE_200_layer = 0;
  double trk_highE_100 = -1e8;
  double trk_highE_200 = -1e8;
  double trk_E_100_tmp = -1e8;
  double trk_E_200_tmp = -1e8;
  double trk_sumE_Total_100 = 0.;
  double trk_sumE_Total_200 = 0.;
  double trk_sumEg0_Total_100 = 0.;
  double trk_sumEg0_Total_200 = 0.;
  for (unsigned int i=0; i<m_layer.size(); i++) { 
    trk_E_100_tmp = trk->auxdata<float>(std::string("CALO_"+m_energyCalib+"_"+m_layer[i]+"_100"))/1e3; 
    trk_E_200_tmp = trk->auxdata<float>(std::string("CALO_"+m_energyCalib+"_"+m_layer[i]+"_200"))/1e3; 
    trk_sumE_Total_100 += trk_E_100_tmp; 
    trk_sumE_Total_200 += trk_E_200_tmp; 
    if (trk_E_100_tmp > 0.) 
      trk_sumEg0_Total_100 += trk_E_100_tmp; 
    if (trk_E_200_tmp > 0.) 
      trk_sumEg0_Total_200 += trk_E_200_tmp; 
    m_trk_E_100_vs_layer -> Fill(i, trk_E_100_tmp, eventWeight); 
    m_trk_E_200_vs_layer -> Fill(i, trk_E_200_tmp, eventWeight); 

    if (trk_E_100_tmp > trk_highE_100) {
      trk_highE_100 = trk_E_100_tmp;
      trk_highE_100_layer = i; 
    }
    if (trk_E_200_tmp > trk_highE_200) {
      trk_highE_200 = trk_E_200_tmp;
      trk_highE_200_layer = i; 
    }
  }

  double trk_TileEfrac_100 = trk_sumE_Tile_100/trk_sumEg0_Total_100;    
  double trk_TileEfrac_200 = trk_sumE_Tile_200/trk_sumEg0_Total_200;    

  int trk_p_i = -1;
  int trk_eta_i = -1;
  if (m_doEtabinsArray && m_doPbinsArray) {
    for (unsigned int i=0; i<nPbinsArray; i++) {
      if (trk_p >= PbinsArray[i] && trk_p < PbinsArray[i+1]) 
        trk_p_i = i;
    }
    for (unsigned int i=0; i<nEtabinsArray; i++) {
      if (trk_etaID >= EtabinsArray[i] && trk_etaID < EtabinsArray[i+1]) 
        trk_eta_i = i; 
    }
  }

  // fill histos with Tile energy fractions
  m_trk_TileEfrac_100 -> Fill(trk_TileEfrac_100, eventWeight);
  m_trk_TileEfrac_100_vs_trk_p -> Fill(trk_p, trk_TileEfrac_100, eventWeight);
  m_trk_TileEfrac_200 -> Fill(trk_TileEfrac_200, eventWeight);
  m_trk_TileEfrac_200_vs_trk_p -> Fill(trk_p, trk_TileEfrac_200, eventWeight);
  // "sanity checks" check so that the energy in all clusters sum to the Total
  m_trk_SumAllLayers_over_Total_100 -> Fill(trk_sumE_Total_100/trk_E_Total_100_all, eventWeight);
  m_trk_SumAllLayers_over_Total_200 -> Fill(trk_sumE_Total_200/trk_E_Total_200_all, eventWeight);

  // fill histos of highest energy layer
  m_trk_E_100_highElayer -> Fill(trk_highE_100_layer, eventWeight);
  m_trk_E_200_highElayer -> Fill(trk_highE_200_layer, eventWeight);
  m_trk_E_100_highElayer_vs_E -> Fill(trk_highE_100, trk_highE_100_layer, eventWeight);
  m_trk_E_200_highElayer_vs_E -> Fill(trk_highE_200, trk_highE_200_layer, eventWeight);

  // fill E/p histograms
  // Total calorimeter energy
  if (m_doCaloTotal) {
    // dR(trk,cluster) < 0.1
    m_trk_E_Total_100 -> Fill(trk_E_Total_100, eventWeight); 
    m_trk_E_Total_100_vs_mu_avg -> Fill(mu_avg, trk_E_Total_100, eventWeight); 
    m_eop_Total_100 -> Fill(trk_E_Total_100/trk_p, eventWeight); 
    m_eop_Total_100_vs_trkP -> Fill(trk_p, trk_E_Total_100/trk_p, eventWeight); 
    m_eop_Total_100_vs_trkEta -> Fill(trk_etaID, trk_E_Total_100/trk_p, eventWeight); 
    m_eop_Total_100_vs_trkPhi -> Fill(trk_phiID, trk_E_Total_100/trk_p, eventWeight); 
    m_eop_Total_100_vs_trkPhiID -> Fill(trk_phiID, trk_E_Total_100/trk_p, eventWeight); 
    m_eop_Total_100_vs_trkPhi_extra -> Fill(trk_phiID, trk_E_Total_100/trk_p, eventWeight); 
    m_eop_Total_100_vs_trkPhi_extra2 -> Fill(trk_phiID, trk_E_Total_100/trk_p, eventWeight); 
    m_eop_Total_100_vs_mu -> Fill(mu, trk_E_Total_100/trk_p, eventWeight); 
    m_eop_Total_100_vs_mu_avg -> Fill(mu_avg, trk_E_Total_100/trk_p, eventWeight); 
    m_eop_Total_100_vs_npv -> Fill(npv, trk_E_Total_100/trk_p, eventWeight); 
    m_eop_Total_100_vs_highE_layer -> Fill(trk_highE_100_layer, trk_E_Total_100/trk_p, eventWeight); 
    if (m_doProfileEta && trk_eta_i >= 0)
      m_eop_Total_100_profileEta[trk_eta_i] -> Fill(trk_p, trk_E_Total_100/trk_p, eventWeight); 
    if (m_doProfileP && trk_p_i >= 0)
      m_eop_Total_100_profileP[trk_p_i] -> Fill(trk_etaID, trk_E_Total_100/trk_p, eventWeight); 
    if (m_doExtraEtaEnergyBinHists && trk_p_i >= 0 && trk_eta_i >= 0)
      m_eop_Total_100_EtaEnergyRanges[trk_p_i][trk_eta_i] -> Fill(trk_E_Total_100/trk_p, eventWeight); 
    // dR(trk,cluster) < 0.2

    // fill zero fraction histograms
    if (m_doPbinsArray && m_doEtabinsArray) {
      m_trk_n_E_200 -> Fill(trk_p, trk_etaID, eventWeight);
      if (trk_E_Total_200 < 0) m_trk_n_E_200_l0 -> Fill(trk_p, trk_etaID, eventWeight);
      if (trk_E_Total_200 == 0) m_trk_n_E_200_eq0 -> Fill(trk_p, trk_etaID, eventWeight);
      if (trk_E_Total_200 <= 0) m_trk_n_E_200_leq0 -> Fill(trk_p, trk_etaID, eventWeight);
    }

    m_trk_E_Total_200 -> Fill(trk_E_Total_200, eventWeight); 
    m_trk_E_Total_200_vs_mu_avg -> Fill(mu_avg, trk_E_Total_200, eventWeight); 
    m_eop_Total_200 -> Fill(trk_E_Total_200/trk_p, eventWeight); 
    m_eop_Total_200_l -> Fill(trk_E_Total_200/trk_p, eventWeight); 
    m_eop_Total_200_vs_trkP -> Fill(trk_p, trk_E_Total_200/trk_p, eventWeight); 
    m_eop_Total_200_vs_trkEta -> Fill(trk_etaID, trk_E_Total_200/trk_p, eventWeight); 
    m_eop_Total_200_vs_trkPhi -> Fill(trk_phiID, trk_E_Total_200/trk_p, eventWeight); 
    m_eop_Total_200_vs_trkPhiID -> Fill(trk_phiID, trk_E_Total_200/trk_p, eventWeight); 
    m_eop_Total_200_vs_trkPhi_extra -> Fill(trk_phiID, trk_E_Total_200/trk_p, eventWeight); 
    m_eop_Total_200_vs_trkPhi_extra2 -> Fill(trk_phiID, trk_E_Total_200/trk_p, eventWeight); 
    m_eop_Total_200_vs_mu -> Fill(mu, trk_E_Total_200/trk_p, eventWeight); 
    m_eop_Total_200_vs_mu_avg -> Fill(mu_avg, trk_E_Total_200/trk_p, eventWeight); 
    m_eop_Total_200_vs_npv -> Fill(npv, trk_E_Total_200/trk_p, eventWeight); 
    m_eop_Total_200_vs_highE_layer -> Fill(trk_highE_200_layer, trk_E_Total_200/trk_p, eventWeight); 
    if (m_doProfileEta && trk_eta_i >= 0)
      m_eop_Total_200_profileEta[trk_eta_i] -> Fill(trk_p, trk_E_Total_200/trk_p, eventWeight); 
    if (m_doProfileP && trk_p_i >= 0)
      m_eop_Total_200_profileP[trk_p_i] -> Fill(trk_etaID, trk_E_Total_200/trk_p, eventWeight); 
    if (m_doExtraEtaEnergyBinHists && trk_p_i >= 0 && trk_eta_i >= 0)
      m_eop_Total_200_EtaEnergyRanges[trk_p_i][trk_eta_i] -> Fill(trk_E_Total_200/trk_p, eventWeight); 
  } // END doCaloTotal

  // EM calorimeter (EMB+EMEC)
  if (m_doCaloEM) {
    // dR(trk,cluster) < 0.1
    m_trk_E_EM_100 -> Fill(trk_E_EM_100, eventWeight); 
    m_eop_EM_100 -> Fill(trk_E_EM_100/trk_p, eventWeight); 
    m_eop_EM_100_vs_trkP -> Fill(trk_p, trk_E_EM_100/trk_p, eventWeight); 
    m_eop_EM_100_vs_trkEta -> Fill(trk_etaID, trk_E_EM_100/trk_p, eventWeight); 
    m_eop_EM_100_vs_trkPhi -> Fill(trk_phiID, trk_E_EM_100/trk_p, eventWeight); 
    m_eop_EM_100_vs_mu -> Fill(mu, trk_E_EM_100/trk_p, eventWeight); 
    m_eop_EM_100_vs_mu_avg -> Fill(mu_avg, trk_E_EM_100/trk_p, eventWeight); 
    m_eop_EM_100_vs_npv -> Fill(npv, trk_E_EM_100/trk_p, eventWeight); 
    m_eop_EM_100_vs_highE_layer -> Fill(trk_highE_100_layer, trk_E_EM_100/trk_p, eventWeight); 
    if (m_doProfileEta && trk_eta_i >= 0)
      m_eop_EM_100_profileEta[trk_eta_i] -> Fill(trk_p, trk_E_EM_100/trk_p, eventWeight); 
    if (m_doProfileP && trk_p_i >= 0)
      m_eop_EM_100_profileP[trk_p_i] -> Fill(trk_etaID, trk_E_EM_100/trk_p, eventWeight); 
    if (m_doExtraEtaEnergyBinHists && trk_p_i >= 0 && trk_eta_i >= 0)
      m_eop_EM_100_EtaEnergyRanges[trk_p_i][trk_eta_i] -> Fill(trk_E_EM_100/trk_p, eventWeight); 
    // dR(trk,cluster) < 0.2
    m_trk_E_EM_200 -> Fill(trk_E_EM_200, eventWeight); 
    m_eop_EM_200 -> Fill(trk_E_EM_200/trk_p, eventWeight); 
    m_eop_EM_200_vs_trkP -> Fill(trk_p, trk_E_EM_200/trk_p, eventWeight); 
    m_eop_EM_200_vs_trkEta -> Fill(trk_etaID, trk_E_EM_200/trk_p, eventWeight); 
    m_eop_EM_200_vs_trkPhi -> Fill(trk_phiID, trk_E_EM_200/trk_p, eventWeight); 
    m_eop_EM_200_vs_mu -> Fill(mu, trk_E_EM_200/trk_p, eventWeight); 
    m_eop_EM_200_vs_mu_avg -> Fill(mu_avg, trk_E_EM_200/trk_p, eventWeight); 
    m_eop_EM_200_vs_npv -> Fill(npv, trk_E_EM_200/trk_p, eventWeight); 
    m_eop_EM_200_vs_highE_layer -> Fill(trk_highE_200_layer, trk_E_EM_200/trk_p, eventWeight); 
    if (m_doProfileEta && trk_eta_i >= 0)
      m_eop_EM_200_profileEta[trk_eta_i] -> Fill(trk_p, trk_E_EM_200/trk_p, eventWeight); 
    if (m_doProfileP && trk_p_i >= 0)
      m_eop_EM_200_profileP[trk_p_i] -> Fill(trk_etaID, trk_E_EM_200/trk_p, eventWeight); 
    if (m_doExtraEtaEnergyBinHists && trk_p_i >= 0 && trk_eta_i >= 0)
      m_eop_EM_200_EtaEnergyRanges[trk_p_i][trk_eta_i] -> Fill(trk_E_EM_200/trk_p, eventWeight); 

    // No clusters in HAD matched to trk 
    if (trk_E_HAD_200 <=0) {
      m_eop_Total_200_noHAD -> Fill(trk_E_Total_200/trk_p, eventWeight);
      m_eop_EM_200_noHAD -> Fill(trk_E_EM_200/trk_p, eventWeight);
      m_eop_HAD_200_noHAD -> Fill(trk_E_HAD_200/trk_p, eventWeight);
    }

  } // END doCaloEM

  // HAD calorimeter (HEC+TileBarrel+TileGap+TileExtBarrel)
  if (m_doCaloHAD) {
    // dR(trk,cluster) < 0.1
    m_trk_E_HAD_100 -> Fill(trk_E_HAD_100, eventWeight); 
    m_eop_HAD_100 -> Fill(trk_E_HAD_100/trk_p, eventWeight); 
    m_eop_HAD_100_vs_trkP -> Fill(trk_p, trk_E_HAD_100/trk_p, eventWeight); 
    m_eop_HAD_100_vs_trkEta -> Fill(trk_etaID, trk_E_HAD_100/trk_p, eventWeight); 
    m_eop_HAD_100_vs_trkPhi -> Fill(trk_phiID, trk_E_HAD_100/trk_p, eventWeight); 
    m_eop_HAD_100_vs_mu -> Fill(mu, trk_E_HAD_100/trk_p, eventWeight); 
    m_eop_HAD_100_vs_mu_avg -> Fill(mu_avg, trk_E_HAD_100/trk_p, eventWeight); 
    m_eop_HAD_100_vs_npv -> Fill(npv, trk_E_HAD_100/trk_p, eventWeight); 
    m_eop_HAD_100_vs_highE_layer -> Fill(trk_highE_100_layer, trk_E_HAD_100/trk_p, eventWeight); 
    if (m_doProfileEta && trk_eta_i >= 0)
      m_eop_HAD_100_profileEta[trk_eta_i] -> Fill(trk_p, trk_E_HAD_100/trk_p, eventWeight); 
    if (m_doProfileP && trk_p_i >= 0)
      m_eop_HAD_100_profileP[trk_p_i] -> Fill(trk_etaID, trk_E_HAD_100/trk_p, eventWeight); 
    if (m_doExtraEtaEnergyBinHists && trk_p_i >= 0 && trk_eta_i >= 0)
      m_eop_HAD_100_EtaEnergyRanges[trk_p_i][trk_eta_i] -> Fill(trk_E_HAD_100/trk_p, eventWeight); 

    // dR(trk,cluster) < 0.2
    m_trk_E_HAD_200 -> Fill(trk_E_HAD_200, eventWeight); 
    m_eop_HAD_200 -> Fill(trk_E_HAD_200/trk_p, eventWeight); 
    m_eop_HAD_200_vs_trkP -> Fill(trk_p, trk_E_HAD_200/trk_p, eventWeight); 
    m_eop_HAD_200_vs_trkEta -> Fill(trk_etaID, trk_E_HAD_200/trk_p, eventWeight); 
    m_eop_HAD_200_vs_trkPhi -> Fill(trk_phiID, trk_E_HAD_200/trk_p, eventWeight); 
    m_eop_HAD_200_vs_mu -> Fill(mu, trk_E_HAD_200/trk_p, eventWeight); 
    m_eop_HAD_200_vs_mu_avg -> Fill(mu_avg, trk_E_HAD_200/trk_p, eventWeight); 
    m_eop_HAD_200_vs_npv -> Fill(npv, trk_E_HAD_200/trk_p, eventWeight); 
    m_eop_HAD_200_vs_highE_layer -> Fill(trk_highE_200_layer, trk_E_HAD_200/trk_p, eventWeight); 
    if (m_doProfileEta && trk_eta_i >= 0)
      m_eop_HAD_200_profileEta[trk_eta_i] -> Fill(trk_p, trk_E_HAD_200/trk_p, eventWeight); 
    if (m_doProfileP && trk_p_i >= 0)
      m_eop_HAD_200_profileP[trk_p_i] -> Fill(trk_etaID, trk_E_HAD_200/trk_p, eventWeight); 
    if (m_doExtraEtaEnergyBinHists && trk_p_i >= 0 && trk_eta_i >= 0)
      m_eop_HAD_200_EtaEnergyRanges[trk_p_i][trk_eta_i] -> Fill(trk_E_HAD_200/trk_p, eventWeight); 

    // MIP requirement for HAD calorimeter
    if (trk_E_EM_100 < 1.1 && trk_E_HAD_100/trk_p > 0.4 && trk_E_HAD_100/trk_p < 0.9) {
      m_eop_Total_200_MIP -> Fill(trk_E_Total_200/trk_p, eventWeight);
      m_eop_EM_200_MIP -> Fill(trk_E_EM_200/trk_p, eventWeight);
      m_eop_HAD_200_MIP -> Fill(trk_E_HAD_200/trk_p, eventWeight);
    }

  } // END doCaloHAD

  // background subtraction, the way it was done in run 1
  if (m_doBgSubtr) {
    if (trk_E_EM_100 < 1.1 && trk_E_HAD_100/trk_p > 0.4 && trk_E_HAD_100/trk_p < 0.9) {
      m_eop_EM_BG -> Fill( (trk_E_EM_200 - trk_E_EM_100)/trk_p, eventWeight); 
      m_eop_EM_BG_vs_trkP -> Fill(trk_p, (trk_E_EM_200 - trk_E_EM_100)/trk_p, eventWeight); 
      m_eop_EM_BG_vs_trkEta -> Fill(trk_etaID, (trk_E_EM_200 - trk_E_EM_100)/trk_p, eventWeight); 
      m_eop_EM_BG_vs_trkPhi -> Fill(trk_phiID, (trk_E_EM_200 - trk_E_EM_100)/trk_p, eventWeight); 
      m_eop_EM_BG_vs_mu -> Fill(mu, (trk_E_EM_200 - trk_E_EM_100)/trk_p, eventWeight); 
      m_eop_EM_BG_vs_mu_avg -> Fill(mu_avg, (trk_E_EM_200 - trk_E_EM_100)/trk_p, eventWeight); 
      m_eop_EM_BG_vs_npv -> Fill(npv, (trk_E_EM_200 - trk_E_EM_100)/trk_p, eventWeight); 
      if (m_doProfileEta && trk_eta_i >= 0)
        m_eop_EM_BG_profileEta[trk_eta_i] -> Fill(trk_p, (trk_E_EM_200 - trk_E_EM_100)/trk_p, eventWeight); 
      if (m_doProfileP && trk_p_i >= 0)
        m_eop_EM_BG_profileP[trk_p_i] -> Fill(trk_etaID, (trk_E_EM_200 - trk_E_EM_100)/trk_p, eventWeight); 
      if (m_doExtraEtaEnergyBinHists && trk_p_i >= 0 && trk_eta_i >= 0)
        m_eop_EM_BG_EtaEnergyRanges[trk_p_i][trk_eta_i] -> Fill((trk_E_EM_200 - trk_E_EM_100)/trk_p, eventWeight); 
    }
  } // END doBgSubtr

  // E/p where E is maximum in either Tile layer A, BC, or D
  if (m_doTileLayer) {
    double highTileLayer_100 = -1e8;
    double highTileLayer_200 = -1e8;
    int highTileLayer_i_100 = -1; 
    int highTileLayer_i_200 = -1; 
    if (trk_highE_100_layer == 12 || trk_highE_100_layer == 18) {
      highTileLayer_100 = trk_highE_100;             
      highTileLayer_i_100 = 1;                         
    } if (trk_highE_100_layer == 13 || trk_highE_100_layer == 19) {
      highTileLayer_100 = trk_highE_100;             
      highTileLayer_i_100 = 2;                         
    } if (trk_highE_100_layer == 14 || trk_highE_100_layer == 20) {
      highTileLayer_100 = trk_highE_100;             
      highTileLayer_i_100 = 3;                         
    } if (trk_highE_200_layer == 12 || trk_highE_200_layer == 18) {
      highTileLayer_200 = trk_highE_200;             
      highTileLayer_i_200 = 1;                         
    } if (trk_highE_200_layer == 13 || trk_highE_200_layer == 19) {
      highTileLayer_200 = trk_highE_200;             
      highTileLayer_i_200 = 2;                         
    } if (trk_highE_200_layer == 14 || trk_highE_200_layer == 20) {
      highTileLayer_200 = trk_highE_200;
      highTileLayer_i_200 = 3;
    }
    // dR(trk,cluster) < 0.1
    // max in Tile Layer A
    if (highTileLayer_i_100 == 1) { 
      m_trk_E_highTileA_100 -> Fill(trk_E_Total_100, eventWeight); 
      m_eop_highTileA_100 -> Fill(trk_E_Total_100/trk_p, eventWeight); 
      m_eop_highTileA_100_vs_trkP -> Fill(trk_p, trk_E_Total_100/trk_p, eventWeight); 
      m_eop_highTileA_100_vs_trkEta -> Fill(trk_etaID, trk_E_Total_100/trk_p, eventWeight); 
      m_eop_highTileA_100_vs_trkPhi -> Fill(trk_phiID, trk_E_Total_100/trk_p, eventWeight); 
      m_eop_highTileA_100_vs_mu -> Fill(mu, trk_E_Total_100/trk_p, eventWeight); 
      m_eop_highTileA_100_vs_mu_avg -> Fill(mu_avg, trk_E_Total_100/trk_p, eventWeight); 
      m_eop_highTileA_100_vs_npv -> Fill(npv, trk_E_Total_100/trk_p, eventWeight); 
      m_eop_highTileA_100_vs_TileA_E -> Fill(highTileLayer_100, trk_E_Total_100/trk_p, eventWeight); 
      if (m_doProfileEta && trk_eta_i >= 0)
        m_eop_highTileA_100_profileEta[trk_eta_i] -> Fill(trk_p, trk_E_Total_100/trk_p, eventWeight); 
      if (m_doProfileP && trk_p_i >= 0)
        m_eop_highTileA_100_profileP[trk_p_i] -> Fill(trk_etaID, trk_E_Total_100/trk_p, eventWeight); 
      if (m_doExtraEtaEnergyBinHists && trk_p_i >= 0 && trk_eta_i >= 0)
        m_eop_highTileA_100_EtaEnergyRanges[trk_p_i][trk_eta_i] -> Fill(trk_E_Total_100/trk_p, eventWeight); 
    }
    // max in Tile Layer BC
    if (highTileLayer_i_100 == 2) { 
      m_trk_E_highTileB_100 -> Fill(trk_E_Total_100, eventWeight); 
      m_eop_highTileB_100 -> Fill(trk_E_Total_100/trk_p, eventWeight); 
      m_eop_highTileB_100_vs_trkP -> Fill(trk_p, trk_E_Total_100/trk_p, eventWeight); 
      m_eop_highTileB_100_vs_trkEta -> Fill(trk_etaID, trk_E_Total_100/trk_p, eventWeight); 
      m_eop_highTileB_100_vs_trkPhi -> Fill(trk_phiID, trk_E_Total_100/trk_p, eventWeight); 
      m_eop_highTileB_100_vs_mu -> Fill(mu, trk_E_Total_100/trk_p, eventWeight); 
      m_eop_highTileB_100_vs_mu_avg -> Fill(mu_avg, trk_E_Total_100/trk_p, eventWeight); 
      m_eop_highTileB_100_vs_npv -> Fill(npv, trk_E_Total_100/trk_p, eventWeight); 
      m_eop_highTileB_100_vs_TileB_E -> Fill(highTileLayer_100, trk_E_Total_100/trk_p, eventWeight); 
      if (m_doProfileEta && trk_eta_i >= 0)
        m_eop_highTileB_100_profileEta[trk_eta_i] -> Fill(trk_p, trk_E_Total_100/trk_p, eventWeight); 
      if (m_doProfileP && trk_p_i >= 0)
        m_eop_highTileB_100_profileP[trk_p_i] -> Fill(trk_etaID, trk_E_Total_100/trk_p, eventWeight); 
      if (m_doExtraEtaEnergyBinHists && trk_p_i >= 0 && trk_eta_i >= 0)
        m_eop_highTileB_100_EtaEnergyRanges[trk_p_i][trk_eta_i] -> Fill(trk_E_Total_100/trk_p, eventWeight); 
    }
    // max in Tile Layer D
    if (highTileLayer_i_100 == 3) { 
      m_trk_E_highTileD_100 -> Fill(trk_E_Total_100, eventWeight); 
      m_eop_highTileD_100 -> Fill(trk_E_Total_100/trk_p, eventWeight); 
      m_eop_highTileD_100_vs_trkP -> Fill(trk_p, trk_E_Total_100/trk_p, eventWeight); 
      m_eop_highTileD_100_vs_trkEta -> Fill(trk_etaID, trk_E_Total_100/trk_p, eventWeight); 
      m_eop_highTileD_100_vs_trkPhi -> Fill(trk_phiID, trk_E_Total_100/trk_p, eventWeight); 
      m_eop_highTileD_100_vs_mu -> Fill(mu, trk_E_Total_100/trk_p, eventWeight); 
      m_eop_highTileD_100_vs_mu_avg -> Fill(mu_avg, trk_E_Total_100/trk_p, eventWeight); 
      m_eop_highTileD_100_vs_npv -> Fill(npv, trk_E_Total_100/trk_p, eventWeight); 
      m_eop_highTileD_100_vs_TileD_E -> Fill(highTileLayer_100, trk_E_Total_100/trk_p, eventWeight); 
      if (m_doProfileEta && trk_eta_i >= 0)
        m_eop_highTileD_100_profileEta[trk_eta_i] -> Fill(trk_p, trk_E_Total_100/trk_p, eventWeight); 
      if (m_doProfileP && trk_p_i >= 0)
        m_eop_highTileD_100_profileP[trk_p_i] -> Fill(trk_etaID, trk_E_Total_100/trk_p, eventWeight); 
      if (m_doExtraEtaEnergyBinHists && trk_p_i >= 0 && trk_eta_i >= 0)
        m_eop_highTileD_100_EtaEnergyRanges[trk_p_i][trk_eta_i] -> Fill(trk_E_Total_100/trk_p, eventWeight); 
    }
    // dR(trk,cluster) < 0.2
    // max in Tile Layer A
    if (highTileLayer_i_200 == 1) { 
      m_trk_E_highTileA_200 -> Fill(trk_E_Total_200, eventWeight); 
      m_eop_highTileA_200 -> Fill(trk_E_Total_200/trk_p, eventWeight); 
      m_eop_highTileA_200_vs_trkP -> Fill(trk_p, trk_E_Total_200/trk_p, eventWeight); 
      m_eop_highTileA_200_vs_trkEta -> Fill(trk_etaID, trk_E_Total_200/trk_p, eventWeight); 
      m_eop_highTileA_200_vs_trkPhi -> Fill(trk_phiID, trk_E_Total_200/trk_p, eventWeight); 
      m_eop_highTileA_200_vs_mu -> Fill(mu, trk_E_Total_200/trk_p, eventWeight); 
      m_eop_highTileA_200_vs_mu_avg -> Fill(mu_avg, trk_E_Total_200/trk_p, eventWeight); 
      m_eop_highTileA_200_vs_npv -> Fill(npv, trk_E_Total_200/trk_p, eventWeight); 
      m_eop_highTileA_200_vs_TileA_E -> Fill(highTileLayer_200, trk_E_Total_200/trk_p, eventWeight); 
      if (m_doProfileEta && trk_eta_i >= 0)
        m_eop_highTileA_200_profileEta[trk_eta_i] -> Fill(trk_p, trk_E_Total_200/trk_p, eventWeight); 
      if (m_doProfileP && trk_p_i >= 0)
        m_eop_highTileA_200_profileP[trk_p_i] -> Fill(trk_etaID, trk_E_Total_200/trk_p, eventWeight); 
      if (m_doExtraEtaEnergyBinHists && trk_p_i >= 0 && trk_eta_i >= 0)
        m_eop_highTileA_200_EtaEnergyRanges[trk_p_i][trk_eta_i] -> Fill(trk_E_Total_200/trk_p, eventWeight); 
    }
    // max in Tile Layer BC
    if (highTileLayer_i_200 == 2) { 
      m_trk_E_highTileB_200 -> Fill(trk_E_Total_200, eventWeight); 
      m_eop_highTileB_200 -> Fill(trk_E_Total_200/trk_p, eventWeight); 
      m_eop_highTileB_200_vs_trkP -> Fill(trk_p, trk_E_Total_200/trk_p, eventWeight); 
      m_eop_highTileB_200_vs_trkEta -> Fill(trk_etaID, trk_E_Total_200/trk_p, eventWeight); 
      m_eop_highTileB_200_vs_trkPhi -> Fill(trk_phiID, trk_E_Total_200/trk_p, eventWeight); 
      m_eop_highTileB_200_vs_mu -> Fill(mu, trk_E_Total_200/trk_p, eventWeight); 
      m_eop_highTileB_200_vs_mu_avg -> Fill(mu_avg, trk_E_Total_200/trk_p, eventWeight); 
      m_eop_highTileB_200_vs_npv -> Fill(npv, trk_E_Total_200/trk_p, eventWeight); 
      m_eop_highTileB_200_vs_TileB_E -> Fill(highTileLayer_200, trk_E_Total_200/trk_p, eventWeight); 
      if (m_doProfileEta && trk_eta_i >= 0)
        m_eop_highTileB_200_profileEta[trk_eta_i] -> Fill(trk_p, trk_E_Total_200/trk_p, eventWeight); 
      if (m_doProfileP && trk_p_i >= 0)
        m_eop_highTileB_200_profileP[trk_p_i] -> Fill(trk_etaID, trk_E_Total_200/trk_p, eventWeight); 
      if (m_doExtraEtaEnergyBinHists && trk_p_i >= 0 && trk_eta_i >= 0)
        m_eop_highTileB_200_EtaEnergyRanges[trk_p_i][trk_eta_i] -> Fill(trk_E_Total_200/trk_p, eventWeight); 
    }
    // max in Tile Layer D
    if (highTileLayer_i_200 == 3) { 
      m_trk_E_highTileD_200 -> Fill(trk_E_Total_200, eventWeight); 
      m_eop_highTileD_200 -> Fill(trk_E_Total_200/trk_p, eventWeight); 
      m_eop_highTileD_200_vs_trkP -> Fill(trk_p, trk_E_Total_200/trk_p, eventWeight); 
      m_eop_highTileD_200_vs_trkEta -> Fill(trk_etaID, trk_E_Total_200/trk_p, eventWeight); 
      m_eop_highTileD_200_vs_trkPhi -> Fill(trk_phiID, trk_E_Total_200/trk_p, eventWeight); 
      m_eop_highTileD_200_vs_mu -> Fill(mu, trk_E_Total_200/trk_p, eventWeight); 
      m_eop_highTileD_200_vs_mu_avg -> Fill(mu_avg, trk_E_Total_200/trk_p, eventWeight); 
      m_eop_highTileD_200_vs_npv -> Fill(npv, trk_E_Total_200/trk_p, eventWeight); 
      m_eop_highTileD_200_vs_TileD_E -> Fill(highTileLayer_200, trk_E_Total_200/trk_p, eventWeight); 
      if (m_doProfileEta && trk_eta_i >= 0)
        m_eop_highTileD_200_profileEta[trk_eta_i] -> Fill(trk_p, trk_E_Total_200/trk_p, eventWeight); 
      if (m_doProfileP && trk_p_i >= 0)
        m_eop_highTileD_200_profileP[trk_p_i] -> Fill(trk_etaID, trk_E_Total_200/trk_p, eventWeight); 
      if (m_doExtraEtaEnergyBinHists && trk_p_i >= 0 && trk_eta_i >= 0)
        m_eop_highTileD_200_EtaEnergyRanges[trk_p_i][trk_eta_i] -> Fill(trk_E_Total_200/trk_p, eventWeight); 
    }
  } // END doTileLayer

  return StatusCode::SUCCESS;
}

std::vector<double> EoverPHists::str2vec(std::string str)
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
