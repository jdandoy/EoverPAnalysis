// E/p analysis for run 2
// Joakim Olsson (joakim.olsson@cern.ch)

#ifndef EoverPAnalysis_EoverPTreeAlgo_H
#define EoverPAnalysis_EoverPTreeAlgo_H

#include "TTree.h"

#include "xAODTracking/VertexContainer.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODTruth/TruthParticle.h"


// algorithm wrapper
#include "xAODAnaHelpers/Algorithm.h"

class EoverPTreeAlgo : public xAH::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
  public:
    std::string m_inTrackContainerName;

    // configuration variables
    std::string m_detailStr;

    // save output tree of key variables
    bool m_fillOutputTree = false;
    // track isolation settings
    // turn on cutflows
    bool m_useCutFlow = false; 

    // energy calibration, either "ClusterEnergy", "ClusterEnergyLCW", or "CellEnergy"
    std::string m_energyCalib = "ClusterEnergy";
    // pileup reweighting
    bool m_doCustomPUreweighting = false;
    std::string m_pileupReweightingFile = "pileup_reweighting.root";
    // pt reweighting
    bool m_doTrkPtReweighting = false;
    std::string m_trkPtReweightingFile = "pt_reweighting.root";

  private:
    // variables to dump to the output ttree
    float trk_etaID;
    float trk_phiID;

    uint8_t m_trkIndex; 
    unsigned long long m_eventNumber; 

    uint8_t trk_nTRT;
    int trk_truthPdgId;

    float trk_truthEnergy;
    float trk_truthP;
    float trk_truthProb;

    TLorentzVector truthPartVec;

    float trk_pt;
    float trk_p;

    float trk_d0;
    float trk_z0sintheta;

    float trk_etaEMB2;
    float trk_phiEMB2;
    
    float trk_etaEME2;
    float trk_phiEME2;

    float trk_nearest_dR;

    float trkWeight;

    float trk_sumEPos_EM_100;
    float trk_sumEPos_EM_200;
    float trk_sumE_EM_200;
    float trk_sumE_EM_100;

    float  trk_sumEPos_HAD_200;
    float  trk_sumEPos_HAD_100;
    float  trk_sumE_HAD_200;
    float  trk_sumE_HAD_100;

    float  trk_sumEPos_Total_200;
    float  trk_sumEPos_Total_100;
    float  trk_sumE_Total_200;
    float  trk_sumE_Total_100;

    float  trk_HADEfrac_200;
    float  trk_HADEfrac_100;

    float trk_E_EM_100;
    float trk_E_EM_nopresampler_100;
    float trk_E_EM_200;
    float trk_E_EM_nopresampler_200;
    float trk_E_HAD_100; 
    float trk_E_HAD_200;
    float trk_E_Total_nopresampler_100;
    float trk_E_Total_nopresampler_200;

    uint8_t trk_NPV_2;
    uint8_t trk_NPV_4;
    char trk_charge;
    char trk_truthFromPileup;
    char trk_truthIsFake;
    float trk_actualmu;
    float trk_averagemu;
    float trk_corrected_averagemu;

    // cutflow
    TH1D* m_cutflowHist; //!
    TH1D* m_cutflowHistW; //!
    int   m_cutflow_bin; //!

    /* object-level cutflow */

    TH1D* m_trk_cutflowHist_1;  //!
    // list of calo layers
    const std::vector<std::string> m_layer = {"PreSamplerB","PreSamplerE", "EMB1", "EMB2", "EMB3", "EME1", "EME2", "EME3", "HEC0", "HEC1", "HEC2", "HEC3", "TileBar0", "TileBar1", "TileBar2", "TileGap1", "TileGap2", "TileGap3", "TileExt0", "TileExt1", "TileExt2"}; //! array of all the calo layers
    const std::vector<std::string> m_layer_EM = {"PreSamplerB","PreSamplerE", "EMB1", "EMB2", "EMB3", "EME1", "EME2", "EME3"};
    const std::vector<std::string> m_layer_HAD = {"TileBar0", "TileBar1", "TileBar2", "TileGap1", "TileGap2", "TileGap3", "TileExt0", "TileExt1", "TileExt2", "HEC0", "HEC1", "HEC2", "HEC3"}; //! array of HAD layers only

    // variables that don't get filled at submission time should be
    // protected from being send from the submission node to the worker
    // node (done by the //!)
    public:
    // Tree *myTree; //!
    // TH1 *myHist; //!

    // tree for saving some key variables
    TTree *m_tree; //!

    // this is a standard constructor
    EoverPTreeAlgo (std::string className = "EoverPTreeAlgo");

    // these are the functions inherited from Algorithm
    virtual EL::StatusCode setupJob (EL::Job& job);
    virtual EL::StatusCode fileExecute ();
    virtual EL::StatusCode histInitialize ();
    virtual EL::StatusCode changeInput (bool firstFile);
    virtual EL::StatusCode initialize ();
    virtual EL::StatusCode execute ();
    virtual EL::StatusCode postExecute ();
    virtual EL::StatusCode finalize ();
    virtual EL::StatusCode histFinalize ();
    const xAOD::TruthParticle* getTruthPtr(const xAOD::TrackParticle* trackParticle); 

    double deltaR (double trk_eta, double trk_phi, double trk2_eta, double trk2_phi);

    std::vector<double> str2vec(std::string str);

    /// @cond
    // this is needed to distribute the algorithm to the workers
    ClassDef(EoverPTreeAlgo, 1);
    /// @endcond
};

#endif
