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
    std::string m_energyCalibList = "ClusterEnergy,CellEnergy,ClusterEnergyLCW";
    // pileup reweighting
    bool m_doCustomPUreweighting = false;
    std::string m_pileupReweightingFile = "pileup_reweighting.root";
    // pt reweighting
    bool m_doTrkPtReweighting = false;
    std::string m_trkPtReweightingFile = "pt_reweighting.root";

  private:

    std::string ClusterEnergy = "ClusterEnergy";
    std::string ClusterEnergyLCW = "ClusterEnergyLCW";
    std::string CellEnergy = "CellEnergy";
    std::string m_energyCalib = "";
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
    char trk_hasTruthParticle;
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
    uint8_t trk_NPV_2;
    uint8_t trk_NPV_4;
    uint8_t trk_nclusters;
    uint8_t trk_nclusters_EM;
    uint8_t trk_nclusters_hadlike;
    uint8_t trk_nclusters_HAD;
    uint8_t trk_nclusters_emlike;
    char trk_charge;
    char trk_truthFromPileup;
    char trk_truthIsFake;
    float trk_actualmu;
    float trk_averagemu;
    float trk_corrected_averagemu;

    float trk_ClusterEnergy_Pos_EM_100;
    float trk_ClusterEnergy_Pos_EM_200;
    float trk_ClusterEnergy_EM_200;
    float trk_ClusterEnergy_EM_100;
    float  trk_ClusterEnergy_Pos_HAD_200;
    float  trk_ClusterEnergy_Pos_HAD_100;
    float  trk_ClusterEnergy_HAD_200;
    float  trk_ClusterEnergy_HAD_100;
    //float  trk_ClusterEnergy_Pos_Total_200;
    //float  trk_ClusterEnergy_Pos_Total_100;
    //float  trk_ClusterEnergy_Total_200;
    //float  trk_ClusterEnergy_Total_100;
    float trk_ClusterEnergy_EM_nopresampler_100;
    float trk_ClusterEnergy_EM_nopresampler_200;
    float trk_ClusterEnergy_Total_nopresampler_100;
    float trk_ClusterEnergy_Total_nopresampler_200;

    float trk_ClusterEnergyLCW_Pos_EM_100;
    float trk_ClusterEnergyLCW_Pos_EM_200;
    float trk_ClusterEnergyLCW_EM_200;
    float trk_ClusterEnergyLCW_EM_100;
    float trk_ClusterEnergyLCW_Pos_HAD_200;
    float trk_ClusterEnergyLCW_Pos_HAD_100;
    float trk_ClusterEnergyLCW_HAD_200;
    float trk_ClusterEnergyLCW_HAD_100;
    //float trk_ClusterEnergyLCW_Pos_Total_200;
    //float trk_ClusterEnergyLCW_Pos_Total_100;
    //float trk_ClusterEnergyLCW_Total_200;
    //float trk_ClusterEnergyLCW_Total_100;
    float trk_ClusterEnergyLCW_EM_nopresampler_100;
    float trk_ClusterEnergyLCW_EM_nopresampler_200;
    float trk_ClusterEnergyLCW_Total_nopresampler_100;
    float trk_ClusterEnergyLCW_Total_nopresampler_200;

    float trk_CellEnergy_Pos_EM_100;
    float trk_CellEnergy_Pos_EM_200;
    float trk_CellEnergy_EM_200;
    float trk_CellEnergy_EM_100;
    float trk_CellEnergy_Pos_HAD_200;
    float trk_CellEnergy_Pos_HAD_100;
    float trk_CellEnergy_HAD_200;
    float trk_CellEnergy_HAD_100;
    //float trk_CellEnergy_Pos_Total_200;
    //float trk_CellEnergy_Pos_Total_100;
    //float trk_CellEnergy_Total_200;
    //float trk_CellEnergy_Total_100;
    float trk_CellEnergy_EM_nopresampler_100;
    float trk_CellEnergy_EM_nopresampler_200;
    float trk_CellEnergy_Total_nopresampler_100;
    float trk_CellEnergy_Total_nopresampler_200;


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
    const std::map<std::string, int> m_layer_to_id = { {"PreSamplerB", 0}, {"PreSamplerE", 1}, {"EMB1", 2} , {"EMB2", 3}, {"EMB3",4}, {"EME1",5}, {"EME2",6}, {"EME3",7} , {"HEC0", 8}, {"HEC1", 9}, {"HEC2", 10}, {"HEC3", 11}, {"TileBar0", 12}, {"TileBar1", 13}, {"TileBar2", 14}, {"TileGap1", 15}, {"TileGap2", 16}, {"TileGap3", 17}, {"TileExt0", 18}, {"TileExt1", 19}, {"TileExt2", 20}};
    const std::map<int, std::string> m_id_to_layer = { {0, "PreSamplerB"}, {1, "PreSamplerE"}, {2, "EMB1"} , {3, "EMB2"}, {4, "EMB3"}, {5, "EME1"}, {6, "EME2"}, {7, "EME3"} , {8, "HEC0"}, {9, "HEC1"}, {10, "HEC2"}, {11, "HEC3"}, {12, "TileBar0"}, {13, "TileBar1"}, {14, "TileBar2"}, {15, "TileGap1"}, {16, "TileGap2"}, {17, "TileGap3"}, {18, "TileExt0"}, {19, "TileExt1"}, {20, "TileExt2"}};

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
