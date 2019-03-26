#ifndef eoverpanalysis_EnergySumHelper_H
#define eoverpanalysis_EnergySumHelper_H

#include <map>
#include <xAODTracking/TrackParticleContainer.h>

namespace EnergySumHelper
{
    const std::vector<std::string> layer = {"PreSamplerB",
					    "PreSamplerE", 
					    "EMB1", 
					    "EMB2", 
					    "EMB3", 
					    "EME1", 
					    "EME2", 
					    "EME3", 
					    "HEC0", 
					    "HEC1", 
					    "HEC2", 
					    "HEC3", 
					    "TileBar0", 
					    "TileBar1", 
					    "TileBar2", 
					    "TileGap1",
					    "TileGap2", 
					    "TileGap3", 
					    "TileExt0",
					    "TileExt1",
					    "TileExt2"}; //! array of all the calo layers
    
    const std::vector<std::string> layer_EM = {"PreSamplerB",
					       "PreSamplerE", 
					       "EMB1",
					       "EMB2",
					       "EMB3", 
					       "EME1",
					       "EME2",
					       "EME3"};

    const std::vector<std::string> layer_HAD = {"TileBar0",
						"TileBar1", 
						"TileBar2",
						"TileGap1",
						"TileGap2",
						"TileGap3",
						"TileExt0",
						"TileExt1",
						"TileExt2", 
						"HEC0",
						"HEC1",
						"HEC2",
						"HEC3"}; //! array of HAD layers only

    const std::map<std::string, unsigned int> layer_to_id = { {"PreSamplerB", 0},
							      {"PreSamplerE", 1}, 
							      {"EMB1", 2},
							      {"EMB2", 3},
							      {"EMB3", 4},
							      {"EME1", 5},
							      {"EME2", 6}, 
							      {"EME3", 7},
							      {"HEC0", 8},
							      {"HEC1", 9},
							      {"HEC2", 10},
							      {"HEC3", 11},
							      {"TileBar0", 12},
							      {"TileBar1", 13},
							      {"TileBar2", 14},
							      {"TileGap1", 15},
							      {"TileGap2", 16},
							      {"TileGap3", 17},
							      {"TileExt0", 18},
							      {"TileExt1", 19},
							      {"TileExt2", 20}};

    const std::map<unsigned int, std::string> id_to_layer = { {0, "PreSamplerB"},
							      {1, "PreSamplerE"},
							      {2, "EMB1"},
							      {3, "EMB2"},
							      {4, "EMB3"},
							      {5, "EME1"},
							      {6, "EME2"},
							      {7, "EME3"},
							      {8, "HEC0"}, 
							      {9, "HEC1"}, 
							      {10, "HEC2"},
							      {11, "HEC3"},
							      {12, "TileBar0"},
							      {13, "TileBar1"},
							      {14, "TileBar2"},
							      {15, "TileGap1"},
							      {16, "TileGap2"},
							      {17, "TileGap3"},
							      {18, "TileExt0"}, 
							      {19, "TileExt1"}, 
							      {20, "TileExt2"}};

    extern std::map<unsigned int, float> getEnergySumInLayers(const xAOD::TrackParticle* trk, 
							      std::string variableToSum,
							      std::string radiusCut);
    
    extern std::map<std::string, float> getEnergySumInCalorimeterRegions(const xAOD::TrackParticle* trk,
									 std::string variableToSum, 
									 std::string radiusCut,
									 bool onlyPositiveEnergy = false);
}

#endif
