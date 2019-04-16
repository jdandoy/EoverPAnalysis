#include "EoverPAnalysis/EnergySumHelper.h"

std::map<unsigned int, float> EnergySumHelper::getEnergySumInLayers(const xAOD::TrackParticle* trk, 
								    std::string variableToSum, 
								    std::string radiusCut)
{
  std::map<unsigned int, float> energySum;

  for (std::string layerName: layer) 
    {
      float trk_E_tmp = trk->auxdata<float>(std::string("CALO_"+variableToSum+"_"+layerName+"_" +radiusCut))/1.e3; 
      unsigned int i = layer_to_id.at(layerName);
      energySum[i] += trk_E_tmp;
    }
  return energySum;
}

std::map<unsigned int, std::map<std::string, std::map<std::string, int > > > EnergySumHelper::getNumberOfClustersInLayers(const xAOD::TrackParticle* trk,
                                     const std::map<std::string,float> map_cut_names_to_values)
{
    std::map<unsigned int, std::map<std::string, std::map< std::string, int > > > to_return;
    for (std::string layerName: layer)
    {
        unsigned int i = layer_to_id.at(layerName);
        for(std::map<std::string, float>::const_iterator it = map_cut_names_to_values.begin(); it != map_cut_names_to_values.end(); ++it)
        {
            std::string cutName = it->first;
            to_return[i][cutName]["NClusters"] = 0;
            to_return[i][cutName]["NClusters_EMLike"] = 0;
            to_return[i][cutName]["NClusters_HADLike"] = 0;

        }
    }
    std::vector<float> cluster_dR_track = trk->auxdata<std::vector<float> >("CALO_ClusterEnergy_dRToTrack");
    std::vector<float> cluster_emProbability = trk->auxdata<std::vector<float> >("CALO_ClusterEnergy_emProbability");
    std::vector<float> cluster_Energy = trk->auxdata<std::vector<float> >("CALO_ClusterEnergy_Energy");
    std::vector<int> cluster_maxEnergyLayer = trk->auxdata<std::vector<int> >("CALO_ClusterEnergy_maxEnergyLayer");

    for (unsigned int clusn = 0; clusn < cluster_maxEnergyLayer.size(); clusn ++)
    {
        unsigned int layer = cluster_maxEnergyLayer.at(clusn);
        for(std::map<std::string,float>::const_iterator it = map_cut_names_to_values.begin(); it != map_cut_names_to_values.end(); ++it)
        {
            std::string cutName = it->first;
            float cutValue = it->second;
            if (cluster_dR_track.at(clusn) < cutValue){
                to_return[layer][cutName]["NClusters"] += 1;
                float em_probability = cluster_emProbability.at(clusn);
                if (em_probability > 0.5) to_return[layer][cutName]["NClusters_EMLike"] += 1;
                else to_return[layer][cutName]["NClusters_HADLike"] += 1;
            }
        }
    }
    return to_return;
}

std::map<std::string, std::map< std::string, std::map<std::string, int> > > EnergySumHelper::getNumberOfClustersInCaloriemterRegions(
                                     const xAOD::TrackParticle* trk,
                                     const std::map<std::string,float> map_cut_names_to_values)
{
    std::map<unsigned int, std::map<std::string, std::map<std::string, int > > > map_layer_to_cutName_to_variable_to_sum =
                                                  getNumberOfClustersInLayers(trk, map_cut_names_to_values);

    std::map<std::string, std::map< std::string, std::map<std::string, int> > > to_return;
    for (std::string layerName: layer_HAD) 
    {
        unsigned int layerId = layer_to_id.at(layerName);
        for(std::map<std::string, std::map<std::string, int > >::iterator it = map_layer_to_cutName_to_variable_to_sum[layerId].begin(); it != map_layer_to_cutName_to_variable_to_sum[layerId].end(); ++it)
        {
            std::string cutName = it->first;
            to_return["HAD"][cutName]["NClusters"] = 0;
            to_return["EM"][cutName]["NClusters"] = 0;

            to_return["HAD"][cutName]["NClusters_EMLike"] = 0;
            to_return["EM"][cutName]["NClusters_EMLike"] = 0;

            to_return["HAD"][cutName]["NClusters_HADLike"] = 0;
            to_return["EM"][cutName]["NClusters_HADLike"] = 0;
        }
    }

    for (std::string layerName: layer_HAD) 
    {
        unsigned int layerId = layer_to_id.at(layerName);
        for(std::map<std::string, std::map<std::string, int > >::iterator it = map_layer_to_cutName_to_variable_to_sum[layerId].begin(); it != map_layer_to_cutName_to_variable_to_sum[layerId].end(); ++it)
        {
            std::string cutName = it->first;

            to_return["HAD"][cutName]["NClusters"] += map_layer_to_cutName_to_variable_to_sum[layerId][cutName]["NClusters"];
            to_return["HAD"][cutName]["NClusters_EMLike"] += map_layer_to_cutName_to_variable_to_sum[layerId][cutName]["NClusters_EMLike"];
            to_return["HAD"][cutName]["NClusters_HADLike"] += map_layer_to_cutName_to_variable_to_sum[layerId][cutName]["NClusters_HADLike"];
        }
    }

    for (std::string layerName: layer_EM) 
    {
        unsigned int layerId = layer_to_id.at(layerName);
        for(std::map<std::string, std::map<std::string, int > >::iterator it = map_layer_to_cutName_to_variable_to_sum[layerId].begin(); it != map_layer_to_cutName_to_variable_to_sum[layerId].end(); ++it)
        {
            std::string cutName = it->first;
            to_return["EM"][cutName]["NClusters"] += map_layer_to_cutName_to_variable_to_sum[layerId][cutName]["NClusters"];
            to_return["EM"][cutName]["NClusters_EMLike"] += map_layer_to_cutName_to_variable_to_sum[layerId][cutName]["NClusters_EMLike"];
            to_return["EM"][cutName]["NClusters_HADLike"] += map_layer_to_cutName_to_variable_to_sum[layerId][cutName]["NClusters_HADLike"];
        }
    }
    return to_return;
}

std::map<std::string, float> EnergySumHelper::getEnergySumInCalorimeterRegions(const xAOD::TrackParticle* trk, 
									       std::string variableToSum, 
									       std::string radiusCut, 
									       bool onlyPositiveEnergy)
{
  std::map<unsigned int, float> energySum =  getEnergySumInLayers(trk, variableToSum, radiusCut);  
  std::map<std::string, float> energySumInCaloRegion;

  energySumInCaloRegion["HAD"] = 0.0;
  energySumInCaloRegion["EM"] = 0.0;

  for (std::string layerName: layer_HAD) 
    {
      float energy = energySum[layer_to_id.at(layerName)];
      if (onlyPositiveEnergy and (energy < 0.0))
	continue;
      energySumInCaloRegion["HAD"] += energy;
    }
  
  for (std::string layerName: layer_EM) 
    {
      float energy = energySum[layer_to_id.at(layerName)];
      if (onlyPositiveEnergy and (energy < 0.0))
	continue;
      energySumInCaloRegion["EM"] += energy;
    }

  return energySumInCaloRegion;
}
