// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Markus Apel, Nora Heese $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/SILACDetector.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/KERNEL/DPeak.h>
#include <OpenMS/KERNEL/MSExperiment.h>
//#include <map>
//#include <cmath>

namespace OpenMS
{

  bool SILACDetector::detectSILAC(MSExperiment experiment)
  {
    MSExperiment experimentMS2;
    for (int i = 0; i < experiment.size(); i++)
    {
      if (2 == experiment[i].getMSLevel())
      {
        experimentMS2.addSpectrum(experiment[i]);
      }
    }
    //int spectrum_number = 0;
    std::vector<int> control_distances = {11, 14, 15, 21, 23, 27};
    std::vector<int> silac_distances = {4, 6, 8, 10};

    std::map<int,int> distance_count = {{4,0},{6,0},{8,0},{10,0},{11,0},{14,0},{15,0},{21,0},{23,0},{27,0}};
    double RT_window = 5;
    int min_index = 0;
    int max_index = 0;
    
    
    // vielleicht von anfang an nur MS2 spektren nehmen

     for (const auto& spectrum : experimentMS2)
    {
       
      //std::cout << "Spectrumnumber: " << spectrum_number << std::endl;
      //spectrum_number++;

      double spectrum_rt = spectrum.getRT();
      while (experimentMS2[min_index].getRT() < spectrum.getRT() - RT_window)
      {
        min_index++;
      } 
      while ((experimentMS2[max_index].getRT() < spectrum.getRT() + RT_window) && max_index < experimentMS2.size())
      {
        max_index++;
      } 

      double spectrum_charge = spectrum.getPrecursors()[0].getCharge();
      for (int i = min_index; i < max_index; i++)
      {
        double experimentMS2_charge = experimentMS2[i].getPrecursors()[0].getCharge();
        if (spectrum_charge == experimentMS2_charge)
        {
          double distance = std::abs((spectrum.getPrecursors()[0].getMZ() - experimentMS2[i].getPrecursors()[0].getMZ()) * spectrum_charge);
          distance += 0.5;
          int rounded_distance = distance;
          //std::cout << distance << std::endl;
          if (distance_count.find(rounded_distance) != distance_count.end())
          {
            distance_count[rounded_distance]++;
          }
  
        }
          
      }
      
    }
    double n = 6;

    double control_mean = (distance_count[11] + distance_count[14] + distance_count[15] + distance_count[21] + distance_count[23] + distance_count[27]) / 6;
    double control_sd = 0;
    double s;
    double sd_sum=0;
    for (int i = 0; i < 6; i++)
    {
      s = distance_count[control_distances[i]]-control_mean;
      s *= s;
      sd_sum += s;

    }
    control_sd = std::sqrt(sd_sum/(n-1));

    bool is_silac = false;
    double z_score;

    for (int i = 4; i <= 10; i += 2)
    {
      z_score = (distance_count[i]-control_mean)/control_sd;
      std::cout << "Distanz " << i << ": Z-score: " << z_score << std::endl;
      is_silac = is_silac || z_score > 2.5; // 2.5 is the cutoff for significance
    }
    
    /*
    std::cout << "Distanz 11: " << distance_count[11] << " Zscore: " << (distance_count[11]-control_mean)/control_sd<<std::endl;
    std::cout << "Distanz 14: " << distance_count[14] << " Zscore: " << (distance_count[14]-control_mean)/control_sd<<std::endl;
    std::cout << "Distanz 15: " << distance_count[15] << " Zscore: " << (distance_count[15]-control_mean)/control_sd<<std::endl;
    std::cout << "Distanz 21: " << distance_count[21] << " Zscore: " << (distance_count[21]-control_mean)/control_sd<<std::endl;
    std::cout << "Distanz 23: " << distance_count[23] << " Zscore: " << (distance_count[23]-control_mean)/control_sd<<std::endl;
    std::cout << "Distanz 27: " << distance_count[27] << " Zscore: " << (distance_count[27]-control_mean)/control_sd<<std::endl;
    */
    return is_silac;

  }
} // namespace OpenMS