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
#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/MATH/STATISTICS/MultipleTesting.h>

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
    std::vector<int> control_distances = {11, 14, 15, 21, 23, 27};
    std::vector<int> silac_distances = {4, 6, 8, 10};

    std::map<int,int> distance_count = {{4,0},{6,0},{8,0},{10,0},{11,0},{14,0},{15,0},{21,0},{23,0},{27,0}};
    double RT_window = 5; // Frage: RT_window = 5 gut oder meh?
    int min_index = 0;
    int max_index = 0;

    for (const auto& spectrum : experimentMS2)
    {
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
          int rounded_distance = distance; // runden auf Int zu grob?
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
    double significance_level = 0.025; // cut-off
    std::vector<bool> significant_distances = {};
    std::vector<String> aminoacids = {"Medium Lysine", "Heavy Lysine(K6) or Medium Arginine", "Heavy Lysine(K8)", "Heavy Arginine"};
    std::vector<double> p_values;

    const double sqrt2 = std::sqrt(2.0);

    for (int i = 4; i <= 10; i += 2)
    {
      z_score = (distance_count[i]-control_mean)/control_sd;
      double tail = 0.5 * std::erfc(z_score / sqrt2);
      p_values.push_back(tail);
      std::cout << "Distance " << i << ": Z-score: " << z_score << " - pValue: " << tail << std::endl;
      //is_silac = is_silac || z_score > 2.5; // 2.5 is the cutoff for significance
      bool is_significant = tail < significance_level;
      is_silac = is_silac || is_significant;
      significant_distances.push_back(is_significant);
    }
    std::cout << std::endl;
        
    if (is_silac)
    {
      std::cout << "Dataset is a SILAC dataset with the following aminoacids:" << std::endl;
      for (auto i = 0; i < significant_distances.size(); i++)
      {
        if (significant_distances[i])
        {
          std::cout << aminoacids[i] << " with p-value of " << std::scientific << p_values[i] << std::endl;
        }
      }
    }
    else
    {
      std::cout << "Dataset unable to be detected as a SILAC dataset" << std::endl;
    }
    return is_silac;
  }
} // namespace OpenMS