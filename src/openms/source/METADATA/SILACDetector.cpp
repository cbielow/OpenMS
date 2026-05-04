// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Markus Apel, Nora Heese $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/KERNEL/DPeak.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/SILACDetector.h>

namespace OpenMS
{

  bool SILACDetector::detectSILAC(MSExperiment experiment)
  {
    if (experiment.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Can not detect a SILAC dataset from an empty MS experiment", "");
    }
    if (!experiment.containsScanOfLevel(2))
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Dataset does not contain any MS2 scans", "");
    }
    if (!experiment.isSorted())
    {
      experiment.sortSpectra();
    }

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
    int experiment_MS2_size = experimentMS2.size();

    for (const auto& spectrum : experimentMS2)
    {
      double spectrum_rt = spectrum.getRT();
      while (experimentMS2[min_index].getRT() < spectrum_rt - RT_window)
      {
        min_index++;
      } 
      while ((experimentMS2[max_index].getRT() < spectrum_rt + RT_window) && max_index < experiment_MS2_size)
      {
        max_index++;
      } 

      Precursor spectrum_precursor = spectrum.getPrecursors()[0];
      double spectrum_charge = spectrum_precursor.getCharge();
      for (int i = min_index; i < max_index; i++)
      {
        if (spectrum_charge == experimentMS2[i].getPrecursors()[0].getCharge())
        {
          double distance = std::abs((spectrum_precursor.getMZ() - experimentMS2[i].getPrecursors()[0].getMZ()) * spectrum_charge);
          distance += 0.5;
          int rounded_distance = distance; // runden auf Int zu grob?
          if (distance_count.find(rounded_distance) != distance_count.end())
          {
            distance_count[rounded_distance]++;
          }
        }
      }
    }
    double n = 6; // size of control distances

    double sum = 0.0;

    for (int i = 0; i < n; i++)
    {
      sum += distance_count[control_distances[i]];
    }

    double control_mean = sum / n;
    double control_sd = 0;
    double sd_sum = 0;
    for (int i = 0; i < 6; i++)
    {
      double s = distance_count[control_distances[i]]-control_mean;
      sd_sum += s * s;
    }
    control_sd = std::sqrt(sd_sum/(n-1));

    if (sd_sum == 0)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No counts for control distances found", 0);
    }

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
      bool is_significant = tail < significance_level;
      is_silac = is_silac || is_significant;
      significant_distances.push_back(is_significant);
    }
    std::cout << std::endl;
        
    if (is_silac)
    {
      std::cout << "Null hypothesis was rejected on significance level " << significance_level << std::endl;
      std::cout << "The following aminoacids have been detected:" << std::endl;
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
      std::cout << "Null hypothesis was not rejected" << std::endl;
    }
    return is_silac;
  }
} // namespace OpenMS