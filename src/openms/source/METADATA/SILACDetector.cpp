/// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
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

    std::vector<MS2Data> MS2experiments;
    int position = 0;
    for (unsigned int i = 0; i < experiment.size(); i++)
    {
      if (2 == experiment[i].getMSLevel())
      {
        position ++;
        MS2Data current_MS2_scan;
        current_MS2_scan.RT = experiment[i].getRT();
        current_MS2_scan.mz = experiment[i].getPrecursors()[0].getMZ();
        current_MS2_scan.charge = experiment[i].getPrecursors()[0].getCharge();
        current_MS2_scan.index = position;
        MS2experiments.push_back(current_MS2_scan);
      }
    }
    double exact_mass_lookup_table[28] = {0, 0, 0, 0, 4.025106983784, 0, 6.020129012016, 0, 8.014198800046, 0, 10.008268588075996, 11, 0, 0, 14, 15, 0 ,0 ,0 ,0 ,0 ,21 ,0 ,23 ,0 ,0 ,0 , 27};
    std::vector<int> control_distances = {11, 14, 15, 21, 23, 27};
    std::vector<int> silac_distances = {4, 6, 8, 10};

    std::map<int,int> distance_count = {{4,0},{6,0},{8,0},{10,0},{11,0},{14,0},{15,0},{21,0},{23,0},{27,0}};
    double RT_window = 5; // Frage: RT_window = 5 gut oder meh?
    int min_index = 0;
    int max_index = 0;
    int experiment_MS2_size = MS2experiments.size();

    for (const auto& spectrum : MS2experiments)
    {
      double spectrum_rt = spectrum.RT;
      min_index = spectrum.index + 1;
      /* while (MS2experiments[min_index].RT < spectrum_rt - RT_window)
      {
        min_index++;
      }  */
      while ((MS2experiments[max_index].RT < spectrum_rt + RT_window) && max_index < experiment_MS2_size)
      {
        max_index++;
      } 

      for (int i = min_index; i < max_index; i++)
      {
        if (spectrum.charge == MS2experiments[i].charge)
        {
          double distance = std::abs((spectrum.mz - MS2experiments[i].mz) * spectrum.charge);
          int rounded_distance = distance + 0.5; // runden auf Int zu grob?
          if (distance_count.find(rounded_distance) != distance_count.end())
          {
            if (distance - exact_mass_lookup_table[rounded_distance] < 0.000005) 
            {
              distance_count[rounded_distance]++;
            }      
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
    std::vector<double> p_values;
    std::vector<double> z_scores;
    const double sqrt2 = std::sqrt(2.0);

    for (int i = 4; i <= 10; i += 2)
    {
      z_score = (distance_count[i]-control_mean)/control_sd;
      double tail = 0.5 * std::erfc(z_score / sqrt2);
      p_values.push_back(tail);
      z_scores.push_back(z_score);
      bool is_significant = tail < significance_level;
      is_silac = is_silac || is_significant;
      significant_distances.push_back(is_significant);
    }

    z_scores_ = z_scores;
    p_values_ = p_values;
    significance_level_ = significance_level;
    is_silac_ = is_silac;
    significant_distances_ = significant_distances;
    distance_count_ = distance_count;
    std::cout << std::endl;

    std::cout << *this;

    std::cout << "Distanz 4: " << distance_count[4] << " Zscore: " << (distance_count[4]-control_mean)/control_sd<< std::endl;
    std::cout << "Distanz 6: " << distance_count[6] << " Zscore: " << (distance_count[6]-control_mean)/control_sd<<std::endl;
    std::cout << "Distanz 8: " << distance_count[8] << " Zscore: " << (distance_count[8]-control_mean)/control_sd<<std::endl;
    std::cout << "Distanz 10: " << distance_count[10] << " Zscore: " << (distance_count[10]-control_mean)/control_sd<<std::endl;
    std::cout << "Distanz 11: " << distance_count[11] << " Zscore: " << (distance_count[11]-control_mean)/control_sd<<std::endl;
    std::cout << "Distanz 14: " << distance_count[14] << " Zscore: " << (distance_count[14]-control_mean)/control_sd<<std::endl;
    std::cout << "Distanz 15: " << distance_count[15] << " Zscore: " << (distance_count[15]-control_mean)/control_sd<<std::endl;
    std::cout << "Distanz 21: " << distance_count[21] << " Zscore: " << (distance_count[21]-control_mean)/control_sd<<std::endl;
    std::cout << "Distanz 23: " << distance_count[23] << " Zscore: " << (distance_count[23]-control_mean)/control_sd<<std::endl;
    std::cout << "Distanz 27: " << distance_count[27] << " Zscore: " << (distance_count[27]-control_mean)/control_sd<<std::endl;
    
    return is_silac;
  }

  std::vector<double> SILACDetector::getZScores() const
  {
    return z_scores_;
  }

  std::vector<double> SILACDetector::getPValues() const
  {
    return p_values_;
  }

  double SILACDetector::getZScoreD4() const
  {
    return z_scores_[0];
  }

  double SILACDetector::getZScoreD6() const
  {
    return z_scores_[1];
  }

  double SILACDetector::getZScoreD8() const
  {
    return z_scores_[2];
  }

  double SILACDetector::getZScoreD10() const
  {
    return z_scores_[3];
  }

  double SILACDetector::getPValueD4() const
  {
    return p_values_[0];
  }

  double SILACDetector::getPValueD6() const
  {
    return p_values_[1];
  }
  
  double SILACDetector::getPValueD8() const
  {
    return p_values_[2];
  }

  double SILACDetector::getPValueD10() const
  {
    return p_values_[3];
  }

  double SILACDetector::getSignificanceLevel() const
  {
    return significance_level_;
  }

  bool SILACDetector::getIsSILAC() const
  {
    return is_silac_;
  }

  std::vector<bool> SILACDetector::getSignificantDistances() const
  {
    return significant_distances_;
  } 

 /*  std::ostream& operator<<(std::ostream& os, const SILACTestStatistics& silac_statistic)
  {
    std::vector<String> aminoacids = {"Medium Lysine", "Heavy Lysine(K6) or Medium Arginine", "Heavy Lysine(K8)", "Heavy Arginine"};
    std::vector<double> p_values = {silac_statistic.d4_p_value, silac_statistic.d6_p_value, silac_statistic.d8_p_value, silac_statistic.d10_p_value};

    os << "Distance 4: Z-score: " << silac_statistic.d4_z_score << " - p-value: " << silac_statistic.d4_p_value << '\n'
       << "Distance 6: Z-score: " << silac_statistic.d6_z_score << " - p-value: " << silac_statistic.d6_p_value << '\n'
       << "Distance 8: Z-score: " << silac_statistic.d8_z_score << " - p-value: " << silac_statistic.d8_p_value << '\n'
       << "Distance 10: Z-score: " << silac_statistic.d10_z_score << " - p-value: " << silac_statistic.d10_p_value << '\n';

    if (silac_statistic.is_silac_dataset)
    {
      os << '\n' << "Null hypothesis was rejected on significance level " << silac_statistic.significance_level << '\n'
         << "The following aminoacids have been detected:" << '\n';
        
      for (auto i = 0; i < 4; i++)
      {
        if(silac_statistic.significant_distances[i])
        {
          os << aminoacids[i] << " with p-value of " << p_values[i] << '\n';
        }
      } 
    }
    else
    {
      os << '\n' << "Null hypothesis was not rejected on significance level " << silac_statistic.significance_level << '\n';
    }
    return os;
  } */

  std::ostream& operator<<(std::ostream& os, const SILACDetector& silac_statistic)
  {
    std::vector<String> aminoacids = {"Medium Lysine", "Heavy Lysine(K6) or Medium Arginine", "Heavy Lysine(K8)", "Heavy Arginine"};
    //std::vector<double> p_values = {silac_statistic.getPValueD4(), silac_statistic.getPValueD6(), silac_statistic.getPValueD8(), silac_statistic.getPValueD10()};

    os << "\nDistance 4: Z-score: " << silac_statistic.getZScoreD4() << " - p-value: " << silac_statistic.getPValueD4() << '\n'
       << "Distance 6: Z-score: " << silac_statistic.getZScoreD6() << " - p-value: " << silac_statistic.getPValueD6() << '\n'
       << "Distance 8: Z-score: " << silac_statistic.getZScoreD8() << " - p-value: " << silac_statistic.getPValueD8() << '\n'
       << "Distance 10: Z-score: " << silac_statistic.getZScoreD10() << " - p-value: " << silac_statistic.getPValueD10() << '\n';

    if (silac_statistic.getIsSILAC())
    {
      os << '\n' << "Null hypothesis was rejected on significance level " << silac_statistic.getSignificanceLevel() << '\n'
         << "The following aminoacids have been detected:" << '\n';
        
      for (auto i = 0; i < 4; i++)
      {
        if(silac_statistic.getSignificantDistances()[i])
        {
          os << aminoacids[i] << " with p-value of " << silac_statistic.getPValues()[i] << '\n';
        }
      } 
    }
    else
    {
      os << '\n' << "Null hypothesis was not rejected on significance level " << silac_statistic.getSignificanceLevel() << '\n';
    }
    return os;
  } 
} // namespace OpenMS