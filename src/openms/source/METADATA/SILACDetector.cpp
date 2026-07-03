/// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Markus Apel, Nora Heese $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/SILACDetector.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <algorithm>
#include <fstream>

namespace OpenMS
{
  bool SILACDetector::detectSILAC(std::vector<MS2Data> MS2Scans)
  {
    if (MS2Scans.empty())
    {
      OPENMS_LOG_WARN << "SILACDetector: Input vector is empty --> SILAC detection will return 'false' (no significant distances found)!" << std::endl;
      return false;
    }
    // sort by RT:
    auto sorted_lambda = [](const MS2Data& a, const MS2Data& b) { return a.RT < b.RT; };
    if (!std::is_sorted(MS2Scans.begin(), MS2Scans.end(), sorted_lambda))
    {
      std::sort(MS2Scans.begin(), MS2Scans.end(), sorted_lambda);
    }

    // Look up table for checking distances                 Lysine (K4)        K6 or R6           Lysine (K8)        Arginine(R10)
    const double exact_mass_lookup_table[28] = {0, 0, 0, 0, 4.025106983784, 0, 6.020129012016, 0, 8.014198800046, 0, 10.008268588075996,
                                                // Control counts
                                                11, 0, 0, 14, 15, 0 ,0 ,0 ,0 ,0 ,21 ,0 ,23 ,0 ,0 ,0 , 27};
    
    const std::array<int, 6> control_distances = {11, 14, 15, 21, 23, 27};
    const std::array<int, 4> silac_distances = {4, 6, 8, 10};

    std::map<int,int> distance_count = {{4,0},{6,0},{8,0},{10,0},{11,0},{14,0},{15,0},{21,0},{23,0},{27,0}};
    const double RT_window = 5; // The size of the window in seconds of the retention time to compare after the current MS2 scan
    const size_t experiment_MS2_size = MS2Scans.size();

    // MS2Scans is sorted by RT (see above), so for each scan we only need to look at the *following* scans
    // (by position) until we leave its RT window. The window walk is purely position-based.
    for (size_t i = 0; i < experiment_MS2_size; ++i)
    {
      const MS2Data& spectrum = MS2Scans[i];
      // Compare the current MS2 scan with all following MS2 scans within the RT window of 5 seconds
      for (size_t idx_window = i + 1; idx_window < experiment_MS2_size && (MS2Scans[idx_window].RT < spectrum.RT + RT_window); ++idx_window)
      {
        if (spectrum.charge == MS2Scans[idx_window].charge)
        {
          // Calculate the mass difference in dalton
          double distance = std::abs((spectrum.mz - MS2Scans[idx_window].mz) * spectrum.charge);
          // Round to int for the count map distance_count
          int rounded_distance = distance + 0.5;
          // Check if the distance is in the map
          if (distance_count.find(rounded_distance) != distance_count.end())
          {
            // Check if the difference of the mass is close enough to the exact mass (5 ppm)
            if (std::abs(distance - exact_mass_lookup_table[rounded_distance]) < 0.005)
            {
              distance_count[rounded_distance]++;
            }
          }
        }
      }
    }
    
    const double n = control_distances.size(); // Size of control distances
    double sum = 0.0;

    // Sum all the counts for the control distances and calculate the mean
    for (int i = 0; i < n; i++)
    {
      sum += distance_count[control_distances[i]];
    }
    double control_mean = sum / n;
    // Calculate the standard deviation
    double sd_sum = 0;
    for (int i = 0; i < n; i++)
    {
      double s = distance_count[control_distances[i]] - control_mean;
      sd_sum += s * s;
    }
    double control_sd = std::sqrt(sd_sum/(n-1));

    // If no control counts have been found
    if (sd_sum == 0)
    {
      OPENMS_LOG_WARN << "SILACDetector: No counts for control distances have been found --> using pseudo counts" << std::endl;
      // add pseudo counts to avoid division by zero
      control_mean = 1.0;
      control_sd = 1.0;
    }

    z_scores_.clear();
    p_values_.clear();
    is_silac_ = false;
    significant_distances_.clear();
    
    const double sqrt2 = std::sqrt(2.0);

    // Calculates z-score and p-value for each SILAC distance, if there are any significant distances then is_silac will be set to true
    for (int i = 4; i <= 10; i += 2)
    {
      double z_score = (distance_count[i] - control_mean) / control_sd;
      double tail = 0.5 * std::erfc(z_score / sqrt2);
      const bool is_significant = tail < significance_level_;
      z_scores_.push_back(z_score);
      p_values_.push_back(tail);
      significant_distances_.push_back(is_significant);
      is_silac_ |= is_significant;
    }

    return is_silac_;
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

  const std::vector<bool>& SILACDetector::getSignificantDistances() const
  {
    return significant_distances_;
  } 

  std::ostream& operator<<(std::ostream& os, const SILACDetector& silac_statistic)
  {
    const std::vector<std::string> aminoacids = {"Medium Lysine(K4)", "Heavy Lysine(K6) or Medium Arginine(R6)", "Heavy Lysine(K8)", "Heavy Arginine(R10)"};

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
      if ((silac_statistic.getSignificantDistances()[0] && silac_statistic.getSignificantDistances()[2])
      || (silac_statistic.getSignificantDistances()[1] && silac_statistic.getSignificantDistances()[3])
      || (silac_statistic.getSignificantDistances()[0] && silac_statistic.getSignificantDistances()[1]))
      {
        os << "\nThis Datasaet is likely a triple SILAC\n";
      }
    }
    else
    {
      os << '\n' << "Null hypothesis was not rejected on significance level " << silac_statistic.getSignificanceLevel() << '\n';
    }
    return os;
  } 

  void SILACDetector::storeMS2Data(const MSExperiment& experiment, const std::string& filename) const
  {
    if (experiment.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "MSExperiment is empty", "");
    }
    if (!experiment.empty() && !experiment.containsScanOfLevel(2))
    { // there are spectra, just no MS2 spectra (e.g. MS1-only run), but we need MS2 spectra for SILAC detection. Omitting MS2 spectra is a user error, so we throw an exception here.
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Dataset does not contain any MS2 scans", "");
    }
    
    std::ofstream out(filename);
    for (unsigned int i = 0; i < experiment.size(); i++)
    {
      if (2 == experiment[i].getMSLevel())
      {
        out << std::setprecision(17);
        out << experiment[i].getRT() << " " << experiment[i].getPrecursors()[0].getMZ() << " " << experiment[i].getPrecursors()[0].getCharge() << "\n";
      }
    }
  }

  std::vector<MS2Data> SILACDetector::msExperimentToMS2Data(const MSExperiment& experiment) const
  {
    if (!experiment.empty() && !experiment.containsScanOfLevel(2))
    { // there are spectra, just no MS2 spectra (e.g. MS1-only run), but we need MS2 spectra for SILAC detection. Omitting MS2 spectra is a user error, so we throw an exception here.
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Dataset does not contain any MS2 scans", "");
    }

    Int32 zero_charge_found {0};
    std::vector<MS2Data> MS2Scans;
    for (const auto& spectrum : experiment)
    {
      if (2 == spectrum.getMSLevel() && !spectrum.getPrecursors().empty())
      {
        MS2Data current_MS2_scan;
        current_MS2_scan.RT = spectrum.getRT();
        current_MS2_scan.mz = spectrum.getPrecursors()[0].getMZ();
        current_MS2_scan.charge = spectrum.getPrecursors()[0].getCharge();
        if (current_MS2_scan.charge == 0)
        { // detection with charge 0 will not work for precursors with true charge 2(or above), since the mass difference computation needs a charge
          // charge=0 happens with old WIFF files + msconvert (even in 03/2026)
          ++zero_charge_found;
          current_MS2_scan.charge = 2; // hack, but we cannot to much else here. See below for user warning message
        }
        MS2Scans.push_back(current_MS2_scan);
      }
    }
    if (zero_charge_found)
    {
      OPENMS_LOG_WARN << "SILACDetector: " << (zero_charge_found == MS2Scans.size() ? "All" : "Some") << 
                         " MS2 scans had charge 0, which is not supported for SILAC detection.These scans were treated as"
                         " charge 2, but the results may be inaccurate. Please check your data and conversion settings."
                      << std::endl;
    }
    return MS2Scans;
  }
} // namespace OpenMS


