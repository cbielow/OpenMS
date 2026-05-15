/// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Markus Apel, Nora Heese $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/KERNEL/DPeak.h>
#include <OpenMS/METADATA/SILACDetector.h>

namespace OpenMS
{
  bool SILACDetector::detectSILAC(const std::vector<MS2Data> MS2Scans)
  {
    if (MS2Scans.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No Input data", "");
    }
    for (Size i = 1; i < MS2Scans.size(); i++)
    {
      if (MS2Scans[i-1].RT > MS2Scans[i].RT)
      {
        throw Exception::NotSorted(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "input data is not sorted by retention time");
      }
    }
    // Look up table for checking distances                 Lysine (K4)        K6 or R6           Lysine (K8)        Arginine(R10)
    const double exact_mass_lookup_table[28] = {0, 0, 0, 0, 4.025106983784, 0, 6.020129012016, 0, 8.014198800046, 0, 10.008268588075996,
                                                // Control counts
                                                11, 0, 0, 14, 15, 0 ,0 ,0 ,0 ,0 ,21 ,0 ,23 ,0 ,0 ,0 , 27};
    
    const std::vector<int> control_distances = {11, 14, 15, 21, 23, 27};
    const std::vector<int> silac_distances = {4, 6, 8, 10};

    std::map<int,int> distance_count = {{4,0},{6,0},{8,0},{10,0},{11,0},{14,0},{15,0},{21,0},{23,0},{27,0}};
    const double RT_window = 5; // The size of the window in seconds of the retention time to compare after the current MS2 scan
    int min_index = 0;
    int max_index = 0;
    const int experiment_MS2_size = MS2Scans.size();

    for (const auto& spectrum : MS2Scans)
    {
      double spectrum_rt = spectrum.RT;
      min_index = spectrum.index + 1;
      // Get the index of the last scan in the window
      while ((MS2Scans[max_index].RT < spectrum_rt + RT_window) && max_index < experiment_MS2_size)
      {
        max_index++;
      } 

      // Compare every scan in the window to the current scan
      for (int i = min_index; i < max_index; i++)
      {
        if (spectrum.charge == MS2Scans[i].charge)
        {
          // Calculate the mass difference in dalton
          double distance = std::abs((spectrum.mz - MS2Scans[i].mz) * spectrum.charge);
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
    
    const double n = 6; // Size of control distances
    double sum = 0.0;

    // Sum all the counts for the control distances and calculate the mean
    for (int i = 0; i < n; i++)
    {
      sum += distance_count[control_distances[i]];
    }
    double control_mean = sum / n;
    // Calculate the standard deviation
    double control_sd = 0;
    double sd_sum = 0;
    for (int i = 0; i < 6; i++)
    {
      double s = distance_count[control_distances[i]]-control_mean;
      sd_sum += s * s;
    }
    control_sd = std::sqrt(sd_sum/(n-1));

    // If no control counts have been found, throw an exception
    if (sd_sum == 0)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No counts for control distances found", "");
    }

    bool is_silac = false;
    double z_score;
    const double significance_level = 0.0125;
    std::vector<bool> significant_distances = {}; // Stores which distances are significant or not (for 4, 6, 8, 10), true is significant, false is not significant
    std::vector<double> p_values;
    std::vector<double> z_scores;
    const double sqrt2 = std::sqrt(2.0);

    // Calculates z-score and p-value for each SILAC distance, if there are any sgnificant distances then is_silac will be set to true
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

    // Stores the results into the object
    z_scores_ = z_scores;
    p_values_ = p_values;
    significance_level_ = significance_level;
    is_silac_ = is_silac;
    significant_distances_ = significant_distances;
    distance_count_ = distance_count;
    std::cout << std::endl;

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

  std::ostream& operator<<(std::ostream& os, const SILACDetector& silac_statistic)
  {
    const std::vector<String> aminoacids = {"Medium Lysine(K4)", "Heavy Lysine(K6) or Medium Arginine(R6)", "Heavy Lysine(K8)", "Heavy Arginine(R10)"};

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

  void SILACDetector::storeMS2Data(MSExperiment experiment, const String filename) const
  {
    if (experiment.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "MSExperiment is empty", "");
    }
    if (!experiment.containsScanOfLevel(2))
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Dataset does not contain any MS2 scans", "");
    }
    if (!experiment.isSorted())
    {
      experiment.sortSpectra();
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

  std::vector<MS2Data> SILACDetector::msExperimentToMS2Data(MSExperiment experiment) const
  {
    if (experiment.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "MSExperiment is empty", "");
    }
    if (!experiment.containsScanOfLevel(2))
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Dataset does not contain any MS2 scans", "");
    }
    if (!experiment.isSorted())
    {
      experiment.sortSpectra();
    }

    std::vector<MS2Data> MS2Scans;
    int position = 0;
    for (unsigned int i = 0; i < experiment.size(); i++)
    {
      if (2 == experiment[i].getMSLevel())
      {
        MS2Data current_MS2_scan;
        current_MS2_scan.RT = experiment[i].getRT();
        current_MS2_scan.mz = experiment[i].getPrecursors()[0].getMZ();
        current_MS2_scan.charge = experiment[i].getPrecursors()[0].getCharge();
        current_MS2_scan.index = position;
        MS2Scans.push_back(current_MS2_scan);
        position ++;
      }
    }
    return MS2Scans;
  }

  std::vector<MS2Data> SILACDetector::txtFileToMS2Data(const std::string file_name) const
  {
    if (file_name.substr(file_name.size() - 4, 4) != ".txt")
    {
      throw Exception::InvalidFileType(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, file_name, "file is not a txt file");
    }
    std::vector<MS2Data> result;
    int index = 0;  // Index of current MS2 scan
    std::ifstream input_file (file_name);   
    std::string current_line;
    if (!input_file.is_open())
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, file_name);
    }
    if (input_file.is_open())
    {
      while (input_file.peek()!=EOF)
      {    
        std::getline (input_file, current_line);
        std::vector<std::string> data_values; // Values of the current line from the file (RT, mz, charge)
        size_t pos = 0; // Current position in the line
        std::string data_value; // Current value of the current line up to position

        // Finds all data values from the file and put them into a vector, the values are separated by a " "
        while ((pos = current_line.find(" ")) != std::string::npos) 
        {
          data_value = current_line.substr(0, pos);
          data_values.push_back(data_value);
          current_line.erase(0, pos + 1); // Removes the current found value
        }
        data_values.push_back(current_line);
        if (data_values.size() !=3 )
        {
          throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, data_values.size(), "File does not have exactly 3 colums");
        }
        MS2Data current_data; 
        try
        {
          current_data.RT = std::stod(data_values[0]);
        }
        catch(const std::exception& e)
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "could not convert value to double", data_values[0]);
        }
        try
        {
          current_data.mz = std::stod(data_values[1]);
        }
        catch(const std::exception& e)
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "could not convert value to double", data_values[1]);
        }
        try
        {
          current_data.charge = std::stoi(data_values[2]);
        }
        catch(const std::exception& e)
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "could not convert value to integer", data_values[2]);
        }
        current_data.index = index;
        result.push_back(current_data);
        index++;
      }
    }
    return result;
  }
} // namespace OpenMS


