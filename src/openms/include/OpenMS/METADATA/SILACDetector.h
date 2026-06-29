// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Markus Apel, Nora Heese $
// --------------------------------------------------------------------------
 
#pragma once
 
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/CONCEPT/LogStream.h>

namespace OpenMS
{

  class MSExperiment;
  /**
    @brief struct which contains the relevant data of a scan for SILAC detection

    See @ref SILACDetector
  */
  struct MS2Data
  {
    double RT;
    double mz;
    int charge;
  };



  /** 
    @ingroup Metadata

    @brief This class is used to detect whether a dataset from a MSexperiment is a SILAC dataset or not
  */  
  class OPENMS_DLLAPI SILACDetector 
  {

public:
  
  /**
    @brief Determines if a dataset is a SILAC dataset

    It counts the distances of 4 aminoacid-isotopes and compares them to the distribution of control distances.
    The SILAC distances are:

    - 4 for medium lysine (K4)
    - 6 for medium arginine (R6) or heavy lysine (K6)
    - 8 for heavy lysine (K8)
    - 10 for heavy arginin (R10)

    Control distances: 11, 14, 15, 21, 23, 27

    The z-scores and p-values will be saved inside the SILACDetector object for each SILAC distance

    The code is based on param-medic (https://github.com/dhmay/param-medic/blob/master/parammedic/mod_inference.py)
    but heavily modified (e.g. mass filtering uses ppm, and window lookup uses RT distance instead of scan count; overall classification performance is a lot better).

    @param MS2Scans A vector of the relevant data of the MS2 scans of an experiment (RT, mz, charge). The vector can be created from an MSExperiment using the msExperimentToMS2Data() function
    @return True if p-value of any SILAC distance is significant on level 0.0125 (1.25%), false if none are significant (or the input vector is empty)
  */
  bool detectSILAC(std::vector<MS2Data> MS2Scans);

  /// Returns the z-scores of the SILAC distances
  std::vector<double> getZScores() const;

  /// Returns the p-values of the SILAC distances
  std::vector<double> getPValues() const;

  /// Returns the z-score of distance 4
  double getZScoreD4() const;

  /// Returns the z-score of distance 6
  double getZScoreD6() const;

  /// Returns the z-score of distance 8
  double getZScoreD8() const;

  /// Returns the z-score of distance 10
  double getZScoreD10() const;

  /// Returns the p-value of distance 4
  double getPValueD4() const;

  /// Returns the p-value of distance 6
  double getPValueD6() const;
  
  /// Returns the p-value of distance 8
  double getPValueD8() const;

  /// Returns the p-value of distance 10
  double getPValueD10() const;

  /// Returns the significance level
  double getSignificanceLevel() const;

  /// Returns if the dataset is a SILAC dataset according to the set siginficance level
  bool getIsSILAC() const;

  /// Returns vector wich contains whether a distance is significant or not(1 if significant, 0 if not sigificant)
  const std::vector<bool>& getSignificantDistances() const;

  /// Ostream iterator to write the statistical data to a stream
  friend OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const SILACDetector& silac_statistic);

  /**
  @brief Takes the relevant data of an MSExperiment for a SILACDetector anlysis and stores it in a txt file

  The data is arranged in three columns. The first column stores the retention time.
  The second column stores the mass to charge ratio. The third column stores the charge.

  @param experiment The experiment from which to store the data
  @param filename The name of the file to store the data in 
  @throw Exception::InvalidValue Throws an exception if the experiment is empty or if the experiment does not contain any MS2 scans
   */
  void storeMS2Data(const MSExperiment& experiment, const std::string& filename) const;

  /**
  @brief Takes the relevant data of an MSExperiment for a SILACDetector anlysis and returns it in a vector
  @param experiment The experiment from which to store the relevant MS2Data for SILACDetector
  @return A vector with the relevant MS2Data (RT, mz, charge) for SILACDetector
  @throw Exception::InvalidValue Throws an exception if the experiment is empty or if the experiment does not contain any MS2 scans
   */
  std::vector<MS2Data> msExperimentToMS2Data(const MSExperiment& experiment) const;

  /**
  @brief Takes a txt file with the relevant MS2Data and stores them into a vector

  The txt-file needs to have 3 columns, seperated by one space:

  - first column: retention time (RT)
  - second column: mass to charge ratio (mz)
  - third column: charge

  @param file_name The name of the input file which contains the MS2Data
  @return A vector of MS2Data for SILACDetector from the file
  @throw Exception::InvalidFileType Throws an exception if the input file is not a txt file
  @throw Exception::FileNotFound Throws an exception if the file can not be found
  @throw Exception::InvalidSize Throws an exception if the file does not contain exactly 3 columns
  @throw Exception::InvalidValue Throws an exception if the data inside the file can not be converted into doubles (RT or mz) or int (charge)
   */
  std::vector<MS2Data> txtFileToMS2Data(const std::string& file_name) const;

private:

  /// Stores the z scores for the distances (4, 6, 8, 10)
  std::vector<double> z_scores_ = {NAN,NAN,NAN,NAN};

  /// Stores the p values for the distances (4, 6, 8, 10)
  std::vector<double> p_values_ = {NAN,NAN,NAN,NAN};

  /// Significance level as a cut off value
  double significance_level_ = 0.0125;

  /// True if the dataset is likely a SILAC dataset, false otherwise
  bool is_silac_ = false;

  /// Stores which distances are significant or not (for 4, 6, 8, 10), true is significant, false is not significant
  std::vector<bool> significant_distances_ = {0, 0, 0, 0};

  /// Map for counting the detected distances for both relevant and control distances
  std::map<int,int> distance_count_ = {{4,0},{6,0},{8,0},{10,0},{11,0},{14,0},{15,0},{21,0},{23,0},{27,0}};

  };
} // namespace OpenMS