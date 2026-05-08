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
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <iostream>

namespace OpenMS
{
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

  - 4 for medium lysine
  - 6 for medium arginine or heavy lysine (K6)
  - 8 for heavy lysine (K8)
  - 10 for heavy arginin

  Control distances: 11, 14, 15, 21, 23, 27

  The z-scores and p-values are being shown on the terminal for each SILAC distance

  The code is based of the code from the param-medic GitHub page: https://github.com/dhmay/param-medic/blob/master/parammedic/mod_inference.py

  @param experiment The MS experiment to check
  @return True if p-value of any SILAC distance is significant on level 0.025 (2.5%), false if none are significant
  @throw Exception::InvalidValue Throws an exception if the experiment is empty, if the experiment does not contain any MS2 scans, and if no counts are counted for the control distances
   
  */
  bool detectSILAC(MSExperiment experiment);

  /// returns the z-scores of the SILAC distances
  std::vector<double> getZScores() const;

  /// returns the p-values of the SILAC distances
  std::vector<double> getPValues() const;

  /// returns the z-score of distance 4
  double getZScoreD4() const;

  /// returns the z-score of distance 6
  double getZScoreD6() const;

  /// returns the z-score of distance 8
  double getZScoreD8() const;

  /// returns the z-score of distance 10
  double getZScoreD10() const;

  /// returns the p-value of distance 4
  double getPValueD4() const;

  /// returns the p-value of distance 6
  double getPValueD6() const;
  
  /// returns the p-value of distance 8
  double getPValueD8() const;

  /// returns the p-value of distance 10
  double getPValueD10() const;

  /// returns the significance level
  double getSignificanceLevel() const;

  /// returns if the dataset is a SILAC dataset according to the set siginficance level
  bool getIsSILAC() const;

  /// returns vector wich contains whether a distance is significant or not(1 if significant, 0 if not sigificant)
  std::vector<bool> getSignificantDistances() const;

  /// ostream iterator to write the statistical data to a stream
  friend OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const SILACDetector& silac_statistic);

private:

  std::vector<double> z_scores_ = {NAN,NAN,NAN,NAN};
  std::vector<double> p_values_ = {NAN,NAN,NAN,NAN};
  double significance_level_ = 0.05;
  bool is_silac_ = false;
  std::vector<bool> significant_distances_ = {0, 0, 0, 0};
  std::map<int,int> distance_count_ = {{4,0},{6,0},{8,0},{10,0},{11,0},{14,0},{15,0},{21,0},{23,0},{27,0}};

  };

  /**
  @brief struct which contains the relevant data of a scan for silac detection
  */
  struct MS2Data
  {
    double RT;
    double mz;
    int charge;
    int index;
  };
} // namespace OpenMS