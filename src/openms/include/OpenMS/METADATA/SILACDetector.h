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

  @param experiment The MS experiment to check
  @return True if p-value of any SILAC distance is significant on level 0.025, false if none are significant
  @throw Exception::InvalidValue Throws an exception if the experiment is empty, if the experiment does not contain any MS2 scans, and if no counts are counted for the control distances
   
  */
  bool detectSILAC(MSExperiment experiment);

  };
} // namespace OpenMS