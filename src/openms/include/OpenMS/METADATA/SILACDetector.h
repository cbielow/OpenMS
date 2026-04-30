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
#include <OpenMS/KERNEL/AreaIterator.h>
#include <functional>
#include <sstream>
//#include <map>

namespace OpenMS
{
  /** 
  @ingroup Kernel

  @brief ...
  */
  
  class OPENMS_DLLAPI SILACDetector 
  {
    std::vector<int> scan_numbers_ = {};

public:
  
  bool detectSILAC(MSExperiment experiment);

  };
} // namespace OpenMS