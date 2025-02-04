// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/SIMULATION/LABELING/BaseLabeler.h>

// include derived classes here
#include <OpenMS/SIMULATION/LABELING/ITRAQLabeler.h>
#include <OpenMS/SIMULATION/LABELING/LabelFreeLabeler.h>
#include <OpenMS/SIMULATION/LABELING/O18Labeler.h>
#include <OpenMS/SIMULATION/LABELING/SILACLabeler.h>
#include <OpenMS/SIMULATION/LABELING/ICPLLabeler.h>

#include <OpenMS/CONCEPT/Factory.h>

namespace OpenMS
{

  void BaseLabeler::registerChildren()
  {
    Factory<BaseLabeler>::registerProduct(LabelFreeLabeler::getProductName(), &LabelFreeLabeler::create);
    Factory<BaseLabeler>::registerProduct(O18Labeler::getProductName(), &O18Labeler::create);
    Factory<BaseLabeler>::registerProduct(ITRAQLabeler::getProductName(), &ITRAQLabeler::create);
    Factory<BaseLabeler>::registerProduct(SILACLabeler::getProductName(), &SILACLabeler::create);
    Factory<BaseLabeler>::registerProduct(ICPLLabeler::getProductName(), &ICPLLabeler::create);
    return;
  }

} // namespace OpenMS

