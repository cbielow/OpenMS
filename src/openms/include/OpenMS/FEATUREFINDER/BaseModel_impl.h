// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FEATUREFINDER/BaseModel.h>

// include derived classes here
#include <OpenMS/FEATUREFINDER/GaussModel.h>
#include <OpenMS/FEATUREFINDER/BiGaussModel.h>
#include <OpenMS/FEATUREFINDER/IsotopeModel.h>
#include <OpenMS/FEATUREFINDER/ExtendedIsotopeModel.h>
#include <OpenMS/FEATUREFINDER/ProductModel.h>
#include <OpenMS/FEATUREFINDER/EmgModel.h>
#include <OpenMS/SIMULATION/EGHModel.h>

#include <OpenMS/CONCEPT/Factory.h>

namespace OpenMS
{

  template <>
  OPENMS_DLLAPI void BaseModel<2>::registerChildren()
  {
    Factory<BaseModel<2> >::registerProduct(ProductModel<2>::getProductName(), &ProductModel<2>::create);
  }

  template <>
  OPENMS_DLLAPI void BaseModel<1>::registerChildren()
  {

    Factory<BaseModel<1> >::registerProduct(GaussModel::getProductName(), &GaussModel::create);
    Factory<BaseModel<1> >::registerProduct(BiGaussModel::getProductName(), &BiGaussModel::create);
    Factory<BaseModel<1> >::registerProduct(IsotopeModel::getProductName(), &IsotopeModel::create);
    Factory<BaseModel<1> >::registerProduct(ExtendedIsotopeModel::getProductName(), &ExtendedIsotopeModel::create);
    Factory<BaseModel<1> >::registerProduct(EmgModel::getProductName(), &EmgModel::create);
    Factory<BaseModel<1> >::registerProduct(EGHModel::getProductName(), &EGHModel::create);

    return;
  }

} // namespace OpenMS

