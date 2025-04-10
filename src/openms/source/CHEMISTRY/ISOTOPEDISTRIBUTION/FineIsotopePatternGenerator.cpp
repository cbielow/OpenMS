// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Rost $
// $Authors: Hannes Rost, Michał Startek, Mateusz Łącki $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/FineIsotopePatternGenerator.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsoSpecWrapper.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>

namespace OpenMS
{

IsotopeDistribution FineIsotopePatternGenerator::run(const EmpiricalFormula& formula) const
{

  if (probabilityMode_ == ProbabilityMode::total_prob)
  {
    IsotopeDistribution result(IsoSpecTotalProbWrapper(formula, 1.0 - stop_condition_, true).run());
    result.sortByMass();
    return result;
  }
  else if (probabilityMode_ == ProbabilityMode::absolute)
  {
    IsotopeDistribution result(IsoSpecThresholdWrapper(formula, stop_condition_, true).run());
    result.sortByMass();
    return result;
  }
  else
  {
    IsotopeDistribution result(IsoSpecThresholdWrapper(formula, stop_condition_, false).run());
    result.sortByMass();
    return result;
  }
}

} // namespace OpenMS
