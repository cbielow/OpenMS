// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FEATUREFINDER/InterpolationModel.h>
#include <OpenMS/MATH/STATISTICS/BasicStatistics.h>

namespace OpenMS
{
  /**
      @brief Exponential-Gaussian hybrid distribution model for elution profiles.

  Lan K, Jorgenson JW.
  A hybrid of exponential and gaussian functions as a simple model of asymmetric chromatographic peaks.
  Journal of Chromatography A. 2001;915(1-2):1-13.
  Available at: http://linkinghub.elsevier.com/retrieve/pii/S0021967301005945

  @htmlinclude OpenMS_EGHModel.parameters

  */
  class OPENMS_DLLAPI EGHModel :
    public InterpolationModel
  {

public:
    typedef InterpolationModel::CoordinateType CoordinateType;
    typedef Math::BasicStatistics<CoordinateType> BasicStatistics;
    typedef LinearInterpolation::container_type ContainerType;

    /// Default constructor
    EGHModel();

    /// copy constructor
    EGHModel(const EGHModel & source);

    /// destructor
    ~EGHModel() override;

    /// assignment operator
    virtual EGHModel & operator=(const EGHModel & source);

    /// create new ElutionModel object (needed by Factory)
    static BaseModel<1> * create()
    {
      return new EGHModel();
    }

    /// name of the model (needed by Factory)
    static const String getProductName()
    {
      return "EGHModel";
    }

    /// set offset without being computing all over and without any discrepancy
    void setOffset(CoordinateType offset) override;

    /// set sample/supporting points of interpolation
    void setSamples() override;

    /// get the center of the Gaussian model i.e. the position of the maximum
    CoordinateType getCenter() const override;

protected:
    CoordinateType  min_;
    CoordinateType  max_;
    BasicStatistics statistics_;
    CoordinateType  height_;     // H in paper
    CoordinateType  apex_rt_;

    CoordinateType  A_;
    CoordinateType  B_;

    CoordinateType  tau_;
    CoordinateType  sigma_square_;
    CoordinateType  sigma_square_2_;


    void updateMembers_() override;

    /// Computes a left & right boundary for the EGH Profile and sets the internal parameters accordingly
    void computeBoundaries_();

    /**
     * @brief Evaluate the EGH function at position rt
     *
     * @param rt        The position where the EGH function should be evaluated. Note that this is the position without the RT offset, meaning that the EGH apex is at position 0
     * @param egh_value The computed value
     */
    inline void evaluateEGH_(CoordinateType & rt, CoordinateType & egh_value) const
    {
      CoordinateType denominator = sigma_square_2_ + tau_ * rt;

      if (denominator > 0)
      {
        // evaluate egh ->
        egh_value = height_ * exp(
          (-1 * rt * rt) / denominator
          );
      }
      else
      {
        egh_value = 0.0;
      }
    }

  };

} // namespace OpenMS

