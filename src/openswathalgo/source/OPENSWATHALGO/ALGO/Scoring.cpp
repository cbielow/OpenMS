// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <iostream>
#include <OpenMS/OPENSWATHALGO/ALGO/Scoring.h>
#include <OpenMS/OPENSWATHALGO/Macros.h>
#include <cmath>
#include <algorithm>

#include <boost/numeric/conversion/cast.hpp>

// Import dependencies from MIToolbox
#include <ArrayOperations.c>
#include <CalculateProbability.c>
#include <Entropy.c>
#include <MutualInformation.c>

namespace OpenSwath::Scoring
{
    void normalize_sum(double x[], unsigned int n)
    { 
      double sumx = std::accumulate(&x[0], &x[0] + n, 0.0);
      if (sumx == 0.0)
      { // avoid divide by zero below
        return;
      }                           
      auto inverse_sum = 1 / sumx; // precompute inverse since division is expensive!
      for (unsigned int i = 0; i < n; ++i)
      {
        x[i] *= inverse_sum;
      }
    }

    double NormalizedManhattanDist(double x[], double y[], int n)
    {
      OPENSWATH_PRECONDITION(n > 0, "Need at least one element");

      double delta_ratio_sum = 0;
      normalize_sum(x, n);
      normalize_sum(y, n);
      for (int i = 0; i < n; i++)
      {
        delta_ratio_sum += std::fabs(x[i] - y[i]);
      }
      return delta_ratio_sum / n;
    }

    double RootMeanSquareDeviation(double x[], double y[], int n)
    {
      OPENSWATH_PRECONDITION(n > 0, "Need at least one element");

      double result = 0;
      for (int i = 0; i < n; i++)
      {

        result += (x[i] - y[i]) * (x[i] - y[i]);
      }
      return std::sqrt(result / n);
    }

    double SpectralAngle(double x[], double y[], int n)
    {
      OPENSWATH_PRECONDITION(n > 0, "Need at least one element");

      double dotprod = 0;
      double x_len = 0;
      double y_len = 0;
      for (int i = 0; i < n; i++)
      {
        dotprod += x[i] * y[i];
        x_len += x[i] * x[i];
        y_len += y[i] * y[i];
      }
      x_len = std::sqrt(x_len);
      y_len = std::sqrt(y_len);

      // normalise, avoiding a divide by zero. See unit tests for what happens
      // when one of the vectors has a length of zero.
      double denominator = x_len * y_len;
      double theta = (denominator == 0) ? 0.0 : dotprod / denominator;

      // clip to range [-1, 1] to save acos blowing up
      theta = std::max(-1.0, std::min(1.0, theta));

      return std::acos(theta);
    }

    XCorrArrayType::const_iterator xcorrArrayGetMaxPeak(const XCorrArrayType& array)
    {
      OPENSWATH_PRECONDITION(array.data.size() > 0, "Cannot get highest apex from empty array.");

      XCorrArrayType::const_iterator max_it = array.begin();
      double max = array.begin()->second;
      for (XCorrArrayType::const_iterator it = array.begin(); it != array.end(); ++it)
      {
        if (it->second > max)
        {
          max = it->second;
          max_it = it;
        }
      }
      return max_it;
    }

    void standardize_data(std::vector<double>& data)
    {
      OPENSWATH_PRECONDITION(data.size() > 0, "Need non-empty array.");

      // subtract the mean and divide by the standard deviation
      double mean = std::accumulate(data.begin(), data.end(), 0.0) / (double) data.size();
      double sqsum = 0;
      for (std::vector<double>::iterator it = data.begin(); it != data.end(); ++it)
      {
        sqsum += (*it - mean) * (*it - mean);
      }
      double stdev = sqrt(sqsum / data.size()); // standard deviation

      if (mean == 0 && stdev == 0)
      {
        return; // all data is zero
      }
      if (stdev == 0)
      {
        stdev = 1; // all data is equal
      }
      stdev = 1/stdev;
      for (std::size_t i = 0; i < data.size(); i++)
      {
        data[i] = (data[i] - mean) * stdev;
      }
    }

    XCorrArrayType normalizedCrossCorrelation(std::vector<double>& data1,
                                              std::vector<double>& data2, int maxdelay, int lag = 1)  //const ref entfernt
    {
      OPENSWATH_PRECONDITION(data1.size() != 0 && data1.size() == data2.size(), "Both data vectors need to have the same length");

      // normalize the data
      standardize_data(data1);
      standardize_data(data2);
      XCorrArrayType result = calculateCrossCorrelation(data1, data2, maxdelay, lag);

      double d = 1.0 / data1.size();
      for(auto& e : result)
      {
        e.second *= d;
      }
      return result;
    }

    XCorrArrayType calculateCrossCorrelation(const std::vector<double>& data1,
                                             const std::vector<double>& data2, int maxdelay, int lag) //const ref entfernt
    {
      OPENSWATH_PRECONDITION(data1.size() != 0 && data1.size() == data2.size(), "Both data vectors need to have the same length");

      XCorrArrayType result;
      result.data.reserve( (size_t)std::ceil((2*maxdelay + 1) / lag));
      int datasize = boost::numeric_cast<int>(data1.size());
      int i, j, delay;

      for (delay = -maxdelay; delay <= maxdelay; delay = delay + lag)
      {
        double sxy = 0;
        for (i = 0; i < datasize; ++i)
        {
          j = i + delay;
          if (j < 0 || j >= datasize)
          {
            continue;
          }
          sxy += (data1[i]) * (data2[j]);
        }
        result.data.push_back(std::make_pair(delay, sxy));
      }
      return result;
    }

    XCorrArrayType calcxcorr_legacy_mquest_(std::vector<double>& data1,
                                            std::vector<double>& data2, bool normalize)
    {
      OPENSWATH_PRECONDITION(!data1.empty() && data1.size() == data2.size(), "Both data vectors need to have the same length");
      int maxdelay = boost::numeric_cast<int>(data1.size());
      int lag = 1;

      double mean1 = std::accumulate(data1.begin(), data1.end(), 0.) / (double)data1.size();
      double mean2 = std::accumulate(data2.begin(), data2.end(), 0.) / (double)data2.size();
      double denominator = 1.0;
      int datasize = boost::numeric_cast<int>(data1.size());
      int i, j, delay;

      // Normalized cross-correlation = subtract the mean and divide by the standard deviation
      if (normalize)
      {
        double sqsum1 = 0;
        double sqsum2 = 0;
        for (std::vector<double>::iterator it = data1.begin(); it != data1.end(); ++it)
        {
          sqsum1 += (*it - mean1) * (*it - mean1);
        }

        for (std::vector<double>::iterator it = data2.begin(); it != data2.end(); ++it)
        {
          sqsum2 += (*it - mean2) * (*it - mean2);
        }
        // sigma_1 * sigma_2 * n
        denominator = sqrt(sqsum1 * sqsum2);
      }
      //avoids division in the for loop
      denominator = 1/denominator;
      XCorrArrayType result;
      result.data.reserve( (size_t)std::ceil((2*maxdelay + 1) / lag));
      int cnt = 0;
      for (delay = -maxdelay; delay <= maxdelay; delay = delay + lag, cnt++)
      {
        double sxy = 0;
        for (i = 0; i < datasize; i++)
        {
          j = i + delay;
          if (j < 0 || j >= datasize)
          {
            continue;
          }
          if (normalize)
          {
            sxy += (data1[i] - mean1) * (data2[j] - mean2);
          }
          else
          {
            sxy += (data1[i]) * (data2[j]);
          }
        }

        if (denominator > 0)
        {
          result.data.emplace_back(delay, sxy*denominator);
        }
        else
        {
          // e.g. if all datapoints are zero
          result.data.emplace_back(delay, 0);
        }
      }
      return result;
    }

    void computeRank(const std::vector<double>& v_temp, std::vector<unsigned int>& ranks_out)
    {
      std::vector<unsigned int> ranks{};
      ranks.resize(v_temp.size());
      std::iota(ranks.begin(), ranks.end(), 0);
      std::sort(ranks.begin(), ranks.end(),
                [&v_temp](unsigned int i, unsigned int j) { return v_temp[i] < v_temp[j]; });
      ranks_out.resize(v_temp.size());
      double x = 0;
      unsigned int y = 0;
      for(unsigned int i = 0; i < ranks.size();++i)
      {
        if(v_temp[ranks[i]] != x)
        {
          x = v_temp[ranks[i]];
          y = i;
        }
        ranks_out[ranks[i]] = y;
      }
    }

    unsigned int maxElem(const std::vector<unsigned int>& arr)
    {
      unsigned int max = arr[0];
      for(auto e : arr)
      {
        if(e > max) max = e;
      }
      return max+1;
    }



    jpstate calcJointProbability(const std::vector<unsigned int>& firstVector,const std::vector<unsigned int>& secondVector,const int& vectorLength)
    {
      jpstate state;
      double length = 1.0 / vectorLength;
      unsigned int firstNumStates = maxElem(firstVector);
      unsigned int secondNumStates = maxElem(secondVector);
      unsigned int jointNumStates = firstNumStates * secondNumStates;

      std::vector<unsigned int> firstStateCounts(firstNumStates, 0);
      std::vector<unsigned int> secondStateCounts(secondNumStates, 0);
      std::vector<unsigned int> jointStateCounts(jointNumStates, 0);
      std::vector<unsigned int> jointPosition(firstNumStates, 0);

      std::vector<double> firstStateProbs(firstNumStates, 0.0);
      std::vector<double> secondStateProbs(secondNumStates, 0.0);
      std::vector<double> jointStateProbs(jointNumStates, 0.0);

      for(int i = 0; i < vectorLength; i++)
      {
        firstStateCounts[firstVector[i]] += 1;
        secondStateCounts[secondVector[i]] += 1;
        jointPosition[i] = secondVector[i] * firstNumStates + firstVector[i];
        jointStateCounts[jointPosition[i]] += 1;
      }

      for (unsigned int i = 0; i < firstNumStates; i++) {
        firstStateProbs[i] = firstStateCounts[i] * length;
      }

      for (unsigned int i = 0; i < secondNumStates; i++) {
        secondStateProbs[i] = secondStateCounts[i] * length;
      }

      for (unsigned int i = 0; i < jointNumStates; i++) {
        jointStateProbs[i] = jointStateCounts[i] * length;
      }

      state.jointPositionVector = jointPosition;
      state.jointProbabilityVector = jointStateProbs;
      state.numJointStates = jointNumStates;
      state.firstProbabilityVector = firstStateProbs;
      state.numFirstStates = firstNumStates;
      state.secondProbabilityVector = secondStateProbs;
      state.numSecondStates = secondNumStates;

      return state;
    }

    double mutualInformation(jpstate& state,const std::vector<unsigned int>& firstVector,const std::vector<unsigned int>& secondVector)
    {
      double mutualInformation = 0.0;
      //int firstIndex,secondIndex;

      /*
      ** I(X;Y) = \sum_x \sum_y p(x,y) * \log (p(x,y)/p(x)p(y))
      */
      for (unsigned int i = 0; i < firstVector.size(); i++)
      {

        int j = state.jointPositionVector[i];
        if(state.jointProbabilityVector[j] != 0)
        {
          /*double division is probably more stable than multiplying two small numbers together
          ** mutualInformation += state.jointProbabilityVector[i] * log(state.jointProbabilityVector[i] / (state.firstProbabilityVector[firstIndex] * state.secondProbabilityVector[secondIndex]));
          */
          mutualInformation += state.jointProbabilityVector[j] *
                               log2(state.jointProbabilityVector[j] / state.firstProbabilityVector[firstVector[i]] /
                                    state.secondProbabilityVector[secondVector[i]]);
          state.jointProbabilityVector[j] = 0;
        }
      }
      return mutualInformation;
    }

    double rankedMutualInformation(std::vector<unsigned int>& data1, std::vector<unsigned int>& data2)
    {
      OPENSWATH_PRECONDITION(data1.size() != 0 && data1.size() == data2.size(), "Both data vectors need to have the same length");

      jpstate state = calcJointProbability(data1, data2, data1.size());

      double result = mutualInformation(state, data1, data2);
      /*
      unsigned int* arr_int_data1 = &data1[0];
      unsigned int* arr_int_data2 = &data2[0];
      double result = calcMutualInformation(arr_int_data1, arr_int_data2, data1.size());
      */
      return result;
    }
}      //namespace OpenMS  // namespace Scoring
