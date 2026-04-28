// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Markus Apel, Nora Heese $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/METADATA/SILACDetector.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(SILACDetector, "$Id$")

START_SECTION(int testfunktion())
{
  MSExperiment experiment;
  SILACDetector test;
  //MSExperiment profile;
  //AASequence seq_1("ACDEF");
  //AASequence seq_2;
  //std::vector<double> vec = profile.computeHydrophobicMoment(seq_1,3,100);
  //TEST_REAL_SIMILAR(vec[0],0.511576803);
  //TEST_REAL_SIMILAR(vec[1],0.435170599);
  //TEST_REAL_SIMILAR(vec[2],0.734926405);
  //TEST_EXCEPTION(Exception::InvalidSize,profile.computeHydrophobicMoment(seq_1,0));
  //TEST_EXCEPTION(Exception::InvalidValue,profile.computeHydrophobicMoment(seq_2,3));
  TEST_EQUAL(test.testfunktion(), 1)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST