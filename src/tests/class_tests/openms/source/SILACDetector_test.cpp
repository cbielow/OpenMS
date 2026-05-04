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

START_SECTION(bool detectSILAC(MSExperiment experiment))
{
  SILACDetector test;
  MzMLFile myfile = MzMLFile();
  MSExperiment experiment = MSExperiment();
  std::cout << "Datei wird geladen..." << std::endl;
  myfile.load("/home/markus/Documents/SoftwarePraktikum/Daten/PXD064675-dualSILAC/J-078-COXIV-HL_3.mzML",experiment);
  std::cout << "Datei wurde geladen" << std::endl;
  TEST_EQUAL(test.detectSILAC(experiment),true)

  /* MzMLFile myfile_2 = MzMLFile();
  MSExperiment experiment_2 = MSExperiment();
  std::cout << "Datei wird geladen..." << std::endl;
  myfile_2.load("/home/markus/Documents/SoftwarePraktikum/Daten/COREAD_MS3_SET1_10NOV14_colonCa-TMT_MS3_33_1st.mzML",experiment_2);
  std::cout << "Datei wurde geladen" << std::endl;
  TEST_EQUAL(test.detectSILAC(experiment_2), false); */
  
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST