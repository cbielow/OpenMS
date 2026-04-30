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
  SILACDetector test;
   MzMLFile myfile = MzMLFile();
  MSExperiment experiment = MSExperiment();
  std::cout << "Datei wird geladen..." << std::endl;
  myfile.load("/home/nora/Documents/SoftwareProjekt/PXD070854-tripleSILAC/20111121_orbi_rk_180_004_7173MS_treated_4.mzML",experiment);
  std::cout << "Datei wurde geladen" << std::endl;
  
  TEST_EQUAL(test.detectSILAC(experiment),true) 
  /* MzMLFile myfile_2 = MzMLFile();
  MSExperiment experiment_2 = MSExperiment();
  std::cout << "Datei wird geladen..." << std::endl;
  myfile_2.load("/home/nora/Documents/SoftwareProjekt/Datensaetze/SILAC/PIAS3_SILAC_REP2_SCX2_1.mzML",experiment_2);
  std::cout << "Datei wurde geladen" << std::endl;
  test.detectSILAC(experiment_2);  */
  
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST