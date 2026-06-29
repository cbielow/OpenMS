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

#include <OpenMS/KERNEL/MSExperiment.h> // SILACDetector.h only forward-declares MSExperiment

using namespace OpenMS;
using namespace std;

START_TEST(SILACDetector, "$Id$")

START_SECTION(std::vector<MS2Data> txtFileToMS2Data(std::string file_name))
{
  SILACDetector test;
  std::vector<MS2Data> example_data = test.txtFileToMS2Data(OPENMS_GET_TEST_DATA_PATH("SILACTestData.txt"));
  TEST_REAL_SIMILAR(example_data[0].mz, 964.2548828125)
  TEST_REAL_SIMILAR(example_data[0].RT, 0.52840948799999998)
  TEST_EQUAL(example_data[0].charge, 2)
  TEST_EQUAL(example_data[0].index, 0) 
  TEST_EXCEPTION(Exception::InvalidFileType, test.txtFileToMS2Data(OPENMS_GET_TEST_DATA_PATH("20171013_HMP_C61_ISO_P1_GA1_UV_VIS_2.mzML")))
  TEST_EXCEPTION(Exception::InvalidSize, test.txtFileToMS2Data(OPENMS_GET_TEST_DATA_PATH("BSpline2d_test_sinus.txt")))
  TEST_EXCEPTION(Exception::InvalidValue, test.txtFileToMS2Data(OPENMS_GET_TEST_DATA_PATH("SILAC_exception_test.txt")))
  TEST_EXCEPTION(Exception::FileNotFound, test.txtFileToMS2Data("fkwjks.txt"))
}
END_SECTION

SILACDetector test;
std::vector<MS2Data> example_data = test.txtFileToMS2Data(OPENMS_GET_TEST_DATA_PATH("SILACTestData.txt"));
SILACDetector test_2;

START_SECTION(bool SILACDetector::detectSILAC(std::vector<MS2Data> MS2Scans))
{
  bool test_result = test.detectSILAC(example_data);
  TEST_EQUAL(test_result, 1)
  // edge cases are handled gracefully (no exceptions). Use a throwaway detector so the example_data
  // results stored in 'test' (checked by the getter sections below) are not overwritten.
  SILACDetector edge;
  // empty input -> no significant distances -> false
  std::vector<MS2Data> test_empty;
  TEST_EQUAL(edge.detectSILAC(test_empty), false)
  MS2Data control_test;
  control_test.charge = 1;
  control_test.mz = 1;
  control_test.RT = 1;
  control_test.index = 0;
  // no counts for any control distance -> pseudo counts are used -> false (no exception)
  std::vector<MS2Data> no_control_count_test = {control_test};
  TEST_EQUAL(edge.detectSILAC(no_control_count_test), false)
  // unsorted input is sorted internally (no exception)
  MS2Data sort_test = control_test;
  sort_test.RT = 0.5;
  no_control_count_test.push_back(sort_test);
  TEST_EQUAL(edge.detectSILAC(no_control_count_test), false)
}
END_SECTION

START_SECTION(std::vector<double> getZScores() const;)
{
  TEST_REAL_SIMILAR(test.getZScores()[0],-0.49148)
}
END_SECTION

START_SECTION(std::vector<double> getPValues() const)
{
  TEST_REAL_SIMILAR(test.getPValues()[0],0.688456)
}
END_SECTION

START_SECTION(double getZScoreD4() const)
{
  TEST_REAL_SIMILAR(test.getZScoreD4(), -0.49148)
}
END_SECTION

START_SECTION(double getZScoreD6() const)
{
  TEST_REAL_SIMILAR(test.getZScoreD6(), -0.277793)
}
END_SECTION

START_SECTION(double getZScoreD8() const)
{
  TEST_REAL_SIMILAR(test.getZScoreD8(), 8.3979)
}
END_SECTION

START_SECTION(double getZScoreD10() const)
{
  TEST_REAL_SIMILAR(test.getZScoreD10(), 6.43198)
}
END_SECTION

START_SECTION(double getPValueD4() const)
{
  TEST_REAL_SIMILAR(test.getPValueD4(), 0.688456)
}
END_SECTION

START_SECTION(double getPValueD6() const)
{
  TEST_REAL_SIMILAR(test.getPValueD6(), 0.609414)
}
END_SECTION

START_SECTION(double getPValueD8() const)
{
  TEST_REAL_SIMILAR(test.getPValueD8(), 2.27275e-17)
}
END_SECTION

START_SECTION(double getPValueD10() const)
{
  TEST_REAL_SIMILAR(test.getPValueD10(), 6.29777e-11)
}
END_SECTION

START_SECTION(double getSignificanceLevel() const)
{
  TEST_REAL_SIMILAR(test.getSignificanceLevel(),0.0125)
}
END_SECTION

START_SECTION(bool getIsSILAC() const)
{
  TEST_EQUAL(test.getIsSILAC(), 1)
}
END_SECTION

START_SECTION(const std::vector<bool>& getSignificantDistances() const)
{
  bool a = test.getSignificantDistances()[0];
  TEST_EQUAL(a, 0)
}
END_SECTION

START_SECTION(friend OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const SILACDetector& silac_statistic))
{
  std::cout << test << std::endl;
  std::cout << test_2 << std::endl;
}
END_SECTION

START_SECTION(void SILACDetector::storeMS2Data(MSExperiment experiment, String filename))
{
  MzMLFile myfile = MzMLFile();
  MSExperiment experiment = MSExperiment();
  TEST_EXCEPTION(Exception::InvalidValue, test.storeMS2Data(experiment, OPENMS_GET_TEST_DATA_PATH("SILACstoreTest.txt"))) // empty experiment -> throws
  myfile.load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"),experiment);
  test.storeMS2Data(experiment, OPENMS_GET_TEST_DATA_PATH("SILACstoreTest.txt"));
  std::vector<MS2Data> output = test.txtFileToMS2Data(OPENMS_GET_TEST_DATA_PATH("SILACstoreTest.txt"));
  TEST_REAL_SIMILAR(output[0].RT, 5.2000000000000002);
}
END_SECTION

START_SECTION(std::vector<MS2Data> msExperimentToMS2Data(MSExperiment experiment))
{
  MzMLFile myfile = MzMLFile();
  MSExperiment experiment = MSExperiment();
  TEST_EQUAL(test.msExperimentToMS2Data(experiment).empty(), true) // empty experiment -> empty result (no exception)
  myfile.load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"),experiment);
  std::vector<MS2Data> output = test.msExperimentToMS2Data(experiment);
  TEST_REAL_SIMILAR(output[0].RT, 5.2)
  TEST_REAL_SIMILAR(output[0].mz, 5.5555)
  TEST_EQUAL(output[0].charge, 2)
  TEST_EQUAL(output[0].index, 0)
  myfile.load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_6_uncompressed.mzML"),experiment);
  TEST_EXCEPTION(Exception::InvalidValue, test.msExperimentToMS2Data(experiment))
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST