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
#include <OpenMS/KERNEL/ConsensusMap.h> // SILACDetector.h only forward-declares ConsensusMap
#include <OpenMS/KERNEL/FeatureMap.h>   // SILACDetector.h only forward-declares FeatureMap
#include <OpenMS/FORMAT/MzMLFile.h>     // used by the storeMS2Data / msExperimentToMS2Data sections

#include <fstream>

using namespace OpenMS;
using namespace std;

namespace
{
  /**
    @brief Test helper: reads a txt file with the relevant MS2Data into a vector.

    The txt-file needs to have 3 columns, separated by one space:
    - first column: retention time (RT)
    - second column: mass to charge ratio (mz)
    - third column: charge

    @throw Exception::InvalidFileType if the input file is not a txt file
    @throw Exception::FileNotFound if the file can not be found
    @throw Exception::InvalidSize if the file does not contain exactly 3 columns
    @throw Exception::InvalidValue if the data inside the file can not be converted into doubles (RT or mz) or int (charge)
  */
  std::vector<MS2Data> txtFileToMS2Data(const std::string& file_name)
  {
    if (!file_name.ends_with(".txt"))
    {
      throw Exception::InvalidFileType(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, file_name, "file is not a txt file");
    }
    std::vector<MS2Data> result;
    std::ifstream input_file (file_name);
    std::string current_line;
    if (!input_file.is_open())
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, file_name);
    }
    while (input_file.peek()!=EOF)
    {
      std::getline (input_file, current_line);
      std::vector<std::string> data_values; // Values of the current line from the file (RT, mz, charge)
      size_t pos = 0; // Current position in the line
      std::string data_value; // Current value of the current line up to position

      // Finds all data values from the file and put them into a vector, the values are separated by a " "
      while ((pos = current_line.find(" ")) != std::string::npos)
      {
        data_value = current_line.substr(0, pos);
        data_values.push_back(data_value);
        current_line.erase(0, pos + 1); // Removes the current found value
      }
      data_values.push_back(current_line);
      if (data_values.size() !=3 )
      {
        throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, data_values.size(), "File does not have exactly 3 colums");
      }
      MS2Data current_data;
      try
      {
        current_data.RT = std::stod(data_values[0]);
      }
      catch(const std::exception& e)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "could not convert value to double", data_values[0]);
      }
      try
      {
        current_data.mz = std::stod(data_values[1]);
      }
      catch(const std::exception& e)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "could not convert value to double", data_values[1]);
      }
      try
      {
        current_data.charge = std::stoi(data_values[2]);
      }
      catch(const std::exception& e)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "could not convert value to integer", data_values[2]);
      }
      result.push_back(current_data);
    }

    return result;
  }
} // namespace

START_TEST(SILACDetector, "$Id$")

START_SECTION([EXTRA] std::vector<MS2Data> txtFileToMS2Data(std::string file_name))
{
  std::vector<MS2Data> example_data = txtFileToMS2Data(OPENMS_GET_TEST_DATA_PATH("SILACTestData.txt"));
  TEST_REAL_SIMILAR(example_data[0].mz, 964.2548828125)
  TEST_REAL_SIMILAR(example_data[0].RT, 0.52840948799999998)
  TEST_EQUAL(example_data[0].charge, 2)
  TEST_EXCEPTION(Exception::InvalidFileType, txtFileToMS2Data(OPENMS_GET_TEST_DATA_PATH("20171013_HMP_C61_ISO_P1_GA1_UV_VIS_2.mzML")))
  TEST_EXCEPTION(Exception::InvalidSize, txtFileToMS2Data(OPENMS_GET_TEST_DATA_PATH("BSpline2d_test_sinus.txt")))
  TEST_EXCEPTION(Exception::InvalidValue, txtFileToMS2Data(OPENMS_GET_TEST_DATA_PATH("SILAC_exception_test.txt")))
  TEST_EXCEPTION(Exception::FileNotFound, txtFileToMS2Data("fkwjks.txt"))
}
END_SECTION

SILACDetector test;
std::vector<MS2Data> example_data = txtFileToMS2Data(OPENMS_GET_TEST_DATA_PATH("SILACTestData.txt"));
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
  std::vector<MS2Data> output = txtFileToMS2Data(OPENMS_GET_TEST_DATA_PATH("SILACstoreTest.txt"));
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
  myfile.load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_6_uncompressed.mzML"),experiment);
  TEST_EXCEPTION(Exception::InvalidValue, test.msExperimentToMS2Data(experiment))
}
END_SECTION

START_SECTION(std::vector<MS2Data> featureMapToMS2Data(const FeatureMap& features) const)
{
  FeatureMap fm;
  TEST_EQUAL(test.featureMapToMS2Data(fm).empty(), true) // empty map -> empty result (no exception)

  Feature f1;
  f1.setRT(100.0);
  f1.setMZ(500.5);
  f1.setCharge(2);
  fm.push_back(f1);

  Feature f2;
  f2.setRT(200.0);
  f2.setMZ(600.25);
  f2.setCharge(0); // charge 0 is remapped to 2 (unsupported for SILAC); warns but does not throw
  fm.push_back(f2);

  std::vector<MS2Data> out = test.featureMapToMS2Data(fm); // one data point per feature
  TEST_EQUAL(out.size(), 2)
  TEST_REAL_SIMILAR(out[0].RT, 100.0)
  TEST_REAL_SIMILAR(out[0].mz, 500.5)
  TEST_EQUAL(out[0].charge, 2)
  TEST_REAL_SIMILAR(out[1].RT, 200.0)
  TEST_REAL_SIMILAR(out[1].mz, 600.25)
  TEST_EQUAL(out[1].charge, 2) // remapped from 0
}
END_SECTION

START_SECTION(std::vector<MS2Data> consensusMapToMS2Data(const ConsensusMap& consensus) const)
{
  ConsensusMap cm;
  TEST_EQUAL(test.consensusMapToMS2Data(cm).empty(), true) // empty map -> empty result (no exception)

  // consensus feature WITH subfeatures -> one data point per subfeature (the consensus feature's own
  // RT/mz/charge are deliberately different to ensure the subfeatures are used, not the parent)
  ConsensusFeature cf_sub;
  cf_sub.setRT(999.0);
  cf_sub.setMZ(999.0);
  cf_sub.setCharge(9);
  Feature sub1;
  sub1.setRT(10.0);
  sub1.setMZ(400.0);
  sub1.setCharge(2);
  Feature sub2;
  sub2.setRT(20.0);
  sub2.setMZ(405.0);
  sub2.setCharge(3);
  cf_sub.insert(0, sub1); // map index 0
  cf_sub.insert(1, sub2); // map index 1
  cm.push_back(cf_sub);

  // consensus feature WITHOUT subfeatures -> the consensus feature itself is one data point
  ConsensusFeature cf_bare;
  cf_bare.setRT(50.0);
  cf_bare.setMZ(700.0);
  cf_bare.setCharge(4);
  cm.push_back(cf_bare);

  std::vector<MS2Data> out = test.consensusMapToMS2Data(cm);
  TEST_EQUAL(out.size(), 3) // 2 subfeatures + 1 bare consensus feature
  // subfeatures are stored in a set ordered by map index, so sub1 (index 0) then sub2 (index 1)
  TEST_REAL_SIMILAR(out[0].RT, 10.0)
  TEST_REAL_SIMILAR(out[0].mz, 400.0)
  TEST_EQUAL(out[0].charge, 2)
  TEST_REAL_SIMILAR(out[1].RT, 20.0)
  TEST_REAL_SIMILAR(out[1].mz, 405.0)
  TEST_EQUAL(out[1].charge, 3)
  // the bare consensus feature contributes its own RT/mz/charge
  TEST_REAL_SIMILAR(out[2].RT, 50.0)
  TEST_REAL_SIMILAR(out[2].mz, 700.0)
  TEST_EQUAL(out[2].charge, 4)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST