#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/SIMULATION/IonMobilitySimulation.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(IonMobilitySimulation, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

IonMobilitySimulation* ptr = nullptr;
IonMobilitySimulation* nullPointer = nullptr;

// Constructor test
START_SECTION(IonMobilitySimulation())
{
  ptr = new IonMobilitySimulation();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

// Destructor test
START_SECTION(~IonMobilitySimulation())
{
  delete ptr;
}
END_SECTION

// Copy constructor test
START_SECTION((IonMobilitySimulation(const IonMobilitySimulation&)))
{
  IonMobilitySimulation source;
  source.setUnit("vssc");
  IonMobilitySimulation copy(source);
  TEST_EQUAL(copy.getParameters(), source.getParameters())
}
END_SECTION

// Assignment operator test
START_SECTION((IonMobilitySimulation & operator=(const IonMobilitySimulation&)))
{
  IonMobilitySimulation ims1;
  ims1.setIM2DeepInputPath("input_path_a");
  ims1.setIM2DeepOutputPath("output_path_a");
  std::map<std::pair<String, int>, double> im_map;
  im_map[{"PEPTIDE", 2}] = 111.1;
  ims1.setIonMobilityMap(im_map);
  ims1.setUnit("vssc");

  IonMobilitySimulation ims2;
  ims2 = ims1; // Assignment

  TEST_EQUAL(ims2.getIM2DeepInputPath(), "input_path_a")
  TEST_EQUAL(ims2.getIM2DeepOutputPath(), "output_path_a")
  TEST_REAL_SIMILAR(ims2.getIonMobilityMap().at({"PEPTIDE", 2}), 111.1)
  TEST_EQUAL(ims2.getUnit(), "vssc")
}
END_SECTION

START_SECTION((bool isIM2DeepAvailable()))
{
  IonMobilitySimulation ims;
  bool available = ims.isIM2DeepAvailable();
  TEST_EQUAL(available, true)
  if (! available) { OPENMS_LOG_WARN << "IM2Deep is not available. Skipping related tests. Install with 'pip install im2deep' to use this function" << std::endl; }
}
END_SECTION

START_SECTION((void createIM2DeepInputCSV(const SimTypes::FeatureMapSim& features)))
{
  IonMobilitySimulation ims;

  String out_path = File::getTemporaryFile();
  ims.setIM2DeepInputPath(out_path);

  SimTypes::FeatureMapSim feature_map;

  // add short sequence (<= 60)
  {
    Feature f;
    PeptideIdentification pid;
    PeptideHit hit;
    hit.setSequence(AASequence::fromString("PEPTIDE"));
    hit.setCharge(2);
    pid.insertHit(hit);
    f.getPeptideIdentifications().push_back(pid);
    feature_map.push_back(f);
  }

  // add long sequence (> 60) that will be split
  {
    Feature f;
    PeptideIdentification pid;
    PeptideHit hit;
    hit.setSequence(AASequence::fromString(String(65, 'A')));
    hit.setCharge(2);
    pid.insertHit(hit);
    f.getPeptideIdentifications().push_back(pid);
    feature_map.push_back(f);
  }

  ims.createIM2DeepInputCSV(feature_map);

  TextFile tf;
  tf.load(out_path);

  Size line_count = 0;
  std::vector<String> lines;
  for (const String& line : tf)
  {
    lines.push_back(line);
    ++line_count;
  }

  TEST_EQUAL(lines.size(), 4)
  TEST_STRING_EQUAL(lines[0], "seq,modifications,charge,CCS")
  TEST_STRING_EQUAL(lines[1], "PEPTIDE,,2,")

  String part1(std::string(32, 'A'));
  String part2(std::string(33, 'A'));

  TEST_EQUAL(lines[2].hasSubstring(part1), true)
  TEST_EQUAL(lines[2].hasSubstring(",,1,"), true)
  TEST_EQUAL(lines[3].hasSubstring(part2), true)
  TEST_EQUAL(lines[3].hasSubstring(",,1,"), true)
}
END_SECTION

// runIM2Deep test
START_SECTION((void runIM2Deep()))
{
  IonMobilitySimulation ims;
  if (! ims.isIM2DeepAvailable()) { OPENMS_LOG_WARN << "Skipping runIM2Deep test (IM2Deep not available)." << std::endl; }
  else
  {
    IonMobilitySimulation ims;

    String input = OPENMS_GET_TEST_DATA_PATH("im2deep_valid_input.csv");
    String output = File::getTemporaryFile();

    ims.setIM2DeepInputPath(input);
    ims.setIM2DeepOutputPath(output);
    ims.setUnit("vssc");

    ims.runIM2Deep();

    TEST_EQUAL(File::exists(output), true)

    TextFile tf;
    tf.load(output);

    Size line_count = 0;
    std::vector<String> lines;
    for (const String& line : tf)
    {
      lines.push_back(line);
      ++line_count;
    }

    TEST_EQUAL(lines[6].hasSubstring("YLQDYGMGPETPLGEPKNK/2,2,1.213545925586197"), true)
  }
}
END_SECTION

START_SECTION((void saveIM2DeepOutput()))
{
  IonMobilitySimulation ims;

  String test_output = OPENMS_GET_TEST_DATA_PATH("IM2Deep_valid_output.csv");
  ims.setIM2DeepCombinedOutputPath(test_output);

  ims.saveIM2DeepOutput();

  const auto& map = ims.getIonMobilityMap();

  TEST_EQUAL(map.size(), 29)
  TEST_REAL_SIMILAR(map.at({"YLQDYGMGPETPLGEPKNK", 2}), 1.213545925586197)
  TEST_REAL_SIMILAR(map.at({"EPIPVRPTAHYTMGGIETDQNCETRIK", 5}), 0.9112299265879632)
}
END_SECTION


// getCCSMap test
START_SECTION((const std::map<String, double>& getCCSMap() const))
{
  IonMobilitySimulation ims;
  std::map<std::pair<String, int>, double> im_map;
  im_map[{"PEPTIDE", 2}] = 123.45;
  ims.setIonMobilityMap(im_map);
  TEST_REAL_SIMILAR(ims.getIonMobilityMap().at({"PEPTIDE", 2}), 123.45)
}
END_SECTION

START_SECTION((float convertCCStoKo(float ccs, float mz, int charge)))
{
  IonMobilitySimulation ims;

  // Test: Normaler Wert
  float ccs = 200.0;
  float mz = 500.0;
  int charge = 2;

  float mass = mz * charge;
  float reduced_mass = (mass * 28.0f) / (mass + 28.0f);
  float expected_k0 = (ccs * std::sqrt(reduced_mass)) / (charge * 1059.62245f);

  TEST_REAL_SIMILAR(ims.convertCCStoKo(ccs, mz, charge), expected_k0)

  // Test: Fehlerhafte Eingaben
  TEST_REAL_SIMILAR(ims.convertCCStoKo(0.0, 500.0, 2), -1.0)
  TEST_REAL_SIMILAR(ims.convertCCStoKo(200.0, 0.0, 2), -1.0)
  TEST_REAL_SIMILAR(ims.convertCCStoKo(200.0, 500.0, 0), -1.0)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
