#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/SIMULATION/IonMobilitySimulation.h>
#include <OpenMS/SYSTEM/ExternalProcess.h>
#include <OpenMS/SYSTEM/File.h>
#include <QDir>
#include <fstream>
#include <iostream>

using namespace std;

namespace OpenMS
{

IonMobilitySimulation::IonMobilitySimulation(): DefaultParamHandler("IonMobilitySimulation")
{
  setDefaultParams_();
  updateMembers_();
}

IonMobilitySimulation::IonMobilitySimulation(const IonMobilitySimulation& source):
    DefaultParamHandler(source),
    im2deep_input_path_(source.im2deep_input_path_),
    im2deep_output_path_(source.im2deep_output_path_),
    ionmobility_map_(source.ionmobility_map_),
    unit_(source.unit_)
{
}

IonMobilitySimulation::~IonMobilitySimulation() = default;

IonMobilitySimulation& IonMobilitySimulation::operator=(const IonMobilitySimulation& source)
{
  if (this != &source)
  {
    DefaultParamHandler::operator=(source);
    im2deep_input_path_ = source.im2deep_input_path_;
    im2deep_output_path_ = source.im2deep_output_path_;
    ionmobility_map_ = source.ionmobility_map_;
    unit_ = source.unit_;
  }
  return *this;
}

void IonMobilitySimulation::setDefaultParams_()
{
  defaults_.setValue("IM_unit", "k0", "Unit of ion mobility. Raw inverse = k0 or collisional cross section = CCS)");
  defaults_.setValidStrings("IM_unit", {"k0", "ccs"});
  defaultsToParam_();
}

void IonMobilitySimulation::updateMembers_()
{
  im2deep_input_path_ = File::getTemporaryFile();
  im2deep_output_path_ = File::getTemporaryFile();
  unit_ = param_.getValue("IM_unit").toString();
}

void IonMobilitySimulation::run(const SimTypes::FeatureMapSim& features)
{
  createIM2DeepInputCSV(features);
  runIM2Deep();
  saveIM2DeepOutput();
}


void IonMobilitySimulation::createIM2DeepInputCSV(const SimTypes::FeatureMapSim& features)
{

  // im2deep_input_path_ = "/buffer/ag_bsc/student_data/mssim/jonnab00/Beispieldaten/MS_IM2Deep/IM2Deep_input.csv";
  std::ofstream file;
  file.open(im2deep_input_path_.c_str(), std::ios::out);

  if (! file.is_open())
  {
    std::cerr << "Fehler: Datei konnte nicht geöffnet werden!" << std::endl;
    return;
  }

  file << "seq,modifications,charge,CCS\n";

  // for (const Feature& feat : *feature_map_)
  for (const Feature& feat : features)
  {
    const PeptideIdentification& pi = feat.getPeptideIdentifications()[0];

    if (! pi.getHits().empty())
    {
      const PeptideHit& hit = pi.getHits()[0];
      const String sequence = hit.getSequence().toString();
      const int charge = hit.getCharge();
      file << sequence << ",," << charge << ",\n";
    }
  }
  file.close();
}

// Callback helpers
void stdoutCallback(const String& output)
{
  std::cout << "stdout: " << output << std::endl;
}

void stderrCallback(const String& output)
{
  std::cerr << "stderr: " << output << std::endl;
}

void IonMobilitySimulation::runIM2Deep()
{
  QString exe = "im2deep";
  QStringList args;

  // im2deep_output_path_ = "/buffer/ag_bsc/student_data/mssim/jonnab00/Beispieldaten/MS_IM2Deep/IM2Deep_output.csv";
  //  Path in Qstring umwandeln, damit als input für ExternalProcess geht
  if (unit_ == "ccs")
  {
    args << QString::fromStdString(im2deep_input_path_) << "-o"
         << QString::fromStdString(
              im2deep_output_path_); //<< "-c" << "/buffer/ag_bsc/student_data/mssim/jonnab00/IM2Deep/im2deep/reference_data/multi_reference_ccs.csv";
  }
  else if (unit_ == "k0")
  {
    args << QString::fromStdString(im2deep_input_path_) << "-o" << QString::fromStdString(im2deep_output_path_)
         << "--ion-mobility"; //<< "-c" << "/buffer/ag_bsc/student_data/mssim/jonnab00/IM2Deep/im2deep/reference_data/multi_reference_ccs.csv";
  }
  QString working_dir = QDir::currentPath();
  String error_msg = "Beim Aufruf von IM2Deep ist etwas schiefgelaufen :(";

  ExternalProcess im2deepCall(stdoutCallback, stderrCallback);
  ExternalProcess::RETURNSTATE result = im2deepCall.run(exe, args, working_dir, true, error_msg);

  if (result == ExternalProcess::RETURNSTATE::SUCCESS) { OPENMS_LOG_INFO << "IM2Deep erfolgreich ausgeführt!\n"; }
  else { std::cerr << "Fehler beim Ausführen von IM2Deep: " << error_msg << std::endl; }

  // löschen von file, weil temporary datei nicht automat. gelöscht wird
  /*if (File::exists("/buffer/ag_bsc/student_data/mssim/jonnab00/Beispieldaten/MS_IM2Deep/IM2Deep_input.csv"))
  {
    File::remove("/buffer/ag_bsc/student_data/mssim/jonnab00/Beispieldaten/MS_IM2Deep/IM2Deep_input.csv");
  }*/
}

void IonMobilitySimulation::saveIM2DeepOutput()
{
  ionmobility_map_.clear();

  std::ifstream input_file(im2deep_output_path_);
  if (! input_file.is_open()) { throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, im2deep_output_path_); }

  String line;
  std::getline(input_file, line); // skip header

  while (std::getline(input_file, line))
  {
    std::vector<String> parts;
    String(line).split(',', parts);
    if (parts.size() != 3) continue;

    String mod_seq = parts[0];
    String seq = mod_seq.hasSubstring("/") ? mod_seq.prefix('/') : mod_seq;

    int charge = parts[1].toInt();
    double ionmobility = parts[2].toDouble();

    ionmobility_map_[{seq, charge}] = ionmobility;
  }

  // JB löschen von file, weil temporary datei nicht automat. gelöscht wird
  /*if (File::exists("/buffer/ag_bsc/student_data/mssim/jonnab00/Beispieldaten/MS_IM2Deep/IM2Deep_output.csv"))
  {
    File::remove("/buffer/ag_bsc/student_data/mssim/jonnab00/Beispieldaten/MS_IM2Deep/IM2Deep_output.csv");
  }*/
}

/*
// JB Funktion von GitHub
void convertVSSCToCCS(MSExperiment& spectra)
{
  OPENMS_LOG_INFO << "Converting 1/k0 to CCS values." << std::endl;
  const double bruker_CCS_coef = 1059.62245; // constant coefficient for Bruker in the Mason-Schamp equation
  const double IM_N2_gas_mass = 28;

  for (auto& s : spectra)
  {
    double IM = s.getDriftTime();
    double mz = s.getPrecursors()[0].getMZ();
    double charge = s.getPrecursors()[0].getCharge();
    double mass = mz * charge;
    double reduced_mass = mass * IM_N2_gas_mass / (mass + IM_N2_gas_mass);
    double CCS = IM * charge * bruker_CCS_coef / std::sqrt(reduced_mass); // Mason-Schamp equation
    s.setDriftTime(CCS);
  }
}
*/

// convert CCS to inverseK0
float IonMobilitySimulation::convertCCStoKo(float ccs, float mz, int charge)
{
  const float bruker_CCS_coef = 1059.62245; // Bruker-spezifischer Faktor
  const float IM_N2_gas_mass = 28.0;

  if (ccs <= 0.0 || mz <= 0.0 || charge == 0) return -1.0;

  float mass = mz * charge;
  float reduced_mass = (mass * IM_N2_gas_mass) / (mass + IM_N2_gas_mass);
  float k0 = (ccs * std::sqrt(reduced_mass)) / (charge * bruker_CCS_coef);
  return k0;
}

} // namespace OpenMS
