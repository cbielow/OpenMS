#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/SIMULATION/IonMobilitySimulation.h>
#include <OpenMS/SYSTEM/ExternalProcess.h>
#include <OpenMS/SYSTEM/File.h>
#include <fstream>
#include <iostream>

using namespace std;

namespace OpenMS
{

IonMobilitySimulation::IonMobilitySimulation(): 
  DefaultParamHandler("IonMobilitySimulation")
{
  setDefaultParams_();
}

IonMobilitySimulation::IonMobilitySimulation(const IonMobilitySimulation& source):
    DefaultParamHandler(source),
    input_path_(source.input_path_),
    output_path_(source.output_path_),
    ccs_map_(source.ccs_map_)
{
}

IonMobilitySimulation::~IonMobilitySimulation() = default;

IonMobilitySimulation& IonMobilitySimulation::operator=(const IonMobilitySimulation& source)
{
  if (this != &source)
  {
    DefaultParamHandler::operator=(source);
    input_path_ = source.input_path_;
    output_path_ = source.output_path_;
    ccs_map_ = source.ccs_map_;
  }
  return *this;
}

void IonMobilitySimulation::setDefaultParams_()
{
  defaults_.setValue("input_path", "/buffer/ag_bsc/student_data/mssim/jonnab00/Beispieldaten/MS_IM2Deep/IM2Deep_input.csv", "Pfad zur Input-CSV-Datei für IM2Deep.");
  defaults_.setValue("output_path", "/buffer/ag_bsc/student_data/mssim/jonnab00/Beispieldaten/MS_IM2Deep/IM2Deep_output.csv", "Pfad zur Output-CSV-Datei für IM2Deep.");
  defaults_.setValue("IM2Deep_working_dir", "/buffer/ag_bsc/student_data/mssim/jonnab00/Beispieldaten/MS_IM2Deep", "Arbeitsverzeichnis für IM2Deep.");
  defaultsToParam_();
}

void IonMobilitySimulation::updateMembers_()
{
  input_path_ = param_.getValue("input_path").toString();
  output_path_ = param_.getValue("output_path").toString();
  im2deep_working_dir_ = param_.getValue("IM2Deep_working_dir").toString();
}

void IonMobilitySimulation::run()
{
  createIM2DeepInputCSV();
  runIM2Deep();
  saveIM2DeepOutput();
}


void IonMobilitySimulation::createIM2DeepInputCSV()
{
  std::ofstream file;
  file.open(input_path_.c_str(), std::ios::out);

  if (! file)
  {
    std::cerr << "Fehler: Datei konnte nicht geöffnet werden!" << std::endl;
    return;
  }

  file << "seq,modifications,charge,CCS\n";

  for (const Feature& feat : *feature_map_)
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
  // Path in Qstring umwandeln, damit als input für ExternalProcess geht
  args << QString::fromStdString(input_path_) << "-o" << QString::fromStdString(output_path_);

  QString working_dir = QString::fromStdString(im2deep_working_dir_);
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
  ccs_map_.clear();

  std::ifstream input_file(output_path_);
  if (! input_file.is_open()) { throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, output_path_); }

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
    double ccs = parts[2].toDouble();

    ccs_map_[{seq, charge}] = ccs;
  }

  // JB löschen von file, weil temporary datei nicht automat. gelöscht wird
  /*if (File::exists("/buffer/ag_bsc/student_data/mssim/jonnab00/Beispieldaten/MS_IM2Deep/IM2Deep_output.csv"))
  {
    File::remove("/buffer/ag_bsc/student_data/mssim/jonnab00/Beispieldaten/MS_IM2Deep/IM2Deep_output.csv");
  }*/
}

}
