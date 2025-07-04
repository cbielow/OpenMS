#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/SIMULATION/IonMobilitySimulation.h>
#include <OpenMS/SIMULATION/RawMSSignalSimulation.h>
#include <OpenMS/SYSTEM/ExternalProcess.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <QDir>
#include <fstream>
#include <iostream>

using namespace std;

// To Do: replace the temporary solution, when im2deep can handle sequences larger than  (search for "temporary solution" in this file AND in the
// Header-file to see what has to be changed) reference: https://github.com/CompOmics/IM2Deep/issues/10

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
    unit_(source.unit_),
    im2deep_combined_output_path_(
      source
        .im2deep_combined_output_path_) // temporary solution for im2deep. To Do: remove this line when im2deep can handle sequences larger than 60
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
    im2deep_combined_output_path_
      = source
          .im2deep_combined_output_path_; // temporary solution for im2deep. To Do: remove this line when im2deep can handle sequences larger than 60
  }
  return *this;
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

bool IonMobilitySimulation::isIM2DeepAvailable()
{
  QString exe = "im2deep";
  QStringList args;
  args << "--help";

  QString working_dir = QDir::currentPath();
  String error_msg = "im2deep was not found.";

  ExternalProcess im2deepCheck;
  ExternalProcess::RETURNSTATE result = im2deepCheck.run(exe, args, working_dir, false, error_msg);

  if (result == ExternalProcess::RETURNSTATE::SUCCESS) { return true; }
  else
  {
    throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  "IM2Deep is not available! Please install it via 'pip install im2deep' to enable Ion Mobility simulation.");
  }
}

void IonMobilitySimulation::setDefaultParams_()
{
  defaults_.setValue("IM_unit", "vssc", "Unit of ion mobility. vssc (= raw inverse reduced ion mobility array) or ccs (= collisional cross section)");
  defaults_.setValidStrings("IM_unit", {"vssc", "ccs"});
  defaultsToParam_();
}

void IonMobilitySimulation::updateMembers_()
{
  im2deep_input_path_ = File::getTemporaryFile();
  im2deep_output_path_ = File::getTemporaryFile();
  unit_ = param_.getValue("IM_unit").toString();
  im2deep_combined_output_path_
    = File::getTemporaryFile(); // temporary solution for im2deep. To Do: remove this line when im2deep can handle sequences larger than 60
}

void IonMobilitySimulation::run(const SimTypes::FeatureMapSim& features)
{
  createIM2DeepInputCSV(features);
  runIM2Deep();
  addsplit_indices(); // temporary solution for im2deep. To Do: remove this line when im2deep can handle sequences larger than 60
  saveIM2DeepOutput();
}


void IonMobilitySimulation::createIM2DeepInputCSV(const SimTypes::FeatureMapSim& features)
{

  im2deep_input_path_ = "/buffer/ag_bsc/student_data/mssim/jonnab00/Beispieldaten/MS_IM2Deep/IM2Deep_input.csv";
  std::ofstream file;
  file.open(im2deep_input_path_.c_str(), std::ios::out);

  if (! file.is_open())
  {
    std::cerr << "Fehler: Datei konnte nicht geöffnet werden!" << std::endl;
    return;
  }

  file << "seq,modifications,charge,CCS\n";
  /// temporary solution for im2deep ///
  int line_index = 1;
  for (const Feature& feat : features)
  {
    const PeptideIdentification& pi = feat.getPeptideIdentifications()[0];

    if (pi.getHits().empty()) continue;

    const PeptideHit& hit = pi.getHits()[0];
    const String sequence = hit.getSequence().toString();
    const int charge = hit.getCharge();
    std::vector<int> split_indices_temp;

    if (sequence.size() > 60)
    {
      int num_parts = (sequence.size() + 59) / 60; // +59 to round up
      int part_length = sequence.size() / num_parts;
      int charge_split = charge / num_parts;
      int remaining_charge = charge % num_parts;

      for (int i = 0; i < num_parts; ++i)
      {
        int start = i * part_length;
        int end = (i == num_parts - 1) ? sequence.size() : start + part_length;
        String subseq = sequence.substr(start, end - start);
        int part_charge = charge_split + (i < remaining_charge ? 1 : 0);

        file << subseq << ",," << part_charge << ",\n";
        split_indices_temp.push_back(line_index);
        ++line_index;
      }
      split_indices_.push_back(split_indices_temp);
    }
    else
    {
      file << sequence << ",," << charge << ",\n";
      ++line_index;
    }
  }
  /// temporary solution for im2deep ///

  /*
  // To Do: replace the "temporary solution for im2deep" above with the following, when im2deep can handle sequences larger than 60:

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
  */
}

void IonMobilitySimulation::runIM2Deep()
{
#ifdef _OPENMP
  Size thread_count = omp_get_max_threads();
#else
  Size thread_count = 1;
#endif

  QString exe = "im2deep";
  QStringList args;

  im2deep_output_path_ = "/buffer/ag_bsc/student_data/mssim/jonnab00/Beispieldaten/MS_IM2Deep/IM2Deep_output.csv";
  //  Path in Qstring umwandeln, damit als input für ExternalProcess geht
  if (unit_ == "ccs")
  {
    args << QString::fromStdString(im2deep_input_path_) << "-o"
         << QString::fromStdString(
              im2deep_output_path_); //<< "-c" << "/buffer/ag_bsc/student_data/mssim/jonnab00/IM2Deep/im2deep/reference_data/multi_reference_ccs.csv";
  }
  else if (unit_ == "vssc")
  {
    args << QString::fromStdString(im2deep_input_path_) << "-o" << QString::fromStdString(im2deep_output_path_) << "--n-jobs" << QString::number(thread_count)
         << "--ion-mobility"; //<< "-c" << "/buffer/ag_bsc/student_data/mssim/jonnab00/IM2Deep/im2deep/reference_data/multi_reference_ccs.csv";
  }
  QString working_dir = QDir::currentPath();
  String error_msg = "Beim Aufruf von IM2Deep ist etwas schiefgelaufen :(";

  ExternalProcess im2deepCall(stdoutCallback, stderrCallback);
  ExternalProcess::RETURNSTATE result = im2deepCall.run(exe, args, working_dir, true, error_msg);

  if (result == ExternalProcess::RETURNSTATE::SUCCESS) { OPENMS_LOG_INFO << "IM2Deep erfolgreich ausgeführt!\n"; }
  else { std::cerr << "Fehler beim Ausführen von IM2Deep: " << error_msg << std::endl; }
}

// temporary solution for im2deep
// To Do: this function can be removed, when im2deep can handle sequences larger than 60
void IonMobilitySimulation::addsplit_indices()
{
  std::ifstream input_file(im2deep_output_path_);
  if (! input_file.is_open()) { throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, im2deep_output_path_); }

  String header;
  std::getline(input_file, header);

  std::vector<String> lines;
  String line;
  while (std::getline(input_file, line))
  {
    lines.push_back(line);
  }

  im2deep_combined_output_path_ = "/buffer/ag_bsc/student_data/mssim/jonnab00/Beispieldaten/MS_IM2Deep/IM2Deep_combined_output.csv";
  std::ofstream combined_file(im2deep_combined_output_path_.c_str());
  combined_file << header << "\n";

  int current_line = 1; // 1-basiert (wegen Header)
  for (const auto& group : split_indices_)
  {
    int start = group.front();
    int end = group.back();

    // Unveränderte Zeilen davor übernehmen
    while (current_line < start && current_line <= (int)lines.size())
    {
      combined_file << lines[current_line - 1] << "\n";
      ++current_line;
    }

    // Zusammensetzen
    String combined_seq;
    int total_charge = 0;
    float total_ccs = 0.0;

    for (int idx : group)
    {
      if (idx < 1 || idx > (int)lines.size()) continue;

      std::vector<String> parts;
      lines[idx - 1].split(',', parts);
      if (parts.size() != 3) continue;

      size_t slash_pos = parts[0].find('/');
      if (slash_pos == String::npos) continue;

      combined_seq += parts[0].substr(0, slash_pos);
      int charge = parts[1].toInt();
      float mz = 0.0;

      // calculated mz from sequence
      String seq = parts[0].substr(0, slash_pos);
      AASequence sequence = AASequence::fromString(seq);
      mz = sequence.getMonoWeight() / charge;


      float im_value = parts[2].toFloat();;

      // convert ccs to vssc to approximate the ion mobility value by adding
      if (unit_ == "vssc")
      {
        im_value = IonMobilitySimulation::convertVSSCToCCS(im_value, mz, charge);
      }
      total_charge += parts[1].toInt();
      total_ccs += im_value;
    }
    if (unit_ == "vssc")
    {
      AASequence final_sequence = AASequence::fromString(combined_seq);
      double mz = final_sequence.getMonoWeight() / total_charge;
      total_ccs = IonMobilitySimulation::convertCCStoKo(total_ccs, mz, total_charge);
    }

    combined_file << combined_seq << "/" << total_charge << "," << total_charge << "," << total_ccs << "\n";
    current_line = end + 1;
  }

  // Restliche Zeilen übernehmen
  while (current_line <= (int)lines.size())
  {
    combined_file << lines[current_line - 1] << "\n";
    ++current_line;
  }

  combined_file.close();
}

void IonMobilitySimulation::saveIM2DeepOutput()
{
  ionmobility_map_.clear();

  // temporary solution for im2deep
  // To Do: use the im2deep_output_path_ instead of the im2deep_combined_ouput_path when im2deep can handle sequences larger than 60

  // std::ifstream input_file(im2deep_output_path_);
  // if (! input_file.is_open()) { throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, im2deep_output_path_); }
  std::ifstream input_file(im2deep_combined_output_path_);
  if (! input_file.is_open()) { throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, im2deep_combined_output_path_); }

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
}

// function based on function from GitHub
// source: https://github.com/OpenMS/OpenMS/issues/6685
float IonMobilitySimulation::convertVSSCToCCS(float vssc, float mz, int charge)
{
  const float bruker_CCS_coef = 1059.62245; // constant coefficient for Bruker in the Mason-Schamp equation
  const float IM_N2_gas_mass = 28;

  float IM = vssc;
  float q = charge;
  float mass = mz * q;
  float reduced_mass = mass * IM_N2_gas_mass / (mass + IM_N2_gas_mass);
  float CCS = IM * q * bruker_CCS_coef / std::sqrt(reduced_mass); // Mason-Schamp equation
  return CCS;
}

// convert CCS to inverseK0
float IonMobilitySimulation::convertCCStoKo(float ccs, float mz, int charge)
{
  const float bruker_CCS_coef = 1059.62245; // Bruker-specific Factor
  const float IM_N2_gas_mass = 28.0;

  if (ccs <= 0.0 || mz <= 0.0 || charge == 0) return -1.0;

  float mass = mz * charge;
  float reduced_mass = (mass * IM_N2_gas_mass) / (mass + IM_N2_gas_mass);
  float k0 = (ccs * std::sqrt(reduced_mass)) / (charge * bruker_CCS_coef);
  return k0;
}

} // namespace OpenMS
