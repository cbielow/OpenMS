#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/SIMULATION/SimTypes.h>

namespace OpenMS
{
/**
  @brief Simulates ion mobility spectra for a given set of peptides

  This class simulates ion mobility spectra for a given set of peptides,
  with charge annotation, given detectabilities, predicted retention times
  and charge values.

  @htmlinclude OpenMS_IonMobilitySimulation.parameters

  @ingroup Simulation
*/
class OPENMS_DLLAPI IonMobilitySimulation : public DefaultParamHandler
{
private:
  String im2deep_input_path_;                                // z. B. "/path/to/im2deep/input.csv"
  String im2deep_output_path_;                               // z. B. "/path/to/im2deep/output.csv"
  std::map<std::pair<String, int>, double> ionmobility_map_; // Map für CCS-Werte
  String unit_;                                              // z. B. "VSSC" oder "CCS"


public:
  /** @name Constructors and Destructors
   */
  //@{
  /// Default constructor
  IonMobilitySimulation();


  /// Copy constructor
  IonMobilitySimulation(const IonMobilitySimulation& source);

  /// Destructor
  ~IonMobilitySimulation() override;
  //@}

  /// Assignment operator
  IonMobilitySimulation& operator=(const IonMobilitySimulation& source);

  /// Set default parameters
  void setDefaultParams_();

  // Save param_ values as members
  void updateMembers_();

  // run IM2Deep and save ccs in map with peptide sequence and charge
  void run(const SimTypes::FeatureMapSim& features);

  /// Simulate ion mobility spectra for a given set of peptides
  void createIM2DeepInputCSV(const SimTypes::FeatureMapSim& features);

  /// Run IM2Deep to predict Ionmobility values
  void runIM2Deep();

  /// Save IM2Deep output to a file
  void saveIM2DeepOutput();

  // IonMobility Map getter
  const std::map<std::pair<String, int>, double>& getIonMobilityMap() const
  {
    return ionmobility_map_;
  }

  // Unit getter
  const String& getUnit() const
  {
    return unit_;
  }

  // Convert CCS to inverse K0
  static float convertCCStoKo(float ccs, float mz, int charge);
};

} // namespace OpenMS