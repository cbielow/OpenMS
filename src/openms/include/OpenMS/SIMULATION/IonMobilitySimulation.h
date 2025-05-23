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
  String im2deep_input_path_;
  String im2deep_output_path_;
  std::map<std::pair<String, int>, double> ionmobility_map_;
  String unit_;
  std::vector<std::vector<int>> split_indices_; // Indices of split sequences in im2deep (temporary till im2deep can handle sequences larger than 60)
  String im2deep_combined_output_path_; // temporary solution for im2deep. To Do: remove this line when im2deep can handle sequences larger than 60

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

  /// temporary solution for im2deep remove this function when im2deep can handle sequences larger than 60
  /// Add split indices together from peptides > 60 for im2deep
  void addsplit_indices();

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