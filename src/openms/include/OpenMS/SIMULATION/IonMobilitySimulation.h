#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

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
  String input_path_;                                // z. B. "/path/to/im2deep/input.csv"
  String output_path_;                               // z. B. "/path/to/im2deep/output.csv"
  std::map<std::pair<String, int>, double> ccs_map_; // Map für CCS-Werte
  const FeatureMap* feature_map_;
  String im2deep_working_dir_; // Arbeitsverzeichnis für IM2Deep


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

  void setFeatureMap(const FeatureMap& fmap)
  {
    feature_map_ = &fmap;
  }

  // run IM2Deep and save ccs in map with peptide sequence and charge
  void run();

  /// Simulate ion mobility spectra for a given set of peptides
  void createIM2DeepInputCSV();

  /// Run IM2Deep to predict CCS values
  void runIM2Deep();

  /// Save IM2Deep output to a file
  void saveIM2DeepOutput();

  // CCS Map getter
  const std::map<std::pair<String, int>, double>& getCCSMap() const
  {
    return ccs_map_;
  }


  /// Add CCS values to the feature map
  /*
  void addCCSToFeature(std::vector<FeatureMap>& feature_maps, const std::map<std::pair<String, int>, double>& ccs_map);
  */
};

} // namespace OpenMS