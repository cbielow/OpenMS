#ifndef OPENMS_CONCEPT_QCPeptideIntensity_H
#define OPENMS_CONCEPT_QCPeptideIntensity_H
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <vector>

class OPENMS_DLLAPI QCPeptideIntensity
{
  std::vector<OpenMS::CsvFile> CsvFilesPeptide;
  public:
    QCPeptideIntensity(std::vector<OpenMS::CsvFile> files):
      CsvFilesPeptide_(files)
      {
      };
    ~QCPeptideIntensity();
    bool PeptideIntensity(OpenMS::MzTab&);

  protected:
    std::vector<OpenMS::CsvFile> CsvFilesPeptide_;
    OpenMS::Size FindPeptideInMzTab_(const OpenMS::String&, const OpenMS::MzTab&);
};
#endif
