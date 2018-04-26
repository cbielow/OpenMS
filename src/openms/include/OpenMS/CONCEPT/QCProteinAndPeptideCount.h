#ifndef OPENMS_CONCEPT_ProteinAndPeptideCount_H
#define OPENMS_CONCEPT_ProteinAndPeptideCount_H
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/CONCEPT/QCMetricMap.h>

class OPENMS_DLLAPI QCProteinAndPeptideCount
{
  public:
    QCProteinAndPeptideCount(std::vector<std::pair<OpenMS::String,OpenMS::CsvFile>> files)://MetricMap& outPep, MetricMap& outProt, vector<pair<string,CsvFile>> files):
      CFile_(files)
      {
      };
    ~QCProteinAndPeptideCount();
    int ProtAndPepCount(MetricMap&,MetricMap&) const;
  protected:
    std::vector<std::pair<OpenMS::String,OpenMS::CsvFile>> CFile_;
};
#endif

