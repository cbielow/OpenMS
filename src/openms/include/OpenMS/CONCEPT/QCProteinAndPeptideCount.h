
#ifndef OPENMS_CONCEPT_QCProteinAndPeptideCount_H
#define OPENMS_CONCEPT_QCProteinAndPeptideCount_H
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <vector>

class OPENMS_DLLAPI QCProteinAndPeptideCount
{
std::vector<OpenMS::CsvFile> CsvFilesPeptide;
std::vector<OpenMS::CsvFile> CsvFilesProtein;
public:
QCProteinAndPeptideCount(std::vector<OpenMS::CsvFile> filespep, std::vector<OpenMS::CsvFile> filesprot):
  CsvFilesPeptide(filespep),
  CsvFilesProtein(filesprot)
      {
      };
    ~QCProteinAndPeptideCount();
    bool ProtAndPepCount(OpenMS::MzTab&);
};
#endif
