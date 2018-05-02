#ifndef OPENMS_CONCEPT_QCMS2IdentificationRate_H
#define OPENMS_CONCEPT_QCMS2IdentificationRate_H
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <vector>
#include <utility>

class OPENMS_DLLAPI QCMS2IdentificationRate
{
std::vector<std::pair<OpenMS::String,std::pair<OpenMS::String,OpenMS::String>>> ivec_;
  public:
    QCMS2IdentificationRate(std::vector<std::pair<OpenMS::String,std::pair<OpenMS::String,OpenMS::String>>> files):
      ivec_(files)
      {
      }
      ~QCMS2IdentificationRate();
      int MS2IDRateidentifier_(OpenMS::MzTabFile&, OpenMS::String);
};
#endif
