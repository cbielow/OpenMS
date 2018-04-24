#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/CONCEPT/QCMetrics.h>
#include <OpenMS/CONCEPT/QCMetricMap.h>
#include <OpenMS/CONCEPT/QCProteinAndPeptideCount.h>
#include <vector>
#include <utility>

using namespace OpenMS;
using namespace std;

Metrics::~Metrics()
 {

 }

//Wenn ihr die Metriken schreibt lasst euch bitte ein Int ausgeben 1 wenn erfolgreich, 0  wenn nicht
//Die Daten die ihr erhaltet am besten als MetricMap, die wichtigen Funktionen dafür stehen oben
void Metrics::runAllMetrics()
{
////////////////Metrik1: Protein And Peptide Count /////////////////////////////////
QCProteinAndPeptideCount papc(CFiles_);
MetricMap PeptideCountData;
MetricMap ProteinCountData;
int a = papc.ProtAndPepCount(PeptideCountData,ProteinCountData);
////////////////Metrik3: ....................../////////////////////////////////////
////////////////Metrik4: ....................../////////////////////////////////////
////////////////Metrik5: ....................../////////////////////////////////////

//MzTab Writer:
MzTabFile MzTabOutputFile;
MzTab mztab;
if(a == 1)
{
    int numOfPeptides = PeptideCountData.size();
    int numOfProteins = ProteinCountData.size();
    vector<String> empty;
    vector<String> allPeptides = numOfPeptides > 0 ? PeptideCountData.getStringsByHead("peptide") : empty;
    vector<String> allProteins = numOfProteins > 0 ? ProteinCountData.getStringsByHead("protein") : empty;
    MzTabPeptideSectionRows PepROWS;
    MzTabProteinSectionRows ProtROWS;
    if (numOfPeptides > 0)
    {
      for( int i = 0; i < numOfPeptides ; i++ )
      {
        MzTabPeptideSectionRow PepROW;
        MzTabString PepSeq;
        PepSeq.set(allPeptides[i]);
        PepROW.sequence = PepSeq;
        PepROWS.push_back(PepROW);
      }
    }
    if( numOfProteins > 0 ){
      for (int i = 0; i < numOfProteins ; i++)
      {
        MzTabProteinSectionRow ProtROW;
        MzTabString ProtSeq;
        ProtSeq.set(allProteins[i]);
        ProtROW.description = ProtSeq;
        ProtROWS.push_back(ProtROW);
      }
    }
    mztab.setPeptideSectionRows(PepROWS);
    mztab.setProteinSectionRows(ProtROWS);
}
MzTabOutputFile.store(out_,mztab);
}
