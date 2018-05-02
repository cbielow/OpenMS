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
#include <OpenMS/CONCEPT/QCProteinAndPeptideCount.h>
#include <OpenMS/CONCEPT/QCMS2IdentificationRate.h>
#include <OpenMS/CONCEPT/QCMBRalignment.h>
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
MzTabFile mzTabOutputFile;
MzTab mzTabOutput;
QCProteinAndPeptideCount ProtAndPepObj(CFiles_);
int papc = ProtAndPepObj.ProtAndPepCount( mzTabOutput);
//QCMS2IdentificationRate MS2IDRate(Idxml_);
//int mid = MS2IDRate.MS2IDRateidentifier_( mzTabOutput);

QCMBRalignment MBRAlign(FeatMapsMBR_);
int mbra = MBRAlign.MBRAlignment( mzTabOutput);
mzTabOutputFile.store(out_,mzTabOutput);
}

////////////////Metrik2: ....................../////////////////////////////////////
////////////////Metrik3: ....................../////////////////////////////////////
////////////////Metrik4: ....................../////////////////////////////////////
////////////////Metrik5: ....................../////////////////////////////////////
