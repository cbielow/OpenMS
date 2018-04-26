#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FILTERING/ID/IDFilter.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/DATASTRUCTURES/String.h>



using namespace OpenMS;
using namespace std;

#include <OpenMS/CONCEPT/QCCollector.h>

class QCCollector:
	public TOPPBase
{
public:
	QCCollector():
	TOPPBase("QCCollector","Will collect several Files from several utils.",false)
	{
	}

protected:
	void registerOptionsAndFlags_()
	{
		registerInputFileList_("in","<files>", StringList(), "Input files",true);
		//registerStringOption_("in_type", "<type>", "", "Input file type -- default: determined from file extension or content\n", false, true); // for TOPPAS
		String formats("idXML,featureXML,csv");
		setValidFormats_("in", ListUtils::create<String>(formats));
	    	//setValidStrings_("in_type", ListUtils::create<String>(formats));
		//registerOutputFile_("MzTAB","<file>", "", "Transformation in MzTab");
		//setValidFormats_("MzTAB", ListUtils::create<String>("MzTAB"));
	}
	ExitCodes main_(int, const char**)
	{
		StringList ins = getStringList_("in");
		typedef vector<PeptideIdentification> VPepID;
		typedef	vector<ProteinIdentification> VProtID;
		vector<VPepID> IDXPeptides;
		vector<VProtID> IDXProteins;
		vector<FeatureMap> VFeatureMaps;
		vector<CsvFile> VecCSV;
		vector<ConsensusMap> CMapVec;
		for (Size i=0;i<ins.size();++i){
			FileTypes::Type in_type = FileHandler::getType(ins[i]);
			if (in_type == FileTypes::IDXML){
				VPepID pep_ids;
				VProtID prot_ids;
				IdXMLFile().load(ins[i], prot_ids, pep_ids);
				IDXPeptides.push_back(pep_ids);
				IDXProteins.push_back(prot_ids);
			}
			else if (in_type == FileTypes::CSV){
				CsvFile fl(ins[i],'	',false,-1);
				VecCSV.push_back(fl);

			}
			else if (in_type == FileTypes::FEATUREXML){
				FeatureMap features;
				FeatureXMLFile().load(ins[i], features);
				VFeatureMaps.push_back(features);
			}
			else if(in_type == FileTypes::CONSENSUSXML){
				ConsensusMap CMap;
				ConsensusXMLFile().load(ins[i],CMap);
				CMapVec.push_back(CMap);
			}
		}
		Metriken metricElemet(VFeatureMaps,IDXPeptides,IDXProteins,VecCSV,CMapVec);
		metricElemet.runAllMetrics();
		return EXECUTION_OK;
	}
};
int main(int argc, const char** argv)
{
QCCollector tool;
return tool.main(argc,argv);
}
