#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/QCProteinAndPeptideCount.h>
#include <boost/regex.hpp>

using namespace OpenMS;
using namespace std;

QCProteinAndPeptideCount::~QCProteinAndPeptideCount(){

}

int QCProteinAndPeptideCount::ProtAndPepCount( MzTab& mztab)//MetricMap& outPep ,MetricMap& outProt) const
{
vector<CsvFile> CsvFiles;
Size control;
for(vector<pair<String,CsvFile>>::const_iterator it = CFile.begin();it!=CFile.end();++it)
{
  if(it->first=="ProteinQuantifier"){
    CsvFiles.push_back(it->second);
  }
}
for(vector<CsvFile>::const_iterator it = CsvFiles.begin(); it!=CsvFiles.end();it++)
{
  String a;
  StringList MetaList;
  MzTabPeptideSectionRows PepROWS;
  MzTabProteinSectionRows ProtROWS;
  StringList DataList;
  bool headfinder = false;
  CsvFile fl = *it;
  Size line = 0;
  StringList CurrentRow;
  Size maxRow = fl.rowCount();
  vector<String> rafiles;
  while(fl.getRow(line,CurrentRow)==false)
  {
    boost::regex rgx("Rawfiles");
    boost::smatch match;
    bool found = boost::regex_search(CurrentRow[0],match,rgx);
    if(found)
    {
      rafiles.push_back(CurrentRow[0]);
    }
    MetaList.push_back(CurrentRow[0]);
    line++;
  }
  while(line<maxRow)
  {
    fl.getRow(line,CurrentRow);
    if(((CurrentRow[0]==("\"peptide\"")||("\"protein\""))&&(headfinder==false)))
    {
      headfinder = true;
      a = CurrentRow[0]=="\"peptide\""?"peptide":"protein";
    }
    else if(headfinder==true)
    {
      DataList.push_back(CurrentRow[0]);
    }
    if(!headfinder)
    {
      line = maxRow;
    }
    line++;
  }
  control = DataList.size();
  for(StringList::const_iterator it=DataList.begin(); it != DataList.end(); ++it)
  {
    if(a=="peptide")
    {
      MzTabPeptideSectionRow PepROW;
      MzTabString PepSeq;
      MzTabString MTDiff;
      MzTabOptionalColumnEntry mTOCE = make_pair("Match_Time_Difference",MTDiff);
      vector<MzTabOptionalColumnEntry> optionals;
      optionals.push_back(mTOCE);
      PepROW.opt_= optionals;
      PepSeq.set(*it);
      PepROW.sequence = PepSeq;
      PepROWS.push_back(PepROW);
    }
    else
    {
      MzTabProteinSectionRow ProtROW;
      MzTabString ProtSeq;
      MzTabString MTDiff;
      MzTabOptionalColumnEntry mTOCE = make_pair("Match_Time_Difference",MTDiff);
      vector<MzTabOptionalColumnEntry> optionals;
      optionals.push_back(mTOCE);
      ProtROW.opt_= optionals;
      ProtSeq.set(*it);
      ProtROW.description = ProtSeq;
      ProtROWS.push_back(ProtROW);
    }
  }
  mztab.setPeptideSectionRows(PepROWS);
  mztab.setProteinSectionRows(ProtROWS);
}
return control!=0 ? 1:0;
}
