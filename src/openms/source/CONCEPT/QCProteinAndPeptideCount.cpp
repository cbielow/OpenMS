#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/QCProteinAndPeptideCount.h>

using namespace OpenMS;
using namespace std;

QCProteinAndPeptideCount::~QCProteinAndPeptideCount()
{

}

int QCProteinAndPeptideCount::ProtAndPepCount(MetricMap& outPep ,MetricMap& outProt) const
{
vector<CsvFile> CsvFiles;
for(vector<pair<String,CsvFile>>::const_iterator it = CFile_.begin();it!=CFile_.end();++it)
{
  if(it->first=="ProteinQuantifier"){
    CsvFiles.push_back(it->second);
  }
}
for(vector<CsvFile>::const_iterator it = CsvFiles.begin(); it!=CsvFiles.end();it++)
{
  String a;
  StringList MetaList;
  StringList DataList;
  bool headfinder = false;
  CsvFile fl = *it;
  Size line = 0;
  StringList CurrentRow;
  Size maxRow = fl.rowCount();
  while(fl.getRow(line,CurrentRow)==false)
  {
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
      for(Size l = 0; l<MetaList.size();l++)
      {
        a=="peptide"?outPep.pushMetaData(MetaList[l]):outProt.pushMetaData(MetaList[l]);
      }
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
  a == "peptide" ? ( headfinder == true ? outPep.pushDataString( a , DataList ) : DataList.clear()) : ((headfinder ==true ? outProt.pushDataString(a,DataList) : DataList.clear()));
  }
  return outPep.isEmpty()&&outProt.isEmpty()?0:1;
}
