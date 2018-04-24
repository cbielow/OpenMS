#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/QCMetricMap.h>
#include <vector>
using namespace OpenMS;
using namespace std;

MetricMap::~MetricMap()
{

}
bool MetricMap::isEmpty()
{
  return !isfilled_;
}  //shows if class Element has been filled

int MetricMap::size()
{
  return Columnlength_;
}

void MetricMap::pushMetaData(String in)
{
  Metadata_.push_back(in);
}  //writes in Metadata

vector<String> MetricMap::getMetaData()
{
  return Metadata_;
}         //returns all saved Metadata

void MetricMap::pushDataString(String Head, vector<String> Data)           //Saves Data with its Head. Example (Proteins,[ProteinA,ProteinB.ProteinC,...])
{
  if ( Columnlength_== 0 || Columnlength_== Data.size())
  {
        DataStrings_.insert(DataStrings_.end(),pair<String,vector<String>>(Head,Data));
        isfilled_ = isfilled_||true;
        Columnlength_ = Data.size();
  }
  else{
    cout << "MetricMap Error: Wrongly sized Datavector" << endl;
  }
}


void MetricMap::pushDataSize(String Head,vector<Size> Data)
{
  if ( Columnlength_==0 || Columnlength_==Data.size() )
  {
    DataStringNum_.insert(DataStringNum_.end(),pair<String,vector<Size>>(Head,Data));
    isfilled_ = isfilled_ || true;
    Columnlength_=Data.size();
  }
  else{
    cout << "MetricMap Error: Wrongly sized Datavector" << '\n';
  }
}

void MetricMap::pushDataFloat(String Head,vector<float> Data)
{
  if(Columnlength_ == 0 || Columnlength_ == Data.size() )
  {
    DataStringFloat_.insert(DataStringFloat_.end(), pair<String,vector<float>>(Head,Data));
    isfilled_=isfilled_ || true;
    Columnlength_ = Data.size();
  }
}
               //Saves Data with its Head (if Data = numbers). Bsp:(Length,[3,5,3,9,5,...])

vector<pair<String,String>> MetricMap::getHeads()//Returns all Heads. With Type of its Data(Sting/Size). I.e:([Proteins ,String],[Length,Size],[NumberOfProteins,Size],..)
{
  vector<pair<String,String>> output;
  for ( map<String,vector<String>>::const_iterator it = DataStrings_.begin() ; it != DataStrings_.end() ; it++ )
  {
      output.push_back(pair<String,String>(it->first,"String"));
  }
  for ( map<String,vector<Size>>::const_iterator it = DataStringNum_.begin() ; it != DataStringNum_.end() ; it++ )
  {
      output.push_back( pair<String,String>(it->first,"Size"));
  } 
  for (map<String,vector<float>>::const_iterator it = DataStringFloat_.begin() ; it != DataStringFloat_.end() ; it ++)
  {
      output.push_back(pair<String,String>(it->first,"Float"));
  }
  return output;
}

vector<Size> MetricMap::getSizesByHead(String WantedHead)
{
  return DataStringNum_[WantedHead];
}//gives one line of Data by its head(if Data is made out of numbers)

vector<String> MetricMap::getStringsByHead(String WantedHead)
{
  return DataStrings_[WantedHead];
}
vector<float> MetricMap::getFloatsByHead(String WantedHead)
{
  return DataStringFloat_[WantedHead];
}
