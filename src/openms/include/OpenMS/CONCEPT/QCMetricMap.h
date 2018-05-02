#ifndef OPENMS_CONCEPT_QCMETRICMAP_H
#define OPENMS_CONCEPT_QCMETRICMAP_H
#include<OpenMS/KERNEL/StandardTypes.h>
#include<vector>

//der Plan für MetricMap ist eine Art Tabelle. Wenn ihr Daten hinein pusht dann wird die Überschrift der Tabellenspalte und die Daten
//    der ganzen Spalte erwartet. In einer MetricMap muessen alle Datenreihen gleich groß sein(Falls das Probleme bringt muss ich
//    vielleicht was aendern).
//alternativ vielleicht schon in MzTab Datenstrukturen speichern,

class OPENMS_DLLAPI MetricMap
{
public:
    MetricMap():
        Metadata_(),
        DataStrings_(),
        DataStringFloat_(),
        DataStringNum_(),
        isfilled_(false),
        Columnlength_(0)
        {
        }
    ~MetricMap();
    //important fuctions for metric map class
    bool isEmpty();  //shows if class Element has been filled

    int size();

    void pushMetaData(OpenMS::String in);  //writes in Metadata

    std::vector<OpenMS::String> getMetaData();         //returns all saved Metadata

    void pushDataString(OpenMS::String, std::vector<OpenMS::String>);            //Saves Data with its Head. Example (Proteins,[ProteinA,ProteinB.ProteinC,...])

    void pushDataSize(OpenMS::String,std::vector<OpenMS::Size>);                 //Saves Data with its Head (if Data = numbers). Bsp:(Length,[3,5,3,9,5,...])

	void pushDataFloat(OpenMS::String,std::vector<float>);

    std::vector<std::pair<OpenMS::String,OpenMS::String>> getHeads();//Returns all Heads. With Type of its Data(Sting/Size). I.e:([Proteins ,String],[Length,Size],[NumberOfProteins,Size],..)

    std::vector<OpenMS::Size> getSizesByHead(OpenMS::String WantedHead);//gives one line of Data by its head(if Data is made out of numbers)
	
	std::vector<float> getFloatsByHead(OpenMS::String WantedHead);
	
    std::vector<OpenMS::String> getStringsByHead(OpenMS::String WantedHead);//gives one line of Data by its head(if Data is made out of Strings)
protected:
    std::vector<OpenMS::String> Metadata_;
    std::map<OpenMS::String,std::vector<OpenMS::String>> DataStrings_;
    std::map<OpenMS::String,std::vector<OpenMS::Size>> DataStringNum_;
    std::map<OpenMS::String,std::vector<float>> DataStringFloat_;
    bool isfilled_;
    OpenMS::Size Columnlength_;
};
#endif
