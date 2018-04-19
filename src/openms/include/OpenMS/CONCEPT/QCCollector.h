// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Anton Haberland, Leo Wurth, Mohammad El-Ismail$
// $Authors: Anton Haberland, Leo Wurth, Mohammad El-Ismail $
// --------------------------------------------------------------------------
#include <vector>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/MzTab.h>

#include <regex>




    //der Plan für MetricMap ist eine Art Tabelle. Wenn ihr Daten hinein pusht dann wird die Überschrift der Tabellenspalte und die Daten
    //    der ganzen Spalte erwartet. In einer MetricMap muessen alle Datenreihen gleich groß sein(Falls das Probleme bringt muss ich
    //    vielleicht was aendern).
    //alternativ vielleicht schon in MzTab Datenstrukturen speichern,
    class OPENMS_DLLAPI MetricMap
    {
    vector<String> Metadata;
    map<String,vector<String>> DataStrings;
    map<String,vector<Size>> DataStringNum;
    bool isfilled;
    Size Columnlength;
    public:
        MetricMap():
            Metadata(),
            DataStrings(),
            DataStringNum(),
            isfilled(false),
            Columnlength(0)
            {
            }
        //important fuctions for metric map class
        bool isEmpty(){return !isfilled;};  //shows if class Element has been filled

        int size(){return Columnlength;};

        void pushMetaData(String in){Metadata.push_back(in);};  //writes in Metadata

        vector<String> getMetaData(){return Metadata;};         //returns all saved Metadata

        void pushDataString(String, vector<String>);            //Saves Data with its Head. Example (Proteins,[ProteinA,ProteinB.ProteinC,...])

        void pushDataSize(String,vector<Size>);                 //Saves Data with its Head (if Data = numbers). Bsp:(Length,[3,5,3,9,5,...])

        vector<pair<String,String>> getHeads();//Returns all Heads. With Type of its Data(Sting/Size). I.e:([Proteins ,String],[Length,Size],[NumberOfProteins,Size],..)

	      vector<Size> getSizesByHead(String WantedHead){return DataStringNum[WantedHead];};//gives one line of Data by its head(if Data is made out of numbers)

        vector<String> getStringsByHead(String WantedHead){return DataStrings[WantedHead];};//gives one line of Data by its head(if Data is made out of Strings)
    };
        void MetricMap::pushDataString(String Head, vector<String> Data){
            if(Columnlength==0||Columnlength==Data.size()){
                  DataStrings.insert(DataStrings.end(),pair<String,vector<String>>(Head,Data));
                  isfilled = isfilled||true;
                  Columnlength = Data.size();
            }
            else{
              cout << "MetricMap Error: Wrongly sized Datavector" << endl;
            }
        }
        void MetricMap::pushDataSize(String Head, vector<Size> Data){
            if(Columnlength==0||Columnlength==Data.size()){
              DataStringNum.insert(DataStringNum.end(),pair<String,vector<Size>>(Head,Data));
              isfilled = isfilled||true;
              Columnlength=Data.size();
            }
            else{
              cout << "MetricMap Error: Wrongly sized Datavector" << '\n';
            }
        }
        vector<pair<String,String>> MetricMap::getHeads(){
            vector<pair<String,String>> output;
            for(map<String,vector<String>>::const_iterator it = DataStrings.begin();it!=DataStrings.end();it++){
                output.push_back(pair<String,String>(it->first,"String"));
            }
            for(map<String,vector<Size>>::const_iterator it = DataStringNum.begin(); it != DataStringNum.end();it++){
                output.push_back(pair<String,String>(it->first,"Size"));
            }
            return output;
        }



    class OPENMS_DLLAPI Metriken
    {
    //das sind die eingelesenen Daten, sie können von den Metriken gelesen aber nicht umgeschrieben werden//
    const vector<FeatureMap> FeatureMaps;			//Alle FeatureXML Datein
    const vector<vector<PeptideIdentification>> PeptidesId;	//Peptide der IDXML's
    const vector<vector<ProteinIdentification>> ProteinsId;	//Proteine der IDXML's
    const vector<CsvFile> CsvFiles;   				//Alle CSV Datein
    const vector<ConsensusMap> ConsensusMaps;    //Alle ConsensusXMLFiles

    public:
        Metriken(const vector<FeatureMap> FM, const vector<vector<PeptideIdentification>> PeI, const vector<vector<ProteinIdentification>> PrI,const vector<CsvFile> CF, const vector<ConsensusMap> CM):
        FeatureMaps(FM),
        PeptidesId(PeI),
        ProteinsId(PrI),
        CsvFiles(CF),
        ConsensusMaps(CM)
        {
        }
        void runAllMetrics();
    protected:
        int ProteinAndPeptideCount_(MetricMap&) const;
	//->Hier Metriken deklarieren<-//
	};



    	//->Hier die Metriken definieren<-//


    Int Metriken::ProteinAndPeptideCount_(MetricMap& outMap) const{
            for(vector<CsvFile>::const_iterator it = CsvFiles.begin(); it!=CsvFiles.end();it++){
              String a;
              StringList MetaList;
              StringList DataList;
              bool headfinder = false;
              CsvFile fl = *it;
              Size line = 0;
              StringList CurrentRow;
              Size maxRow = fl.rowCount();
              while(fl.getRow(line,CurrentRow)==false){
                MetaList.push_back(CurrentRow[0]);
                line++;
              }
            while(line<maxRow){
                fl.getRow(line,CurrentRow);
                if(((CurrentRow[0]=="\"peptide\"")||(CurrentRow[0]=="\"protein\""))&&(headfinder==false)){
                    headfinder = true;
                    for(int l = 0; l<MetaList.size();l++){
                        outMap.pushMetaData(MetaList[l]);
                    }
                    a = CurrentRow[0]=="\"peptide\""?"peptide":"protein";
                }
                else if(headfinder==true){
                    DataList.push_back(CurrentRow[0]);
                }
                if(!headfinder){
                    line = maxRow;
                }
                line++;
            }
            headfinder==true ? outMap.pushDataString(a,DataList): DataList.clear();
        }
        return outMap.isEmpty()?0:1;
    }



        //Wenn ihr die Metriken schreibt lasst euch bitte ein Int ausgeben 1 wenn erfolgreich, 0  wenn nicht
        //Die Daten die ihr erhaltet am besten als MetricMap, die wichtigen Funktionen dafür stehen oben
    void Metriken::runAllMetrics(){
        ////////////////Metrik1: Protein And Peptide Count /////////////////////////////////
        MetricMap ProteinAndPeptideCountData;
        int a = this->ProteinAndPeptideCount_(ProteinAndPeptideCountData);
	      ////////////////Metrik2: ....................../////////////////////////////////////
        ////////////////Metrik3: ....................../////////////////////////////////////
        ////////////////Metrik4: ....................../////////////////////////////////////
        ////////////////Metrik5: ....................../////////////////////////////////////

        //MzTab Writer:
        if(a == 1){
          MzTab MzTabAusgabe;
          int numOfPeptides = ProteinAndPeptideCountData.size();
          vector<String> allPeptides = ProteinAndPeptideCountData.getStringsByHead("peptide");
          MzTabPeptideSectionRows ROWS;
          for(int i = 0; i< numOfPeptides;i++){
            MzTabPeptideSectionRow ROW;
            MzTabString Seq;
            Seq.set(allPeptides[i]);
            ROW.sequence = Seq;
            ROWS.push_back(ROW);
          }
          MzTabAusgabe.setPeptideSectionRows(ROWS);
          MzTabFile().store("ausgabe1",MzTabAusgabe);
        }




    }
