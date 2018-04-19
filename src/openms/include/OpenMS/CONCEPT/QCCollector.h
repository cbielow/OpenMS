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
#include <OpenMS/KERNEL/FeatureMap.h>


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
     //  die Vectoren enthalten alle eingelesenen Datein in den entsperechenden Formaten. Der erste Wert in jedem Paar
    //  gibt an aus welchem TOPPTOOL die Datei kommt.
    const vector<pair<string,FeatureMap>> FeatureMaps;			//Alle FeatureXML Datein
    const vector<pair<string,string>> Idxml;								//Peptide der IDXML's
    const vector<pair<string,CsvFile>> CFiles;							//Proteine der IDXML's;
    const vector<pair<string,ConsensusMap>> ConsensusMaps;  //Alle ConsensusXMLFiles
    
		public:
        Metriken(const vector<pair<string,FeatureMap>> fvec, const vector<pair<string,string>> ivec, const vector<pair<string,CsvFile>> cvec,const vector<pair<string,ConsensusMap>> CMapVec):
        FeatureMaps(fvec),
        Idxml(ivec),
        CFiles(cvec),
        ConsensusMaps(CMapVec)
        {
        }
        void runAllMetrics();
    protected:
        int ProteinAndPeptideCount_(MetricMap&,MetricMap&) const;
//---------------------------------------------------------------
				int MassCalibration(MetricMap&, MetricMap&, MetricMap&, MetricMap&, MetricMap&) const;
//---------------------------------------------------------------
	//->Hier Metriken deklarieren<-//
	};



    	//->Hier die Metriken definieren<-//


    Int Metriken::ProteinAndPeptideCount_(MetricMap& outPep,MetricMap& outProt) const{
      vector<CsvFile> CsvFiles;
      for(vector<pair<string,CsvFile>>::const_iterator it = CFiles.begin();it!=CFiles.end();++it){
        if(it->first=="ProteinQuantifier"){CsvFiles.push_back(it->second);}
      }
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
                  a = CurrentRow[0]=="\"peptide\""?"peptide":"protein";
                  for(int l = 0; l<MetaList.size();l++){
                      a=="peptide"?outPep.pushMetaData(MetaList[l]):outProt.pushMetaData(MetaList[l]);
                  }
              }
              else if(headfinder==true){
                  DataList.push_back(CurrentRow[0]);
              }
              if(!headfinder){
                  line = maxRow;
              }
              line++;
          }
          headfinder==true ? a=="peptide" ? outPep.pushDataString(a,DataList) : outProt.pushDataString(a,DataList) : DataList.clear();
      }
      return outPep.isEmpty()&&outProt.isEmpty()?0:1;
    }

//--------------------------------------------------------------------------------------------
		Int Metriken::MassCalibration(MetricMap& outRT, MetricMap& outMZT, MetricMap& outMZO, MetricMap& outPPMB, MetricMap& outPPMA) const
		{
			vector<CsvFile> CsvFiles;
			for (vector<pair<string,CsvFile>>::const_iterator it = CFiles.begin(); it != CFiles.end(); it++)
			{

				if (it->first == "InternalCalibration")
				{
					CsvFiles.push_back(it->second);
				}
			}

				for (vector<CsvFile>::const_iterator it = CsvFiles.begin(); it != CsvFiles.end(); it++)
				{
					
					String a;
					StringList MetaList;
					StringList DataList;
					bool headfinder = false;
					CsvFile fl = *it;
					Size line = 0;
					StringList CurrentRow;
					Size maxRow = fl.rowCount();
					
					while (fl.getRow(line,CurrentRow) == false)
					{
						MetaList.push_back(CurrentRow[0]);
						line++;
					}
		
					while (line < maxRow)
					{
						fl.getRow(line,CurrentRow);
					
						if (((CurrentRow[0] == "\"RT\"") || (CurrentRow[0] == "\"mz ref\"") || (CurrentRow[0] == "\"mz before\"") || (CurrentRow[0] == "\"ppm before\"") || (CurrentRow[0] == "\"ppm after\"")) && (headfinder == false))
						{
							headfinder = true;

							if (CurrentRow[0] == "\"RT\"")
							{
								a = "RT";
							}
						
							else if (CurrentRow[0] == "\"mz ref\"")
							{
								a = "mz ref";
							}

							else if (CurrentRow[0] == "\"mz before\"")
							{
								a = "mz before";
							}

							else if (CurrentRow[0] == "\"ppm before\"")
							{
								a = "ppm before";
							}

							else if (CurrentRow[0] == "\"ppm after\"")
							{
								a = "ppm after";
							}
						
							for (int l = 0; l < MetaList.size(); l++)
							{
								if (a == "RT")
								{
									outRT.pushMetaData(MetaList[l]);
								}
								if (a == "mz ref")
								{
									outMZT.pushMetaData(MetaList[l]);
								}
								if (a == "mz after")
								{
									outMZO.pushMetaData(MetaList[l]);
								}
								if (a == "ppm before")
								{
									outPPMB.pushMetaData(MetaList[l]);
								}
								if (a == "ppm after")
								{
									outPPMA.pushMetaData(MetaList[l]);
								}
							}
						}

						else if (headfinder == true)
						{
							DataList.push_back(CurrentRow[0]);
						}
					
						if (!headfinder)
						{
							line = maxRow;
						}
						line++;
					}
			
					if (headfinder == true)
					{
				
						if (a  == "RT")
						{
							outRT.pushDataString(a, DataList);	
						}
						else if (a  == "mz ref")
						{
							outMZT.pushDataString(a, DataList);	
						}
						else if (a  == "mz before")
						{
							outMZO.pushDataString(a, DataList);	
						}
						else if (a  == "ppm before")
						{
							outPPMB.pushDataString(a, DataList);	
						}
						else if (a  == "ppm after")
						{
							outPPMA.pushDataString(a, DataList);	
						}
					}
					else
					{
						DataList.clear();
					}
				std::sort(DataList.begin(), DataList.end());
				
				std::cout << DataList << std::endl;
				}
			
				return outRT.isEmpty() && outMZT.isEmpty() && outMZO.isEmpty() && outPPMB.isEmpty() && outPPMA.isEmpty() ? 0:1;	
			
		}
		
		/*vector<string> Idxmls;
		for (vector<pair<string, string>>::const_iterator it = Idxml.begin(); it! != Idxml.end(); it++)
		{

			if (it->first == "XTandemAdapter")
			{
				Idxml.push_back(it->second);
			}
		}
		
		for (vector<string>::const_iterator it = 
		*/
		 
		
//------------------------------------------------------------------------------------------------

        //Wenn ihr die Metriken schreibt lasst euch bitte ein Int ausgeben 1 wenn erfolgreich, 0  wenn nicht
        //Die Daten die ihr erhaltet am besten als MetricMap, die wichtigen Funktionen dafür stehen oben
    void Metriken::runAllMetrics(){
        ////////////////Metrik1: Protein And Peptide Count /////////////////////////////////
			MetricMap PeptideCountData;
			MetricMap ProteinCountData;
			int papc = this->ProteinAndPeptideCount_(PeptideCountData,ProteinCountData);
	    	////////////////Metrik2: Mass Calibration ----------------------////////////////////
			MetricMap RTData;
			MetricMap MZTData;
			MetricMap MZOData;
			MetricMap PPMBData;
			MetricMap PPMAData;

			int cal = this->MassCalibration(RTData, MZTData, MZOData, PPMBData, PPMAData);
        ////////////////Metrik3: ....................../////////////////////////////////////
        ////////////////Metrik4: ....................../////////////////////////////////////
        ////////////////Metrik5: ....................../////////////////////////////////////

        //MzTab Writer:
        MzTabFile MzTabOutputFile;
        MzTab MzTabAusgabe;
        if(papc == 1){
            int numOfPeptides = PeptideCountData.size();
            int numOfProteins = ProteinCountData.size();
            vector<String> empty;
            vector<String> allPeptides = numOfPeptides>0?PeptideCountData.getStringsByHead("peptide"):empty;
            vector<String> allProteins = numOfProteins>0?ProteinCountData.getStringsByHead("protein"):empty;
            MzTabPeptideSectionRows PepROWS;
            MzTabProteinSectionRows ProtROWS;
            if(numOfPeptides>0){
              for(int i = 0; i< numOfPeptides;i++){
                MzTabPeptideSectionRow PepROW;
                MzTabString PepSeq;
                PepSeq.set(allPeptides[i]);
                PepROW.sequence = PepSeq;
                PepROWS.push_back(PepROW);
              }
            }
            if(numOfProteins>0){
              for(int i = 0; i< numOfProteins;i++){
                MzTabProteinSectionRow ProtROW;
                MzTabString ProtSeq;
                ProtSeq.set(allProteins[i]);
                ProtROW.description = ProtSeq;
                ProtROWS.push_back(ProtROW);
              }
            }
            MzTabAusgabe.setPeptideSectionRows(PepROWS);
            MzTabAusgabe.setProteinSectionRows(ProtROWS);
        }			
				MzTabOutputFile.store("output",MzTabAusgabe); 
				
				

    }

