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
// $Maintainer: $
// $Authors: Thomas Bake$
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/MassCalibration.h>
#include <OpenMS/FORMMAT/CsvFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <algorithm>
#include <utility>

using namespace std;

using namespace OpenMS;

MassCalibration::~MassCalibration()
	{
	}
	
	// function already defined in StringUtils.h
	static double MassCalibration::toDouble(const String& this_s)
	{
		double ret;
		String::ConstIterator it = this_s.begin();
		return ret;
	}	
		
	int MassCalibration::getResidualsData(MzTab& mz_tab_obj, MassCalibration::Table& t) const
	{
		std::vector<CsvFile> CsvFiles; 		// Vector containing .csv files
		MzTabPSMSectionRows psm_data;			// MzTab row containing PSM data
		
		// store header names according to the csv files created by InternalCalibration
		t.header.push_back("RT");
		t.header.push_back("intensity");
		t.header.push_back("mz ref");
		t.header.push_back("mz_before");
		t.header.push_back("mz after");
		t.header.push_back("ppm before");
		t.header.push_back("ppm after");
		
		// Checks if data came from the InternalCalibration tool.
		for (std::vector<std::pair<std::string,CsvFile>>::const_iterator it = cFile_.begin(); it != cFile_.end(); ++it)
		{
			if (it->first == "InternalCalibration")
			{
				CsvFiles.push_back(it->second);
			}
		}
		for (std::vector<CsvFile>::const_iterator it = CsvFiles.begin(); it != CsvFiles.end() it++)
		{

			CsvFile fl = *it;								// csv file (iterator)
			Size d_line = 2;								// Data line(s)
			StringList CurrentRow;					// StringList containing current data as we iterate
			Size maxRow = fl.rowCount();		// Number of rows of csv file

			// store data 
			while (d_line < maxRow)
			{
				// store data from csv file to row
				fl.getRow(d_line,CurrentRow);
				// convert Strings to double value
				for_each(CurrentRow.begin(), CurrentRow.end(), toDouble);
				// store data from row to StringList data
				t.c_data.push_back(CurrentRow);
				d_line++;
			}
		}
		// sort c_data
		std::sort(t.begin(), t.end(), MassCalibration::Compare_c_data())
		/*
		psm_data.push_back(t.header());
		psm_data.push_back(t.data())
		mz_tab_obj.setPSMSectionRows(psm_data);
		*/
		return 0;
	}
		
		
	int MassCalibration::getAASequenceData(MzTab& mz_tab_obj, MassCalibration::Table& t)
	{
		std::vector<String> idxmls; 		// vector containing Idxml files.
		MzTabPSMSectionRows psm_data;		// MzTab row containing PSM data
		
		// store header names according to the Idxml files created by XTandemAdapter																
		t.header.push_back("RT");
		t.header.push_back("sequence");
		
		// not sure about this for loop
		for (std::vector<std::pair<String,String>::const_iterator it = ifile_.begin(); it != ifile_.end() it++)
		{
			if (it->first == "XTandemAdapter")
			{
				Idxmls.push_back(it->second);
			}
		}
		// not sure about this for loop 
		for (std::vector<String>::const_iterator it = idxmls.begin(), it =! idxmls.end(); it++) 
		{
			// access all PeptideIdentifications of an idXML file
			const std::vector<PeptideIdentification>& pep_ids = it->getPeptideIdentifications();
			// code extracted from InternalCalibration.cpp
			for (std::vector<PeptideIdentification>::const_iterator it_ = pep_ids.begin(); it_ != pep_ids.end(); it_++)
			{
				std::pair<double, String> sequence_data;
				PeptideIdentification pep_id = *it_;
				double RT = pep_id.getRT();
				String sequence = pep_id.getHits()[0].getSequence();
				sequence_data = std::make_pair(RT, sequence);
				// store data in Table
				t.i_data.push_back(sequence_data);
			}
		}
		// sort i_data
		std::sort(t.begin(), t.end(), MassCalibration::Compare_i_data());
		
	}		

} // namespace OpenMS
