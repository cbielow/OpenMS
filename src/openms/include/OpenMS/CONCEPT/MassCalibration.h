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
// $Maintainer: Thomas Bake$
// $Authors: Thomas Bake$
// --------------------------------------------------------------------------

#ifndef OPENMS_CONCEPT_MASSCALIBRATION_H
#define OPENMS_CONCEPT_MASSCALIBRATION_H

//includes

#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <OpenMS/FORMAT/MzTabFile.h>

#include <vector>
#include <utility>

namespace OpenMS;
{
	class CsvFile;
	class String;
	class PeptideIdentification;
	class FeatureMap;
	
	class OPENMS_DLLAPI MassCalibration
	{
	
		public:
		
		// Constructor
		MassCalibration(std::vector<std::pair<OpenMS::String, OpenMS::CsvFile>> files, std::vector<std::pair<OpenMS::String, OpenMS::String>> ifiles):
			cFile_(files),
			idxml_(ifiles)
			{
			};
		// Destructor
		~MassCalibration();
		
		struct Table
		{
			StringList header;
			std::vector<StringList> c_data;
			std::vector<std::pair<double, String>> i_data;
			
			Table(StringList h_, std::vector<StringList> c_d_, std::vector<std::pair<double, String>> i_d_): 
				header(h_),
				c_data(c_d_),
				i_data(i_d_) 
				{}
				
				
				double getcRT(Size i) const
				{
					return c_data[i][0];
				}
				
				double getiRT(Size i) const
				{
					return i_data[i].first;
				}
		};
		
		struct Compare_c_data :
		std::binary_function<Table, Table, bool>
		{
			bool operator()(const Table& lhs, const Table& rhs) const
			{
				return lhs.getcRT() < rhs.getcRT();
			}
		};
		
		struct Compare_i_data : 
		std::binary_function<Table, Table, bool>
		{
			bool operator()(const Table& lhs, const Table& rhs) const
			{
				return lhs.getiRT() < rhs.getiRT();
			}
		};
	}
	
			static double toDouble(const String& this_s);
			
			int getResidualsData(MzTab& mz_tab_obj, MassCalibration::Table& t) const;
		
			int getAASequenceData(MzTab&) const;

		protected:	
		 	std::vector<std::pair<OpenMS::String, OpenMS::CsvFile>> cFile_;
		 	std::vector<std::pair<OpenMS::String, OpenMS::String>> idxml_;
			
	}; // class MassCalibration
		
}// namespace OpenMS

#endif
