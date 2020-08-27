// --------------------------------------------------------------------------
//           OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//  notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//  notice, this list of conditions and the following disclaimer in the
//  documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//  may be used to endorse or promote products derived from this software
//  without specific prior written permission.
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
// $Maintainer: Chris Bielow $
// $Authors: Patricia Scheil $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/TextFile.h>

#include <OpenMS/ANALYSIS/ID/AhoCorasickDA.h>
#include <iostream>
#include <algorithm>
#include <string>

using namespace OpenMS;


START_TEST(AhoCorasickDA, "$Id$")


  std::ifstream in(OPENMS_GET_TEST_DATA_PATH("sortedSeq.txt"));
  String str;
  std::vector<String> vec_seq;
  while (std::getline(in, str))
  {
    if(!str.empty())
    {
      vec_seq.push_back(str);
    }
  }

  std::vector<String> vec_short;
  vec_short.insert(vec_short.begin(), vec_seq.begin(), vec_seq.end());

  AhoCorasickDA ac_da(vec_short);

  ac_da.printDA(false);

  std::vector<String> keys = {"VAR", "AAR", "AA", "DVAR", "DKL", "ARK", "ARK"};
  std::vector<String> keysAmbigious = {"OK", "WITHB"};
  std::vector<String> empty = {};
  //std::vector<String> large = ...;

  AhoCorasickDA ac_keys(keys);

  ac_keys.printDA(true);





  START_SECTION(bool findNext(Size& pos_in_protein, Size& peptide_index))
  {

    TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidParameter, AhoCorasickDA ac_empty(empty), "No sequences given");
    TEST_EXCEPTION_WITH_MESSAGE(Exception::IllegalArgument, AhoCorasickDA ac_ambigious(keysAmbigious), "Input peptide must NOT contain ambiguous amino acids (B/J/Z/X) or characters other than amino acids!");
    //TEST_EXCEPTION(Exception::InvalidSize, AhoCorasickDA ac_large(large))

    Size pos = 0;
    Size idx = 0;

    TEST_EXCEPTION_WITH_MESSAGE(Exception::MissingInformation, ac_keys.findNext(pos, idx), "No protein for retrieval. Use function 'setProtein()'");

    ac_keys.setProtein("AARDVARK");

    std::vector<Size>position{};
    std::vector<Size>index{};
    while(ac_keys.findNext(pos, idx))
    {
      position.push_back(pos);
      index.push_back(idx);
    }

    TEST_EQUAL(position[0], 0);
    TEST_EQUAL(position[1], 0);
    TEST_EQUAL(position[2], 3);
    TEST_EQUAL(position[3], 4);
    TEST_EQUAL(position[4], 5);
    TEST_EQUAL(position[5], 5);

    TEST_EQUAL(index[0], 2);
    TEST_EQUAL(index[1], 1);
    TEST_EQUAL(index[2], 3);
    TEST_EQUAL(index[3], 0);
    TEST_EQUAL(index[4], 5);
    TEST_EQUAL(index[5], 6);



    ac_keys.setProtein("VAR");

    ac_keys.findNext(pos, idx);

    TEST_EQUAL(pos, 0);
    TEST_EQUAL(idx, 0);


  }
  END_SECTION

END_TEST
