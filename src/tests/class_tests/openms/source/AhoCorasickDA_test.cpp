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

  std::vector<String> keys = {"VARK", "AAR", "AA", "DLAR", "DKL", "ARK", "ARK", "A", "KVAR", "KVAK"};
  std::vector<String> ambigious_tail = {"WITHB"};
  std::vector<String> ambigious_bc = {"ABCDE", "ABCFG"};

  std::vector<String> empty = {};
  //std::vector<String> large = ...;

  AhoCorasickDA ac_keys(keys);

  ac_keys.printDA(true);

  START_SECTION(bool findNext(Size& pos_in_protein, Size& peptide_index))
  {


    TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidParameter, AhoCorasickDA ac_empty(empty), "No sequences given");
    TEST_EXCEPTION_WITH_MESSAGE(Exception::IllegalArgument, AhoCorasickDA ac_ambigious_t(ambigious_tail ), "Input peptide must NOT contain ambiguous amino acids (B/J/Z/X) or characters other than amino acids!");
    TEST_EXCEPTION_WITH_MESSAGE(Exception::IllegalArgument, AhoCorasickDA ac_ambigious_b(ambigious_bc ), "Input peptide must NOT contain ambiguous amino acids (B/J/Z/X) or characters other than amino acids!");
    //TEST_EXCEPTION(Exception::InvalidSize, AhoCorasickDA ac_large(large))

    Size pos = 0;
    Size idx = 0;

    //TEST_EXCEPTION_WITH_MESSAGE(Exception::MissingInformation, ac_keys.findNext(pos, idx), "No protein for retrieval. Use function 'setProtein()'");

    std::vector<String> observed{};

    ac_keys.setProtein("AARDVARK", 0);
    while(ac_keys.findNext(pos, idx))
    {
      observed.push_back(keys[idx] + String(idx) + " @ " + String(pos));
      //std::cout << keys[idx] + String(idx) + " @ " + String(pos) << std::endl;
      //std::cout << std::endl;
      //std::cout << std::endl;
    }

    for (auto i: observed)
    {
      std::cout << i << " ,  ";
    }
    std::cout << std::endl;

    TEST_EQUAL(observed.size(), 8)
    TEST_EQUAL(observed[0], "A7 @ 0");
    TEST_EQUAL(observed[1], "AA2 @ 0");
    TEST_EQUAL(observed[2], "A7 @ 1")
    TEST_EQUAL(observed[3], "AAR1 @ 0");
    TEST_EQUAL(observed[4], "A7 @ 5")
    TEST_EQUAL(observed[5], "VARK0 @ 4");
    TEST_EQUAL(observed[6], "ARK6 @ 5");
    TEST_EQUAL(observed[7], "ARK5 @ 5");



    ac_keys.setProtein("DKL", 0);

    ac_keys.findNext(pos, idx);

    TEST_EQUAL(keys[idx] + String(idx) + " @ " + String(pos), "DKL4 @ 0");


    ac_keys.setProtein("CBKJ", 3);

    ac_keys.findNext(pos, idx);

    TEST_EQUAL(keys[idx] + String(idx) + " @ " + String(pos), "DKL4 @ 1");

    observed.clear();
    ac_keys.setProtein("AXX", 3);
    while(ac_keys.findNext(pos, idx))
    {
      observed.push_back(keys[idx] + String(idx) + " @ " + String(pos));
      //std::cout << keys[idx] + String(idx) + " @ " + String(pos) << std::endl;
      //std::cout << std::endl;
      //std::cout << std::endl;
    }

    for (auto i: observed)
    {
      std::cout << i << " ,  ";
    }
    std::cout << std::endl;

    TEST_EQUAL(observed.size(), 8)
    TEST_EQUAL(observed[0], "A7 @ 0");
    TEST_EQUAL(observed[1], "A7 @ 2");
    TEST_EQUAL(observed[2], "ARK5 @ 0");
    TEST_EQUAL(observed[3], "ARK6 @ 0");
    TEST_EQUAL(observed[4], "AA2 @ 0");
    TEST_EQUAL(observed[5], "A7 @ 1");
    TEST_EQUAL(observed[6], "AAR1 @ 0");
    TEST_EQUAL(observed[7], "AA2 @ 1");



    observed.clear();
    ac_keys.setProtein("XXX", 3);
    while(ac_keys.findNext(pos, idx))
    {
      observed.push_back(keys[idx] + String(idx) + " @ " + String(pos));
    }

    for (auto i: observed)
    {
      std::cout << i << " ,  ";
    }
    std::cout << std::endl;

    TEST_EQUAL(observed.size(), 9)
    TEST_EQUAL(observed[0], "A7 @ 2");
    TEST_EQUAL(observed[1], "A7 @ 1");
    TEST_EQUAL(observed[2], "AA2 @ 1");
    TEST_EQUAL(observed[3], "DKL4 @ 0");
    TEST_EQUAL(observed[4], "A7 @ 0");
    TEST_EQUAL(observed[5], "ARK5 @ 0");
    TEST_EQUAL(observed[6], "ARK6 @ 0");
    TEST_EQUAL(observed[7], "AA2 @ 0");
    TEST_EQUAL(observed[8], "AAR1 @ 0");

    observed.clear();
    ac_keys.setProtein("XVAXX", 3);
    while(ac_keys.findNext(pos, idx))
    {
      observed.push_back(keys[idx] + String(idx) + " @ " + String(pos));
      //std::cout << keys[idx] + String(idx) + " @ " + String(pos) << std::endl;
      //std::cout << std::endl;
      //std::cout << std::endl;
    }

    for (auto i: observed)
    {
      std::cout << i << " ,  ";
    }
    std::cout << std::endl;

    TEST_EQUAL(observed.size(), 12)
    TEST_EQUAL(observed[0], "A7 @ 2");
    TEST_EQUAL(observed[1], "A7 @ 4");
    TEST_EQUAL(observed[2], "VARK0 @ 1");
    TEST_EQUAL(observed[3], "ARK6 @ 2");
    TEST_EQUAL(observed[4], "ARK5 @ 2");
    TEST_EQUAL(observed[5], "AA2 @ 2");
    TEST_EQUAL(observed[6], "A7 @ 3");
    TEST_EQUAL(observed[7], "AAR1 @ 2");
    TEST_EQUAL(observed[8], "AA2 @ 3");
    TEST_EQUAL(observed[9], "KVAR8 @ 0");
    TEST_EQUAL(observed[10], "KVAK9 @ 0");
    TEST_EQUAL(observed[11], "A7 @ 0");


    observed.clear();
    ac_keys.setProtein("VARXKJ", 2);
    while(ac_keys.findNext(pos, idx))
    {
      observed.push_back(keys[idx] + String(idx) + " @ " + String(pos));
    }
    for (auto i: observed)
    {
      std::cout << i << " ,  ";
    }
    std::cout << std::endl;

    TEST_EQUAL(observed.size(), 6)
    TEST_EQUAL(observed[0], "A7 @ 1");
    TEST_EQUAL(observed[1], "DKL4 @ 3");
    TEST_EQUAL(observed[2], "VARK0 @ 0");
    TEST_EQUAL(observed[3], "ARK6 @ 1");
    TEST_EQUAL(observed[4], "ARK5 @ 1");
    TEST_EQUAL(observed[5], "A7 @ 3");


    observed.clear();
    ac_keys.setProtein("XXAX");
    while(ac_keys.findNext(pos, idx))
    {
      observed.push_back(keys[idx] + String(idx) + " @ " + String(pos));
    }

    for (auto i: observed)
    {
      std::cout << i << " ,  ";
    }
    std::cout << std::endl;

    TEST_EQUAL(observed.size(), 11)
    TEST_EQUAL(observed[0], "A7 @ 2");
    TEST_EQUAL(observed[1], "AA2 @ 2");
    TEST_EQUAL(observed[2], "A7 @ 3");
    TEST_EQUAL(observed[3], "A7 @ 1");
    TEST_EQUAL(observed[4], "AA2 @ 1");
    TEST_EQUAL(observed[5], "AAR1 @ 1");
    TEST_EQUAL(observed[6], "DLAR3 @ 0");
    TEST_EQUAL(observed[7], "KVAR8 @ 0");
    TEST_EQUAL(observed[8], "KVAK9 @ 0");
    TEST_EQUAL(observed[9], "A7 @ 0");
    TEST_EQUAL(observed[10], "AA2 @ 0");



    observed.clear();
    ac_keys.setProtein("XXXX", 2);
    while(ac_keys.findNext(pos, idx))
    {
      observed.push_back(keys[idx] + String(idx) + " @ " + String(pos));
    }

    for (auto i: observed)
    {
      std::cout << i << " ,  ";
    }
    std::cout << std::endl;

    TEST_EQUAL(observed.size(), 7)
    TEST_EQUAL(observed[0], "A7 @ 3");
    TEST_EQUAL(observed[1], "A7 @ 2");
    TEST_EQUAL(observed[2], "AA2 @ 2");
    TEST_EQUAL(observed[3], "A7 @ 1");
    TEST_EQUAL(observed[4], "AA2 @ 1");
    TEST_EQUAL(observed[5], "A7 @ 0");
    TEST_EQUAL(observed[6], "AA2 @ 0");





    observed.clear();
    ac_keys.setProtein("VAXXXXKL", 1);
    while(ac_keys.findNext(pos, idx))
    {
      observed.push_back(keys[idx] + String(idx) + " @ " + String(pos));

    }
    for (auto i: observed)
    {
      std::cout << i << " ,  ";
    }
    std::cout << std::endl;

    TEST_EQUAL(observed.size(), 7)
    TEST_EQUAL(observed[0], "A7 @ 1")
    TEST_EQUAL(observed[1], "DKL4 @ 5");
    TEST_EQUAL(observed[2], "A7 @ 5");
    TEST_EQUAL(observed[3], "A7 @ 4");
    TEST_EQUAL(observed[4], "A7 @ 3")
    TEST_EQUAL(observed[5], "AA2 @ 1");
    TEST_EQUAL(observed[6], "A7 @ 2");


    observed.clear();
    ac_keys.setProtein("DKLBXJX", 1);
    while(ac_keys.findNext(pos, idx))
    {
      observed.push_back(keys[idx] + String(idx) + " @ " + String(pos));
      //std::cout << keys[idx] + String(idx) + " @ " + String(pos) << std::endl;
      //std::cout << std::endl;
      //std::cout << std::endl;
    }
    for (auto i: observed)
    {
      std::cout << i << " ,  ";
    }
    std::cout << std::endl;

    TEST_EQUAL(observed.size(), 3)
    TEST_EQUAL(observed[0], "DKL4 @ 0");
    TEST_EQUAL(observed[1], "A7 @ 6")
    TEST_EQUAL(observed[2], "A7 @ 4");


    observed.clear();
    ac_keys.setProtein("DKLBXJX", 3);
    while(ac_keys.findNext(pos, idx))
    {
      observed.push_back(keys[idx] + String(idx) + " @ " + String(pos));

    }
    for (auto i: observed)
    {
      std::cout << i << " ,  ";
    }
    std::cout << std::endl;

    TEST_EQUAL(observed.size(), 4)
    TEST_EQUAL(observed[0], "DKL4 @ 0");
    TEST_EQUAL(observed[1], "A7 @ 6")
    TEST_EQUAL(observed[2], "A7 @ 4");
    TEST_EQUAL(observed[3], "DKL4 @ 3");


    observed.clear();
    ac_keys.setProtein("aard", 0);
    while(ac_keys.findNext(pos, idx))
    {
      observed.push_back(keys[idx] + String(idx) + " @ " + String(pos));
      //std::cout << keys[idx] + String(idx) + " @ " + String(pos) << std::endl;
      //std::cout << std::endl;
      //std::cout << std::endl;
    }

    for (auto i: observed)
    {
      std::cout << i << " ,  ";
    }
    std::cout << std::endl;

    TEST_EQUAL(observed.size(), 4)
    TEST_EQUAL(observed[0], "A7 @ 0");
    TEST_EQUAL(observed[1], "AA2 @ 0");
    TEST_EQUAL(observed[2], "A7 @ 1")
    TEST_EQUAL(observed[3], "AAR1 @ 0");

    TEST_FILE_EQUAL(OPENMS_GET_TEST_DATA_PATH("PeptideIndexerAmbiguous.idXML"), OPENMS_GET_TEST_DATA_PATH("PeptideIndexerDA.idXML"))
  }
  END_SECTION

END_TEST
