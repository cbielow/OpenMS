// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Chris Bielow$
// $Authors: Patricia Scheil$
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/Types.h>
#include <vector>
#include <stdint.h> // for empty_head


namespace OpenMS
{

  class OPENMS_DLLAPI AhoCorasickDA
  {

  public:

    /// Default constructor
    AhoCorasickDA() = default;

    /// Destructor
    virtual ~AhoCorasickDA() = default;


    //needed for the construction
    //contains informations about used nodes, used base values and empty elements
    struct BuildInformation
    {
      UInt32 next        :31 ;    //next empty node in DA
      UInt32 used_flag   :1  ;    //node is used
      UInt32 prev        :31 ;    //previous empty node in the DA
      UInt32 base_flag   :1  ;    //this value is used as base
      UInt32 check           ;    //parent node
    };

    struct BC
    {
      UInt32 base       :22 ;    // base value of DA
      UInt32 lcheck     :8  ;    // arc label from parent node
      UInt32 leaf_flag  :1  ;    // node is leaf
      UInt32 term_flag  :1  ;    // string end

    };

    std::vector<BuildInformation> build_info{};
    std::vector<BC> bc{};

    //first empty element
    UInt32 empty_head = UINT32_MAX;



    //Code Table, only capital letter
    //                      A, B,  C,  D, E, F,  G, H,  I, J,  K, L, M,  N,  O,  P,  Q,  R, S, T,  U,  V, W,  X,  Y,  Z
    UInt code_tabel [26] = {1, 22, 18, 8, 4, 14, 2, 17, 6, 24, 7, 0, 16, 12, 20, 11, 13, 9, 5, 10, 21, 3, 19, 25, 15, 23};

    // tail
    std::vector<char> tail{};

    // strings to build the Aho Corasick automaton
    std::vector<String> pattern{};

    // stores all edges from a node
    std::vector<char> edges{};

    void construct(std::vector<String>& seq);

    bool retrieval(const char *str);



  private:


    // returns the index of child of the node pos with the transition label
    UInt32 child_(const UInt32& pos, const char& label);

    // expands the BC array and initializes linkage
    void expand_();

    void buildTrie_(Size begin, const Size& end, const Size& depth, const UInt32& node_pos);

    // returns base value
    UInt32 xCheck_();

    // tests the base value for all labels in edges
    bool checkAllEdges_(const UInt32& base);

    // delete a used empty element form linkage
    void arrangeEmpty_(const UInt32& node_pos);

    void storeInTail_(const String& suffix);


    //temp workaround
    UInt32 mini_();

    // returns code of label
    UInt32 code_(const char& label);


    // isAmbiguous



    //void getStatistic();
    ///stat??
    //getEmpty
    //getSize
    //getNumStr

  };


} // namespace OpenMS