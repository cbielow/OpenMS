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
#include <unordered_map>


namespace OpenMS
{

  class OPENMS_DLLAPI AhoCorasickDA
  {

  public:

    /// constructor
    explicit AhoCorasickDA(const std::vector<String>& seq) : sequences_(seq)
    {
      construct_();
    }

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
      UInt32 check           ;    //parent node TODO not necessary
    };

    struct BC
    {
      UInt32 base       :22 ;    // base value of DA
      UInt32 lcheck     :8  ;    // arc label from parent node
      UInt32 leaf_flag  :1  ;    // node is leaf
      UInt32 term_flag  :1  ;
    };

    struct Spawn
    {
      Int32 node = 0;
      UInt32 prot_pos = 0;
      std::vector<UInt> a_acids{};
    };





    bool findNext(Size& pos_in_protein, Size& peptide_index);

    void setProtein(const char* prot);

    // for debugging
    bool retrievalCDA();

    // for debugging
    void printDA(bool arrays);



  private:

    // given sequences for trie
    std::vector<String> sequences_{};

    // indices of the sorted sequences
    std::vector<UInt32> sorted_idx_{};

    // indices of the sorted sequences, duplicate sequences removed
    std::vector<UInt32> unique_idx_{};

    // stores information needed for the construction
    std::vector<BuildInformation> build_info_{};

    // BC array
    std::vector<BC> bc_{};

    // tail
    std::vector<char> tail_ = {'\0'};

    // failure function
    std::vector<Int32> failure_tail_;
    std::vector<Int32> failure_bc_;

    // key: node position in DA, value: index of Sequence in sequences_
    std::unordered_map<Int32, UInt32> output_{};

    //first empty element in DA
    UInt32 empty_head_ = 0;

    // stores all edges from a node
    std::vector<char> edges_{};

    // protein sequence
    const char* protein_ = "";

    // position in protein
    UInt prot_pos_ = 0;

    // position in DA
    Int32 node_pos_= 0;

    std::vector<Spawn> spawns_{};


    // code table
    UInt code_tabel_ [256] =
        { // char --> amino acid (unsigned char with values from 0..26)
            25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25, //0
            25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25, //1
            25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  26,  25,  25,  25,  25,  25, //2
            25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25, //3

        //    ,   A,   B,   C,   D,   E,   F,   G,   H,   I,   J,   K,   L,   M,   N,   O,
            25,   1,  22,   18,  8,   4,  14,   2,  17,   6,  24,   7,   0,  16,  12,  20, //4

        //   P,   Q,   R,   S,   T,   U,   V,   W,   X,   Y,   Z,    ,    ,    ,    ,    ,
            11,  13,   9,   5,  10,  21,   3,  19,  25,  15,  23,  25,  25,  25,  25,  25, //5

        //    ,   a,   b,   c,   d,   e,   f,   g,   h,   i,   j,   k,   l,   m,   n,   o,
            25,   1,  22,  18,   8,   4,  14,   2,  17,   6,  24,   7,   0,  16,  12,  20, //6

        //   p,   q,   r,   s,   t,   u,   v,   w,   x,   y,   z,    ,    ,    ,    ,    ,
            11,  13,   9,   5,  10,  21,   3,  19,  25,  15,  23,  25,  25,  25,  25,      //7

            25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25, //8
            25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25, //9
            25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25, //10
            25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25, //11
            25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25, //12
            25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25, //13
            25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25, //14
            25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25,  25  //15
        };



    // returns code of label
    UInt32 code_(const char& label);

    // returns the index of child of the node with base_value with the transition label uses XOR
    UInt32 child_(const UInt32& base_value, const char& label);

    // calls functions to construct the Aho Corasick
    void construct_();

    // recursive function to create the trie by subtries
    void buildTrie_(Size begin, const Size& end, const Size& depth, const UInt32& node_pos);


    bool retrieval_();

    //expands the BC array and initializes linkage
    void expand_();


    // returns base value
    UInt32 xCheck_();
    UInt32 plusCheck_();


    // tests the base value for all labels in edges
    bool checkAllEdges_(const UInt32& base);

    // delete a used empty element form linkage
    void arrangeEmpty_(const UInt32& node_pos);

    // store the suffix string in tail as chars
    void storeInTail_(const String& suffix);

    // tests if a aa is ambiguous (B, X, J, Z)
    bool isAmbiguous_(const char& s);

    bool testAmbiguous(const char aa);

    //Construct the failure function
    void constructFailure_();

    /**
     * @brief set the given supply link for node. If node negative in failure_tail, if positiv in failure_bc
     * @param node the node number in DA if positiv, index in tail if negative
     * @param sl the computed supply link
     */
    void setSupplyLink_(const Int32& node, const Int32& sl);

    // returns the supply link for node either from failure_tail if node is negative or from failure_bc if positive
    Int32 getSupplyLink_(const Int32& node);

    /**
     * @brief compute child node for failure function. If node is postive and not a leaf node child node is computed with child_ function.
     * if node is a leaf node or negative it returns negative tail index if label is correct
     * @param node
     * @param label
     * @return either the child node or 0 if child doesn't exist
     */
    Int32 getNode_(const Int32& node, const UInt& label);


  };


} // namespace OpenMS
