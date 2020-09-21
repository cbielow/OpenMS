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
#include <stdint.h>


namespace OpenMS
{

  /**
   * @brief Extended Aho-Corasick algorithm capable of matching ambiguous amino acids in the pattern (i.e. proteins).
   * The trie of the Aho Corasick automaton is based on a compact double-array (DA) :
   *
   * Susumu  Yata  Masaki  Oono  Toru  Sumitomo  Kazuhiro  Morita,  Masao  Fuketa  and  Jun’ichiAoe.
   * A  compact  static  double-array  keeping  character  codes.
   * Information Processing and Managemnet : 237-247 (2007)
   */
  class OPENMS_DLLAPI AhoCorasickDA
  {

  public:

    /// constructor support DNA, RNA and Protein sequences (lower and capital letter) except ambiguous aa (B/J/X/Z)
    explicit AhoCorasickDA(const std::vector<String>& seq) : sequences_(seq)
    {
      construct_();
    }

    /// Destructor
    ~AhoCorasickDA() = default;


    /**
     * @brief 
     * @param pos_in_protein position in the set protein
     * @param peptide_index index of sequence in given vector
     * @return  False if end of protein is reached. True if a hit is found.
     */
    bool findNext(Size& pos_in_protein, Size& peptide_index);

    /// Set new protein sequence. All previous data is forgotten.
    void setProtein(const char* prot, const uint8_t amb_max);

    // for debugging
    bool retrievalCDA();

    // for debugging
    void printDA(bool arrays);



  private:

    /// Needed for the construction of the BC array. Contains informations about used nodes, used base values and empty elements
    struct BuildInformation_
    {
      explicit BuildInformation_(const UInt32 index) : next(index+1), used_flag(0), prev(index-1),  base_flag(0), check(0){};
      UInt32 next        :31 ;    ///< next empty node in DA
      UInt32 used_flag   :1  ;    ///< node is used
      UInt32 prev        :31 ;    ///< previous empty node in the DA
      UInt32 base_flag   :1  ;    ///< this value is used as base
      UInt32 check           ;    ///< parent node TODO not necessary
    };

    /// The properties of a node in the BC array
    struct BC_
    {
      UInt32 base       :22 ;    ///< base value of DA
      UInt32 lcheck     :8  ;    ///< arc label from parent node. Only code of the aa is stored
      UInt32 leaf_flag  :1  ;    ///< node is leaf
      UInt32 term_flag  :1  ;    ///< marks the end of a sequence if the sequence is a substring of another sequence
    };

    /// The properties of node in Tail
    struct Tail_
    {
      explicit Tail_(const uint8_t code) : label(code), term_flag(0){};
      uint8_t label      :7  ;  ///< code of aa
      uint8_t term_flag  :1  ;  ///< marks the end of a sequence if the sequence is a substring of another sequence
    };

    /// properties of Spawn
    struct Spawn_
    {
      Int32 node;
      UInt32 prot_pos;
      uint8_t counter;
      uint16_t max_depth_decrease;
    };


    struct SupplyLink_
    {
      Int32 link      :24 ;   ///< node to which the supply link leads
      Int32 depth     :8  ;   ///< depth of node
    };

    /// given sequences for trie
    const std::vector<String>& sequences_{};

    /// indices of the sorted sequences
    std::vector<UInt32> sorted_idx_{};

    /// indices of the sorted sequences, duplicate sequences removed
    std::vector<UInt32> unique_idx_{};

    /// stores information needed for the construction
    std::vector<BuildInformation_> build_info_{};

    /// BC array. Contains all nodes which have at least two children.
    std::vector<BC_> bc_{};

    /// tail
    std::vector<Tail_> tail_ = {Tail_(0)};

    /// failure function
    std::vector<SupplyLink_> failure_tail_;
    std::vector<SupplyLink_> failure_bc_;

    /// key: node position in DA, value: index of Sequence in sequences_
    std::unordered_map<Int32, UInt32> output_{};

    /// first empty element in DA
    UInt32 first_empty_elem_ = 0;

    /// stores all edges from a node
    std::vector<char> edges_{};

    /// protein sequence
    const char* protein_ = "";

    /// position in protein
    UInt32 prot_pos_ = 0;

    /// largest index of already read amb aa
    UInt32 pos_max_ = 0;

    /// depth of the first read aa
    uint16_t max_depth_decrease_ = 0;

    /// position in DA
    Int32 node_pos_= 0;

    /// contains all spawns that have not yet been processed
    std::vector<Spawn_> spawns_{};

    /// contains indices of sequences if more than one hit was found at a node
    std::vector<UInt32> multi_hit_;

    /// Maximum permitted number of amb aa
    uint8_t amb_max_ = 0;

    /// current count of amb aa
    uint8_t amb_count_ = 0;


    /**
     * @brief code table
     * All characters not representing an amino acid are treated as X.
     * Code 0 is used in the BC array for label of unused nodes and in the tail as end of a sequence.
     */
    uint8_t code_table_ [256] =
        { // char --> amino acid (unsigned char with values from 1..26)
            0 ,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26, //0
            26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26, //1
            26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26, //2
            26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26, //3

        //    ,   A,   B,   C,   D,   E,   F,   G,   H,   I,   J,   K,   L,   M,   N,   O,
            26,   3,  23,   19,  10,  8,  15,   4,  18,   2,  25,   7,   1,  17,  11,  21, //4

        //   P,   Q,   R,   S,   T,   U,   V,   W,   X,   Y,   Z,    ,    ,    ,    ,    ,
            14,   9,   12,  6,  13,  22,   5,  20,  26,  16,  24,  26,  26,  26,  26,  26, //5

        //    ,   a,   b,   c,   d,   e,   f,   g,   h,   i,   j,   k,   l,   m,   n,   o,
            26,   3,  23,   19,  10,  8,  15,   4,  18,   2,  25,   7,   1,  17,  11,  21,//6

        //   p,   q,   r,   s,   t,   u,   v,   w,   x,   y,   z,    ,    ,    ,    ,    ,
            14,   9,   12,  6,  13,  22,   5,  20,  26,  16,  24,  26,  26,  26,  26,  26, //7

            26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26, //8
            26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26, //9
            26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26, //10
            26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26, //11
            26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26, //12
            26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26, //13
            26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26, //14
            26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26,  26  //15
        };

    /// Number of non ambiguous aa
    static constexpr const uint8_t AA_COUNT_ = 22;

    //-------------------------------------------------------------------------------------------------------------------------------------
    //  CDA - functions
    //-------------------------------------------------------------------------------------------------------------------------------------

    /// returns code of label
    uint8_t code_(const char label);

    /// returns the index of child of the node with base_value with the transition label uses XOR
    UInt32 child_(const UInt32 base_value, const char label);
    UInt32 child_(const UInt32 base_value, const uint8_t label);

    /// tests if a aa is ambiguous (B, X, J, Z)
    bool isAmbiguous_(const char s);


    //-------------------------------------------------------------------------------------------------------------------------------------
    //  CDA - retrieval
    //-------------------------------------------------------------------------------------------------------------------------------------

    //-----------------------------------------------------------------------------------------------------------------
    // TODO - Question:
    //  The functions 'failure_' and 'consumeAmbiguousAA_' call 'retrieval_' if a valid transition or spawn was
    //  found and return false otherwise. If 'retrieval_' returns false, it means that nothing else was found.
    //  If 'retrieval_' is called in one of the two functions, the result is returned accordingly. But because
    //  'retrieval_' calls these two functions itself, it is very interlaced. Should the functions only return
    //  true and in 'retrieval_' even then retrieval_() should be called at true?
    //  'useSpawn_' too.
    //-----------------------------------------------------------------------------------------------------------------

    /**
     * @brief Follows the supply links of a node marked with term_flag to the root and stores all substrings found in multi_seq.
     */
    void findSubstrOfHit_();

    /**
     * @brief checks if a sequence occurs more than once in the given sequences and stores all indices in multi_seq_
     * @param pep_idx index of found sequence
     */
    void multipleOccurrence_(const Size pep_idx );

    /**
     * @brief
     * @return returns true if hit was found and false if there is nothing more to found
     */
    bool retrieval_();

    /**
     * @brief compares the peptide with the end of the sequences stored in Tail
     * @param tail_str pointer to a given position in Tail
     * @return true if a sequences was found, otherwise failure is called and returned
     */
    bool compareTail_(const Tail_* tail_str);

    /**
     * @brief set new node_pos if retrieval fails
     * @return false if no transition was found with the help of the supply link,
     * if no transition with the next aa was found and no spawn exists.
     * Otherwise retrieval is called with the new node and peptide position.
     */
    bool failure_();

    /**
     * @brief returns the supply link if this does not eliminate ambiguous aa. Don't do anything special with 0. Must be considered when calling
     * @param node node from which the supply link should be taken
     * @param output_node set to the linked node
     * @param amb_depth depth of the first ambiguous aa in trie
     * @return false, if a supply link would cause the first read aa to be dropped. Otherwise true
     */
    bool followSupplyLink_(Int32 node, Int32& output_node, uint16_t& amb_depth);

    /**
     * @brief processes ambiguous aa
     * @param aa the ambiguous aa
     * @return false if no spawn was found and no transition with the next aa was found.
     * Otherwise retrieval is called and returned
     */
    bool consumeAmbiguousAA_(const char aa);


    /**
     * @brief sets to the smallest and largest code value of the given aa
     * @param c Character of an ambiguous aa
     * @param idxFirst given by reference and set to the smallest code value
     * @param idxLast given by reference and set to the largest code value
     */
    void getSpawnRange_(const char c, uint8_t & idxFirst, uint8_t & idxLast);

    /**
     * @brief if spawns exist, the next spawn is set
     * @return false if there are no spwans, otherwise retrieval
     */
    bool useSpawn_();



    //TODO macht exakt das gleiche wie 'getNode_' mit true/false.
    /**
     * @brief returns the child node for each node depending on whether it is a node in trie, tail or leaf node
     * @param node current node
     * @param label code of label
     * @param output_node child node
     * @return true if a child node was found
     */
    bool getChildNode_(const Int32 node, const uint8_t label,  Int32& output_node);

    /**
     * @brief calls getChildNode_. If false, the supply links are followed to 0 and the transition is tested for each node.
     * @param node current node
     * @param label code of label for transition
     * @param output_node child node if found
     * @param amb_depth depth of the first ambiguous aa
     * @return true if child node is found
     */
    bool getNextNode_(const Int32 node, const uint8_t label,  Int32& output_node, uint16_t & amb_depth);



    /**
     * @brief tests if a sequence ends in node
     * @param node node to be tested
     * @return true if hit else false
     */
    bool isHit_(const Int32 node);


    /// returns depth of the given node
    uint8_t getDepth_(const Int32 node);

    //-------------------------------------------------------------------------------------------------------------------------------------
    //  CDA - construction
    //-------------------------------------------------------------------------------------------------------------------------------------

    /// calls functions to construct the Aho Corasick
    void construct_();

    /// expands the BC array and initializes linkage
    void expand_();

    /// recursive function to create the trie by subtries
    void buildTrie_(Size begin, const Size end, const Size depth, const UInt32 node_pos);

    /// returns base value
    UInt32 xCheck_();
    UInt32 plusCheck_();


    /// tests the base value for all labels in edges
    bool checkAllEdges_(const UInt32 base);

    /// deletes a used empty element form linkage
    void arrangeEmpty_(const UInt32 node_pos);

    /// stores the suffix string in tail as chars
    void appendToTail_(const String& suffix);

    /// Construct the failure function
    void constructFailure_();


    /**
     * @brief compute child node for failure function. If node is postive and not a leaf node child node is computed with child_ function.
     * if node is a leaf node or negative it returns negative tail index if label is correct
     * @param node node
     * @param label code of label
     * @return either the child node or 0 if child doesn't exist
     */
    Int32 getNode_(const Int32 node, const uint8_t label);

    /**
     * @brief set the given supply link for node. If node negative in failure_tail, if positiv in failure_bc
     * @param node the node number in DA if positiv, index in tail if negative
     * @param sl the computed supply link
     */
    void setSupplyLink_(const Int32 node, Int32 sl, const uint8_t depth);

    /// returns the supply link for node either from failure_tail if node is negative or from failure_bc if positive
    Int32 getSupplyLink_(const Int32 node);


    /// tests if a supply link leads to a hit and marks the node with a term_flag
    void isEndOfSeq_(const Int32 node);

    /// tests if a sequence ends in node. During construction node+1 in tail must be used
    bool isHitConstr_(const Int32 node);
  };


} // namespace OpenMS
