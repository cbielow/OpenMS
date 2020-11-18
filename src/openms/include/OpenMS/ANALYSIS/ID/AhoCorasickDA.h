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

    /// Constructor support DNA, RNA and Protein sequences (accept lower and capital letter except ambiguous aa (B/J/X/Z))
    explicit AhoCorasickDA(const std::vector<String>& seq) : sequences_(seq)
    {
      construct_();
    }


    /// Destructor
    ~AhoCorasickDA() = default;


    /// Set new protein sequence. All previous data is forgotten.
    void setProtein(const char* prot, const uint8_t amb_max = 3);


    /**
     * @brief search for and numerate all occurrences of the peptides in the protein
     * @param[out] pos_in_protein Position (zero-based) in the protein
     * @param[out] peptide_index index of sequence in peptide database
     * @return False if end of protein is reached. True if a hit is found.
     */
    bool findNext(Size& pos_in_protein, Size& peptide_index);


    // for debugging
    bool retrievalCDA();

    // for debugging
    void printDA(bool arrays);

    UInt32 stat();


  private:

    /// Node informations needed for the construction of the BC array. Contains informations about used nodes, used base values and empty elements.
    struct BuildInformation_
    {
      explicit BuildInformation_(const UInt32 index) : next(index+1), used_flag(0), prev(index-1),  base_flag(0), check(0){};
      UInt32 next        :31 ;    ///< next empty node in DA
      UInt32 used_flag   :1  ;    ///< node is used
      UInt32 prev        :31 ;    ///< previous empty node in the DA
      UInt32 base_flag   :1  ;    ///< this value is used as base
      UInt32 check           ;    ///< parent node TODO not necessary
    };

    /// The properties of a node in the base-check-array (BC array)
    struct BCNode_
    {
      UInt32 base       :22 ;    ///< base value of the node
      UInt32 lcheck     :8  ;    ///< arc label from parent node. Only code of the aa is stored
      UInt32 leaf_flag  :1  ;    ///< marks if a node is a leaf
      UInt32 term_flag  :1  ;    ///< marks the end of a sequence if the sequence is a substring of another sequence
      //Int32 link        :24 ;    ///< node to which the supply link leads (longest suffix)
      //Int32 depth       :8  ;    ///< depth of node
    };

    /// The properties of node in Tail
    struct TailNode_
    {
      explicit TailNode_(const uint8_t code) : label(code), term_flag(0){};
      uint8_t label      :7  ;  ///< code of aa, 0 marks the end of a sequence
      uint8_t term_flag  :1  ;  ///< marks the end of a sequence if the sequence is a substring of another sequence
      //Int32 link         :24 ;  ///< node to which the supply link leads (longest suffix)
      //Int32 depth        :8  ;  ///< depth of node
    };

    /// The properties of spawn. A spawn is created when an ambiguous aa is read and the possible aa form a valid transition
    struct Spawn_
    {
      Int32 node;                   ///< Node after the ambiguous aa transition
      UInt32 prot_pos;              ///< Protein position after the ambiguous aa
      uint8_t counter;              ///< Count of all read ambiguous aa
      uint16_t max_depth_decrease;  ///< Depth of first occurred aa
    };


    struct SupplyLink_
    {
      Int32 link      :24 ;   ///< node to which the supply link leads (longest suffix)
      Int32 depth     :8  ;   ///< depth of node
    };

    /// Given peptides for query
    const std::vector<String>& sequences_{};

    /// Indices of the sorted peptides in sequences_
    std::vector<UInt32> sorted_idx_{};

    UInt32 aa_count_ = 0;

    /// Indices of the sorted sequences, duplicate sequences are removed
    std::vector<UInt32> unique_idx_{};

    /// contains informations about duplicated sequences. Key: Peptide index, Value: position in sorted_idx_ vector, count of occurence
    std::unordered_map<UInt32, std::pair<UInt32, UInt32>> multi_idx_{};


    /// Stores information needed for the construction
    std::vector<BuildInformation_> build_info_{};

    /// BC array. Contains all nodes which have at least two children.
    std::vector<BCNode_> bc_{};

    /// Tail array. Stores only the suffixes of the sequences, if all following nodes have exactly one child.
    std::vector<TailNode_> tail_ = {TailNode_(0)};

    /// failure function
    std::vector<SupplyLink_> failure_tail_;
    std::vector<SupplyLink_> failure_bc_;

    /// Output function: contains all nodes, which are terminal nodes and the index of the Peptid.  Key: node position in Trie, value: original index of peptide in sequences_
    std::unordered_map<Int32, UInt32> output_{};

    /// First empty element in BC array
    UInt32 first_empty_elem_ = 0;

    /// Stores all edges from a node during trie construction
    std::vector<char> edges_{};

    /// Protein sequence for query
    const char* protein_ = "";

    /// Current position in protein
    UInt32 prot_pos_ = 0;



    /// Depth of the first read aa
    uint16_t max_depth_decrease_ = 0;

    /// Current position in trie
    Int32 node_pos_= 0;

    /// Contains all spawns that have not yet been processed
    std::vector<Spawn_> spawns_{};

    /// Contains indices of sequences if more than one hit was found at a node
    std::vector<UInt32> multi_hit_;

    /// Maximum allowed ambiguous characters in the matching protein sequence
    uint8_t amb_max_ = 0;

    /// Current count of amb aa
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

    /// Marks that supply links of the node leads to a hit
    static const UInt32 SUBSTR_IS_HIT = UINT32_MAX;

    //-------------------------------------------------------------------------------------------------------------------------------------
    //  CDA - functions
    //-------------------------------------------------------------------------------------------------------------------------------------

    /// returns code of a aa
    uint8_t code_(const char label);


    /**
     * @brief compute the child node. Uses XOR for computation
     * @param base_value base value of current node
     * @param label aa for transition
     * @return returns the node index of child
     */
    UInt32 child_(const UInt32 base_value, const char label);
    UInt32 child_(const UInt32 base_value, const uint8_t label);

    /// tests if a aa is ambiguous (B, X, J, Z)
    bool isAmbiguous_(const char s);


    //-------------------------------------------------------------------------------------------------------------------------------------
    //  CDA - retrieval
    //-------------------------------------------------------------------------------------------------------------------------------------


    /**
     * @brief Follows the supply links of a node marked with term_flag to the root and stores the index of all substrings found in multi_seq_.
     */
    void findSubstrOfHit_();

    /**
     * @brief Checks if a sequence occurs more than once in the given peptide set and stores all indices in multi_seq_
     * @param pep_idx Index of found peptide
     */
    void multipleOccurrence_(const Size pep_idx );

    /**
     * @brief Do the retrieval of the given protein in BC array. If node is negative it calls compareTail_. If no transition is found it calls failure.
     * Calls itself with the new node and peptide position given by the failure_ or consumeAmbiguousAA_.
     * @return Returns true if hit was found and false if there is nothing more to found
     */
    bool retrieval_();

    /**
     * @brief Compares the peptide with the end of the sequences stored in Tail. If there is a mismatch failure and retrieval are called
     * @param tail_str Pointer to a given position in tail
     * @return True if a sequences was found and false if there is nothing more to found
     */
    bool compareTail_(const TailNode_* tail_str);

    /**
     * @brief Set new node_pos if there is no transition in trie
     * @param node pointer to supply link of the current node
     * @return false if no transition was found with the help of the supply link,
     * if no transition with the next aa was found and no spawn exists.
     * Otherwise true
     */
    bool failure_(const SupplyLink_* node);


    /**
     * @brief Returns the supply link if this does not eliminate ambiguous aa. Don't do anything special with 0. Must be considered when calling
     * @param node pointer to supply link
     * @param[out] output_node Is set to the linked node
     * @param amb_depth Depth of the first read ambiguous aa
     * @return False if a supply link would cause the first read aa to be dropped. Otherwise true
     */
    bool followSupplyLink_(const SupplyLink_* node, Int32& output_node, uint16_t& amb_depth);

    /**
     * @brief Processes ambiguous aa
     * @param aa The ambiguous aa
     * @return False if no spawn was found and no transition with the next aa was found.
     * Otherwise true
     */
    bool consumeAmbiguousAA_(const char aa);


    /**
     * @brief Returns the smallest and largest code value of amino acids that can stand for the given ambiguous aa
     * @param c Character of an ambiguous aa
     * @param[out] idxFirst Given by reference and set to the smallest code value
     * @param[out] idxLast Given by reference and set to the largest code value
     */
    void getSpawnRange_(const char c, uint8_t & idxFirst, uint8_t & idxLast);

    /**
     * @brief If spawns exist, the next spawn is set
     * @return False if there are no spwans, otherwise true
     */
    bool useSpawn_();


     /**
      * @brief Compute the child node for each node depending on whether it is a node in trie, tail or leaf node
      * @param node Current node
      * @param label Code of transition label
      * @param[out] output_node Child node
      * @param[out] sl_pos supply link of the given node
      * @return True if a child node was found else false
      */
    bool getChildNode_(const Int32 node, const uint8_t label,  Int32& output_node, SupplyLink_*& sl_pos);

    /**
     * @brief Calls getChildNode_. If false, the supply links are followed to 0 and the transition is tested for each node.
     * @param node Current node
     * @param label Code of label for transition
     * @param[out] output_node Child node if found
     * @param amb_depth Depth of the first read ambiguous aa
     * @return True if child node is found else false
     */
    bool getNextNode_(const Int32 node, const uint8_t label,  Int32& output_node, uint16_t & amb_depth);



    /**
     * @brief Tests if a Peptide ends in node
     * @param node node to be tested
     * @return True if node is a hit else false
     */
    bool isHit_(const Int32 node);


    /// Returns depth of the given node
    uint8_t getDepth_(const Int32 node);

    AhoCorasickDA::SupplyLink_* getSupplyLink_(Int32 node);

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

    /// deletes a used empty element from linkage
    void arrangeEmpty_(const UInt32 node_pos);

    /// stores the suffix string in tail as chars
    void appendToTail_(const char* suffix);

    /// Construct the failure function
    void constructFailure_();


    /**
     * @brief compute child node for failure function. If node is postive and not a leaf node child node is computed with child_ function.
     * if node is a leaf node or negative it returns negative tail index if label is correct
     * @param node node
     * @param label code of transition label
     * @return Either the child node or 0 if child doesn't exist
     */
    Int32 getNode_(const Int32 node, const uint8_t label);

    /**
     * @brief Set the given supply link for node. If node negative in failure_tail, if positiv in failure_bc
     * @param node the node number in bc array if positiv, index in tail if negative
     * @param sl The computed supply link node
     */
    void setSupplyLink_(const Int32 node, Int32 sl, const uint8_t depth);

    /// returns the supply link for node either from failure_tail if node is negative or from failure_bc if positive
    Int32 getSupplyLinkConstr_(const Int32 node);


    /// tests if a supply link leads to a hit and marks the node with a term_flag
    void isEndOfSeq_(const Int32 node);

    /// tests if a sequence ends in node. During construction node+1 in tail must be used
    bool isHitConstr_(const Int32 node);
  };


} // namespace OpenMS
