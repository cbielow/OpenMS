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

#include <OpenMS/ANALYSIS/ID/AhoCorasickDA.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Types.h>

#include <algorithm>
#include <math.h>
#include <iostream>
#include <numeric>
#include <iterator>


namespace OpenMS
{

  //-------------------------------------------------------------------------------------------------------------------------------------
  //  Debugging
  //-------------------------------------------------------------------------------------------------------------------------------------

  void AhoCorasickDA::printDA(bool arrays)
  {
    std::cout << std::endl;
    std::cout << "BC size: " << bc_.size() << std::endl;
    std::cout << "TAIL size: " << tail_.size() << std::endl;

    for (Size i = 0; i < tail_.size(); i++)
    {
      if ( i == 189706 || i == 189705|| i == 189704)
      {
        std::cout << i << " " << tail_[i].label+0 << " " << failure_tail_[i].depth<< std::endl;
      }

    }

    for (Size i = 0; i < bc_.size(); i++)
    {
      if ( i == 58509 || i == 55151  || i == 59629 || i == 60855)
      {
        std::cout << i << " " << bc_[i].lcheck+0  << " " << failure_bc_[i].depth << std::endl;
      }

    }



    if (arrays)
    {
      String base;
      String check;
      String leaf;
      String term;
      String termTail;
      String tail;

      for (auto i : bc_)
      {
        base += std::to_string(i.base);
        base += " ";

        leaf += std::to_string(i.leaf_flag);
        term += std::to_string(i.term_flag);
        check += std::to_string(i.lcheck);
        check += " ";
      }

      std::cout << "BASE " << base << std::endl;
      std::cout << "CHECK " << check << std::endl;
      std::cout << "LEAF " << leaf << std::endl;
      std::cout << "TERM " << term << std::endl;



      for (auto t: tail_ )
      {
        tail += std::to_string(t.label);
        tail += " ";
        termTail += std::to_string(t.term_flag);
        termTail += " ";
      }

      std::cout << "Tail: " << tail << std::endl;
      std::cout << "Tail term:" << termTail << std::endl;


      std::cout << "Parent nodes: ";
      for (auto k: build_info_ )
      {
        std::cout << k.check <<" ";
      }
      std::cout << std::endl;

      std::cout << "Failure BC: ";
      for (auto l : failure_bc_)
      {
        std::cout << l.link <<" ";
      }
      std::cout << std::endl;

      std::cout << "Failure TAIL: ";
      for (auto m : failure_tail_)
      {
        std::cout << m.link << " ";
      }
      std::cout << std::endl;

      std::cout << "depth BC: ";
      for (Size i = 0; i < failure_bc_.size(); i++)
      {
        std::cout << i << "-" << failure_bc_[i].depth <<" ";
      }
      std::cout << std::endl;

      std::cout << "depth TAIL: ";
      for (Size i = 0; i < failure_tail_.size(); i++)
      {
        std::cout << i << "-" << failure_tail_[i].depth << " ";
      }
      std::cout << std::endl;

      for( const auto& n : output_ )
      {
        std::cout << "Key:[" << n.first << "] Value:[" << n.second << "]\n";
      }
    }




    UInt32 count = 0;
    UInt32 out = 0;
    for (const String& i: sequences_)
    {
      const char* c = i.c_str();

      setProtein(c, 0);

      if (retrievalCDA())
      {
        ++count;
        if (i == sequences_[(output_[node_pos_])])
        {
          ++out;
        }
        else
        {
          std::cout << i << std::endl;
          std::cout << sequences_[(output_[node_pos_])] << "   " << "  node: " << node_pos_ << std::endl;
        }
      }
      else
      {
        std::cout << i << std::endl;
      }
    }


    UInt32 num_emps_ = 0;

    if (first_empty_elem_ != 0) {
      auto pos = first_empty_elem_ ;
      do
      {
        //std::cout << pos << std::endl;
        pos = build_info_[pos].next;
        ++num_emps_;
      } while (first_empty_elem_  != pos);
    }

    UInt32 sub_str = 0;
    for (auto l : bc_)
    {
      if (l.leaf_flag)
      {
        sub_str += l.term_flag;
      }

    }
    std::cout << "num_emps_ " << num_emps_ << std::endl;
    std::cout << "sub_str " << sub_str << std::endl;
    std::cout << "# Strings: " << sequences_.size() << "  -  Found: " << count << " Out: " << out << std::endl;

  }

  bool AhoCorasickDA::retrievalCDA()
  {
    while (!bc_[node_pos_].leaf_flag)
    {
      if (*protein_ == '\0')
      {
        return bc_[node_pos_].term_flag == 1;
      }
      node_pos_ = child_(bc_[node_pos_].base, *protein_);
      if(bc_[node_pos_].lcheck != code_(*protein_))
      {
        return false;
      }
      ++protein_;
    }
    if (*protein_ == '\0')
    {
      ++protein_;
      return bc_[node_pos_].term_flag == 1;
    }
    auto tail_str = &tail_[bc_[node_pos_].base];
    node_pos_ = -bc_[node_pos_].base;
    while (tail_str->label != 0 && tail_str->label == code_(*protein_))
    {
      ++tail_str;
      ++protein_;
      --node_pos_;
    }
    return tail_str->label == code_(*protein_);
  }





  //-------------------------------------------------------------------------------------------------------------------------------------
  //  CDA
  //-------------------------------------------------------------------------------------------------------------------------------------




  UInt32 AhoCorasickDA::child_(const UInt32 base_value, const char label)
  {
    return child_(base_value, code_(label));
  }

  UInt32  AhoCorasickDA::child_(const UInt32 base_value, const uint8_t label)
  {
    return base_value ^ label;
    //return base_value + label;
  }

  uint8_t AhoCorasickDA::code_(const char label)
  {
    return code_table_[static_cast<uint8_t>(label)];
  }

  bool AhoCorasickDA::isAmbiguous_(const char aa)
  {
    //between B (23) and X (26) (B,Z,J,X)
    return (23 <= code_(aa) && code_(aa) <= 26 );
  }

  //void AhoCorasickDA::setProtein(const String& prot, const uint8_t amb_max)
  void AhoCorasickDA::setProtein(const char* prot, const uint8_t amb_max)
  {
    protein_ = prot;
    prot_pos_ = 0;
    node_pos_ = 0;
    amb_max_ = amb_max;
    amb_count_ = 0;
    max_depth_decrease_ = 0;
    pos_max_ = 0;
  }

  bool AhoCorasickDA::findNext(Size& pos_in_protein, Size& peptide_index)
  {
    // Exception if no protein has been set
    if (strcmp(protein_, "") == 0 && prot_pos_ == 0)
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No protein for retrieval. Use function 'setProtein()'");
    }

    // Returns hits if more than one hit was found at the node position
    if (!multi_hit_.empty())
    {
      peptide_index = multi_hit_.back();
      pos_in_protein = prot_pos_  - sequences_[peptide_index].size();
      multi_hit_.pop_back();
      return true;
    }


    // A transition must be performed without testing for a hit, because findNext is either called with node_pos_ 0,
    // which is never a hit, or is still in state because a hit was found and otherwise always the same hit is
    // reported.
    // makes only the next valid transition, no test for hit

    if (node_pos_ != 0)
    {
      if (isAmbiguous_(*protein_))
      {
        if (!consumeAmbiguousAA_(*protein_))
        {
          return false;
        }
      }
      else
      {
        if (getNextNode_(node_pos_, code_(*protein_), node_pos_, max_depth_decrease_))
        {
          if (*protein_ != '\0')
          {
            ++prot_pos_;
            ++protein_;
          }
        }
        else
        {
          if (!useSpawn_())
          {
            return false;
          }
        }
      }
    }


    while (retrieval_())
    {

      // Node itself is a hit
      if (output_[node_pos_] != SUBSTR_IS_HIT)
      {
        peptide_index = output_[node_pos_];
        pos_in_protein = prot_pos_ - sequences_[peptide_index].size();
        multipleOccurrence_(peptide_index);
        findSubstrOfHit_();
      }

      // Node itself is not a hit but a supply link leads to a hit
      else
      {
        findSubstrOfHit_();

        // Returns the last hit, all other hits are returned when calling findNext again
        if (!multi_hit_.empty())
        {
          peptide_index = multi_hit_.back();
          pos_in_protein = prot_pos_ - sequences_[peptide_index].size();
          multi_hit_.pop_back();
        }

        // In the case of ambiguous aa, the supply link leading to a hit may drop the first ambiguous aa and nothing is returned
        else
        {
          return findNext(pos_in_protein, peptide_index);
        }
      }
      return true;
    }
    return false;
  }


  //-------------------------------------------------------------------------------------------------------------------------------------
  //  CDA - retrieval
  //-------------------------------------------------------------------------------------------------------------------------------------

  void AhoCorasickDA::findSubstrOfHit_()
  {

    Int32 sl_node = 0;
    // Global variable should not be changed
    uint16_t amb_depth = max_depth_decrease_;
    Int32 node = node_pos_;

    // Follows the supply links to root
    while ((node != 0) && followSupplyLink_(node, sl_node, amb_depth))
    {
      if (isHit_(sl_node))
      {
        // Test if node is an output node or just marked that a substr is a hit
        if (output_[sl_node] != SUBSTR_IS_HIT)
        {
          multi_hit_.push_back(output_[sl_node]);
          multipleOccurrence_(output_[sl_node]);
        }
      }
      node = sl_node;
    }
  }

  void AhoCorasickDA::multipleOccurrence_(const Size pep_idx )
  {
    UInt32 idx = static_cast<UInt32>(std::distance(sorted_idx_.begin(), std::find(sorted_idx_.begin(), sorted_idx_.end(), pep_idx)));
    while (idx+1 < sequences_.size() && sequences_[sorted_idx_[idx]] == sequences_[sorted_idx_[idx+1]])
    {
      multi_hit_.push_back(sorted_idx_[idx+1]);
      ++idx;
    }
  }


  bool AhoCorasickDA::retrieval_()
  {
    TailNode_ * tail_str;

    // Node is in BC array
    if (node_pos_ >= 0)
    {
      while (!bc_[node_pos_].leaf_flag)
      {
        if (bc_[node_pos_].term_flag == 1)
        {
          return true;
        }

        if (*protein_ == '\0')
        {
          if (failure_())
          {
            return retrieval_();
          }
          return false;
        }

        if (isAmbiguous_(*protein_))
        {
          if (consumeAmbiguousAA_(*protein_))
          {
            return retrieval_();
          }
          return false;
        }

        // Compute child node
        UInt32 next_node = child_(bc_[node_pos_].base, *protein_);
        if(bc_[next_node].lcheck != code_(*protein_))
        {
          if (failure_())
          {
            return retrieval_();
          }
          return false;
        }

        node_pos_ = next_node;
        ++protein_;
        ++prot_pos_;
      }

      // Node is a leaf (following node is in Tail)
      if (bc_[node_pos_].term_flag == 1)
      {
        return true;
      }

      if (isAmbiguous_(*protein_))
      {
        if (consumeAmbiguousAA_(*protein_))
        {
          return retrieval_();
        }
        return false;
      }

      // Find tail position
      if (tail_[bc_[node_pos_].base].label == code_(*protein_))
      {
        tail_str = &tail_[bc_[node_pos_].base];
        node_pos_ = -bc_[node_pos_].base;
      }
      else
      {
        if (failure_())
        {
          return retrieval_();
        }
        return false;
      }
    }

    // Node_pos_ is negative which means node is in Tail
    // This is called after either a new spawn is consumed or after a hit
    else
    {
      // The current node position was determined by getChildNode_ and is a valid node.
      // Before the next node you have to test if it is a terminal node. The protein position is already one position further.
      tail_str = &tail_[-node_pos_];
      // If it is a terminal node, then the node position must be one position further
      // Because of the query in Tail, the node is always one position further.
      --node_pos_;
      if (tail_str->term_flag == 1)
      {
        return true;
      }
      ++tail_str;
    }

    return compareTail_(tail_str);
  }


  bool AhoCorasickDA::compareTail_(const TailNode_* tail_str)
  {

    // As long as the aa are the same, comparing and testing if a substr ends
    while (tail_str->label != 0 && tail_str->label == code_(*protein_))
    {
      // If it is a terminal node, then the node position must be one position further, because because of the query in tail, the node is always one position further.
      ++protein_;
      ++prot_pos_;
      --node_pos_;

      if (tail_str->term_flag == 1)
      {
        return true;
      }

      ++tail_str;

    }

    // 0 marks end of sequence
    if (tail_str->label == 0)
    {
      return true;
    }


    // Ambiguous aa always result in a mismatch, because no amb aa must occur in the patterns
    if (isAmbiguous_(*protein_))
    {
      if (consumeAmbiguousAA_(*protein_))
      {
        return retrieval_();
      }
      return false;
    }

    if (failure_())
    {
      return retrieval_();
    }
    return false;
  }


  bool AhoCorasickDA::failure_()
  {

    if (*protein_ == '\0')
    {
      return useSpawn_();
    }

    Int32 cur_node = node_pos_;
    Int32 sl_node;
    node_pos_ = 0;

    if (isAmbiguous_(*protein_))
    {
      return consumeAmbiguousAA_(*protein_);
    }

    // Tests for each supply link if there is a vaild trasition, if it reachs root, try next aa
    while (followSupplyLink_(cur_node, sl_node, max_depth_decrease_))
    {
      // Test valid transition from supply link node
      if (getChildNode_(sl_node, code_(*protein_), cur_node))
      {
        ++protein_;
        ++prot_pos_;
        node_pos_ = cur_node;
        return true;
      }

      // No futher supply links
      if (sl_node == 0 )
      {
        // As long as there is a next aa use next aa
        if (*protein_ != '\0')
        {
          ++protein_;
          ++prot_pos_;
          if (isAmbiguous_(*protein_))
          {
            return consumeAmbiguousAA_(*protein_);
          }
        }
        // No transition found in the remaining protein
        else
        {
          break;
        }
      }
      cur_node = sl_node;
    }

    return useSpawn_();
  }

  bool AhoCorasickDA::followSupplyLink_(Int32 node, Int32& output_node, uint16_t& amb_depth)
  {
    // Get the node the supply link leads to
    Int32 sl = 0;
    if (node < 0)
    {
      sl = failure_tail_[-node].link;
    }
    else
    {
      sl = failure_bc_[node].link;
    }

    // If the supply link would cause the first read ambiguous aa to be dropped, the supply link is not valid.
    // Else the number of aa before the first ambiguous aa is updated (Depth)
    if (amb_count_ > 0)
    {


      uint8_t up_count = getDepth_(node) - getDepth_(sl);
      if (up_count > amb_depth || ((amb_depth == 0 ) && (up_count == 0)))
      {
        output_node = 0;
        return false;
      }
      amb_depth -= up_count;
    }
    output_node = sl;
    return true;
  }

  bool AhoCorasickDA::consumeAmbiguousAA_(const char aa)
  {

    Spawn_ cur_spawn;

    // Maximum permitted count of ambiguous aa is reached
    if (amb_count_ == amb_max_)
    {

      // If the current aa has not yet been read, it forgets everything before and starts again from the root.
      // If it has been read, the spawn dies here and the next one is used
      if (pos_max_ < prot_pos_ )
      {
        pos_max_ = prot_pos_;
        node_pos_ = 0;
        amb_count_ = 0;
        max_depth_decrease_ = 0;
        return true;
      }

      return useSpawn_();
    }
    // Depth is only set for the first amb aa
    if (amb_count_ == 0)
    {
       max_depth_decrease_ = getDepth_(node_pos_);
    }
    ++amb_count_;
    uint8_t first = 0;
    uint8_t last = 0;
    getSpawnRange_(aa, first, last);

    //-----------------------------------------------------------------------------------------------------------------
    // TODO:
    //  At the moment the query for BC, Leaf or Tail in 'getChildNode_' in 'getNextNode_'
    //  is done every time for the same node. Same problem with getNode_ in failure array construction
    //-----------------------------------------------------------------------------------------------------------------

    // Search child of node_pos_ with all aa for which the given amb aa can represent
    for (uint8_t a = first; a <= last; ++a)
    {
      uint16_t amb_depth = max_depth_decrease_;
      Int32 child_pos;
      if (getNextNode_(node_pos_, a, child_pos, amb_depth))
      {
        cur_spawn.node = child_pos;
        cur_spawn.prot_pos = prot_pos_+1;
        cur_spawn.max_depth_decrease = amb_depth;
        cur_spawn.counter = amb_count_;
        spawns_.push_back(cur_spawn);
      }
    }

    // When the amb aa is seen to the first amn, the following aa is started from the root. Otherwise this has already been done and a spawn is used.
    if (pos_max_ < prot_pos_ || (pos_max_ == 0 && prot_pos_ == 0))
    {
      pos_max_ = prot_pos_;

      node_pos_ = 0;
      ++protein_;
      ++prot_pos_;
      amb_count_ = 0;
      max_depth_decrease_ = 0;
      return true;
    }
    else
    {
      return useSpawn_();
    }

  }


  void AhoCorasickDA::getSpawnRange_(const char c, uint8_t & idxFirst, uint8_t & idxLast)
  {
    static const uint8_t jump[4][2] = { { code_('D'), code_('N') },   // B = D,N
                                        { code_('E'), code_('Q') },   // Z = E,Q
                                        { code_('L'), code_('I') },   // J = L,I
                                        { 1         , AA_COUNT_  } }; // X = A..V
    static const uint8_t ord_b = code_('B');
    idxFirst = jump[code_(c) - ord_b][0];
    idxLast = jump[code_(c) - ord_b][1];
  }

  bool AhoCorasickDA::useSpawn_()
  {
    if (!spawns_.empty())
    {
      auto &spawn = spawns_.back();
      node_pos_ = spawn.node;
      protein_ = &protein_[static_cast<Int32>(spawn.prot_pos - prot_pos_)];
      prot_pos_ = spawn.prot_pos;
      amb_count_ = spawn.counter;
      max_depth_decrease_  = spawn.max_depth_decrease;
      spawns_.pop_back();
      return true;
    }
    return false;
  }


  bool AhoCorasickDA::getChildNode_(const Int32 node, const uint8_t label,  Int32& output_node)
  {
    // Node is in tail, child is next index if label match
    if (node < 0)
    {
      // test for 0, because then the end of a sequence has already been reached
      if (tail_[-node].label != 0 && tail_[-node].label == label)
      {
        output_node = node ;
        return true;
      }
      output_node = 0;
      return false;
    }

    // Child is tail[base] if label match, returns negative base value because child node is in tail
    if (bc_[node].leaf_flag)
    {
      if (tail_[bc_[node].base].label == label)
      {
        output_node = -bc_[node].base;
        return true;
      }
      output_node = 0;
      return false;

    }

    // Child is in bc
    else
    {
      if (label != 0)
      {
        UInt32 c_node = child_(bc_[node].base, label);
        if (bc_[c_node].lcheck == label)
        {
          output_node = c_node;
          return true;
        }
      }
      output_node = 0;
      return false;
    }
  }


  bool AhoCorasickDA::getNextNode_(const Int32 node, const uint8_t label, Int32& output_node, uint16_t& amb_depth)
  {
    // Test for valid transition of current node position
    if (getChildNode_(node, label, output_node))
    {
      return true;
    }
    Int32 cur_node = node;

    // If child was not found, the supply links are followed to root and the transition is tested for each node.
    while (cur_node != 0 && followSupplyLink_(cur_node, output_node, amb_depth))
    {
      if (getChildNode_(output_node, label, cur_node))
      {
        output_node = cur_node;
        return true;
      }
      cur_node = output_node;
    }
    return false;
  }



  bool AhoCorasickDA::isHit_(const Int32 node)
  {
    if (node < 0)
    {
      return tail_[-node].label == 0;
    }
    return bc_[node].term_flag == 1;
  }

  uint8_t AhoCorasickDA::getDepth_(const Int32 node)
  {
    if ( node < 0)
    {
      return failure_tail_[-node].depth;
    }
    return failure_bc_[node].depth;
  }


  //-------------------------------------------------------------------------------------------------------------------------------------
  //  CDA - construction
  //-------------------------------------------------------------------------------------------------------------------------------------




  void AhoCorasickDA::construct_()
  {
    if (sequences_.empty())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No sequences given");
    }

    // Sort strings

    // Index vector of sorted pattern
    sorted_idx_.resize(sequences_.size());
    std::iota(sorted_idx_.begin(), sorted_idx_.end(), 0);
    std::sort(sorted_idx_.begin(), sorted_idx_.end(), [this](UInt32 idx1, UInt32 idx2){return sequences_[idx1].compare(sequences_[idx2]) < 0;});

    // unique vector
    std::unique_copy(sorted_idx_.begin(), sorted_idx_.end(), std::back_inserter(unique_idx_), [this](UInt32 idx1, UInt32 idx2){return sequences_[idx1]==sequences_[idx2];});

    // Double-array construction
    Size init_capa = 1;
    while (init_capa < unique_idx_.size())
    {
      init_capa <<= 1;
    }

    bc_.reserve(init_capa);
    build_info_.reserve(init_capa);
    edges_.reserve(22); // 23 with node method


    expand_();
    arrangeEmpty_(0);
    buildTrie_(0, unique_idx_.size(), 0, 0);
    constructFailure_();
  }

  void AhoCorasickDA::expand_()
  {
    auto old_size = static_cast<uint32_t>(bc_.size());
    auto new_size = old_size + 256;


    // Exception size of DA becomes to large to store
    if (new_size >= 1<<22)
    {
      throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, new_size);
    }


    // Add new empty nodes
    bc_.resize(new_size);


    // Set build informations and linkage
    for (auto i = old_size; i < new_size; ++i)
    {
      build_info_.emplace_back(i);
    }


    // Add the new nodes to linkage
    if (first_empty_elem_ == 0)
    {
      build_info_[old_size].prev = new_size - 1;
      build_info_[new_size - 1].next = old_size;
      first_empty_elem_ = old_size;
    }


    else
    {
      auto emp_tail = build_info_[first_empty_elem_].prev;
      build_info_[old_size].prev = emp_tail;
      build_info_[emp_tail].next = old_size;
      build_info_[first_empty_elem_].prev = new_size - 1;
      build_info_.back().next = first_empty_elem_;
    }
  }

  void AhoCorasickDA::buildTrie_(Size begin, const Size end, const Size depth, const UInt32 node_pos)
  {
    // Reaches end of a string
    if (sequences_[unique_idx_[begin]].size() == depth)
    {
      bc_[node_pos].term_flag = 1;
      output_[node_pos] = unique_idx_[begin];
      ++begin;

      // String is not substring of another string
      if (begin == end)
      {
        bc_[node_pos].base = 0;
        bc_[node_pos].leaf_flag = 1;
        return;
      }
    }

    // Only one string in subtrie -> store in tail
    if (begin+1 == end) // && !bc_[node_pos].term_flag)
    {
      bc_[node_pos].leaf_flag = 1;
      bc_[node_pos].base = static_cast<UInt32>(tail_.size());
      appendToTail_(sequences_[unique_idx_[begin]].substr(depth));
      output_[-(tail_.size()-1)] = unique_idx_[begin];
      return;
    }

    // Creates vector containing all transition labels
    edges_.clear();

    char current_label = sequences_[unique_idx_[begin]][depth];

    for (auto i = begin; i < end; ++i)
    {
      char next_label = sequences_[unique_idx_[i]][depth];
      if (current_label != next_label)
      {
        if (isAmbiguous_(current_label))
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Input peptide must NOT contain ambiguous amino acids (B/J/Z/X) or characters other than amino acids!");
        }
        edges_.push_back(current_label);
        current_label = next_label;
      }
    }
    if (isAmbiguous_(current_label))
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Input peptide must NOT contain ambiguous amino acids (B/J/Z/X) or characters other than amino acids!");
    }
    edges_.push_back(current_label);


    // Return a base value which is valid for all labels in edges
     auto base_value = xCheck_();
    //auto base_value = plusCheck_();

    // Set values
    bc_[node_pos].base = base_value;

    build_info_[base_value].base_flag = 1;


    for (auto const label : edges_)
    {
      UInt32 child_pos = child_(base_value, label);
      // More nodes needed
      if (child_pos >= bc_.size())
      {
        expand_();
      }
      arrangeEmpty_(child_pos);
      bc_[child_pos].lcheck = code_(label);
      build_info_[child_pos].check = node_pos;
    }

    // Calls buildTrie_ for each subtrie
    Size begin_sub = begin;
    auto cur_label = sequences_[unique_idx_[begin]][depth];
    for (auto end_sub = begin +1; end_sub < end; ++end_sub)
    {
      auto n_label = sequences_[unique_idx_[end_sub]][depth];
      if (cur_label != n_label)
      {
        buildTrie_(begin_sub, end_sub, depth+1, child_(base_value, cur_label));
        cur_label = n_label;
        begin_sub = end_sub;
      }
    }
    buildTrie_(begin_sub, end, depth+1, child_(base_value, cur_label));
  }

  UInt32 AhoCorasickDA::xCheck_()
  {
    // If empty head doesn't exist the nodes are attached at the end of the bc array
    if (first_empty_elem_ == 0)
    {
      expand_();
      auto base_value = static_cast<UInt32>(bc_.size());
      while (!checkAllEdges_(base_value))
      {
        ++base_value;
      }
      return base_value;
    }
    // Try every empty element in DA
    UInt32 pos = first_empty_elem_;
    do
    {
      UInt32 base_value = child_(pos, edges_[0]);
      if (checkAllEdges_(base_value))
      {
        return base_value;
      }
      pos = build_info_[pos].next;
    } while (first_empty_elem_ != pos);

    // If the empty elements can't be used the nodes are attached at the end of the bc array
    auto base_value = static_cast<UInt32>(bc_.size());
    expand_();
    while (!checkAllEdges_(base_value))
    {
      ++base_value;
    }
    return base_value;
  }


  bool AhoCorasickDA::checkAllEdges_(const UInt32 base)
  {
    // Tests if base value was already used
    if (build_info_[base].base_flag)
    {
      return false;
    }
    // Tests for each child in the subtrie if the calculated node position is available
    for (auto label : edges_)
    {
      UInt32 pos = child_(base, label);
      // In case child_node is greater than the current array size, used_flag does not need to be tested, because the nodes could not be used yet
      if (build_info_[pos].used_flag && pos < bc_.size())
      {
        return false;
      }
    }
    return true;
  }


  void AhoCorasickDA::arrangeEmpty_(const UInt32 node_pos)
  {
    build_info_[node_pos].used_flag = 1;
    // Removes the node from the link by linking the previous empty node and the following empty node
    auto prev = build_info_[node_pos].prev;
    auto next = build_info_[node_pos].next;
    build_info_[next].prev = prev;
    build_info_[prev].next = next;

    // If the empty node was the empty head, then empty head is either set to the next empty node or to 0 (there is no more empty node in the bc array)
    if (node_pos == first_empty_elem_)
    {
      if (node_pos == next)
      {
        first_empty_elem_ = 0;
      }
      else
      {
        first_empty_elem_ = next;
      }
    }
  }

  void AhoCorasickDA::appendToTail_(const String& suffix)
  {
     const char* str = suffix.c_str();
     while (*str != '\0')
     {
       if (isAmbiguous_(*str))
       {
         throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Input peptide must NOT contain ambiguous amino acids (B/J/Z/X) or characters other than amino acids!");
       }
      tail_.emplace_back(TailNode_(code_(*str++)));
     }
     tail_.emplace_back(TailNode_(0));
  }


  void AhoCorasickDA::constructFailure_()
  {
    // Create empty failure vectors of size bc and tail
    failure_bc_.resize(bc_.size());
    failure_tail_.resize(tail_.size());

    // All nodes of depth d-1 which have already been processed
    std::vector<Int32> nodes;
    // Child nodes of nodes (depth d)
    std::vector<Int32> next_nodes;

    uint8_t depth = 1;


    // Create vector of nodes for depth 1 and set supply link and depth
    for (uint8_t label = 1; label <= AA_COUNT_; ++label)
    {
      Int32 node_pos = getNode_(0, label);

      if (node_pos != 0)
      {
        nodes.push_back(node_pos);
        setSupplyLink_(node_pos, 0, 1);
      }
    }

    // For all nodes in each depth of the trie
    while(!nodes.empty())
    {
     ++depth;

     // Depth is only used for the restriction of ambiguous aa
     if (depth > UINT8_MAX && amb_max_ > 0)
     {
       throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Depth of the tree cannot be stored. Set the maximum number of ambiguous amino acids allowed to 0 to continue with this data");
     }

     // For every child node of node the supply link is computed
     for (Int32 node : nodes)
     {
       // Compute child nodes of current node
       for (uint8_t label = 1; label <= AA_COUNT_; ++label)
       {
         Int32 child_pos = getNode_(node, label);


         if (child_pos != 0)
         {
           // Stores all child nodes as parent for next iteration
           next_nodes.push_back(child_pos);

           // Supply link of parent of child_pos
           Int32 down = getSupplyLink_(node);


           // Follow the supply links (longest suffix) and test transition with label of child_pos
           while ((down != 0) && (getNode_(down, label) == 0))
           {
             down = getSupplyLink_(down);
           }

           // Set supply link if transition found else supply link points to root
           setSupplyLink_(child_pos, getNode_(down, label), depth);

           // Test if this node is an end of a substr for retrieval and mark them so it just has to be done once for each node
           isEndOfSeq_(child_pos);
         }
       }
     }
     // Processed nodes become new parent nodes
     nodes = next_nodes;
     next_nodes.clear();
    }

  }

  Int32 AhoCorasickDA::getNode_(const Int32 node, const uint8_t label)
  {
    // Node is in tail, child is next index if label match
    if (node < 0)
    {
      if (tail_[-node+1].label == label )
      {
        return node -1;
      }
      return 0;
    }

    // Child is tail[base] if label match, returns negative base value because child node is in tail
    if (bc_[node].leaf_flag)
    {
      if (tail_[bc_[node].base].label == label)
      {
        return -bc_[node].base;
      }
      return 0;
    }

    // Child is in bc
    else
    {
      UInt32 c_node = child_(bc_[node].base, label);
      if (bc_[c_node].lcheck == label)
      {
        return c_node;
      }
      return 0;
    }
  }

  void AhoCorasickDA::setSupplyLink_(const Int32 node, Int32 sl, const uint8_t depth)
  {
    // In the tail query, the stored sequence is compared with the protein sequence letter by letter (no transtion anymore).

    // The supply link leads to the node up to which the sequence still matches. In tail the query must start with the following characters.
    if (sl < 0)
    {
      --sl;
    }

    // supply link for node in tail, absolute value as index
    if (node < 0)
    {
      // If a mismatch occures, the node position is already the position of the mismatch. The supply links are calculated so that they go from the last valid node
      // to the node with the longest matching suffix. In the tail, a mismatch would always require going back one node to find the correct supply link.
      // Therefore the supply link for each node is stored one position further (at the node the missmatch occures).
      failure_tail_[-node+1].link = sl;

      // If a label is read that forms the transition between a leaf node and a tail node, the depth of this label would be the depth of the leaf node and not
      // the depth of the tail node. However, the depth is only determined after the transition in the tail node.
      // Therefore depth-1 is calculated for each tail node.
      failure_tail_[-node].depth = depth-1;
      // Because the depth is shifted to the 'left' and therefore the leaf and tail node have the same depth, the end node must also have a depth
      if (tail_[-node+1].label == 0)
      {
        failure_tail_[-node+1].depth = depth;
      }
    }

    // BC array
    else
    {
      failure_bc_[node].link = sl;
      failure_bc_[node].depth = depth;
    }
  }

  Int32 AhoCorasickDA::getSupplyLink_(const Int32 node)
  {
    // Because of setting in setSupplyLink_
    Int32 sl;

    if ( node < 0)
    {
      sl =  failure_tail_[-node+1].link;
    }
    else
    {
      sl = failure_bc_[node].link;
    }
    if (sl < 0)
    {
      ++sl;
    }
    return sl;
  }

  void AhoCorasickDA::isEndOfSeq_(const Int32 node)
  {
    // Is a hit, the retrieval already searches for substrings
    if (isHitConstr_(node))
    {
      return;
    }

    Int32 sl = getSupplyLink_(node);

    // Follow supply links to root a check for a peptide end
    while (sl != 0)
    {
      if (isHitConstr_(sl))
      {
        if (node < 0)
        {
          tail_[-node].term_flag = 1;
          // Because of the query in tail, the node is always one position further.
          output_[node-1] = SUBSTR_IS_HIT;
        }
        else
        {
          bc_[node].term_flag = 1;
          output_[node] = SUBSTR_IS_HIT;
        }
        return;
      }
      sl = getSupplyLink_(sl);
    }
  }

  bool AhoCorasickDA::isHitConstr_(const Int32 node)
  {
    if (node < 0)
    {
      // 0 marks the end of a sequence. If the next node is 0, the current aa is the last of the peptide
      return tail_[-node+1].label == 0;
    }
    return bc_[node].term_flag==1;
  }


  // Only test difference between + and XOR
  UInt32 AhoCorasickDA::plusCheck_()
  {
    // if no empty head
    if (first_empty_elem_ == 0)
    {
      return static_cast<UInt32>(bc_.size());
    }
    UInt32 pos = first_empty_elem_;
    do
    {
      if ((static_cast<Int32 >(pos) - static_cast<Int32 >(code_(edges_[0]))) >= 0)
      {
        UInt32 base_value = pos - code_(edges_[0]);
        if (checkAllEdges_(base_value))
        {
          return base_value;
        }
      }
      pos = build_info_[pos].next;
    } while (first_empty_elem_ != pos);

    return static_cast<UInt32>(bc_.size());
  }

} //namespace OpenMS
