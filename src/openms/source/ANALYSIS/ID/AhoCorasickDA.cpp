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

    if (arrays)
    {
      String t(tail_.begin(), tail_.end());
      String base;
      String check;
      String leaf;
      String term;

      for (auto i : bc_)
      {
        base += std::to_string(i.base);
        base += " ";

        leaf += std::to_string(i.leaf_flag);
        term += std::to_string(i.term_flag);
        check += static_cast<char>(i.lcheck);
      }

      std::cout << "BASE " << base << std::endl;
      std::cout << "CHECK " << check << std::endl;
      std::cout << "LEAF " << leaf << std::endl;
      std::cout << "TERM " << term << std::endl;
      std::cout << "TAIL " << t << std::endl;


      std::cout << "Parent nodes: ";
      for (auto k: build_info_ )
      {
        std::cout << k.check <<" ";
      }
      std::cout << std::endl;

      std::cout << "Failure BC: ";
      for (auto l : failure_bc_)
      {
        std::cout << l <<" ";
      }
      std::cout << std::endl;

      std::cout << "Failure TAIL: ";
      for (auto m : failure_tail_)
      {
        std::cout << m << " ";
      }
      std::cout << std::endl;
    }


    /*for( const auto& n : ac_da.output_ ) {
        std::cout << "Key:[" << n.first << "] Value:[" << n.second << "]\n";
      }*/

    UInt32 count = 0;
    UInt32 out = 0;
    for (String& i: sequences_)
    {
      const char* c = i.c_str();
      setProtein(c);

      if (retrievalCDA())
      {
        count++;
        if (i == sequences_[(output_[node_pos_])])
        {
          out++;
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

    if (empty_head_ != 0) {
      auto pos = empty_head_ ;
      do
      {
        pos = build_info_[pos].next;
        ++num_emps_;
      } while (empty_head_  != pos);
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
      if(bc_[node_pos_].lcheck != *protein_)
      {
        return false;
      }
      protein_++;
    }
    if (*protein_ == '\0')
    {
      protein_++;
      return bc_[node_pos_].term_flag == 1;
    }
    auto tail_str = &tail_[bc_[node_pos_].base];
    node_pos_ = -bc_[node_pos_].base;
    while (*tail_str != '\0' && *tail_str == *protein_)
    {
      ++tail_str;
      ++protein_;
      --node_pos_;
    }
    ++node_pos_;
    return *tail_str == *protein_;
  }





  //-------------------------------------------------------------------------------------------------------------------------------------
  //  CDA
  //-------------------------------------------------------------------------------------------------------------------------------------




  UInt32 AhoCorasickDA::child_(const UInt32& base_value, const char& label)
  {
    return base_value ^ code_(label);
    //return base_value + code_(label);
  }

  UInt32 AhoCorasickDA::code_(const char& label)
  {
    return code_tabel_[static_cast<UInt>(label)];
  }


  void AhoCorasickDA::setProtein(const char* prot)
  {
    protein_ = prot;
    prot_pos_ = 0;
    node_pos_ = 0;
  }

  bool AhoCorasickDA::testAmbiguous(const char aa)
  {
    if (!isAmbiguous_(aa))
    {
      return true;
    }

    Spawn cur_spawn;
    
    // B = D, N
    if (code_(aa) == 22)
    {
      cur_spawn.a_acids.push_back(code_('D'));
      cur_spawn.a_acids.push_back(code_('N'));
    }

    // J = I, L
    if (code_(aa) == 24)
    {
      cur_spawn.a_acids.push_back(code_('I'));
      cur_spawn.a_acids.push_back(code_('L'));
    }

    // Z = E, Q
    if (code_(aa) == 23)
    {
      cur_spawn.a_acids.push_back(code_('E'));
      cur_spawn.a_acids.push_back(code_('Q'));
    }

    if (code_(aa) == 25)
    {
      cur_spawn.a_acids.resize(22);
      std::iota(cur_spawn.a_acids.begin(), cur_spawn.a_acids.end(), 0);
    }

    while (!cur_spawn.a_acids.empty())
    {
      Int32 child_pos = getNode_(node_pos_, cur_spawn.a_acids[0]);
      cur_spawn.a_acids.erase(cur_spawn.a_acids.begin());
      if (child_pos != 0)
      {
        if (!cur_spawn.a_acids.empty())
        {
          cur_spawn.node = node_pos_;
          cur_spawn.prot_pos = prot_pos_;
          spawns_.push_back(cur_spawn);
        }
        node_pos_ = child_pos;
        protein_++;
        prot_pos_++;
        return true;
      }
    }
    return false;
  }

  bool AhoCorasickDA::retrieval_()
  {
    char* tail_str;
    if (node_pos_ >= 0)
    {
      while (!bc_[node_pos_].leaf_flag)
      {
        if (bc_[node_pos_].term_flag == 1)
        {
          return true;
        }
        UInt32 next_node = child_(bc_[node_pos_].base, *protein_);
        if(bc_[next_node].lcheck != *protein_)
        {
          protein_++;
          prot_pos_++;
          return false;
        }
        node_pos_ = next_node;
        protein_++;
        prot_pos_++;
      }

      if (bc_[node_pos_].term_flag == 1)
      {
        protein_--;
        prot_pos_--;
        return true;
      }
      tail_str = &tail_[bc_[node_pos_].base];
      node_pos_ = -bc_[node_pos_].base;
    }
    else
    {
      tail_str = &tail_[-node_pos_];
    }

    while (*tail_str != '\0' && *tail_str == *protein_)
    {
      ++tail_str;
      ++protein_;
      --node_pos_;
      ++prot_pos_;
    }
    if (*tail_str == '\0')
    {
      ++node_pos_;
      --protein_;
      --prot_pos_;
      return true;
    }
    return false;
  }



  bool AhoCorasickDA::findNext(Size& pos_in_protein, Size& peptide_index)
  {
    if (strcmp(protein_, "") == 0 && prot_pos_ == 0)
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No protein for retrieval. Use function 'setProtein()'");
    }

    // after a hit is found, it is tested if the sequence occurs more than once
    if (prot_pos_ != 0)
    {
      UInt32 idx = static_cast<UInt32>(std::distance(sorted_idx_.begin(), std::find(sorted_idx_.begin(), sorted_idx_.end(), peptide_index)));
      if (idx+1 < sequences_.size() && sequences_[sorted_idx_[idx]] == sequences_[sorted_idx_[idx+1]])
      {
        peptide_index = sorted_idx_[idx+1];
        pos_in_protein = pos_in_protein;
        return true;
      }
    }

    while (*protein_ != '\0')
    {
      if (retrieval_())
      {
        peptide_index = output_[node_pos_];
        pos_in_protein = prot_pos_ + 1 - sequences_[peptide_index].size();
        protein_++;
        prot_pos_++;
        Int32 next_node = getNode_(node_pos_, code_(*protein_));
        if (next_node == 0)
        {
          node_pos_ = getSupplyLink_(node_pos_);
          // TODO for the retrieval it must be the next node in Tail
          if (node_pos_ < 0 )
          {
            node_pos_ = node_pos_-1;
          }
        }
        else
        {
          node_pos_ = next_node;
        }
        return true;
      }
      else
      {
        node_pos_ = getSupplyLink_(node_pos_);
        // TODO for the retrieval it must be the next node in Tail
        if (node_pos_ < 0 )
        {
          node_pos_ = node_pos_-1;
        }
      }
    }
    return false;
  }


  void AhoCorasickDA::construct_()
  {
    if (sequences_.empty())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No sequences given");
    }

    // sort strings

    // index vector of sorted pattern
    sorted_idx_.resize(sequences_.size());
    std::iota(sorted_idx_.begin(), sorted_idx_.end(), 0);
    std::sort(sorted_idx_.begin(), sorted_idx_.end(), [this](UInt32 idx1, UInt32 idx2){return sequences_[idx1].compare(sequences_[idx2]) < 0;});

    // unique vector
    std::unique_copy(sorted_idx_.begin(), sorted_idx_.end(), std::back_inserter(unique_idx_), [this](UInt32 idx1, UInt32 idx2){return sequences_[idx1]==sequences_[idx2];});

    //DA

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


    // add new empty nodes
    for (uint32_t i = 0; i < 256; ++i)
    {
      bc_.push_back(BC{});
      build_info_.push_back(BuildInformation{});
    }


    // set linkage
    for (auto i = old_size; i < new_size; i++)
    {
      build_info_[i].prev = i - 1;
      build_info_[i].next = i + 1;
    }


    // add new node to bc and set empty head
    if (empty_head_ == 0)
    {
      build_info_[old_size].prev = new_size - 1;
      build_info_[new_size - 1].next = old_size;
      empty_head_ = old_size;
    }

    // add new nodes to empty linkage
    else
    {
      auto emp_tail = build_info_[empty_head_].prev;
      build_info_[old_size].prev = emp_tail;
      build_info_[emp_tail].next = old_size;
      build_info_[empty_head_].prev = new_size - 1;
      build_info_[new_size - 1].next = empty_head_;
    }
  }

  void AhoCorasickDA::buildTrie_(Size begin, const Size& end, const Size& depth, const UInt32& node_pos)
  {
    // reaches end of a string
    if (sequences_[unique_idx_[begin]].size() == depth)
    {
      bc_[node_pos].term_flag = 1;
      output_[node_pos] = unique_idx_[begin];
      ++begin;

      // string is not substring of another string
      if (begin == end)
      {
        bc_[node_pos].base = 0;
        bc_[node_pos].leaf_flag = 1;
        return;
      }
    }

    // only one string in subtrie -> store in tail
    if (begin+1 == end) // && !bc_[node_pos].term_flag)
    {
      bc_[node_pos].leaf_flag = 1;
      bc_[node_pos].base = static_cast<UInt32>(tail_.size());
      storeInTail_(sequences_[unique_idx_[begin]].substr(depth));
      output_[-(tail_.size()-2)] = unique_idx_[begin];
      return;
    }

    // creates vector containing all transition labels
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
    edges_.push_back(current_label);


    // return a base value which is vaild for all labels in edges
     auto base_value = xCheck_();
    //auto base_value = plusCheck_();

    // set values
    bc_[node_pos].base = base_value;

    build_info_[base_value].base_flag = 1;


    for (auto label : edges_)
    {
      UInt32 child_pos = child_(base_value, label);
      // need more nodes  //TODO other possibility ?
      if (child_pos >= bc_.size())
      {
        expand_();
      }
      arrangeEmpty_(child_pos);
      bc_[child_pos].lcheck = label;
      build_info_[child_pos].check = node_pos;
    }

    // calls buildTrie_ for each subtrie
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
    // if no empty head the nodes are attached at the end of the DA
    if (empty_head_ == 0)
    {
      expand_();
      auto base_value = static_cast<UInt32>(bc_.size());
      while (!checkAllEdges_(base_value))
      {
        base_value++;
      }
      return base_value;
    }
    // try every empty element in DA
    UInt32 pos = empty_head_;
    do
    {
      UInt32 base_value = child_(pos, edges_[0]);
      if (checkAllEdges_(base_value))
      {
        return base_value;
      }
      pos = build_info_[pos].next;
    } while (empty_head_ != pos);

    // if no empty element can be used the nodes are attached at the end of the DA
    auto base_value = static_cast<UInt32>(bc_.size());
    expand_();
    while (!checkAllEdges_(base_value))
    {
      base_value++;
    }
    return base_value;
  }


  bool AhoCorasickDA::checkAllEdges_(const UInt32& base)
  {
    // tests if base value was already used
    if (build_info_[base].base_flag)
    {
      return false;
    }
    // tests for each child in the subtrie if the calculated node position is available
    for (auto label : edges_)
    {
      UInt32 pos = child_(base, label);
      if (build_info_[pos].used_flag && pos < bc_.size())
      {
        return false;
      }
    }
    return true;
  }


  void AhoCorasickDA::arrangeEmpty_(const UInt32& node_pos)
  {
    build_info_[node_pos].used_flag = 1;
    // removes the node from the link by linking the previous empty node and the following empty node
    auto prev = build_info_[node_pos].prev;
    auto next = build_info_[node_pos].next;
    build_info_[next].prev = prev;
    build_info_[prev].next = next;
    // if the empty node was the empty head, then empty head is either set to the next empty node or to 0 (there is no more empty node in the DA)
    if (node_pos == empty_head_)
    {
      if (node_pos == next)
      {
        empty_head_ = 0;
      }
      else
      {
        empty_head_ = next;
      }
    }
  }

  void AhoCorasickDA::storeInTail_(const String& suffix)
  {
     const char* str = suffix.c_str();
     while (*str != '\0')
     {
       if (isAmbiguous_(*str))
       {
         throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Input peptide must NOT contain ambiguous amino acids (B/J/Z/X) or characters other than amino acids!");
       }
      tail_.push_back(*str++);
     }
     tail_.push_back('\0');
  }

  bool AhoCorasickDA::isAmbiguous_(const char& aa)
  {
    //between B (22) and X (25) (B,Z,J,X)
    return (22 <= code_(aa) && code_(aa) <= 25);
  }

  void AhoCorasickDA::constructFailure_()
  {
    // create empty failure vectors of size bc and tail
    for (Size i = 0; i < bc_.size(); i++)
    {
      failure_bc_.push_back(0);
    }

    for (Size j = 0; j < tail_.size(); j++)
    {
      failure_tail_.push_back(0);
    }

    std::vector<Int32>nodes;
    std::vector<Int32>next_nodes;

    //TODO need label (getNode)
    //edges_ = {'A', 'C', 'D', 'E', 'F', 'G', 'H','I', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'Y'};


    // compute child nodes of current node_pos
    for (UInt label = 0; label <= 21; label++)
    {
      Int32 node_pos = getNode_(0, label);

      if (node_pos != 0)
      {
        nodes.push_back(node_pos);
      }
    }

    while(!nodes.empty())
    {

     for (Int32 node : nodes)
     {
       for (UInt label = 0; label <= 21; label++)
       {
         Int32 node_pos = getNode_(node, label);

         if (node_pos != 0)
         {
           next_nodes.push_back(node_pos);
           // supply link of parent of node_pos
           Int32 down = getSupplyLink_(node);

           while ((down != 0) && (getNode_(down, label) == 0)) //
           {
             down = getSupplyLink_(down);
           }
           // if (down != 0)
           //{
           setSupplyLink_(node_pos,getNode_(down, label));
           //}
           //else
           //{
           // no suffix exists: point suffix link of current node to root
           //setSupplyLink_(node_pos, 0);
           //}
         }
       }
     }
     nodes = next_nodes;
     next_nodes.clear();
    }
  }

  Int32 AhoCorasickDA::getNode_(const Int32& node, const UInt& label)
  {
    // node is in tail, child is next index if label match
    if (node < 0)
    {
      if (code_(tail_[abs(node)+1]) == label)
      {
        return node -1;
      }
      return 0;
    }
    // child is tail[base] if label match, returns negative base value because child node is in tail
    if (bc_[node].leaf_flag)
    {
      if (code_(tail_[bc_[node].base]) == label)
      {
        return -bc_[node].base;
      }
      return 0;
    }

    // child is in bc
    else
    {
      UInt32 c_node = bc_[node].base ^ label;
      if (code_(bc_[c_node].lcheck) == label)
      {
        return c_node;
      }
      else
      {
        return 0;
      }
    }
  }

  void AhoCorasickDA::setSupplyLink_(const Int32& node, const Int32& sl)
  {
    // supply link for node in tail, absolute value as index
    if (node < 0)
    {
      failure_tail_[abs(node)] = sl;
    }
    else
    {
      failure_bc_[node] = sl;
    }
  }

  Int32 AhoCorasickDA::getSupplyLink_(const Int32& node)
  {
    // supply link of node in tail, absolute value as index
    if ( node < 0)
    {
      return failure_tail_[-node];
    }
    else
    {
      return failure_bc_[node];
    }
  }


  UInt32 AhoCorasickDA::plusCheck_()
  {
    // if no empty head
    if (empty_head_ == 0)
    {
      return static_cast<UInt32>(bc_.size());
    }
    UInt32 pos = empty_head_;
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
    } while (empty_head_ != pos);

    return static_cast<UInt32>(bc_.size());
  }

} //namespace OpenMS
