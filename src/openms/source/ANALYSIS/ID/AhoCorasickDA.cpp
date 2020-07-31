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


namespace OpenMS
{


  UInt32 AhoCorasickDA::child_(const UInt32& base_value, const char& label)
  {
    return base_value ^ code_(label);
  }

  UInt32 AhoCorasickDA::code_(const char& label)
  {
    return code_tabel[static_cast<UInt>(label) - 65];
  }

  bool AhoCorasickDA::retrieval(const char *str)
  {
    UInt32 node_pos = 0;
    while (!bc[node_pos].leaf_flag)
    {
      if (*str == '\0')
      {
        if (!bc[node_pos].term_flag)
        {
          return true;
        }
        return false;
      }
      node_pos = child_(bc[node_pos].base, *str);
      if(bc[node_pos].lcheck != *str)
      {
        return false;
      }
      str++;
    }
    auto tail_str = &tail[bc[node_pos].base];
    while (*tail_str != '\0' && *tail_str == *str)
    {
      ++tail_str;
      ++str;
    }
    return *tail_str == *str;
  }


  void AhoCorasickDA::construct(std::vector<String>& seq)
  {
    if (seq.empty())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No sequences given");
    }

    // sort strings
    std::sort(seq.begin(), seq.end());
    pattern = seq;

    //DA

    // reasonable ?
    Size init_capa = 1;
    while (init_capa < pattern.size())
    {
      init_capa <<= 1;
    }

    bc.reserve(init_capa);
    build_info.reserve(init_capa);
    edges.reserve(22); // 23 with node method


    expand_();
    build_info[0].used_flag = 1;
    buildTrie_(0, pattern.size(), 0, 0);

  }

  void AhoCorasickDA::expand_()
  {
    auto old_size = static_cast<uint32_t>(bc.size());
    auto new_size = old_size + 22;


    // Exception size() of DA becomes to large to store?
    if (new_size >= 4194304)
    {
      throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, new_size);
    }

    // add new empty nodes
    for (uint32_t i = 0; i < 22; ++i)
    {
      bc.push_back(BC{});
      build_info.push_back(BuildInformation{});
    }

    // set linkage
    for (auto i = old_size; i < new_size; ++i)
    {
      build_info[i].prev = i - 1;
      build_info[i].next = i + 1;
    }

    // add new node to bc and set empty head
    if (empty_head == UINT32_MAX)
    {
      build_info[old_size].prev = new_size - 1;
      build_info[new_size - 1].next = old_size;
      empty_head = old_size;
    }

    //add new nodes to empty linkage
    else
    {
      auto emp_tail = build_info[empty_head].prev;
      build_info[old_size].prev = emp_tail;
      build_info[emp_tail].next = old_size;
      build_info[empty_head].prev = new_size - 1;
      build_info[new_size - 1].next = empty_head;
    }
  }

  void AhoCorasickDA::buildTrie_(Size begin, const Size& end, const Size& depth, const UInt32& node_pos)
  {
    // reaches end of a string

    if (pattern[begin].size() == depth)
    {
      ++begin;
      bc[node_pos].term_flag = 1;

 /*     if (begin +1 == end)
      {
        bc[node_pos].leaf_flag = 1;
        bc[node_pos].base = static_cast<UInt32>(tail.size());
        storeInTail_(pattern[begin].substr(depth));
        return;
      } */
    }


    if (begin +1 == end ) //&& !bc[node_pos].term_flag
    {
      bc[node_pos].leaf_flag = 1;
      bc[node_pos].base = static_cast<UInt32>(tail.size());
      storeInTail_(pattern[begin].substr(depth));
      return;
    }

    // creates vector containing all transition labels
    edges.clear();

    char current_label = pattern[begin][depth];

    for (auto i = begin; i < end; ++i)
    {
      char next_label = pattern[i][depth];
      if (current_label != next_label)
      {
        edges.push_back(current_label);
        current_label = next_label;
      }
    }
    edges.push_back(current_label);


    // return base for all labels in edges
    auto base_value = xCheck_();

    if (bc.size() <= base_value)
    {
      expand_();
    }


    // set values
    bc[node_pos].base = base_value;
    build_info[base_value].base_flag = 1;

    for (auto label : edges)
    {
      UInt32 child_pos = child_(base_value, label);
      arrangeEmpty_(child_pos);
      bc[child_pos].lcheck = label;
      build_info[child_pos].check = node_pos;
    }


    // calls buildTrie_ for each subtree
    Size begin_sub = begin;
    auto cur_label = pattern[begin][depth];
    for (auto end_sub = begin +1; end_sub < end; ++end_sub)
    {
      auto n_label = pattern[end_sub][depth];
      if (cur_label != n_label)
      {
        buildTrie_(begin_sub, end_sub, depth+1, child_(base_value, cur_label));
        cur_label = n_label;
        begin_sub = end_sub;
      }
    }
    buildTrie_(begin_sub, end, depth+1, child_(base_value, cur_label));
  }

  UInt32 AhoCorasickDA::mini_()
  {
    UInt32 min = child_(bc.size(), edges[0]);
    for (Size i = 1; i < edges.size(); i++)
    {
      if ( child_(bc.size(), edges[i])< min)
      {
        min = child_(bc.size(), edges[i]);
      }
    }
    return min;
  }

  UInt32 AhoCorasickDA::xCheck_()
  {
    // if no empty head
    if (empty_head == UINT32_MAX)
    {
      UInt32 min_index = mini_();
      return static_cast<UInt32>(bc.size()) - min_index;
    }
    UInt32 pos = empty_head;
   do
    {
      UInt32 base_value = child_(pos, edges[0]);
      if (checkAllEdges_(base_value))
      {
        return base_value;
      }
      pos = build_info[pos].next;
    } while (empty_head != pos);
    UInt32 min_index = mini_();
    return static_cast<UInt32>(bc.size()) - min_index;
  }


  bool AhoCorasickDA::checkAllEdges_(const UInt32& base)
  {
    if (build_info[base].base_flag)
    {
      return false;
    }
    for (auto label : edges)
    {
      UInt32 pos = child_(base, label);
      if (build_info[pos].used_flag)
      {
        return false;
      }
    }
    return true;
  }


  void AhoCorasickDA::arrangeEmpty_(const UInt32& node_pos)
  {
    build_info[node_pos].used_flag = 1;
    auto prev = build_info[node_pos].prev;
    auto next = build_info[node_pos].next;
    build_info[next].prev = prev;
    build_info[prev].next = next;
    if (node_pos == empty_head)
    {
      if (node_pos == next)
      {
        empty_head = UINT32_MAX;
      }
      else
      {
        empty_head = next;
      }
    }
  }

  void AhoCorasickDA::storeInTail_(const String& suffix)
  {
     const char* str = suffix.c_str();
     while (*str != '\0')
     {
      tail.push_back(*str++);
     }
     tail.push_back('\0');
  }




} //namespace OpenMS