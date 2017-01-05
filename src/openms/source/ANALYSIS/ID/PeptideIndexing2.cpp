// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Chris Bielow, Knut Reinert $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/PeptideIndexing2.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/DATASTRUCTURES/SeqanIncludeWrapper.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/SYSTEM/StopWatch.h>
#include <OpenMS/METADATA/PeptideEvidence.h>
#include <OpenMS/CHEMISTRY/EnzymesDB.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <fstream>
#include <iostream>
#include <algorithm>

using namespace OpenMS;
using namespace std;


struct PeptideProteinMatchInformation {
    /// index of the protein the peptide is contained in
    OpenMS::Size protein_index;

    /// index Peptide was found in
    OpenMS::Size indexType;

    /// the amino acid after the peptide in the protein
    char AABefore;

    /// the amino acid before the peptide in the protein
    char AAAfter;

    /// the position of the peptide in the protein
    OpenMS::Int position;

    bool operator<(const PeptideProteinMatchInformation &other) const {
        if (protein_index != other.protein_index) {
            return protein_index < other.protein_index;
        }
        else if (position != other.position) {
            return position < other.position;
        }
        else if (AABefore != other.AABefore) {
            return AABefore < other.AABefore;
        }
        else if (AAAfter != other.AAAfter) {
            return AAAfter < other.AAAfter;
        }
        else if (indexType != other.indexType){
            return indexType < other.indexType;
        }
        return false;
    }

    bool operator==(const PeptideProteinMatchInformation &other) const {
        return protein_index == other.protein_index &&
               position == other.position &&
               AABefore == other.AABefore &&
               AAAfter == other.AAAfter &&
               indexType == other.indexType;
    }

};

/* */
struct FMind {
    FMind() { }
};

/* */
struct SAind {
    SAind() { }
};

/* */
struct WOTDind {
    WOTDind() { }
};


namespace seqan
{

    struct FoundProteinFunctor {
    public:
        typedef OpenMS::Map<OpenMS::Size, std::set<PeptideProteinMatchInformation> > MapType;

        /// peptide index --> protein indices
        MapType pep_to_prot;

        /// number of accepted hits (passing addHit() constraints)
        OpenMS::Size filter_passed;

        /// number of rejected hits (not passing addHit())
        OpenMS::Size filter_rejected;

    private:
        EnzymaticDigestion enzyme_;

    public:
        explicit FoundProteinFunctor(const EnzymaticDigestion &enzyme) :
                pep_to_prot(), filter_passed(0), filter_rejected(0), enzyme_(enzyme) {
        }

        template<typename TIter1, typename TIter2>
        void operator()(const TIter1 &iter_pep, const TIter2 &iter_prot, OpenMS::Size indexType) {
            // the peptide sequence (will not change)
            const OpenMS::String tmp_pep(begin(representative(iter_pep)),
                                         end(representative(iter_pep)));

            // remember mapping of proteins to peptides and vice versa
            const OpenMS::Size count_occ = countOccurrences(iter_pep);
            for (OpenMS::Size i_pep = 0; i_pep < count_occ; ++i_pep) {
                const OpenMS::Size idx_pep = getOccurrences(iter_pep)[i_pep].i1;
                const OpenMS::Size count_occ_prot = countOccurrences(iter_prot);
                for (OpenMS::Size i_prot = 0; i_prot < count_occ_prot; ++i_prot) {
                    const seqan::Pair<int> prot_occ = getOccurrences(iter_prot)[i_prot];
                    // the protein sequence (will change for every Occurrence -- hitting
                    // multiple proteins)
                    const OpenMS::String tmp_prot(
                            begin(indexText(container(iter_prot))[getSeqNo(prot_occ)]),
                            end(indexText(container(iter_prot))[getSeqNo(prot_occ)]));
                    // check if hit is valid and add (if valid)
                    addHit(idx_pep, prot_occ.i1, tmp_pep, tmp_prot,
                           getSeqOffset(prot_occ),indexType);
                }
            }
        }

        bool checkAmbigous_(const OpenMS::String &tmp_pep, const OpenMS::String &tmp_prot, long &beg, const long unsigned int &max_err){
            OpenMS::Size error = 0;
            for (unsigned i = 0; i < tmp_pep.length(); i++){
                // check for ambigous AA
                if (tmp_pep[i] == 'B' || tmp_pep[i] == 'Z' || tmp_pep[i] == 'X') {
                    // check if position in Protein is different to found AAA
                    if (tmp_prot.at(beg + i) != tmp_pep[i]) {
                        // this is a real error
                        error++;
                    }
                }
            }
            if (error > max_err){
                return false;
            }
            return true;
        }

        void addHit(OpenMS::Size idx_pep,
                    OpenMS::Size idx_prot,
                    const OpenMS::String &seq_pep,
                    const OpenMS::String &protein,
                    OpenMS::Size position,
                    OpenMS::Size indexType) {
            if (enzyme_.isValidProduct(AASequence::fromString(protein), position,
                                       seq_pep.length(), true)) {
                PeptideProteinMatchInformation match;
                match.protein_index = idx_prot;
                match.position = position;
                match.AABefore = (position == 0) ? PeptideEvidence::N_TERMINAL_AA : protein[position - 1];
                match.AAAfter = (position + seq_pep.length() >= protein.size()) ? PeptideEvidence::C_TERMINAL_AA
                                                                                : protein[position + seq_pep.length()];
                match.indexType = indexType;
                pep_to_prot[idx_pep].insert(match);
                //std:: cout <<"Treffer hinzugefügt. \n \t Index des Proteins in der DB: "<<  idx_prot<< std::endl;
                //std:: cout <<"\t Index des PeptideHit aus der idXML: " << idx_pep<< std::endl;
                ++filter_passed;
            }
            else {
                //std::cerr << "REJECTED Peptide " << seq_pep << " with hit to protein "
                //  << protein << " at position " << position << std::endl;
                   ++filter_rejected;
            }
        }

        bool operator==(const FoundProteinFunctor &rhs) const {
            if (pep_to_prot.size() != rhs.pep_to_prot.size()) {
                LOG_ERROR << "Size " << pep_to_prot.size() << " "
                << rhs.pep_to_prot.size() << std::endl;
                return false;
            }

            MapType::const_iterator it1 = pep_to_prot.begin();
            MapType::const_iterator it2 = rhs.pep_to_prot.begin();
            while (it1 != pep_to_prot.end()) {
                if (it1->first != it2->first) {
                    LOG_ERROR << "Index of " << it1->first << " " << it2->first
                    << std::endl;
                    return false;
                }
                if (it1->second.size() != it2->second.size()) {
                    LOG_ERROR << "Size of " << it1->first << " " << it1->second.size()
                    << "--" << it2->second.size() << std::endl;
                    return false;
                }
                if (!equal(it1->second.begin(), it1->second.end(), it2->second.begin())) {
                    LOG_ERROR << "not equal set for " << it1->first << std::endl;
                    return false;
                }
                ++it1;
                ++it2;
            }
            return true;
        }

        bool operator!=(const FoundProteinFunctor &rhs) const {
            return !(*this == rhs);
        }

    };


    // saving some memory for the SA
    template<>
    struct SAValue<Index<StringSet<Peptide>, IndexWotd<> > > {
        typedef Pair<unsigned> Type;
    };

    template<typename T = void>
    struct EquivalenceClassAA_ {
        static unsigned const VALUE[24];
    };
    template<typename T>
    unsigned const EquivalenceClassAA_<T>::VALUE[24] =
            {
                    1, // 0 Ala Alanine (A)
                    2, // 1 Arg Arginine (R)
                    4, // 2 Asn Asparagine (N)
                    8, // 3 Asp Aspartic Acid (D)
                    16, // 4 Cys Cystine (C)
                    32, // 5 Gln Glutamine (Q)
                    64, // 6 Glu Glutamic Acid (E)
                    128, // 7 Gly Glycine (G)
                    256, // 8 His Histidine (H)
                    512, // 9 Ile Isoleucine (I)
                    1024, // 10 Leu Leucine (L)
                    2048, // 11 Lys Lysine (K)
                    4096, // 12 Met Methionine (M)
                    8192, // 13 Phe Phenylalanine (F)
                    16384, // 14 Pro Proline (P)
                    32768, // 15 Ser Serine (S)
                    65536, // 16 Thr Threonine (T)
                    131072, // 17 Trp Tryptophan (W)
                    262144, // 18 Tyr Tyrosine (Y)
                    524288, // 19 Val Valine (V)
                    // ambiguous AA's
                    4 + 8, //  Aspartic Acid (D), Asparagine(N) == (B)
                    32 + 64, // Glutamic Acid(E), Glutamine(Q) == (Z)
                    static_cast<unsigned>(-1),  // 22 Unknown (matches ALL)
                    static_cast<unsigned>(-1), // 23 Terminator (dummy)
            };
    template <bool enumerateA, bool enumerateB, typename TOnFoundFunctor,
            typename TTreeIteratorA, typename TIterPosA,
            typename TTreeIteratorB, typename TIterPosB, typename TErrors>
    inline void _approximateAminoAcidTreeSearch(TOnFoundFunctor& onFoundFunctor,
                                                TTreeIteratorA iterA,
                                                TIterPosA iterPosA,
                                                TTreeIteratorB iterB_,
                                                TIterPosB iterPosB,
                                                TErrors errorsLeft, // usually 0 for our case
                                                TErrors classErrorsLeft) // ambiguous AA's allowed
    {
        if (enumerateA && !goDown(iterA)) return;

        if (enumerateB && !goDown(iterB_)) return;

        do
        {
            //std::cout << "B: " << repLength(iterB_) << "..." <<  representative(iterB_) << std::endl;
            TTreeIteratorB iterB = iterB_;
            do
            {
                TErrors e = errorsLeft;
                TErrors ec = classErrorsLeft;
                TIterPosA ipA = iterPosA;
                TIterPosB ipB = iterPosB;

                while (true)
                {
                    if (ipA == repLength(iterA))
                    {
                        if (isLeaf(iterA))
                        {
                            onFoundFunctor(iterA, iterB);
                        }
                        else if (ipB == repLength(iterB) && !isLeaf(iterB))
                        {
                            _approximateAminoAcidTreeSearch<true, true>(
                                    onFoundFunctor, iterA, ipA, iterB, ipB, e, ec);
                        }
                        else
                        {
                            _approximateAminoAcidTreeSearch<true, false>(
                                    onFoundFunctor, iterA, ipA, iterB, ipB, e, ec);
                        }
                        break;
                    }
                    else
                    {
                        if (ipB == repLength(iterB))
                        {
                            if (!isLeaf(iterB))
                            {
                                _approximateAminoAcidTreeSearch<false, true>(
                                        onFoundFunctor, iterA, ipA, iterB, ipB, e, ec);
                            }
                            break;
                        }
                    }

                    if (_charComparator(representative(iterA)[ipA],
                                        representative(iterB)[ipB],
                                        EquivalenceClassAA_<char>::VALUE))
                    {
                        // matched (including character classes) - look at ambiguous AA in
                        // PROTEIN tree (peptide tree is not considered!)
                        const char x_prot = representative(iterB)[ipB];
                        if ((x_prot == 'X') || (x_prot == 'B') || (x_prot == 'Z'))
                        {
                            if (ec == 0) break;
                            --ec; // decrease class error tokens
                        }

                        // dealing with 'X' in peptide sequence: only match exactly 'X' in
                        // proteinDB, not just any representative of 'X' (this is how X!
                        // Tandem would report results)
                        const char x_pep = representative(iterA)[ipA];
                        if ((x_pep == 'X') || (x_pep == 'B') || (x_pep == 'Z'))
                        {
                            if (x_pep != x_prot) break;
                        }
                    }
                    else // real mismatches (e.g., when checking for SNP's)
                    {
                        if (e == 0) break;
                        --e;
                    }

                    ++ipA;
                    ++ipB;
                }
            }
            while (enumerateB && goRight(iterB));
        }
        while (enumerateA && goRight(iterA));
    }

    template <typename TEquivalenceTable>
    inline bool _charComparator(AminoAcid charA, AminoAcid charB,
                                TEquivalenceTable equivalence)
    {
        const unsigned a_index = ordValue(charA);
        const unsigned b_index = ordValue(charB);
        return (equivalence[a_index] & equivalence[b_index]) != 0;
    }

    template <typename TContainer, typename TFunctor, typename TIterTag, typename TParallelTag>
    inline void iterate(TContainer & c, TFunctor f, Tag<TIterTag> const & iterTag, Tag<TParallelTag> const &)
    {
        typedef Tag<TIterTag> const                                     TIterSpec;
        typedef typename Iterator<TContainer, TIterSpec>::Type          TIter;

        for (TIter it = begin(c, iterTag); !atEnd(it, c); ++it)
            f(it);
    }

    template <typename TIndexIt, typename TNeedle, typename TNeedleIt,
              typename TThreshold, typename TDelegate, typename TDistance>
    inline void _findBacktracking(TIndexIt indexIt,
                      TNeedle const & needle,
                      TNeedleIt needleIt,
                      TThreshold errors,
                      TThreshold threshold,
                      TDelegate && delegate,
                      TDistance)
    {
        // Exact case.
        if (errors == threshold)
        {
            //std::cout << suffix(needle, position(needleIt, needle)) << std::endl;
            if (goDown(indexIt, suffix(needle, position(needleIt, needle)))){

                //std::cout << "rep: " << representative(indexIt) << std::endl;
                delegate(indexIt, errors);
            }

        }
        // Approximate case.
        else if (errors < threshold)
        {
            // Base case.
            if (atEnd(needleIt, needle))
            {
                delegate(indexIt, errors);
            }
            // Recursive case.
            else
            {
                // Insertion.
                if (IsSameType<TDistance, EditDistance>::VALUE)
                {
                    _findBacktracking(indexIt, needle, needleIt + 1,
                                      static_cast<TThreshold>(errors + 1), threshold, delegate, TDistance());
                }

                if (goDown(indexIt))
                {
                    do
                    {
                        //std::cout << "rep: " << representative(indexIt) << std::endl;
                        // Mismatch.
                        TThreshold delta = !ordEqual(parentEdgeLabel(indexIt), value(needleIt));
                        _findBacktracking(indexIt, needle, needleIt + 1,
                                          static_cast<TThreshold>(errors + delta), threshold, delegate, TDistance());

                        // Deletion.
                        if (IsSameType<TDistance, EditDistance>::VALUE)
                        {
                            _findBacktracking(indexIt, needle, needleIt,
                                              static_cast<TThreshold>(errors + 1), threshold, delegate, TDistance());
                        }
                    }
                    while (goRight(indexIt));
                }
            }
        }
    }
    template <typename TState, typename TIndex, typename TNeedle,
              typename TThreshold, typename TDelegate, typename TDistance, typename TSpec>
    //SEQAN_FUNC_ENABLE_IF(And<Is<StringTrieConcept<TIndex> >, IsSequence<TNeedle> >, void)
    void _findImpl(TState & /* indexIt */,
              TIndex & index,
              TNeedle const & needle,
              TThreshold threshold,
              TDelegate && delegate,
              Backtracking<TDistance, TSpec>)
    {
        typedef typename Iterator<TIndex, TopDown<> >::Type       TIndexIt;
        typedef typename Iterator<TNeedle const, Standard>::Type  TNeedleIt;

        TIndexIt indexIt(index);
        TNeedleIt needleIt = begin(needle, Standard());
        TThreshold errors = 0;

        _findBacktracking(indexIt, needle, needleIt, errors, threshold, delegate, TDistance());
    }


    template <typename TText, typename TPattern>
    struct DefaultFind
    {
        typedef void Type;
    };

    // ----------------------------------------------------------------------------
    // Metafunction FindState_<>
    // ----------------------------------------------------------------------------

    template <typename TText, typename TPattern, typename TAlgorithm>
    struct FindState_
    {
        typedef Nothing Type;
    };

    // ----------------------------------------------------------------------------
    // Metafunction HasStatesPool_<>
    // ----------------------------------------------------------------------------

    template <typename TState, typename TThreading>
    struct HasStatesPool_
    {
        typedef typename IsSameType<TThreading, Parallel>::Type         IsParallel;
        typedef typename Not<IsSameType<TState, Nothing> >::Type        IsStateful;
        typedef typename And<IsParallel, IsStateful>::Type              Type;
    };

    // ----------------------------------------------------------------------------
    // Metafunction StatesPool_<>
    // ----------------------------------------------------------------------------

    template <typename TState, typename TThreading>
    struct StatesPool_
    {
        typedef typename HasStatesPool_<TState, TThreading>::Type           HasStatesPool;
        typedef typename If<HasStatesPool, String<TState>, TState>::Type    Type;
    };

    // ============================================================================
    // Classes
    // ============================================================================

    // ----------------------------------------------------------------------------
    // Class FindDelegator_
    // ----------------------------------------------------------------------------

    template <typename TDelegate, typename TNeedlesIt>
    struct FindDelegator_
    {
        TDelegate & delegate;
        TNeedlesIt & needlesIt;

        FindDelegator_(TDelegate & delegate, TNeedlesIt & needlesIt) :
                delegate(delegate),
                needlesIt(needlesIt)
        {}

        template <typename TTextIt, typename TScore>
        void operator()(TTextIt & textIt, TScore score)
        {
            delegate(textIt, needlesIt, score);
        }
    };

    // ----------------------------------------------------------------------------
    // Function find(text, pattern, errors, [](...){}, Algorithm());
    // ----------------------------------------------------------------------------

    template <typename TText, typename TPattern, typename TThreshold, typename TDelegate, typename TAlgorithm>
    inline void find(TText & text, TPattern const & pattern, TThreshold threshold, TDelegate && delegate, TAlgorithm)
    {
        typedef typename FindState_<TText, TPattern const, TAlgorithm>::Type    TState;

        TState state;
        _findStateInit(state, text, pattern, threshold, TAlgorithm());
        _findImpl(state, text, pattern, threshold, delegate, TAlgorithm());
    }

    // ----------------------------------------------------------------------------
    // Function find(text, pattern, [](...){});
    // ----------------------------------------------------------------------------

    template <typename TText, typename TPattern, typename TDelegate>
    inline void find(TText & text, TPattern const & pattern, TDelegate && delegate)
    {
        typedef typename DefaultFind<TText, TPattern const>::Type   TAlgorithm;

        find(text, pattern, unsigned(), delegate, TAlgorithm());
    }

    // ----------------------------------------------------------------------------
    // Function _findStateInit()
    // ----------------------------------------------------------------------------
    // Stateless algorithm.

    template <typename TState, typename TText, typename TPattern, typename TThreshold, typename TAlgorithm>
    inline void
    _findStateInit(TState & /* s */, TText const & /* t */, TPattern const & /* p */, TThreshold /* t */, TAlgorithm) {}

    // ----------------------------------------------------------------------------
    // Function _findStatesPoolInit()
    // ----------------------------------------------------------------------------

    template <typename TState>
    inline void _findStatesPoolInit(TState & /* pool */, False) {}

    template <typename TStates>
    inline void _findStatesPoolInit(TStates & pool, True)
    {
        resize(pool, omp_get_num_threads(), Exact());
    }

    // ----------------------------------------------------------------------------
    // Function _findPickState()
    // ----------------------------------------------------------------------------

    template <typename TState>
    inline TState &
    _findPickState(TState & state, False)
    {
        return state;
    }

    template <typename TStates>
    inline typename Value<TStates>::Type &
    _findPickState(TStates & pool, True)
    {
        return pool[omp_get_thread_num()];
    }

    // ----------------------------------------------------------------------------
    // Function find(text, needles, errors, [](...){}, Algorithm(), Parallel());
    // ----------------------------------------------------------------------------

    template <typename TText, typename TNeedle_, typename TSSetSpec,
            typename TThreshold, typename TDelegate, typename TAlgorithm, typename TThreading>
    inline void find(TText & text,
                     StringSet<TNeedle_, TSSetSpec> const & needles,
                     TThreshold threshold,
                     TDelegate && delegate,
                     TAlgorithm,
                     TThreading)
    {
        typedef StringSet<TNeedle_, TSSetSpec> const                    TNeedles;
        typedef typename Value<TNeedles>::Type                          TNeedle;
        typedef typename Iterator<TNeedles, Rooted>::Type               TNeedlesIt;
        typedef typename Reference<TNeedlesIt>::Type                    TNeedleRef;

        typedef typename FindState_<TText, TNeedle, TAlgorithm>::Type   TState;
        typedef typename StatesPool_<TState, TThreading>::Type          TStatesPool;
        typedef typename HasStatesPool_<TState, TThreading>::Type       HasStatesPool;

        typedef FindDelegator_<TDelegate, TNeedlesIt const>             TDelegator;

        TStatesPool pool;
        _findStatesPoolInit(pool, HasStatesPool());

        iterate(needles, [&](TNeedlesIt const & needlesIt)
                {
                    TState & state = _findPickState(pool, HasStatesPool());
                    TNeedleRef needle = value(needlesIt);
                    TDelegator delegator(delegate, needlesIt);

                    _findStateInit(state, text, needle, threshold, TAlgorithm());
                    _findImpl(state, text, needle, threshold, delegator, TAlgorithm());
                },
                Rooted(), TThreading());
    }

    // ----------------------------------------------------------------------------
    // Function find(text, needles, errors, [](...){}, Algorithm());
    // ----------------------------------------------------------------------------

    template <typename TText, typename TNeedle, typename TSSetSpec,
            typename TThreshold, typename TDelegate, typename TAlgorithm>
    inline void find(TText & text,
                     StringSet<TNeedle, TSSetSpec> const & needles,
                     TThreshold threshold,
                     TDelegate && delegate,
                     TAlgorithm)
    {
        find(text, needles, threshold, delegate, TAlgorithm(), Serial());
    }




} // end namespace seqan

PeptideIndexing2::PeptideIndexing2() :
        DefaultParamHandler("PeptideIndexing2") {

    defaults_.setValue("decoy_string", "DECOY_", "String that was appended (or prefixed - see 'decoy_string_position' flag below) to the accessions in the protein database to indicate decoy proteins.");

    defaults_.setValue("decoy_string_position", "prefix", "Should the 'decoy_string' be prepended (prefix) or appended (suffix) to the protein accession?");
    defaults_.setValidStrings("decoy_string_position", ListUtils::create<String>("prefix,suffix"));

    defaults_.setValue("missing_decoy_action", "error", "Action to take if NO peptide was assigned to a decoy protein (which indicates wrong database or decoy string): 'error' (exit with error, no output), 'warn' (exit with success, warning message)");
    defaults_.setValidStrings("missing_decoy_action", ListUtils::create<String>("error,warn"));

    defaults_.setValue("enzyme:name", "Trypsin", "Enzyme which determines valid cleavage sites - e.g. trypsin cleaves after lysine (K) or arginine (R), but not before proline (P).");

    StringList enzymes;
    EnzymesDB::getInstance()->getAllNames(enzymes);
    defaults_.setValidStrings("enzyme:name", enzymes);

    defaults_.setValue("enzyme:specificity", EnzymaticDigestion::NamesOfSpecificity[0], "Specificity of the enzyme."
                                                                                                "\n  '" + EnzymaticDigestion::NamesOfSpecificity[0] + "': both internal cleavage sites must match."
                                                                                                "\n  '" + EnzymaticDigestion::NamesOfSpecificity[1] + "': one of two internal cleavage sites must match."
                                                                                                "\n  '" + EnzymaticDigestion::NamesOfSpecificity[2] + "': allow all peptide hits no matter their context. Therefore, the enzyme chosen does not play a role here");

    StringList spec;
    spec.assign(EnzymaticDigestion::NamesOfSpecificity, EnzymaticDigestion::NamesOfSpecificity + EnzymaticDigestion::SIZE_OF_SPECIFICITY);
    defaults_.setValidStrings("enzyme:specificity", spec);

    defaults_.setValue("write_protein_sequence", "false", "If set, the protein sequences are stored as well.");
    defaults_.setValidStrings("write_protein_sequence", ListUtils::create<String>("true,false"));

    defaults_.setValue("write_protein_description", "false", "If set, the protein description is stored as well.");
    defaults_.setValidStrings("write_protein_description", ListUtils::create<String>("true,false"));

    defaults_.setValue("keep_unreferenced_proteins", "false", "If set, protein hits which are not referenced by any peptide are kept.");
    defaults_.setValidStrings("keep_unreferenced_proteins", ListUtils::create<String>("true,false"));

    defaults_.setValue("allow_unmatched", "false", "If set, unmatched peptide sequences are allowed. By default (i.e. if this flag is not set) the program terminates with an error on unmatched peptides.");
    defaults_.setValidStrings("allow_unmatched", ListUtils::create<String>("true,false"));

    defaults_.setValue("aaa_max", 4, "Maximal number of ambiguous amino acids (AAAs) allowed when matching to a protein database with AAAs. AAAs are 'B', 'Z' and 'X'");
    defaults_.setMinInt("aaa_max", 0);

    defaults_.setValue("mismatches_max", 0, " Maximal number of real mismatches (will be used after checking for ambiguous AA's (see 'aaa_max' option). In general this param should only be changed if you want to look for other potential origins of a peptide which might have unknown SNPs or the like.");
    defaults_.setMinInt("mismatches_max", 0);

    defaults_.setValue("IL_equivalent", "false", "Treat the isobaric amino acids isoleucine ('I') and leucine ('L') as equivalent (indistinguishable)");
    defaults_.setValidStrings("IL_equivalent", ListUtils::create<String>("true,false"));

    defaults_.setValue("suffix_array", "false", "The index specified under -index is a suffix array");
    defaults_.setValidStrings("suffix_array", ListUtils::create<String>("true,false"));

    defaults_.setValue("FM_index", "false", "The index specified under -index is an FM-index");
    defaults_.setValidStrings("FM_index", ListUtils::create<String>("true,false"));

    defaults_.setValue("log", "", "Name of log file (created only when specified)");
    defaults_.setValue("debug", 0, "Sets the debug level");

    defaultsToParam_();
}

PeptideIndexing2::~PeptideIndexing2() {

}

void PeptideIndexing2::writeLog_(const String &text) const {
    LOG_INFO << text << endl;
    if (!log_file_.empty()) {
        log_ << QDateTime::currentDateTime().toString("yyyy-MM-dd hh:mm:ss").toStdString() << ": " << text << endl;
    }
}

void PeptideIndexing2::writeDebug_(const String &text, const Size min_level) const {
    if (debug_ >= min_level) {
        writeLog_(text);
    }
}

void PeptideIndexing2::updateMembers_() {
    decoy_string_ = static_cast<String>(param_.getValue("decoy_string"));
    prefix_ = (param_.getValue("decoy_string_position") == "prefix" ? true : false);
    missing_decoy_action_ = static_cast<String>(param_.getValue("missing_decoy_action"));
    enzyme_name_ = static_cast<String>(param_.getValue("enzyme:name"));
    enzyme_specificity_ = static_cast<String>(param_.getValue("enzyme:specificity"));

    write_protein_sequence_ = param_.getValue("write_protein_sequence").toBool();

    write_protein_description_ = param_.getValue("write_protein_description").toBool();
    keep_unreferenced_proteins_ = param_.getValue("keep_unreferenced_proteins").toBool();
    allow_unmatched_ = param_.getValue("allow_unmatched").toBool();
    IL_equivalent_ = param_.getValue("IL_equivalent").toBool();

    aaa_max_ = static_cast<Size>(param_.getValue("aaa_max"));
    mismatches_max_ = static_cast<Size>(param_.getValue("mismatches_max"));

    suffix_array_ = param_.getValue("suffix_array").toBool();
    FM_index_ = param_.getValue("FM_index").toBool();

    log_file_ = param_.getValue("log");
    debug_ = static_cast<Size>(param_.getValue("debug")) > 0;
}

PeptideIndexing2::ExitCodes PeptideIndexing2::checkUserInput_(
                std::vector<ProteinIdentification> &prot_ids,
        std::vector<PeptideIdentification> &pep_ids) {


    // check if Peptides to search for are empty
    if (pep_ids.empty()) // Aho-Corasick requires non-empty input
    {
        LOG_WARN << "Warning: An empty set of peptide identifications was provided. Output will be empty as well." <<
        std::endl;
        if (!keep_unreferenced_proteins_) {
            // delete only protein hits, not whole ID runs incl. meta data:
            for (vector<ProteinIdentification>::iterator it = prot_ids.begin();
                 it != prot_ids.end(); ++it) {
                it->getHits().clear();
            }
        }
        return PEPTIDE_IDS_EMPTY;
    }


    // check for Parameter Error
    if (mismatches_max_ > 4) // this would talke to much time!
    {
        LOG_ERROR <<
        "Fatal error: Searching for Peptides allowing > 4 mistakes would blow up search time.\n"
        << "Please adapt your settings." << endl;
        return ILLEGAL_PARAMETERS;
    }

    if (suffix_array_ && FM_index_)
    {
        LOG_ERROR <<
                  "Fatal error: Cannot use more than one format to load index.\n"
                  << "Please adapt your settings." << endl;
        return ILLEGAL_PARAMETERS;
    }

    if (!search_for_aaa_proteins_ && !search_for_normal_proteins_){
        LOG_ERROR <<
                  "Fatal error: Cannot search for nothing.\n"
                  << "Please adapt your settings." << endl;
        return ILLEGAL_PARAMETERS;
    }


    return CHECKPOINT_OK;
}


PeptideIndexing2::ExitCodes PeptideIndexing2::buildPepDB_(seqan::StringSet<seqan::Peptide> &pep_DB,
                                                          std::vector<PeptideIdentification> &pep_ids) {
    for (vector<PeptideIdentification>::const_iterator it1 = pep_ids.begin(); it1 != pep_ids.end(); ++it1) {
        //String run_id = it1->getIdentifier();
        vector<PeptideHit> hits = it1->getHits();
        for (vector<PeptideHit>::iterator it2 = hits.begin(); it2 != hits.end(); ++it2) {
            String seq = it2->getSequence().toUnmodifiedString().remove('*');
            if (IL_equivalent_) // convert  L to I; warning: do not use 'J', since Seqan does not know about it and will convert 'J' to 'X'
            {
                seq.substitute('L', 'I');
            }
            appendValue(pep_DB, seq.c_str());
        }
    }
    return CHECKPOINT_OK;
}

PeptideIndexing2::ExitCodes PeptideIndexing2::buildReversePepDB_(seqan::StringSet<seqan::Peptide> &pep_DB, std::vector<PeptideIdentification> &pep_ids) {
    for (vector<PeptideIdentification>::const_iterator it1 = pep_ids.begin(); it1 != pep_ids.end(); ++it1) {
        //String run_id = it1->getIdentifier();
        vector<PeptideHit> hits = it1->getHits();
        for (vector<PeptideHit>::iterator it2 = hits.begin(); it2 != hits.end(); ++it2) {
            String seq = it2->getSequence().toUnmodifiedString().reverse().remove('*');
            if (IL_equivalent_) // convert  L to I; warning: do not use 'J', since Seqan does not know about it and will convert 'J' to 'X'
            {
                seq.substitute('L', 'I');
            }
            appendValue(pep_DB, seq.c_str());
        }
    }
    return CHECKPOINT_OK;
}


PeptideIndexing2::ExitCodes PeptideIndexing2::readAcc_to_prot_(Map<String, Size> &acc_to_prot, String path){
    ifstream infile;
    infile.open(path.c_str(), std::ios::in);
    if (!infile.is_open()){
        return INPUT_ERROR;
    }
    String line;
    while(getline(infile, line)) // To get you all the lines.
    {
        //getline(infile, line);
        std::vector<String> substr;
        line.split(';',substr);
        // catch Format error
        if(substr.size() != 2){
            return INPUT_ERROR;
        }
        acc_to_prot[substr.at(0)] = (Size)substr.at(1).toInt();
    }
    infile.close();
    return CHECKPOINT_OK;
}

PeptideIndexing2::ExitCodes PeptideIndexing2::readProteins_(std::vector<FASTAFile::FASTAEntry> &proteins, String path){
    // read Proteins
    fstream infile;
    infile.open((path + "_proteins").c_str(), std::ios::in);
    if (!infile.is_open()){
        return INPUT_ERROR;
    }
    String line;
    while(getline(infile, line)) // To get you all the lines.
    {
        //getline(infile, line);
        std::vector<String> substr;
        line.split(';',substr);
        if(substr.size() < 2){
            writeLog_(String("Format Error in proteins file. Please rebuild index."));
            return INPUT_ERROR;
        }
        FASTAFile::FASTAEntry fastaEnt;
        fastaEnt.identifier = substr.at(0);
        fastaEnt.sequence = substr.at(1);
        if(substr.size() == 3){
            fastaEnt.description = substr.at(2);
        } else {
            fastaEnt.description = "";
        }
        proteins.push_back(fastaEnt);
    }
    infile.close();

    return CHECKPOINT_OK;
}

PeptideIndexing2::ExitCodes PeptideIndexing2::loadInfo_(std::vector<FASTAFile::FASTAEntry> &proteins,
                                                        std::vector<FASTAFile::FASTAEntry> &proteinsAAA,
                                                        Map<String, Size> &acc_to_prot,
                                                        Map<String, Size> &acc_to_AAAprot,
                                                        String &path,
                                                        String &pathAAA) {
    // read Acc_to_prot and proteins only if we search for normal proteins
    if (search_for_normal_proteins_) {
        if (readAcc_to_prot_(acc_to_prot, (path + "_acc_to_prot")) != CHECKPOINT_OK) {
            return INPUT_ERROR;
        };
        if (readProteins_(proteins, path) != CHECKPOINT_OK) {
            return INPUT_ERROR;
        };
    }
    // read _acc_to_AAAprot only if we search for AAA proteins
    if (search_for_aaa_proteins_) {
        if (readAcc_to_prot_(acc_to_AAAprot, (pathAAA + "_acc_to_prot")) != CHECKPOINT_OK) {
            return INPUT_ERROR;
        };
        if (readProteins_(proteinsAAA, pathAAA) != CHECKPOINT_OK) {
            return INPUT_ERROR;
        };
    }
    return CHECKPOINT_OK;
}

inline void PeptideIndexing2::searchWrapper_(seqan::FoundProteinFunctor &func_SA,
                                             seqan::Index<seqan::StringSet<seqan::Peptide>, seqan::IndexSa<> > &prot_Index,
                                             seqan::StringSet<seqan::Peptide> &pep_DB,
                                             int mm,
                                             OpenMS::Size max_aaa,
                                             OpenMS::Size indexType){
    typedef seqan::Index<seqan::StringSet<seqan::Peptide>, seqan::IndexSa<> > TIndex;
    typedef typename seqan::Iterator<seqan::StringSet<seqan::Peptide> const, seqan::Rooted>::Type TPatternsIt;
    typedef typename seqan::Iterator<TIndex, seqan::TopDown<seqan::PreorderEmptyEdges> >::Type TIndexIt;
    auto delegate = [&func_SA, &prot_Index, max_aaa, indexType](TIndexIt const &it, TPatternsIt const &patternsIt, unsigned /*score*/) {
        //auto pattern_len = length(*patternsIt);
        //std::cout << "laenge: " << pattern_len << std::endl;
        long indexOfPeptideHit = position(patternsIt);
        //std::cout << "Mindestens ein Treffer in PeptideHit mit der Nummer: " << indexOfPeptideHit << std::endl;
        // loop through occurences of pattern in text
        for (auto &occ: getOccurrences(it)){

            long entryNr = seqan::getSeqNo(occ);
            //std::cout << "Treffer Nummer: " << countHit << " in der DB an Position: " << entryNr << std::endl;
            long beg = seqan::getSeqOffset(occ);
            //long end_ = beg + pattern_len;

            const OpenMS::String tmp_pep(begin(*patternsIt),end(*patternsIt));
            //std::cout << "Peptid: " << tmp_pep << std::endl;

            const OpenMS::String tmp_prot(begin(seqan::indexText(prot_Index)[entryNr]),
                                          end(seqan::indexText(prot_Index)[entryNr]));


            if(func_SA.checkAmbigous_(tmp_pep, tmp_prot, beg, max_aaa)){
               func_SA.addHit(indexOfPeptideHit,entryNr,tmp_pep,tmp_prot,beg,indexType);
            }
        }
    };
    find(prot_Index, pep_DB, mm, delegate, seqan::Backtracking<seqan::EditDistance>(), seqan::Parallel());
}

inline void PeptideIndexing2::searchWrapper_(seqan::FoundProteinFunctor &func_SA,
                                             seqan::Index<seqan::StringSet<seqan::Peptide>, seqan::FMIndex<> > &prot_Index,
                                             seqan::StringSet<seqan::Peptide> &pep_DB,
                                             int mm,
                                             Size max_aaa,
                                             OpenMS::Size indexType){
    typedef seqan::Index<seqan::StringSet<seqan::Peptide>, seqan::FMIndex<> > TIndex;
    typedef typename seqan::Iterator<seqan::StringSet<seqan::Peptide> const, seqan::Rooted>::Type TPatternsIt;
    typedef typename seqan::Iterator<TIndex, seqan::TopDown<seqan::PreorderEmptyEdges> >::Type TIndexIt;
    auto delegate = [&func_SA, &prot_Index, max_aaa, indexType](TIndexIt const &it, TPatternsIt const &patternsIt, unsigned /*score*/) {
        long indexOfPeptideHit = position(patternsIt);
        for (auto occ: getOccurrences(it)){
            long entryNr = seqan::getSeqNo(occ);
            long beg = seqan::getSeqOffset(occ);
            const OpenMS::String tmp_pep(begin(*patternsIt),end(*patternsIt));
            const OpenMS::String tmp_prot(begin(seqan::indexText(prot_Index)[entryNr]),
                                          end(seqan::indexText(prot_Index)[entryNr]));
            if(func_SA.checkAmbigous_(tmp_pep, tmp_prot, beg, max_aaa)){
                func_SA.addHit(indexOfPeptideHit,entryNr,tmp_pep,tmp_prot,beg, indexType);
            }
        }
    };
    find(prot_Index, pep_DB, mm, delegate, seqan::Backtracking<seqan::EditDistance>(), seqan::Parallel());
}


PeptideIndexing2::ExitCodes PeptideIndexing2::mappingPepToProt_(std::vector<FASTAFile::FASTAEntry> &proteins,
                                                                std::vector<FASTAFile::FASTAEntry> &proteinsAAA,
                                                                std::vector<ProteinIdentification> &prot_ids,
                                                                std::vector<PeptideIdentification> &pep_ids,
                                                                Map<String, bool> &protein_is_decoy,
                                                                Map<Size, std::set<Size> > &runidx_to_protidx,
                                                                Map<Size, std::set<Size> > &runidx_to_protidxAAA,
                                                                Size &stats_unmatched,
                                                                seqan::FoundProteinFunctor &func){
    /* do mapping */
    writeDebug_("Reindexing peptide/protein matches...", 1);

    /// index existing proteins
    Map<String, Size> runid_to_runidx; // identifier to index
    // mapping von identifier des <ProteinIdentification> zu Position in idXML file
    // hat die idXML 2 <ProteinIdentification>...</ProteinIdentification> so gibt es 2 Einträge in runid_to_runidx
    for (Size run_idx = 0; run_idx < prot_ids.size(); ++run_idx)
    {
        //std::cout << prot_ids[run_idx].getIdentifier() << std::endl;
        runid_to_runidx[prot_ids[run_idx].getIdentifier()] = run_idx;
    }

    /// for peptides --> proteins
    Size stats_matched_unique(0);
    Size stats_matched_multi(0);
    //Size stats_unmatched(0);
    Size stats_count_m_t(0);
    Size stats_count_m_d(0);
    Size stats_count_m_td(0);
    //Map<Size, set<Size> > runidx_to_protidx; // in which protID do appear which proteins (according to mapped peptides)

    Size pep_idx(0);
    //unsigned peptidIdentificationCounter = 0;
    // iteriere durch alle <PeptidIdentification>...</PeptidIdentification> der idXML
    for (vector<PeptideIdentification>::iterator it1 = pep_ids.begin(); it1 != pep_ids.end(); ++it1)
    {
        //std::cout << "PeptidIdentification Nummer: " << peptidIdentificationCounter << std::endl;
        //peptidIdentificationCounter ++;
        // which ProteinIdentification does the peptide belong to?
        Size run_idx = runid_to_runidx[it1->getIdentifier()];
        // findet heraus ob der Identifier des <PeptidIdentification>...</PeptidIdentification> zum
        // 1. <ProteinIdentification>...</ProteinIdentification> in der idXML gehört oder zu, 2. 3. etc. (sofern vorhanden)
        // anhand des identifiers


        vector<PeptideHit>& hits = it1->getHits();

        // Iteriere durch alle <PeptideHit> ...</PeptideHit> der <PeptidIdentification>... </PeptidIdentification>
        for (vector<PeptideHit>::iterator it2 = hits.begin(); it2 != hits.end(); ++it2)
        {

            //std::cout << "PeptideHit Nummer: " << pep_idx << std::endl;
            // clear protein accessions
            it2->setPeptideEvidences(vector<PeptideEvidence>());

            // Map<OpenMS::Size, std::set<PeptideProteinMatchInformation> > ist der typ des pep_to_prot
            // in pep_to_prot sind nach <PeptidHit> NUMMER!!! sortiert alle hits die gefunden sind gespeichert.
            // gibt es also in der Pep_ids aus der idXML 10 <PeptideIdentification> mit jeweils 2 Hits so hatt
            // pep_to_prot den maximalen Index der Größe 20 (genau dann wenn alle Hits gefunden werden)

            // pep_idx wird bei jedem Hit der in allen PeptidIdentification existeirt hochgezählt...
            // zählt also die <PeptideHit> in allen <PeptidIdentification>...
            // um auf das obrige Beispiel zurückzukommen also genau 20!

            // add new protein references
            for (set<PeptideProteinMatchInformation>::const_iterator it_i = func.pep_to_prot[pep_idx].begin();
                 it_i != func.pep_to_prot[pep_idx].end(); ++it_i) {

                // check in wich index hit was found
                if (!it_i->indexType) {
                    // für alle Treffer zu dem Hit hinter der Nummer pep_idx wurd nun der Identifier eingetragen.
                    // wenn also der Hit mit der nummer pep_idx aus der idXML in 5 Einträgen aus der DB gefunden wird dann werden
                    // hier 5 Einträge gemacht.
                    // in proteins stehen die Proteine in denen gesucht werden geordnet nach der ursprünglichen FASTA file
                    const String &accession = proteins[it_i->protein_index].identifier;
                    //Index des Proteins in der DB (index): protein_index
                    PeptideEvidence pe;
                    pe.setProteinAccession(accession);
                    pe.setStart(it_i->position);
                    pe.setEnd(it_i->position + it2->getSequence().size() - 1);
                    pe.setAABefore(it_i->AABefore);
                    pe.setAAAfter(it_i->AAAfter);
                    it2->addPeptideEvidence(pe);

                    // runidx_to_protidx ist am anfang leer und vom typ Map<Size, set<Size> >
                    // run_idx ist der counter der bestimmt in welchem <ProteinIdentification> das Peptid zu finden ist
                    // (meistens 0)
                    // hier wird also die Position festgehalten bei welchem Protein in der Datenbank der Hit gefunden wurde
                    // Spaeter koennen diese nun geupdated werden in der idXML
                    runidx_to_protidx[run_idx].insert(it_i->protein_index); // fill protein hits
                    //std::cout << "Index des Proteins in der Datenbank: "<< it_i->protein_index << std::endl;
                    // der rest ist statistik!
                    if (!protein_is_decoy.has(accession)) {
                        protein_is_decoy[accession] = (prefix_ && accession.hasPrefix(decoy_string_)) ||
                                                      (!prefix_ && accession.hasSuffix(decoy_string_));
                    }
                }else{
                    // hier der fall der hit kommt aus AAA db müsste in runidx_to_protidxAAA kommen
                    const String &accession = proteinsAAA[it_i->protein_index].identifier;
                    //Index des Proteins in der DB (index): protein_index
                    PeptideEvidence pe;
                    pe.setProteinAccession(accession);
                    pe.setStart(it_i->position);
                    pe.setEnd(it_i->position + it2->getSequence().size() - 1);
                    pe.setAABefore(it_i->AABefore);
                    pe.setAAAfter(it_i->AAAfter);
                    it2->addPeptideEvidence(pe);

                    // runidx_to_protidx ist am anfang leer und vom typ Map<Size, set<Size> >
                    // run_idx ist der counter der bestimmt in welchem <ProteinIdentification> das Peptid zu finden ist
                    // (meistens 0)
                    // hier wird also die Position festgehalten bei welchem Protein in der Datenbank der Hit gefunden wurde
                    // Spaeter koennen diese nun geupdated werden in der idXML
                    runidx_to_protidxAAA[run_idx].insert(it_i->protein_index); // fill protein hits
                    //std::cout << "Index des Proteins in der Datenbank: "<< it_i->protein_index << std::endl;
                    // der rest ist statistik!
                    if (!protein_is_decoy.has(accession)) {
                        protein_is_decoy[accession] = (prefix_ && accession.hasPrefix(decoy_string_)) ||
                                                      (!prefix_ && accession.hasSuffix(decoy_string_));
                    }
                }
            }

            ///
            /// is this a decoy hit?
            ///
            bool matches_target(false);
            bool matches_decoy(false);

            set<String> protein_accessions = it2->extractProteinAccessions();
            for (set<String>::const_iterator it = protein_accessions.begin(); it != protein_accessions.end(); ++it)
            {
                if (protein_is_decoy[*it])
                {
                    matches_decoy = true;
                }
                else
                {
                    matches_target = true;
                }
                // this is rare in practice, so the test may not really save time:
                // if (matches_decoy && matches_target)
                // {
                //   break; // no need to check remaining accessions
                // }
            }
            String target_decoy;
            if (matches_decoy && matches_target)
            {
                target_decoy = "target+decoy";
                ++stats_count_m_td;
            }
            else if (matches_target)
            {
                target_decoy = "target";
                ++stats_count_m_t;
            }
            else if (matches_decoy)
            {
                target_decoy = "decoy";
                ++stats_count_m_d;
            }
            it2->setMetaValue("target_decoy", target_decoy);

            if (protein_accessions.size() == 1)
            {
                it2->setMetaValue("protein_references", "unique");
                ++stats_matched_unique;
            }
            else if (protein_accessions.size() > 1)
            {
                it2->setMetaValue("protein_references", "non-unique");
                ++stats_matched_multi;
            }
            else
            {
                it2->setMetaValue("protein_references", "unmatched");
                ++stats_unmatched;
                if (stats_unmatched < 15) LOG_INFO << "Unmatched peptide: " << it2->getSequence() << "\n";
                else if (stats_unmatched == 15) LOG_INFO << "Unmatched peptide: ...\n";
            }

            ++pep_idx; // next hit
        }
    }
    LOG_INFO << "-----------------------------------\n";
    LOG_INFO << "Peptides statistics\n";
    LOG_INFO << "\n";
    LOG_INFO << "  target/decoy:\n";
    LOG_INFO << "    match to target DB only: " << stats_count_m_t << "\n";
    LOG_INFO << "    match to decoy DB only : " << stats_count_m_d << "\n";
    LOG_INFO << "    match to both          : " << stats_count_m_td << "\n";
    LOG_INFO << "\n";
    LOG_INFO << "  mapping to proteins:\n";
    LOG_INFO << "    no match (to 0 protein)         : " << stats_unmatched << "\n";
    LOG_INFO << "    unique match (to 1 protein)     : " << stats_matched_unique << "\n";
    LOG_INFO << "    non-unique match (to >1 protein): " << stats_matched_multi << std::endl;


    /// exit if no peptides were matched to decoy
    if ((stats_count_m_d + stats_count_m_td) == 0)
    {
        String msg("No peptides were matched to the decoy portion of the database! Did you provide the correct concatenated database? Are your 'decoy_string' (=" + String(decoy_string_) + ") and 'decoy_string_position' (=" + String(param_.getValue("decoy_string_position")) + ") settings correct?");
        if (missing_decoy_action_== "error")
        {
            LOG_ERROR << "Error: " << msg << "\nSet 'missing_decoy_action' to 'warn' if you are sure this is ok!\nAborting ..." << std::endl;
            return UNEXPECTED_RESULT;
        }
        else
        {
            LOG_WARN << "Warn: " << msg << "\nSet 'missing_decoy_action' to 'error' if you want to elevate this to an error!" << std::endl;
        }
    }
    return CHECKPOINT_OK;
}

PeptideIndexing2::ExitCodes PeptideIndexing2::updateProtHit_(std::vector<FASTAFile::FASTAEntry> &proteins,
                                                             std::vector<FASTAFile::FASTAEntry> &proteinsAAA,
                                                             std::vector<ProteinIdentification> &prot_ids,
                                                             Map<String, Size> &acc_to_prot,
                                                             Map<String, Size> &acc_to_AAAprot,
                                                             Map<String, bool> &protein_is_decoy,
                                                             Map<Size, std::set<Size> > &runidx_to_protidx,
                                                             Map<Size, std::set<Size> > &runidx_to_protidxAAA,
                                                             Size &stats_unmatched){
    Int stats_new_proteins(0);
    Int stats_orphaned_proteins(0);

    // all peptides contain the correct protein hit references, now update the protein hits
    for (Size run_idx = 0; run_idx < prot_ids.size(); ++run_idx)
    {
        // bei dieser for-schleife wird durch die <ProteinIdentification></ProteinIdentification> der idXML iteriert

        // es wird die Information in masterset gespeichert wo festgehalten ist zu welchem
        // Protein in der DB Treffer gefunden wurden.
        set<Size> masterset = runidx_to_protidx[run_idx]; // all found protein matches
        set<Size> mastersetAAA = runidx_to_protidxAAA[run_idx]; // all found protein AAA matches

        vector<ProteinHit> new_protein_hits;
        // go through existing hits and update (do not create from anew, as there might be other information (score, rank, etc.) which
        // we want to preserve
        for (vector<ProteinHit>::iterator p_hit = prot_ids[run_idx].getHits().begin(); p_hit != prot_ids[run_idx].getHits().end(); ++p_hit)
        {
            // for schleife durch alle <ProteinHit></ProteinHit> der idXML
            const String& acc = p_hit->getAccession();
            if (acc_to_prot.has(acc) // accession needs to exist in new FASTA file
                && masterset.find(acc_to_prot[acc]) != masterset.end())
            { // this accession was there already
                String seq;
                if (write_protein_sequence_)
                {
                    seq = proteins[acc_to_prot[acc]].sequence;
                }
                p_hit->setSequence(seq);

                if (write_protein_description_)
                {
                    const String& description = proteins[acc_to_prot[acc]].description;
                    //std::cout << "Description = " << description << "\n";
                    p_hit->setDescription(description);
                }

                new_protein_hits.push_back(*p_hit);
                masterset.erase(acc_to_prot[acc]); // remove from master (at the end only new proteins remain)
            } else if (acc_to_AAAprot.has(acc) && mastersetAAA.find(acc_to_AAAprot[acc]) != mastersetAAA.end()){
                String seq;
                if (write_protein_sequence_)
                {
                    seq = proteinsAAA[acc_to_AAAprot[acc]].sequence;
                }
                p_hit->setSequence(seq);

                if (write_protein_description_)
                {
                    const String& description = proteinsAAA[acc_to_AAAprot[acc]].description;
                    //std::cout << "Description = " << description << "\n";
                    p_hit->setDescription(description);
                }

                new_protein_hits.push_back(*p_hit);
                mastersetAAA.erase(acc_to_AAAprot[acc]); // remove from master (at the end only new proteins remain)
            }

            else // old hit is orphaned
            {
                ++stats_orphaned_proteins;
                if (keep_unreferenced_proteins_) new_protein_hits.push_back(*p_hit);
            }
        }

        // add remaining new hits
        for (set<Size>::const_iterator it = masterset.begin();
             it != masterset.end(); ++it)
        {
            // schleife durch alle Indexe von Proteinen der DB die noch  nicht vorhanden waren in der idXML
            // da ich hier eventuell in einer größeren DB mit mehr Proteinen suche!

            // in *it ist also der index des Proteins aus der db wo ein treffer gefunden wurde.
            ProteinHit hit;
            hit.setAccession(proteins[*it].identifier);
            if (write_protein_sequence_)
            {
                //std::cout << "sequence = " << proteins[*it].sequence << "\n";
                hit.setSequence(proteins[*it].sequence);
            }

            if (write_protein_description_)
            {
                //std::cout << "Description = " << proteins[*it].description << "\n";
                hit.setDescription(proteins[*it].description);
            }

            new_protein_hits.push_back(hit);
            ++stats_new_proteins;
        }
        // add remaining AAA new hits schaue aber hier in proteinsAAA nach!
        for (set<Size>::const_iterator it = mastersetAAA.begin();
             it != mastersetAAA.end(); ++it)
        {
            // schleife durch alle Indexe von Proteinen der DB die noch  nicht vorhanden waren in der idXML
            // da ich hier eventuell in einer größeren DB mit mehr Proteinen suche!

            // in *it ist also der index des Proteins aus der db wo ein treffer gefunden wurde.
            ProteinHit hit;
            hit.setAccession(proteinsAAA[*it].identifier);
            if (write_protein_sequence_)
            {
                //std::cout << "sequence = " << proteins[*it].sequence << "\n";
                hit.setSequence(proteinsAAA[*it].sequence);
            }

            if (write_protein_description_)
            {
                //std::cout << "Description = " << proteins[*it].description << "\n";
                hit.setDescription(proteinsAAA[*it].description);
            }

            new_protein_hits.push_back(hit);
            ++stats_new_proteins;
        }

        prot_ids[run_idx].setHits(new_protein_hits);
    }
    // annotate target/decoy status of proteins:
    for (vector<ProteinIdentification>::iterator id_it = prot_ids.begin(); id_it != prot_ids.end(); ++id_it)
    {
        for (vector<ProteinHit>::iterator hit_it = id_it->getHits().begin(); hit_it != id_it->getHits().end(); ++hit_it)
        {
            hit_it->setMetaValue("target_decoy", (protein_is_decoy[hit_it->getAccession()] ? "decoy" : "target"));
        }
    }

    LOG_INFO << "-----------------------------------\n";
    LOG_INFO << "Protein statistics\n";
    LOG_INFO << "\n";
    LOG_INFO << "  new proteins: " << stats_new_proteins << "\n";
    LOG_INFO << "  orphaned proteins: " << stats_orphaned_proteins << (keep_unreferenced_proteins_ ? " (all kept)" : " (all removed)") << "\n";

    writeDebug_("Reindexing finished!", 1);

    if ((!allow_unmatched_) && (stats_unmatched > 0))
    {
        LOG_WARN << "PeptideIndexer found unmatched peptides, which could not be associated to a protein.\n"
        << "Potential solutions:\n"
        << "   - check your FASTA database for completeness\n"
        << "   - set 'enzyme:specificity' to match the identification parameters of the search engine\n"
        << "   - some engines (e.g. X! Tandem) employ loose cutting rules generating non-tryptic peptides;\n"
        << "     if you trust them, disable enzyme specificity\n"
        << "   - increase 'aaa_max' to allow more ambiguous amino acids\n"
        << "   - as a last resort: use the 'allow_unmatched' option to accept unmatched peptides\n"
        << "     (note that unmatched peptides cannot be used for FDR calculation or quantification)\n";

        LOG_WARN << "Result files were written, but PeptideIndexer will exit with an error code." << std::endl;
        return UNEXPECTED_RESULT;
    }
    return CHECKPOINT_OK;
}

PeptideIndexing2::ExitCodes PeptideIndexing2::processMap_(EnzymaticDigestion enzyme,
                                                          String path,
                                                          String pathAAA,
                                                          std::vector<ProteinIdentification> &prot_ids,
                                                          std::vector<PeptideIdentification> &pep_ids,
                                                          FMind /**/) {
    seqan::FoundProteinFunctor func(enzyme); // stores the matches
    seqan::Index<seqan::StringSet<seqan::Peptide>, seqan::FMIndex<> > index;
    if (search_for_normal_proteins_) {
        if (!seqan::open(index, path.c_str())) {
            return INPUT_ERROR;
        }
    }
    seqan::Index<seqan::StringSet<seqan::Peptide>, seqan::FMIndex<> > indexAAA;
    if (search_for_aaa_proteins_) {
        if (!seqan::open(indexAAA, pathAAA.c_str())) {
            return INPUT_ERROR;
        }
    }

    std::vector<FASTAFile::FASTAEntry> proteins;
    std::vector<FASTAFile::FASTAEntry> proteinsAAA;
    Map<String, Size> acc_to_prot;
    Map<String, Size> acc_to_AAAprot;
    if (loadInfo_(proteins, proteinsAAA, acc_to_prot, acc_to_AAAprot, path , pathAAA) != CHECKPOINT_OK){
        return INPUT_ERROR;
    };

    /* BUILD PEPTIDE DB */
    writeLog_(String("Build Peptide DB"));
    seqan::StringSet<seqan::Peptide> pep_DB;
    if (buildReversePepDB_(pep_DB, pep_ids)!= CHECKPOINT_OK){
        return UNEXPECTED_RESULT;
    }

    /* SEARCH WITH FMINDEX */
    writeLog_(String("Searching ..."));
    // check options
    if (search_for_normal_proteins_ && search_for_aaa_proteins_) {
        // search normal
        searchWrapper_(func, index, pep_DB, mismatches_max_, aaa_max_,0);
        std::cout << "suche 1 " << func.pep_to_prot.size() << std::endl;
        // search AAA
        searchWrapper_(func, indexAAA, pep_DB, mismatches_max_, aaa_max_,1);
        std::cout << "suche 2 " << func.pep_to_prot.size() << std::endl;
    } else if (search_for_normal_proteins_ && !search_for_aaa_proteins_ ){
        // search normal
        searchWrapper_(func, index, pep_DB, mismatches_max_, aaa_max_,0);
    } else {
        // search AAA
        searchWrapper_(func, indexAAA, pep_DB, mismatches_max_, aaa_max_,1);
    }

    // write some stats
    LOG_INFO << "Peptide hits passing enzyme filter: " << func.filter_passed << "\n"
             << "     ... rejected by enzyme filter: " << func.filter_rejected << std::endl;

    /* MAPPING */
    writeLog_(String("Mapping"));
    Map<Size, set<Size> > runidx_to_protidx;
    Map<Size, set<Size> > runidx_to_protidxAAA;
    /// store target/decoy status of proteins
    Map<String, bool> protein_is_decoy; // accession -> is decoy?
    Size stats_unmatched(0);
    if (mappingPepToProt_(proteins, proteinsAAA, prot_ids, pep_ids, protein_is_decoy, runidx_to_protidx, runidx_to_protidxAAA, stats_unmatched, func)!= CHECKPOINT_OK) {
        return UNEXPECTED_RESULT;
    }

    /* UPDATE HITS */
    writeLog_(String("Updating"));
    if (updateProtHit_(proteins, proteinsAAA, prot_ids, acc_to_prot,acc_to_AAAprot, protein_is_decoy, runidx_to_protidx, runidx_to_protidxAAA, stats_unmatched)!= CHECKPOINT_OK){
        return UNEXPECTED_RESULT;
    }
    return CHECKPOINT_OK;
}


PeptideIndexing2::ExitCodes PeptideIndexing2::processMap_(EnzymaticDigestion enzyme,
                                                          String path,
                                                          String pathAAA,
                                                          std::vector<ProteinIdentification> &prot_ids,
                                                          std::vector<PeptideIdentification> &pep_ids,
                                                          SAind /**/) {
    writeLog_(String("Read Input"));
    seqan::FoundProteinFunctor func(enzyme); // stores the matches
    seqan::Index<seqan::StringSet<seqan::Peptide>, seqan::IndexSa<> > index;
    if (search_for_normal_proteins_) {
        if (!seqan::open(index, path.c_str())) {
            writeLog_(String("ERROR: Could not open Index!"));
            return INPUT_ERROR;
        }
    }
    seqan::Index<seqan::StringSet<seqan::Peptide>, seqan::IndexSa<> > indexAAA;
    if (search_for_aaa_proteins_) {
        if (!seqan::open(indexAAA, pathAAA.c_str())) {
            writeLog_(String("ERROR: Could not open Index for ambiguous amino acid!"));
            return INPUT_ERROR;
        }
    }
    std::vector<FASTAFile::FASTAEntry> proteins;
    std::vector<FASTAFile::FASTAEntry> proteinsAAA;
    Map<String, Size> acc_to_prot;
    Map<String, Size> acc_to_AAAprot;
    if (loadInfo_(proteins, proteinsAAA, acc_to_prot, acc_to_AAAprot, path, pathAAA) != CHECKPOINT_OK){
        return INPUT_ERROR;
    };
    writeLog_(String("Build Peptide DB"));
    /* BUILD PEPTIDE DB */
    seqan::StringSet<seqan::Peptide> pep_DB;
    if (buildPepDB_(pep_DB, pep_ids)!= CHECKPOINT_OK){
        return UNEXPECTED_RESULT;
    }

    /* SEARCH WITH SUFFIX ARRAY */

    writeLog_(String("Searching ..."));
    // check options
    if (search_for_normal_proteins_ && search_for_aaa_proteins_) {
        // search normal
        searchWrapper_(func, index, pep_DB, mismatches_max_, aaa_max_,0);
        std::cout << "suche 1 " << func.pep_to_prot.size() << std::endl;
        // search AAA
        searchWrapper_(func, indexAAA, pep_DB, mismatches_max_, aaa_max_,1);
        std::cout << "suche 2 " << func.pep_to_prot.size() << std::endl;
    } else if (search_for_normal_proteins_ && !search_for_aaa_proteins_ ){
        // search normal
        searchWrapper_(func, index, pep_DB, mismatches_max_, aaa_max_,0);
    } else {
        // search AAA
        searchWrapper_(func, indexAAA, pep_DB, mismatches_max_, aaa_max_,1);
    }

    // write some stats

    LOG_INFO << "Peptide hits passing enzyme filter: " << func.filter_passed << "\n"
             << "     ... rejected by enzyme filter: " << func.filter_rejected << std::endl;

    /* MAPPING */
    writeLog_(String("Mapping" ));
    Map<Size, set<Size> > runidx_to_protidx;
    Map<Size, set<Size> > runidx_to_protidxAAA;
    /// store target/decoy status of proteins
    Map<String, bool> protein_is_decoy; // accession -> is decoy?
    Size stats_unmatched(0);
    if (mappingPepToProt_(proteins, proteinsAAA, prot_ids, pep_ids, protein_is_decoy, runidx_to_protidx, runidx_to_protidxAAA, stats_unmatched, func)!= CHECKPOINT_OK) {
        return UNEXPECTED_RESULT;
    }

    /* UPDATE HITS */
    writeLog_(String("Updating"));
    if (updateProtHit_(proteins, proteinsAAA, prot_ids, acc_to_prot, acc_to_AAAprot, protein_is_decoy, runidx_to_protidx, runidx_to_protidxAAA,stats_unmatched)!= CHECKPOINT_OK){
        return UNEXPECTED_RESULT;
    }
    return CHECKPOINT_OK;
}

PeptideIndexing2::ExitCodes PeptideIndexing2::run(String &index,
                                                  String &indexAAA,
                                std::vector<ProteinIdentification> &prot_ids,
                                std::vector<PeptideIdentification> &pep_ids) {
    // proteins = UniprotKB -> FASTA  -> prot_DB aka IndexStruktur
    // prot_ids = IdXML -> User INput, zunächst leer, zu jdeem Peptide wird nun gemappt welche PROTEIN ID aus "DB" passt
    // pep_ids = IdXML -> User INput, gefüllt mit kurzen Peptiden -> haben id, sequence, etc. -> pep_DB
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    EnzymaticDigestion enzyme;
    enzyme.setEnzyme(enzyme_name_);
    enzyme.setSpecificity(enzyme.getSpecificityByName(enzyme_specificity_));

    if (!log_file_.empty()) {
        log_.open(log_file_.c_str());
    }


    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    search_for_aaa_proteins_ = true;
    search_for_normal_proteins_ = true;
    // Error checking
    if (indexAAA.empty()){
        search_for_aaa_proteins_ = false;
    }
    if (index.empty()){
        search_for_normal_proteins_ = false;
    }

    ExitCodes exCheck = checkUserInput_(prot_ids, pep_ids);
    if (exCheck != CHECKPOINT_OK) {
        return exCheck;
    };


//    writeDebug_("Collecting peptides...", 1);
//    std::cout << "Process start with option ";
    if (suffix_array_) {
        ExitCodes exCheck = processMap_(enzyme, index, indexAAA, prot_ids, pep_ids, SAind());
        if (exCheck != CHECKPOINT_OK) {
            return exCheck;
        };
    } else{
        ExitCodes exCheck = processMap_(enzyme, index, indexAAA, prot_ids, pep_ids, FMind());
        if (exCheck != CHECKPOINT_OK) {
            return exCheck;
        };
    }
//    std::cout << "Process end" << std::endl;
    if (!log_file_.empty()) {
        log_.close();
    }
    sleep(1);
    return EXECUTION_OK;
}
/// @endcond

