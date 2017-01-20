#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/FORMAT/FASTAFile.h>
//#include <OpenMS/DATASTRUCTURES/SeqanIncludeWrapper.h>
#include <seqan2/index.h>
#include <seqan2/find.h>
#include <seqan2/basic.h>
#include <seqan2/sequence.h>
#include <seqan2/stream.h>

#include <fstream>


struct FMind;
struct SAind;

namespace OpenMS {

    class OPENMS_DLLAPI BuildIndexing :
            public DefaultParamHandler {
    public:

        /// Exit codes
        enum ExitCodes {
            EXECUTION_OK,
            DATABASE_EMPTY,
            DATABASE_CONTAINS_MULTIPLES,
            ILLEGAL_PARAMETERS,
            UNEXPECTED_RESULT,
            OUTPUT_ERROR,
            CHECKPOINT_OK,
            CHECKPOINT_DONE
        };

        /// Default constructor
        BuildIndexing();

        /// Default destructor
        virtual ~BuildIndexing();

        /// main method of PeptideIndexerBuilder
        ExitCodes run(const String &in, const String &out);

    protected:

        virtual void updateMembers_();

        void writeLog_(const String &text) const;

        void writeDebug_(const String &text, const Size min_level) const;

        /// Output stream for log/debug info
        String log_file_;
        mutable std::ofstream log_;
        /// debug flag
        bool debug_;
        /// program options
        bool suffix_array_;
        bool FM_index_;
        bool IL_equivalent_;
        bool has_aaa_(const String seq);

        /// function to build protein database
        ExitCodes buildProtDB_(const String &in,
                std::vector<FASTAFile::FASTAEntry>& proteins,
                               std::vector<FASTAFile::FASTAEntry>& proteinsAAA,
                               seqan2::StringSet<seqan2::CharString>& ids,
                                seqan2::StringSet<seqan2::Peptide>& seqs,
                               Map<String, Size> &acc_to_prot,
                               seqan2::StringSet<seqan2::Peptide> &prot_DB,
                               Map<String, Size> &acc_to_AAAprot,
                               seqan2::StringSet<seqan2::Peptide> &prot_DB_AAA,
                               std::vector<String> &duplicate_accessions);

        /// function to check for duplicate fasta entrys
        bool check_duplicate_(const String seq,
                           std::vector<String> &duplicate_accessions,
                           const Map<String, Size> &acc_to_prot,
                           const String &acc,
                           const seqan2::StringSet<seqan2::Peptide> &prot_DB);

        /// function to safe additional information on disk
        ExitCodes saveOutput_(const Map<String, Size> &acc_to_prot,
                              const Map<String, Size> &acc_to_AAAprot,
                              const std::vector<FASTAFile::FASTAEntry>& proteins,
                              const std::vector<FASTAFile::FASTAEntry>& proteinsAAA,
                              const String &out);

        /// function to build suffix array index
        ExitCodes build_index_(const seqan2::StringSet<seqan2::Peptide> &prot_DB,
                               const seqan2::StringSet<seqan2::Peptide> &prot_DB_AAA,
                               const String out,
                               SAind /**/);

        /// function to build FM index
        ExitCodes build_index_(seqan2::StringSet<seqan2::Peptide> &prot_DB,
                               seqan2::StringSet<seqan2::Peptide> &prot_DB_AAA,
                               const String out,
                               FMind /**/);

        /// function to check user input
        ExitCodes checkUserInput_();

    };
}
