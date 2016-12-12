#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/DATASTRUCTURES/SeqanIncludeWrapper.h>


#include <fstream>


struct FMind;
struct SAind;
struct WOTDind;

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

                /// main method of PeptideIndexing
        ExitCodes run(std::vector<FASTAFile::FASTAEntry> &proteins, String &out);

    protected:

        virtual void updateMembers_();

        void writeLog_(const String &text) const;

        void writeDebug_(const String &text, const Size min_level) const;

        ExitCodes buildProtDB_(std::vector<FASTAFile::FASTAEntry>& proteins,
                               Map<String, Size> &acc_to_prot,
                               seqan::StringSet<seqan::Peptide> &prot_DB,
                               Map<String, Size> &acc_to_AAAprot,
                               seqan::StringSet<seqan::Peptide> &prot_DB_AAA,
                               std::vector<String> &duplicate_accessions);

        ExitCodes check_duplicate_(std::vector<FASTAFile::FASTAEntry>& proteins,
                                   String seq,
                                   std::vector<String> &duplicate_accessions,
                                   Map<String, Size> &acc_to_prot,
                                   String &acc,
                                   seqan::StringSet<seqan::Peptide>
                                   &prot_DB, Size protIndex);

        /// Output stream for log/debug info
        String log_file_;
        mutable std::ofstream log_;
        /// debug flag
        bool debug_;
        bool WOTD_;
        bool suffix_array_;
        bool FM_index_;
        bool IL_equivalent_;

        bool has_aaa_(String seq);

        ExitCodes saveOutput_(Map<String, Size> &acc_to_prot,
                              Map<String, Size> &acc_to_AAAprot,
                              std::vector<FASTAFile::FASTAEntry>& proteins,
                              String &out);

        ExitCodes build_index_(seqan::StringSet<seqan::Peptide> &prot_DB,
                               seqan::StringSet<seqan::Peptide> &prot_DB_AAA,
                               String out,
                               SAind /**/);

        ExitCodes build_index_(seqan::StringSet<seqan::Peptide> &prot_DB,
                               seqan::StringSet<seqan::Peptide> &prot_DB_AAA,
                               String out,
                               FMind /**/);

        ExitCodes checkUserInput_(std::vector<FASTAFile::FASTAEntry> &proteins);
    };
}
