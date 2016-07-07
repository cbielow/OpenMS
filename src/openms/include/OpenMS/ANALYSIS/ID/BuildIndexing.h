#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/FORMAT/FASTAFile.h>


#include <fstream>

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
            OUTPUT_ERROR
        };

        /// Default constructor
        BuildIndexing();

        /// Default destructor
        virtual ~BuildIndexing();

        /// main method of PeptideIndexing
        ExitCodes run(std::vector<FASTAFile::FASTAEntry> &proteins, String &out);

    protected:

        void writeLog_(const String &text) const;

        void writeDebug_(const String &text, const Size min_level) const;

        /// Output stream for log/debug info
        String log_file_;
        mutable std::ofstream log_;
        /// debug flag
        bool debug_;

        String decoy_string_;
        bool prefix_;
        String missing_decoy_action_;
        String enzyme_name_;
        String enzyme_specificity_;

        bool write_protein_sequence_;
        bool write_protein_description_;
        bool keep_unreferenced_proteins_;
        bool allow_unmatched_;
        bool full_tolerant_search_;
        bool IL_equivalent_;

        Size aaa_max_;
        UInt mismatches_max_;
        bool filter_aaa_proteins_;

    };
}
