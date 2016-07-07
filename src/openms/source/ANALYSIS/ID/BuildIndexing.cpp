#include <OpenMS/ANALYSIS/ID/BuildIndexing.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/DATASTRUCTURES/SeqanIncludeWrapper.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/SYSTEM/StopWatch.h>
#include <OpenMS/METADATA/PeptideEvidence.h>
#include <OpenMS/CHEMISTRY/EnzymesDB.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <algorithm>

using namespace OpenMS;
using namespace std;

BuildIndexing::BuildIndexing() :
        DefaultParamHandler("BuildIndexing") {
    defaults_.setValue("SuffixArray", "true", "Specify the Format of the Output File.");
    defaults_.setValidStrings("SuffixArray", ListUtils::create<String>("true,false"));

    defaults_.setValue("FMIndex", "false", "Specify the Format of the Output File.");
    defaults_.setValidStrings("SuffixArray", ListUtils::create<String>("true,false"));

//    defaults_.setValue("decoy_string", "_rev", "String that was appended (or prefixed - see 'prefix' flag below) to the accessions in the protein database to indicate decoy proteins.");
//
//    defaults_.setValue("prefix", "false", "If set, protein accessions in the database contain 'decoy_string' as prefix.");
//    defaults_.setValidStrings("prefix", ListUtils::create<String>("true,false"));
//
//    defaults_.setValue("missing_decoy_action", "error", "Action to take if NO peptide was assigned to a decoy protein (which indicates wrong database or decoy string): 'error' (exit with error, no output), 'warn' (exit with success, warning message)");
//    defaults_.setValidStrings("missing_decoy_action", ListUtils::create<String>("error,warn"));
//
//    defaults_.setValue("write_protein_sequence", "false", "If set, the protein sequences are stored as well.");
//    defaults_.setValidStrings("write_protein_sequence", ListUtils::create<String>("true,false"));
//
//    defaults_.setValue("write_protein_description", "false", "If set, the protein description is stored as well.");
//    defaults_.setValidStrings("write_protein_description", ListUtils::create<String>("true,false"));
//
//    defaults_.setValue("keep_unreferenced_proteins", "false", "If set, protein hits which are not referenced by any peptide are kept.");
//    defaults_.setValidStrings("keep_unreferenced_proteins", ListUtils::create<String>("true,false"));
//
//    defaults_.setValue("allow_unmatched", "false", "If set, unmatched peptide sequences are allowed. By default (i.e. if this flag is not set) the program terminates with an error on unmatched peptides.");
//    defaults_.setValidStrings("allow_unmatched", ListUtils::create<String>("true,false"));
//
//    defaults_.setValue("full_tolerant_search", "false", "If set, all peptide sequences are matched using tolerant search. Thus potentially more proteins (containing ambiguous amino acids) are associated. This is much slower!");
//    defaults_.setValidStrings("full_tolerant_search", ListUtils::create<String>("true,false"));
//
//    defaults_.setValue("aaa_max", 4, "[tolerant search only] Maximal number of ambiguous amino acids (AAAs) allowed when matching to a protein database with AAAs. AAAs are 'B', 'Z' and 'X'");
//    defaults_.setMinInt("aaa_max", 0);
//
//    defaults_.setValue("mismatches_max", 0, "[tolerant search only] Maximal number of real mismatches (will be used after checking for ambiguous AA's (see 'aaa_max' option). In general this param should only be changed if you want to look for other potential origins of a peptide which might have unknown SNPs or the like.");
//    defaults_.setMinInt("mismatches_max", 0);
//
    defaults_.setValue("IL_equivalent", "false",
                       "Treat the isobaric amino acids isoleucine ('I') and leucine ('L') as equivalent (indistinguishable)");
    defaults_.setValidStrings("IL_equivalent", ListUtils::create<String>("true,false"));
//
//    defaults_.setValue("filter_aaa_proteins", "false", "In the tolerant search for matches to proteins with ambiguous amino acids (AAAs), rebuild the search database to only consider proteins with AAAs. This may save time if most proteins don't contain AAAs and if there is a significant number of peptides that enter the tolerant search.");
//    defaults_.setValidStrings("filter_aaa_proteins", ListUtils::create<String>("true,false"));
//
    defaults_.setValue("log", "", "Name of log file (created only when specified)");
    defaults_.setValue("debug", 0, "Sets the debug level");

    defaultsToParam_();
}

BuildIndexing::~BuildIndexing() {

}

void BuildIndexing::writeLog_(const String &text) const {
    LOG_INFO << text << endl;
    if (!log_file_.empty()) {
        log_ << QDateTime::currentDateTime().toString("yyyy-MM-dd hh:mm:ss").toStdString() << ": " << text << endl;
    }
}

void BuildIndexing::writeDebug_(const String &text, const Size min_level) const {
    if (debug_ >= min_level) {
        writeLog_(text);
    }
}

BuildIndexing::ExitCodes BuildIndexing::run(std::vector<FASTAFile::FASTAEntry>& proteins, String &out) {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    if (!log_file_.empty()) {
        log_.open(log_file_.c_str());
    }


    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    if (proteins.empty()) // we do not allow an empty database
    {
        LOG_ERROR << "Error: An empty database was provided. Mapping makes no sense. Aborting..." << std::endl;
        return DATABASE_EMPTY;
    }
    writeDebug_("Collecting peptides...", 1);

    /**
       BUILD Protein DB
    */
    Map<String, Size> acc_to_prot; // build map: accessions to FASTA protein index
    seqan::StringSet<seqan::Peptide> prot_DB;
    vector<String> duplicate_accessions;

    for (Size i = 0; i != proteins.size(); ++i) {
        String seq = proteins[i].sequence.remove('*');
        if (IL_equivalent_) // convert  L to I; warning: do not use 'J', since Seqan does not know about it and will convert 'J' to 'X'
        {
            seq.substitute('L', 'I');
        }

        String acc = proteins[i].identifier;
        // check for duplicate proteins
        if (acc_to_prot.has(acc)) {
            duplicate_accessions.push_back(acc);
            // check if sequence is identical
            const seqan::Peptide &tmp_prot = prot_DB[acc_to_prot[acc]];
            if (String(begin(tmp_prot), end(tmp_prot)) != seq) {
                LOG_ERROR << "Fatal error: Protein identifier '" << acc <<
                "' found multiple times with different sequences" << (IL_equivalent_ ? " (I/L substituted)" : "") <<
                ":\n"
                << tmp_prot << "\nvs.\n" << seq << "\nPlease fix the database and run PeptideIndexer again." <<
                std::endl;
                return DATABASE_CONTAINS_MULTIPLES;
            }
            // Remove duplicate entry from 'proteins', since 'prot_DB' and 'proteins' need to correspond 1:1 (later indexing depends on it)
            // The other option would be to allow two identical entries, but later on, only the last one will be reported (making the first protein an orphan; implementation details below)
            // Thus, the only safe option is to remove the duplicate from 'proteins' and not to add it to 'prot_DB'
            proteins.erase(proteins.begin() + i);
            --i;  // try this index again
        } else {
            // extend protein DB
            seqan::appendValue(prot_DB, seq.c_str());
            acc_to_prot[acc] = i;
        }

    }


    /**
    BUILD Protein DB Index
    */

    seqan::Index<seqan::StringSet<seqan::Peptide>, seqan::FMIndex<> > index(prot_DB);
    seqan::indexRequire(index, seqan::FibreSaLfTable());


    /**
    SAVE Protein DB Index
    */

    if (!seqan::save(index, out.c_str() )) {
        return OUTPUT_ERROR;
    }


    return EXECUTION_OK;
}