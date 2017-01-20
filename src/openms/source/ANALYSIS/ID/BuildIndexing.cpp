#include <OpenMS/ANALYSIS/ID/BuildIndexing.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
//#include <OpenMS/DATASTRUCTURES/SeqanIncludeWrapper.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/SYSTEM/StopWatch.h>
#include <OpenMS/METADATA/PeptideEvidence.h>
#include <OpenMS/CHEMISTRY/EnzymesDB.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <seqan2/seq_io.h>
#include <seqan2/index.h>
#include <seqan2/seeds.h>
#include <seqan2/arg_parse.h>
#include <seqan2/bam_io.h>
#include <seqan2/find.h>
#include <seqan2/basic.h>
#include <seqan2/stream.h>
#include <fstream>
#include <iostream>
#include <algorithm>

#include <OpenMS/FORMAT/FASTAFile.h>

using namespace OpenMS;
using namespace std;

/* */
struct FMind {
    FMind() { }
};

/* */
struct SAind {
    SAind() { }
};

//namespace seqan2 {
//    template<>
//    struct SAValue<StringSet<Peptide> > {
//    // anzahl strings - laenge strings
//        typedef Pair<unsigned, unsigned> Type;
//    };
//}

BuildIndexing::BuildIndexing() :
        DefaultParamHandler("BuildIndexing") {
    defaults_.setValue("suffix_array", "false", "Using the Suffix Array as index");
    defaults_.setValidStrings("suffix_array", ListUtils::create<String>("true,false"));

    defaults_.setValue("FM_index", "false", "Using the FM-Index as index");
    defaults_.setValidStrings("FM_index", ListUtils::create<String>("true,false"));

    defaults_.setValue("IL_equivalent", "false",
                       "Treat the isobaric amino acids isoleucine ('I') and leucine ('L') as equivalent (indistinguishable)");
    defaults_.setValidStrings("IL_equivalent", ListUtils::create<String>("true,false"));
//
    defaults_.setValue("log", "", "Name of log file (created only when specified)");
    defaults_.setValue("debug", 0, "Sets the debug level");

    defaultsToParam_();
}

BuildIndexing::~BuildIndexing() {

}

BuildIndexing::ExitCodes BuildIndexing::saveOutput_(const Map<String, Size> &acc_to_prot,
                                                    const Map<String, Size> &acc_to_AAAprot,
                                                    const std::vector<FASTAFile::FASTAEntry>& proteins,
                                                    const std::vector<FASTAFile::FASTAEntry>& proteinsAAA,
                                                    const String &out){
    // save them on disk
    // jede accession nummer wird die Position in dem Suffix-Array gespeichert
    // Eintrag: CAQ30518;1383
    // accession Nummern sind alphabetisch geordet
    std::ofstream acc_to_prot_out((out + "_acc_to_prot").c_str());
    for (map<OpenMS::String, unsigned long>::const_iterator i = acc_to_prot.begin(); i != acc_to_prot.end(); i++){
        acc_to_prot_out << (*i).first << ";" << (*i).second << "\n";
    }
    acc_to_prot_out.close();

    std::ofstream acc_to_AAAprot_out((out + "_AAA_acc_to_prot").c_str());
    for (map<OpenMS::String, unsigned long>::const_iterator i = acc_to_AAAprot.begin(); i != acc_to_AAAprot.end(); i++){
        acc_to_AAAprot_out << (*i).first << ";" << (*i).second << "\n";
    }
    acc_to_AAAprot_out.close();

    // speichert quasi die Fasta erneut.
    std::ofstream proteins_out((out + "_proteins").c_str());
    for (vector<OpenMS::FASTAFile::FASTAEntry>::const_iterator i = proteins.begin(); i != proteins.end(); i++){
        proteins_out << (*i).identifier << ";" << (*i).sequence << ";" << (*i).description << "\n";
    }
    proteins_out.close();

    std::ofstream proteinsAAA_out((out + "_AAA_proteins").c_str());
    for (vector<OpenMS::FASTAFile::FASTAEntry>::const_iterator i = proteinsAAA.begin(); i != proteinsAAA.end(); i++){
        proteinsAAA_out << (*i).identifier << ";" << (*i).sequence << ";" << (*i).description << "\n";
    }
    proteinsAAA_out.close();

    return CHECKPOINT_OK;

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

BuildIndexing::ExitCodes BuildIndexing::run(const String &in, const String &out) {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    if (!log_file_.empty()) {
        log_.open(log_file_.c_str());
    }

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    ExitCodes erg = checkUserInput_();
    if (erg != CHECKPOINT_OK){
        return erg;
    }

    /**
       BUILD Protein DB
    */
    writeLog_(String("Building protein database..."));
    seqan2::StringSet<seqan2::CharString> ids;
    seqan2::StringSet<seqan2::Peptide> seqs;
    seqan2::StringSet<seqan2::CharString> quals;
    std::vector<FASTAFile::FASTAEntry> proteinsAAA;
    std::vector<FASTAFile::FASTAEntry> proteins;
    Map<String, Size> acc_to_prot; // build map: accessions to FASTA protein index
    Map<String, Size> acc_to_AAAprot;
    seqan2::StringSet<seqan2::Peptide> prot_DB;
    seqan2::StringSet<seqan2::Peptide> prot_DB_AAA;
    vector<String> duplicate_accessions;
    // check for duplicated accessions & build prot DB
    if (buildProtDB_(in,proteins, proteinsAAA, ids, seqs, acc_to_prot, prot_DB, acc_to_AAAprot, prot_DB_AAA, duplicate_accessions) != CHECKPOINT_OK ) {
        return DATABASE_CONTAINS_MULTIPLES;
    }

    /**
    BUILD Index & Save
    */
    writeLog_(String("Building protein database index..."));
    if (suffix_array_){
        if (build_index_(prot_DB,prot_DB_AAA,out,SAind()) != CHECKPOINT_OK){
            return OUTPUT_ERROR;
        }
    }else {
        if (build_index_(prot_DB,prot_DB_AAA,out,FMind()) != CHECKPOINT_OK){
            return OUTPUT_ERROR;
        }
    }

    writeLog_(String("saving protein database index..."));
    saveOutput_(acc_to_prot,acc_to_AAAprot, proteins, proteinsAAA, out);


    return EXECUTION_OK;
}

BuildIndexing::ExitCodes BuildIndexing::buildProtDB_(const String &in,
                                                    std::vector<FASTAFile::FASTAEntry>& proteins,
                                                     std::vector<FASTAFile::FASTAEntry>& proteinsAAA,
                                                     seqan2::StringSet<seqan2::CharString> & ids,
                                                     seqan2::StringSet<seqan2::Peptide> & seqs,
                                                     Map<String, Size> &acc_to_prot,
                                                     seqan2::StringSet<seqan2::Peptide> &prot_DB,
                                                     Map<String, Size> &acc_to_AAAprot,
                                                     seqan2::StringSet<seqan2::Peptide> &prot_DB_AAA,
                                                     std::vector<String> &duplicate_accessions) {
    seqan2::SeqFileIn seqFileIn;
    if (!seqan2::open(seqFileIn, seqan2::toCString(in)))
    {
        std::cerr << "ERROR: Could not open the file.\n";
        return UNEXPECTED_RESULT;
    }
    seqan2::StringSet<seqan2::CharString> quals;
    seqan2::readRecords(ids, seqs, quals, seqFileIn);
    if (seqan2::length(ids) == 0) // we do not allow an empty database
    {
        LOG_ERROR << "Error: An empty database was provided. Preparing Index makes no sense. Aborting..." << std::endl;
        return DATABASE_EMPTY;
    }

    unsigned counterAcc_to_prot=0;
    unsigned counterAcc_to_protAAA=0;

    String::size_type position = String::npos;
    for (unsigned i = 0; i < seqan2::length(ids); ++i) {
        // new FASTA entry
        FASTAFile::FASTAEntry newEntry;
        // get ID
        String id = seqan2::toCString(ids[i]);
        // handle ID
        id = id.trim();
        position = id.find_first_of(" \v\t");
        if (position == String::npos) {
            newEntry.identifier = id;
            newEntry.description = "";
        } else {
            newEntry.identifier = id.substr(0, position);
            newEntry.description = id.suffix(id.size() - position - 1);
        }
        // get sequence
        if (IL_equivalent_) {
            for (unsigned j = 0;j<seqan2::length(seqs[i]);++j){
                if (seqs[i][j] == 'L') seqs[i][j] = 'I';
            }
        }
        seqan2::CharString castedSeq(seqs[i]);
        newEntry.sequence = seqan2::toCString(castedSeq);
        // add to newEntry to proteins and add sequence to ProteinDB separated normal and AAA
        if (has_aaa_(newEntry.sequence)) { // case AAA Protein
            // check if already added
            if (!check_duplicate_(newEntry.sequence, duplicate_accessions, acc_to_AAAprot, newEntry.identifier, prot_DB_AAA)) {
                // add FastaEntry to proteins
                proteinsAAA.push_back(newEntry);
                // add sequence to DB
                seqan2::appendValue(prot_DB_AAA, seqs[i]);
                // link ID to Position in DB
                acc_to_AAAprot[newEntry.identifier] = counterAcc_to_protAAA;
                ++counterAcc_to_protAAA;
            }
        } else {
            // check if already added
            if (!check_duplicate_(newEntry.sequence, duplicate_accessions, acc_to_prot, newEntry.identifier, prot_DB)) {
                // add FastaEntry to proteins
                proteins.push_back(newEntry);
                // add sequence to DB
                seqan2::appendValue(prot_DB, seqs[i]);
                // link ID to Position in DB
                acc_to_prot[newEntry.identifier] = counterAcc_to_prot;
                ++counterAcc_to_prot;
            }
        }
    }
    return BuildIndexing::CHECKPOINT_OK;
}

BuildIndexing::ExitCodes BuildIndexing::build_index_(const seqan2::StringSet<seqan2::Peptide> &prot_DB,
                                                     const seqan2::StringSet<seqan2::Peptide> &prot_DB_AAA,
                                                     const String out,
                                                     SAind /**/){
    if (!seqan2::empty(prot_DB_AAA)){
        seqan2::Index<seqan2::StringSet<seqan2::Peptide>, seqan2::IndexSa<> > indexAAA(prot_DB_AAA);
        seqan2::indexRequire(indexAAA, seqan2::FibreSA());
        if (!seqan2::save(indexAAA, (out + "_AAA").c_str() )) {
            LOG_ERROR << "Error: Could not save output to disk. Check for correct name, available space and privileges..." << std::endl;
            return OUTPUT_ERROR;
        }
    }
    if (!seqan2::empty(prot_DB)){
        seqan2::Index<seqan2::StringSet<seqan2::Peptide>, seqan2::IndexSa<> > index(prot_DB);
        //seqan2::Iter<seqan2::StringSet<seqan2::Peptide>, seqan2::IndexSa<> > iter(index);
        seqan2::indexRequire(index, seqan2::FibreSA());
        if (!seqan2::save(index, out.c_str() )) {
            LOG_ERROR << "Error: Could not save output to disk. Check for correct name, available space and privileges..." << std::endl;
            return OUTPUT_ERROR;
        }
    }
    return CHECKPOINT_OK;
}

BuildIndexing::ExitCodes BuildIndexing::build_index_(seqan2::StringSet<seqan2::Peptide> &prot_DB,
                                                     seqan2::StringSet<seqan2::Peptide> &prot_DB_AAA,
                                                     const String out,
                                                     FMind /**/){
    //typedef seqan2::FastFMIndexConfig<void,std::uint32_t,2,1> TFastFMConfig;
    if (!seqan2::empty(prot_DB_AAA)){
        //seqan2::Index<seqan2::StringSet<seqan2::Peptide>, seqan2::FMIndex<void, TFastFMConfig> > indexAAA(prot_DB_AAA);
        seqan2::Index<seqan2::StringSet<seqan2::Peptide>, seqan2::FMIndex<> > indexAAA(prot_DB_AAA);
        seqan2::indexRequire(indexAAA, seqan2::FibreSALF());
        if (!seqan2::save(indexAAA, (out + "_AAA").c_str() )) {
            LOG_ERROR << "Error: Could not save output to disk. Check for correct name, available space and privileges..." << std::endl;
            return OUTPUT_ERROR;
        }
    }
    if (!seqan2::empty(prot_DB)){
        //seqan2::Index<seqan2::StringSet<seqan2::Peptide>, seqan2::FMIndex<void, TFastFMConfig> > index(prot_DB);
        seqan2::Index<seqan2::StringSet<seqan2::Peptide>, seqan2::FMIndex<> > index(prot_DB);
        seqan2::indexRequire(index, seqan2::FibreSALF());
        if (!seqan2::save(index, out.c_str() )) {
            LOG_ERROR << "Error: Could not save output to disk. Check for correct name, available space and privileges..." << std::endl;
            return OUTPUT_ERROR;
        }
    }
    return CHECKPOINT_OK;
}

bool BuildIndexing::check_duplicate_(const String seq,
                                     std::vector<String> &duplicate_accessions,
                                     const Map<String, Size> &acc_to_prot,
                                     const String &acc,
                                     const seqan2::StringSet<seqan2::Peptide>  &prot_DB){
    if (acc_to_prot.has(acc)) {
        duplicate_accessions.push_back(acc);
        // check if sequence is identical
        const seqan2::Peptide &tmp_prot = prot_DB[acc_to_prot[acc]];
        if (String(begin(tmp_prot), end(tmp_prot)) != seq) {
            LOG_ERROR << "Fatal error: Protein identifier '" << acc <<
                      "' found multiple times with different sequences" << (IL_equivalent_ ? " (I/L substituted)" : "") <<
                      ":\n"
                      << tmp_prot << "\nvs.\n" << seq << "\nPlease fix the database and build index again." <<
                      std::endl;
            return BuildIndexing::DATABASE_CONTAINS_MULTIPLES;
        }
        return true;
    }
    return false;
}

bool BuildIndexing::has_aaa_(const String seq){
    for (Size j = 0; j < seq.length(); ++j) {
        if ((seq[j] == 'B') || (seq[j] == 'X') ||
            (seq[j] == 'Z')) {
            return true;
        }
    }
    return false;
}

void BuildIndexing::updateMembers_() {
    suffix_array_ = param_.getValue("suffix_array").toBool();
    FM_index_ = param_.getValue("FM_index").toBool();
    IL_equivalent_ = param_.getValue("IL_equivalent").toBool();

    log_file_ = param_.getValue("log");
    debug_ = static_cast<Size>(param_.getValue("debug")) > 0;
}

BuildIndexing::ExitCodes BuildIndexing::checkUserInput_() {
    if (!suffix_array_ && !FM_index_){
        LOG_ERROR <<
                  "Fatal error: Specify a format to build index.\n"
                  << "Please adapt your settings." << endl;
        return ILLEGAL_PARAMETERS;
    }
    if (suffix_array_ && FM_index_){
        LOG_ERROR <<
                  "Fatal error: Specify only one format to build index.\n"
                  << "Please adapt your settings." << endl;
        return ILLEGAL_PARAMETERS;
    }
    return CHECKPOINT_OK;
}