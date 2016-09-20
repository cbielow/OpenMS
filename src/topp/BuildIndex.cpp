#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/ANALYSIS/ID/BuildIndexing.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/CHEMISTRY/EnzymesDB.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/SYSTEM/File.h>


using namespace OpenMS;

class TOPPBuildIndex :
  public TOPPBase
{
public:
  TOPPBuildIndex() :
    TOPPBase("BuildIndex",
             "Build Index out of FASTA Input")
  {
  }

protected:
  void registerOptionsAndFlags_()
  {
    registerInputFile_("fasta", "<file>", "", "Input sequence database in FASTA format. Non-existing relative filenames are looked up via 'OpenMS.ini:id_db_dir'", true, false, ListUtils::create<String>("skipexists"));
    setValidFormats_("fasta", ListUtils::create<String>("fasta"));
    registerOutputFile_("out", "<file>", "", "Output Index file.");
    setValidFormats_("out", ListUtils::create<String>("txt"));

    Param temp = BuildIndexing().getParameters();
    registerFullParam_(temp);
   }

  ExitCodes main_(int, const char**)
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String in = getStringOption_("fasta");
    String out = getStringOption_("out");
//
    BuildIndexing builder;
    Param param = getParam_().copy("", true);
    Param param_pi = builder.getParameters();
    param_pi.update(param, false, Log_debug); // suppress param. update message
    builder.setParameters(param_pi);

    String db_name = getStringOption_("fasta");
    if (!File::readable(db_name))
    {
      String full_db_name;
      try
      {
        full_db_name = File::findDatabase(db_name);
      }
      catch (...)
      {
        printUsage_();
        return ILLEGAL_PARAMETERS;
      }
      db_name = full_db_name;
    }


    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------

    // we stream the Fasta file
    std::vector<FASTAFile::FASTAEntry> proteins;
    FASTAFile().load(db_name, proteins);

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    BuildIndexing::ExitCodes builder_exit = builder.run(proteins, out);
    if ((builder_exit != BuildIndexing::EXECUTION_OK))
    {
      if (builder_exit == BuildIndexing::DATABASE_EMPTY)
      {
        return INPUT_FILE_EMPTY;
      }
      else if (builder_exit == BuildIndexing::UNEXPECTED_RESULT)
      {
        return UNEXPECTED_RESULT;
      }
      else if (builder_exit == BuildIndexing::OUTPUT_ERROR)
      {
        return CANNOT_WRITE_OUTPUT_FILE;
      }
      else
      {
        return UNKNOWN_ERROR;
      }
    }

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPBuildIndex tool;
  return tool.main(argc, argv);
}
