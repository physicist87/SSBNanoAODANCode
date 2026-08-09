/////////////////////////////////////////////////
//                                             //
//                                             //
//  Author: Seungkyu Ha @ seungkyu.ha@cern.ch  //
//                                             //
//                                             //
/////////////////////////////////////////////////

#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <ctime>

#include <TROOT.h>
#include <TUnixSystem.h>
#include <TChain.h>
#include <TStyle.h>
#include <TApplication.h>
#include <TString.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>

#include "./interface/Analysis.h"
#include "./interface/DebugTools.h"

using namespace std;

TROOT root ("Plots", "Program for CMS Analysis");

//argc: # of arguments, argv:array for arguments
int main(int argc, char **argv)
{
    Logger logger;

    logger.Info() << "The number of options is: " << argc - 1 << std::endl;

    if (argc < 2)
    {
       logger.Error() << "At least, you have to set 1, 2" << std::endl;
       logger.Error() << "1. Input filelist" << std::endl;
       logger.Error() << "2. Output file" << std::endl;
       logger.Error() << "3. Config file name" << std::endl;
       logger.Error() << "4. SE directory output name" << std::endl;
       logger.Error() << "5. RunPeriod" << std::endl;
       logger.Error() << "6. Max Number Events" << std::endl;
       logger.Error() << "7. Brach List" << std::endl;
       return 1;
    }

    for (int iopt=0; iopt<argc; iopt++)
    {
       logger.Info() << "Option " << iopt << " = " << argv[iopt] << std::endl;
    }

    char *flist = argv[1];
    logger.Info() << "Input filelist = " << flist << std::endl;

    char *outname = argv[2];
    logger.Info() << "Output file name = " << outname << std::endl;

    char *confname = argv[3];
    logger.Info() << "Config file name = " << confname << std::endl;

    char *sedirname = argv[4];
    logger.Info() << "SE Dir name = " << sedirname << std::endl;

    char *runPeriod = argv[5];
    logger.Info() << "Run period = " << runPeriod << std::endl;

    int maxEvt = std::stoi(argv[6]);
    logger.Info() << "Max Events = " << maxEvt << std::endl;

    char *brachList = argv[7];
    logger.Info() << "Brach List File = " << brachList << std::endl;

    // Merge input files into a single TChain
    FILE *filelist;
    char filename[1000];
    string filelistDir, filelistName, filelistPath;

    filelistDir = "./input/";

    filelistName = argv[1];
    filelistPath = filelistDir + filelistName;
    filelist = fopen(filelistPath.c_str(),"r");

    std::vector<double> genentries_pertree;
    std::vector<double> entries_pertree;

    if (filelist == NULL)
    {
        logger.Error() << "no filelist " << filelistPath << std::endl;
        return 0;
    }

    auto ch = std::make_unique<TChain>("Events");

    logger.Info() << "start merging file(s)" << std::endl;

    while (fscanf(filelist, "%s", filename) != EOF)
    {
       logger.Info() << "adding: " << filename << std::endl;
       ch->Add(filename, 0);
       entries_pertree.push_back(ch->GetEntries());
    }
    fclose(filelist);
    logger.Info() << "Total number of events after merging root files: " << ch->GetEntries() << std::endl;

    // Analysis only observes the TChain (see interface/Analysis.h) - ch.get()
    // keeps ownership here, in main, for the lifetime of the job.
    Analysis analysis(ch.get(), filelistName, sedirname, outname, Form("./branchlist/%s", brachList), confname, maxEvt);

    analysis.SetVariables();
    analysis.Loop();

   return 0;
}

