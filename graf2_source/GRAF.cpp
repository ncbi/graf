/*  $Id$
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including all
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 */
#include "Graf.h"
#include "XchrFpSnps.h"
#include "BimFileSnps.h"
#include "FamSubjects.h"
#include "SexCheck.h"

#include <limits.h>
#include <libgen.h>
#include <unistd.h>

#if defined(__sun)
#define PROC_SELF_EXE "/proc/self/path/a.out"
#else
#define PROC_SELF_EXE "/proc/self/exe"
#endif

int minBimFpSnps = 64;          /// Minimum number of FP SNPs in bim file to read genotypes from file
int minSmpFpSnps = 200;         /// Minimum number of FP SNPs for a sample to be checked
int numThreads = 1;             /// Number of threads to be launched
SampleGeno* popSmpGeno = NULL;  /// Sample genotypes for populations
SampleGeno* relSmpGeno = NULL;  /// Sample genotypes for subject relationships
Pedigree* ped = NULL;           /// Pedigree info read from the input file
vector<int> rowCompSmps;        /// Number of samples to comapre at each row (of the pairwise comparison matrix) for each thread
RelParams params;               /// User inputed parameters

vector< vector<char*> > threadRelLines; /// Related subject pairs detected by different threads
bool debug = false;

int main(int argc, char* argv[])
{
    string usage = "Usage: graf [options]\n"
    "    -plink  PLINK set root:  File root of PLINK .bed, .bim and .fam files\n"
    "    -geno   fpg file:        Specify GRAF .fgp file\n"
    "    -exfp   PLINK set list:  Extract fingerprint genotypes from a list of PLINK sets (file roots) separated by commas\n"
    "    -pop    output file:     Check subject populations and save results to the output file\n"
    "    -sex    output file:     Do sex check using fingerprint genotypes on X-chromosome and save results to the output file\n"
    "    -out    output file:     Output file to save the results\n"
    "    -appd   DS No.:          Append extracted fingerprinting genotypes to the output file.  The integer is dataset No. of the first PLINK set\n"
    "    -ssrs   SS-RS mapping:   Specify SS# to RS# mapping file (Two columns: SS# and RS# without column headers)\n"
    "    -ped    pedigree file:   Specify pedigree file of subject IDs (with column headers)\n"
    "    -ssm    SSM file:        Specify dbGaP subject sample mapping file\n"
    "    -maxhm  max HGMR value:  Specify maximum HGMR values for a pair of subjects to be reported by GRAF\n"
    "    -xpmr   type:            Specify how expected HGMR and AGMR values are calculated for each type of relationship\n"
    "                             1: Use input dataset to calculate the expected HGMR and AGMR values\n"
    "                             2: Use average HGMR and AGMR values in dbGaP database for the expected values\n"
    "    -type   relation_type:   Specify relation type.  Acceptable values are 1, 2, 3, or 4 (default 3)\n"
    "                             1: Find all duplicates and PO pairs\n"
    "                             2: Find all duplicates, PO and FS pairs\n"
    "                             3: Find all duplicates, PO, FS and second degree relatives\n"
    "                             4: Check all the 10,000 SNPs to find all the related subjects\n"
    "NOTE:\n"
    "    1. Exactly one of the following two options should be selected: -plink or -geno.\n"
    "    2. When option -exfp is selected, -out must also be selected and output file should have .fpg extension.\n"
    "    3. When multiple PLINK sets are used, each dataset will be assigned an integer dataset ID starting with 1.\n"
    "    4. The above PLINK set starting index can be specified using option -appd.\n"
    "       When -appd is selected, the out file should be an existing GRAF .fpg file.\n"
    "    5. Multiple datasets can be combined into a single .fpg file using the -exfp and -appd options.\n"
    "    6. When multiple datasets are used, the program does pairwise comparisons to find related samples both within and across datasets.\n"
    "\n";

    string disclaimer =
    "\n *==========================================================================="
    "\n *  GRAF 2.4"
    "\n *  Yumi (Jimmy) Jin, PhD"
    "\n *  jinyu@ncbi.nlm.nih.gov"
    "\n *  Date: 04/24/2020"
    "\n *"
    "\n *                            PUBLIC DOMAIN NOTICE"
    "\n *               National Center for Biotechnology Information"
    "\n *"
    "\n *  This software/database is a \"United States Government Work\" under the"
    "\n *  terms of the United States Copyright Act.  It was written as part of"
    "\n *  the author's official duties as a United States Government employee and"
    "\n *  thus cannot be copyrighted.  This software/database is freely available"
    "\n *  to the public for use. The National Library of Medicine and the U.S."
    "\n *  Government have not placed any restriction on its use or reproduction."
    "\n *"
    "\n *  Although all reasonable efforts have been taken to ensure the accuracy"
    "\n *  and reliability of the software and data, the NLM and the U.S."
    "\n *  Government do not and cannot warrant the performance or results that"
    "\n *  may be obtained by using this software or data. The NLM and the U.S."
    "\n *  Government disclaim all warranties, express or implied, including"
    "\n *  warranties of performance, merchantability or fitness for any particular"
    "\n *  purpose."
    "\n *"
    "\n *  Please cite the author in any work or product based on this material."
    "\n *"
    "\n *===========================================================================";

    if (argc < 3 || argc % 2 == 0) {
        cout << disclaimer << "\n\n";
        cout << usage;
        exit(1);
    }

    struct timeval t1, t2;
    gettimeofday(&t1, NULL);

    int i, j, k;

    /// Make sure long integer has 64 bits.
    int numLongBytes = sizeof(unsigned long);
    if (numLongBytes != 8) {
        cout << "This program only runs on 64-bit machine with 8-byte long integers.\n";
        exit(1);
    }

    char* ssmFile    = NULL;
    char* outFile    = NULL;
    char* outPopFile = NULL;
    char* outSexFile = NULL;
    char* svSnpFile  = NULL;
    char* baseFile   = NULL;
    char* ssRsFile   = NULL;
    char* genoFile   = NULL;
    char* pedigreeFile = NULL;
    char* plinkList  = NULL;
    char* fpSnpFile  = NULL;

    int checkRelType = 3;
    int chkWordNum = 10;
    int minHomoDiff = 58;
    double maxHgmr = 18;
    bool isUserHgmr = false;;
    int stDsNo = 0;
    int argPairs = (argc - 1) / 2;
    int hasSbjId = 0;
    bool argValid = true;
    int expMrCalc = 0;
    bool ssmIsGap = false;

    bool checkPop = false;
    bool checkSex = false;

    for (i = 0; i < argPairs; i++) {
        char* option = argv[i*2+1];
        char* val = argv[i*2+2];

        if (option[0] == '-') {
            if (strcmp(option, "-ped") == 0) {
                pedigreeFile = val;
            }
            else if (strcmp(option, "-ssm") == 0) {
                ssmFile = val;
            }
            else if (strcmp(option, "-gapssm") == 0) {
	        ssmIsGap = true;
            }
            else if (strcmp(option, "-fpsnps") == 0) {
                fpSnpFile = val;
            }
            else if (strcmp(option, "-ssrs") == 0) {
                ssRsFile = val;
            }
            else if (strcmp(option, "-exfp") == 0) {
                plinkList = val;
            }
            else if (strcmp(option, "-maxhm") == 0) {
                maxHgmr = atof(val);
                isUserHgmr = true;
            }
            else if (strcmp(option, "-sbj") == 0) {
                hasSbjId = 1;
            }
            else if (strcmp(option, "-plink") == 0) {
                baseFile = val;
                int fileLen = strlen(baseFile);
                if (fileLen > 4 && baseFile[fileLen-4] == '.' &&
                    ( (baseFile[fileLen-3] == 'b' && baseFile[fileLen-2] == 'e' && baseFile[fileLen-1] == 'd') ||
                      (baseFile[fileLen-3] == 'b' && baseFile[fileLen-2] == 'i' && baseFile[fileLen-1] == 'm') ||
                      (baseFile[fileLen-3] == 'f' && baseFile[fileLen-2] == 'a' && baseFile[fileLen-1] == 'm') )  ) {
                    baseFile[fileLen-4] = 0;
                }
            }
            else if (strcmp(option, "-geno") == 0) {
                if (strstr(val, ".fpg")) {
                    genoFile = val;
                }
                else {
                    cout << "\n********** ERROR: Geno file should have .fpg extension!\n";
                    argValid = false;
                    break;
                }
            }
            else if (strcmp(option, "-out") == 0) {
                outFile = val;
            }
            else if (strcmp(option, "-pop") == 0) {
                checkPop = true;
                outPopFile = val;
            }
            else if (strcmp(option, "-sex") == 0) {
                checkSex = true;
                outSexFile = val;
            }
            else if (strcmp(option, "-appd") == 0) {
                stDsNo = atoi(val);
                stDsNo = stDsNo - 1;
            }
            else if (strcmp(option, "-type") == 0) {
                checkRelType = atoi(val);
                if (checkRelType < 1 || checkRelType > 4) {
                    cout << "\n********** ERROR: Invalid relation type " << val << ". Should be 1, 2, 3, or 4!\n";
                    argValid = false;
                    break;
                }
            }
            else if (strcmp(option, "-m") == 0) {
                maxHgmr = atof(val);
                isUserHgmr = true;
            }
            else if (strcmp(option, "-s") == 0) {
                minSmpFpSnps = atof(val);
            }
            else if (strcmp(option, "-xpmr") == 0) {
                expMrCalc = atoi(val);
                if (expMrCalc < 1 || expMrCalc > 2) {
                    cout << "\n********** ERROR: Invalid xpmr value " << val << ". Should be 1 or 2!\n";
                    argValid = false;
                    break;
                }
            }
            else if (strcmp(option, "-debug") == 0) {
                debug = true;
            }
        }
        else {
            argValid = false;
            break;
        }
    }

    if (checkRelType == 1) {
        chkWordNum = 3;
        minHomoDiff = 5;
    }
    else if (checkRelType == 2) {
        chkWordNum = 6;
        minHomoDiff = 23;
    }
    else if (checkRelType == 3) {
        chkWordNum = 10;
        minHomoDiff = 58;
    }
    else if (checkRelType > 3) {
        chkWordNum = 1000;
        minHomoDiff = 1000;
    }

    if (debug) cout << "Minimum number of FP SNPs required: " << minSmpFpSnps << "\n";

    numThreads = thread::hardware_concurrency();
    numThreads--;

    if (plinkList) {
        if (!outFile) {
            cout << "\n********** ERROR: No output file specified for extracting fingerprinting genotypes option!\n";
        }
        else {
            int nameLen = strlen(outFile);
            if (outFile[nameLen-4] != '.' || outFile[nameLen-3] != 'f' || outFile[nameLen-2] != 'p' || outFile[nameLen-1] != 'g') {
                cout << "\n********** ERROR: Invalid output file " << outFile << ". With -exfp option, output file should have .fpg extension.\n" << "\n";
                argValid = false;
            }
        }
    }
    else {
        if (!baseFile && !genoFile) {
            cout << "\n********** ERROR: Either option -plink or option -geno must be selected!\n";
            argValid = false;
        }

        if (FileExists(outFile)) {
            cout << "\n********** ERROR: Output file " << outFile << " already exists! Please use a different output file or delete the old file.\n";
            argValid = false;
        }
    }

    if (!argValid) {
        cout << "\n" << usage << "\n";
        exit(1);
    }

    /// Read SS to RS mapping if mapping file exists
    map<int, int> ssToRs;
    if (ssRsFile) {
        char ssLine[256];
        int numSs = 0;
        FILE *ssRsPtr = fopen(ssRsFile, "r");
        if(ssRsPtr) {
            int ss, rs;
            while (fgets(ssLine, 256, ssRsPtr) != NULL) {
                sscanf(ssLine, "%d %d", &ss, &rs);
                if (ss && rs) {
                    ssToRs[ss] = rs;
                    numSs++;
                }
            }
            fclose(ssRsPtr);
            cout << "Read " << numSs << " SS numbers from " << ssRsFile << ".\n";
        }
    }

    params.plinkBase    = baseFile;
    params.genoFile     = genoFile;
    params.outFile      = outFile;
    params.outPopFile   = outPopFile;
    params.outSexFile   = outSexFile;
    params.ssmFile      = ssmFile;
    params.fpSnpFile    = fpSnpFile;
    params.pedigreeFile = pedigreeFile;
    params.checkRelType = checkRelType;
    params.chkWordNum   = chkWordNum;
    params.minHomoDiff  = minHomoDiff;
    params.maxHgmr      = maxHgmr;
    params.expMrCalc    = expMrCalc;
    params.isUserHgmr   = isUserHgmr;
    params.ssmIsGap     = ssmIsGap;

    bool isPlink = false;
    if (baseFile) isPlink = true;

    vector<string> logLines;
    string logLine = "";
    for (i = 0; i < argc; i++) logLine = logLine + argv[i] + " ";
    logLines.push_back(logLine);

    if (plinkList) {
        int totSmps = CombinePlinkSets(plinkList, outFile, ssToRs, stDsNo);
        char* logLine = (char*)malloc(sizeof(char) * 256);
        sprintf(logLine, "Extrated FP genotypes for %d samples.", totSmps);
        logLines.push_back(logLine);
    }
    else if (checkPop) {
        if (!FileWriteable(outPopFile)) {
            cout << "********** ERROR: Couldn't open " << outPopFile << " for writing!\nPlease specify an output file to save the population results.\n\n"
            << usage << "\n";
            exit(1);
        }

        int numPopSmps = CheckSubjectPopulations(params, isPlink, ssToRs, hasSbjId);
    }
    else if (checkSex) {
        if (!FileWriteable(outSexFile)) {
            cout << "********** ERROR: Couldn't open " << outSexFile << " for writing!\nPlease specify an output file to save the sex check results.\n\n"
            << usage << "\n";
            exit(1);
        }

        int numSexSmps = CheckSubjectSexes(params, hasSbjId);
    }
    else {
        if (!outFile) {
            time_t t = time(NULL);
            struct tm tm = *localtime(&t);
            outFile = (char*)malloc(sizeof(char) * 128);
            sprintf(outFile, "graf_rel_%04d%02d%02d_%02d%02d.txt", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min);
        }

        params.outFile = outFile;

        if (!FileWriteable(outFile)) {
            cout << "********** ERROR: Couldn't open " << outFile << " for writing!\nPlease specify an output file to save the results.\n\n"
            << usage << "\n";
            exit(1);
        }

        for (int thNo = 0; thNo < numThreads; thNo++) {
            vector<char*> thisThIds;
            threadRelLines.push_back(thisThIds);
        }

        Relatives* relatives = FindRelatedSubjects(params, isPlink, ssToRs, hasSbjId);

        if (relatives) {
            char* logLine = (char*)malloc(sizeof(char) * 256);
            sprintf(logLine, "Found %d pairs of related samples.", relatives->numPairs);
            logLines.push_back(logLine);
        }
    }

    gettimeofday(&t2, NULL);
    ShowTimeDiff(t1, t2);

    return 0;
}

int CheckSubjectPopulations(RelParams params, bool isPlink, map<int, int> ssToRs, int hasSbjId)
{
    char* baseFile     = params.plinkBase;
    char* genoFile     = params.genoFile;
    char* outFile      = params.outFile;
    char* outPopFile   = params.outPopFile;
    char* ssmFile      = params.ssmFile;
    char* pedigreeFile = params.pedigreeFile;
    int   checkRelType = params.checkRelType;
    int   chkWordNum   = params.chkWordNum;
    int   minHomoDiff  = params.minHomoDiff;

    minSmpFpSnps = 50;

    int numSmps = 0;
    int numStds = 1;
    int i, j, k;

    char bimFile[256];
    char bedFile[256];
    char famFile[256];

    if (isPlink) {
        sprintf(bimFile, "%s.bim", baseFile);
        sprintf(bedFile, "%s.bed", baseFile);
        sprintf(famFile, "%s.fam", baseFile);

        if (!FileExists(bimFile)) {
            cout << "ERROR: file " << bimFile << " doesn't exist.\n\n";
            return 0;
        }
        if (!FileExists(bedFile)) {
            cout << "ERROR: file " << bedFile << " doesn't exist.\n\n";
            return 0;
        }
        if (!FileExists(famFile)) {
            cout << "ERROR: file " << famFile << " doesn't exist.\n\n";
            return 0;
        }
        if (ssmFile && !FileExists(ssmFile)) {
            cout << "ERROR: file " << ssmFile << " doesn't exist.\n\n";
            return 0;
        }
        if (pedigreeFile && !FileExists(pedigreeFile)) {
            cout << "ERROR: file " << pedigreeFile << " doesn't exist.\n\n";
            return 0;
        }

        popSmpGeno = new SampleGeno(famFile, 0);            /// Read genotyped sample list from fam file
        popSmpGeno->SetMaxSnpMiss(minSmpFpSnps);
        numSmps = popSmpGeno->numSmps;
        cout << "\nFound " << numSmps << " samples in " << famFile << "\n";

        bool hasErr = popSmpGeno->ReadGenotypesFromBedFile(bimFile, bedFile, minBimFpSnps, ssToRs);
        if (hasErr) {
            cout << "ERROR: Failed to read genotypes from PLINK set " << bedFile << "\n";
            return 0;
        }

        cout << "\tRead genotypes of " << numSmps << " samples from " << bedFile << ".\n";
    }
    else {
        popSmpGeno = new SampleGeno(genoFile);
        popSmpGeno->SetMaxSnpMiss(minSmpFpSnps);
        numSmps = popSmpGeno->numSmps;
        numStds = popSmpGeno->GetNumberStudies();
        cout << "Found " << numSmps << " samples from " << numStds << " datasets in " << genoFile << "\n";

        int genoSmps = popSmpGeno->ReadGenotypesFromHexGenoFile(genoFile, hasSbjId);
        if (!genoSmps) {
            cout << "ERROR: failed to read genotypes from " << genoFile << "\n";
            return 0;
        }

        numSmps = popSmpGeno->numSmps;
        cout << "\tRead genotypes of " << genoSmps << " samples from geno file.\n";
        if (debug) cout << "Samples with enough genotypes = " << numSmps << "\n";
    }

    popSmpGeno->InitPopPvalues();

    cout << "\n\tLaunching " << numThreads << " threads to calculate genetic distances.\n";
    popSmpGeno->numThreads = numThreads;

    mutex iomutex;
    vector<thread> threads(numThreads);

    for (unsigned i = 0; i < numThreads; ++i) {
        threads[i] = thread([&iomutex, i] {
            {
                lock_guard<mutex> iolock(iomutex);
            }

            popSmpGeno->SetPopulationPvalues(i);
        });
    }

    for (auto& t : threads) {
        t.join();
    }

    cout << "Done. Calculated genetic distances for " << numSmps << " samples.\n\n";

    popSmpGeno->SavePopulationResults(outPopFile);
    return numSmps;
}

int CheckSubjectSexes(RelParams params, int hasSbjId)
{
    char* baseFile     = params.plinkBase;
    char* outSexFile   = params.outSexFile;
    char* ssmFile      = params.ssmFile;
    char* pedigreeFile = params.pedigreeFile;

    string bimFile = baseFile + string(".bim");
    string bedFile = baseFile + string(".bed");
    string famFile = baseFile + string(".fam");

    if (!FileExists(bimFile.c_str())) {
        cout << "ERROR: file " << bimFile << " doesn't exist.\n\n";
        return 0;
    }
    if (!FileExists(bedFile.c_str())) {
        cout << "ERROR: file " << bedFile << " doesn't exist.\n\n";
        return 0;
    }
    if (!FileExists(famFile.c_str())) {
        cout << "ERROR: file " << famFile << " doesn't exist.\n\n";
        return 0;
    }

    XchrFpSnps *sexFpSnps = new XchrFpSnps();

    string grafDir = GetExecutablePath();
    string fpSnpFile = grafDir + "/data/FP_SNPs.txt";
    if (params.fpSnpFile) {
        fpSnpFile = string(params.fpSnpFile);
        if (!FileExists(fpSnpFile.c_str())) {
            cout << "ERROR: file " << fpSnpFile << " doesn't exist.\n\n";
            return 0;
        }
    }
    else {
        if (!FileExists(fpSnpFile.c_str())) {
            cout << "ERROR: file " << fpSnpFile << " doesn't exist. Please use -fpsnps option to specify the file including FP SNPs.\n\n";
            return 0;
        }
    }
    int numXfpSnps = sexFpSnps->ReadXchrFpSnpsFromFile(fpSnpFile);

    BimFileSnps *bimSnps = new BimFileSnps(30);
    int numBimSnps = bimSnps->ReadXchrFpSnpsFromFile(bimFile, sexFpSnps);
    int numSmps = 0;

    if (bimSnps->HasEnoughXchrFpSnps()) {
        bimSnps->ShowSummary();

        FamSubjects *famSbjs = new FamSubjects(famFile);
        int numFamSbjs = famSbjs->ReadSubjectsFromFile();
        famSbjs->ShowSummary();

        SexCheck *sexCheck = new SexCheck(bedFile, sexFpSnps, bimSnps, famSbjs);
        bool hasPlinkErr = sexCheck->CalculateHetroRates();
        if (hasPlinkErr) {
            cout << "Output file " << outSexFile << " is not created due to the errors.\n\n";
        }
        else {
            sexCheck->ShowSummary();
            int numSmps = sexCheck->SaveResults(outSexFile);

            cout << "Saved sex check results of " << numSmps << " samples to " << outSexFile << "\n";
        }
    }
    else {
        int numFpSnps = bimSnps->GetXchrFpSnps();
        int numGoodFpSnps = bimSnps->GetGoodFpSnps();

        cout << "\nWARNING: Found " << numFpSnps << " X-chromosome FP SNPs in bim file. "
        << numGoodFpSnps << " SNPs have expected alleles, not enough to do sex checking.\n"
        << "If the SNP IDs in the bim file are not RS IDs, please run SetRsIdsInBimFile.pl to convert them into RS IDs.\n\n";
        numSmps = 0;
    }

    return numSmps;
}

Relatives* FindRelatedSubjects(RelParams params, bool isPlink, map<int, int> ssToRs, int hasSbjId)
{
    char* baseFile     = params.plinkBase;
    char* genoFile     = params.genoFile;
    char* outFile      = params.outFile;
    char* ssmFile      = params.ssmFile;
    char* pedigreeFile = params.pedigreeFile;
    int   checkRelType = params.checkRelType;
    int   chkWordNum   = params.chkWordNum;
    int   minHomoDiff  = params.minHomoDiff;
    double maxHgmr     = params.maxHgmr;
    bool isUserHgmr    = params.isUserHgmr;

    int numSmps = 0;
    int numStds = 1;
    int i, j, k;

    char bimFile[256];
    char bedFile[256];
    char famFile[256];

    if (isPlink) {
        sprintf(bimFile, "%s.bim", baseFile);
        sprintf(bedFile, "%s.bed", baseFile);
        sprintf(famFile, "%s.fam", baseFile);

        if (!FileExists(bimFile)) {
            cout << "ERROR: file " << bimFile << " doesn't exist.\n\n";
            return NULL;
        }
        if (!FileExists(bedFile)) {
            cout << "ERROR: file " << bedFile << " doesn't exist.\n\n";
            return NULL;
        }
        if (!FileExists(famFile)) {
            cout << "ERROR: file " << famFile << " doesn't exist.\n\n";
            return NULL;
        }
        if (ssmFile && !FileExists(ssmFile)) {
            cout << "ERROR: file " << ssmFile << " doesn't exist.\n\n";
            return NULL;
        }
        if (pedigreeFile && !FileExists(pedigreeFile)) {
            cout << "ERROR: file " << pedigreeFile << " doesn't exist.\n\n";
            return NULL;
        }

        relSmpGeno = new SampleGeno(famFile, 0);            /// Read genotyped sample list from fam file
        numSmps = relSmpGeno->numSmps;
        cout << "Found " << numSmps << " samples in " << famFile << "\n";

        if (numSmps == 1) {
            cout << "\nNOTE: there is only one sample in " << famFile << ". No relationship check can be done.\n";
            return 0;
        }
    }
    else {
        relSmpGeno = new SampleGeno(genoFile);
        numSmps = relSmpGeno->numSmps;
        numStds = relSmpGeno->GetNumberStudies();
        cout << "Found " << numSmps << " samples from " << numStds << " datasets in " << genoFile << "\n";
    }

    relSmpGeno->SetMaxSnpMiss(minSmpFpSnps);

    ped = new Pedigree(relSmpGeno->samples, numSmps);     /// Initialize sample mapping to itself
    ped->SetSsmIsGapFormat(params.ssmIsGap);
    if (ssmFile) {
        bool readSsm = ped->ReadSmpSbjMapping(ssmFile);     /// Set sample mapping
        if (!readSsm) return NULL;
    }
    relSmpGeno->SetSampleSubjectIds(ped->ssm);

    if (pedigreeFile) {
        ped->ReadPedigree(pedigreeFile, false);
    }

    if (ssmFile || pedigreeFile) {
        ped->GetPedigreeRelatives();
        ped->ShowSummary(0);
        int numSexUpds = relSmpGeno->SetSampleGenders(ped->sbjSexes);
        if (numSexUpds > 0) cout << "Set genders for " << numSexUpds << " samples.\n";
    }

    if (isPlink) {
        bool hasErr = relSmpGeno->ReadGenotypesFromBedFile(bimFile, bedFile, minBimFpSnps, ssToRs);
        if (hasErr) {

            return NULL;
        }

        cout << "\tRead genotypes of " << numSmps << " samples from bed file.\n";
    }
    else {
        int genoSmps = relSmpGeno->ReadGenotypesFromHexGenoFile(genoFile, hasSbjId);
        if (!genoSmps) {
            cout << "ERROR: failed to read genotypes from " << genoFile << "\n";
            return NULL;
        }

        numSmps = relSmpGeno->numSmps;
        cout << "\tRead genotypes of " << genoSmps << " samples from geno file.\n";
        if (debug) cout << "Samples with enough genotypes = " << numSmps << "\n";
    }

    /// In case we get too many related pairs for large datasets
    if (!isUserHgmr) {
        if (numSmps >  10000) params.maxHgmr = 17;
        if (numSmps >  30000) params.maxHgmr = 16;
        if (numSmps > 100000) params.maxHgmr = 15;
        if (numSmps > 300000) params.maxHgmr = 14;
    }

    cout << "\nChecking genotypes of " << numSmps << " samples to find related subjects.\n";

    relSmpGeno->ReorderSnps();

    double hgmrAdj = relSmpGeno->hgmrAdj;
    int minSmpsForExpMr = 100;
    for (int i = 0; i < numSmps; i++) rowCompSmps.push_back(numSmps - i - 1);

    Relatives* allRelatives = new Relatives(hgmrAdj);
    /// expMrCalc = 0:
    ///    #samples >= 100: using genotypes in dataset to calculat expected MRs
    ///    #samples <  100: using average MRs from dbGaP database
    /// expMrCalc = 1:      using genotypes in dataset to calculat expected MRs
    /// expMrCalc = 2:      using average MRs from dbGaP database
    if ((params.expMrCalc == 0 && numSmps >= minSmpsForExpMr) || params.expMrCalc == 1) {
        allRelatives->SetExpectedMismatchRates(relSmpGeno->poExpAgmr, relSmpGeno->fsExpAgmr, relSmpGeno->d2ExpAgmr, relSmpGeno->d3ExpAgmr,
        relSmpGeno->fsExpHgmr, relSmpGeno->d2ExpHgmr, relSmpGeno->d3ExpHgmr);
    }

    cout << "\tLaunching " << numThreads << " threads to pairwisely compare subject genotypes.\n";
    mutex iomutex;
    vector<thread> threads(numThreads);

    for (unsigned i = 0; i < numThreads; ++i) {
        threads[i] = thread([&iomutex, i] {
            {
                lock_guard<mutex> iolock(iomutex);
            }
            FindRelatedSubjectsForThread(i);
        });
    }

    for (auto& t : threads) {
        t.join();
    }

    int numRelPairs = 0;
    for (int i = 0; i < numThreads; i++) {
        int thNumPairs = threadRelLines[i].size();
        numRelPairs += thNumPairs;


        for (int j = 0; j < thNumPairs; j++) {
            int phs1, phs2, sex1, sex2, allSame, allDiff, homoSame, homoDiff, pedRelNo;
            double hgmr, agmr;
            char smpName1[80], smpName2[80], sbjName1[80], sbjName2[80];
            sscanf(threadRelLines[i][j], "%d\t%s\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%lf\t%d\t%d\t%lf\t%d",
            &phs1, smpName1, sbjName1, &sex1, &phs2, smpName2, sbjName2, &sex2, &allSame, &allDiff, &agmr, &homoSame, &homoDiff, &hgmr, &pedRelNo);
            SampleInfo sample1(smpName1, sbjName1, sex1, 0, 0, phs1);
            SampleInfo sample2(smpName2, sbjName2, sex2, 0, 0, phs2);
            allRelatives->AddOnePair(sample1, sample2, homoSame, homoDiff, hgmr, allSame, allDiff, agmr, pedRelNo);

            delete threadRelLines[i][j];
        }
        threadRelLines[i].clear();
    }
    threadRelLines.clear();

    allRelatives->ShowSummary();
    allRelatives->SaveRelativesToFile(outFile, numStds);

    return allRelatives;
}

void FindRelatedSubjectsForThread(int thNo)
{
    int   checkRelType = params.checkRelType;
    int   chkWordNum   = params.chkWordNum;
    int   minHomoDiff  = params.minHomoDiff;
    double maxHgmr     = params.maxHgmr;

    unsigned long totComps = 0;
    int numRelPairs = 0;
    unsigned long pairNo = 0;

    int numSmps = relSmpGeno->numSmps;
    SampleInfo* samples = relSmpGeno->samples;
    double hgmrAdj = relSmpGeno->hgmrAdj;

    int numValidSnps = relSmpGeno->numValidSnps;
    int numGenoInts = relSmpGeno->numFpGenoInts;

    unsigned long totPairs = (unsigned long) numSmps * (numSmps-1) / 2;

    unsigned long meanPairs = totPairs / numThreads;
    int rmPairs = totPairs % numThreads;
    unsigned long chkPairs = meanPairs;
    if (thNo+1 <= rmPairs) chkPairs = meanPairs + 1;

    unsigned long stPair = thNo * chkPairs;
    if (thNo >= rmPairs) stPair = thNo * chkPairs + rmPairs;
    unsigned long edPair = stPair + chkPairs - 1;

    int stRow, edRow, stCol, edCol;

    unsigned long nPairs = 0;  /// Number of pairs compared for each thread
    for (int row = 0; row < numSmps; row++) {
        nPairs += rowCompSmps[row];

        if (nPairs > stPair) {
            stRow = row;
            stCol = stPair - (nPairs - rowCompSmps[row]);
            break;
        }
    }

    nPairs = 0;
    for (int row = 0; row < numSmps; row++) {
        nPairs += rowCompSmps[row];

        if (nPairs > edPair) {
            edRow = row;
            edCol = edPair - (nPairs - rowCompSmps[row]);
            break;
        }
    }

    if (debug)
        cout << "No. " << thNo << ": " << stPair << " -> " << edPair << " pairs " << chkPairs << " tot " << totPairs
        << " Start: " << stRow << " : " << stCol
        << " End: " << edRow << " : " << edCol
        << "\n" ;

        unsigned long iG1, iG2, jG1, jG2;
        unsigned long isOi, isOj, isXi, isXj;
        unsigned long homoSameBits, homoDiffBits, hetroSameBits, ibs1Bits;
        int homoSame, homoDiff, allSame, allDiff, hetroSame, ibs1, pedRelNo, isRel;
        double hgmr, agmr;
        int smpId1, sbjId1, sex1, phs1, smpId2, sbjId2, sex2, phs2;
        char *smpName1, *sbjName1, *smpName2, *sbjName2;
        unsigned long smp1Val;
        map<string, int>::iterator it1, it2;
        string pair1, pair2;
        int totMPairs = chkPairs/1000000;
        if (thNo == 0) cout << "\tThread No. 1: comparing total " << chkPairs << " pairs of samples.\n";
        if (debug) cout << "HGMR cutoff = " << maxHgmr << "\n";
        int j, k;

        int row = stRow;
        while(row <= edRow) {
        int lastCol = rowCompSmps[row] - 1;
        if (row == edRow) lastCol = edCol;
        int firstCol = 0;
        if (row == stRow) firstCol = stCol;

        smpId1   = samples[row].smpId;
        smpName1 = samples[row].smpName;
        sbjName1 = samples[row].sbjName;
        sex1     = samples[row].sex;
        phs1     = samples[row].phs;
        smp1Val = (unsigned long)smpId1 * 10000000;

        for (j = firstCol; j <= lastCol; j++) {
            int col = row + j + 1;

            smpId2   = samples[col].smpId;
            smpName2 = samples[col].smpName;
            sbjName2 = samples[col].sbjName;
            sex2     = samples[col].sex;
            phs2     = samples[col].phs;
            homoSame = 0;
            homoDiff = 0;
            allSame = 0;
            allDiff = 0;
            hetroSame = 0;
            ibs1 = 0;
            pedRelNo = 0;

            pair1 = (string)smpName1 + "\t" + (string)smpName2;
            pair2 = (string)smpName2 + "\t" + (string)smpName1;
            it1 = ped->smpRels.find(pair1);
            if (it1 == ped->smpRels.end()) {
                it2 = ped->smpRels.find(pair2);
                if (it2 != ped->smpRels.end()) pedRelNo = ped->smpRels[pair2];
            }
            else {
                pedRelNo = ped->smpRels[pair1];
            }

            isRel = 1;
            for (k = 0; k < numGenoInts; k++) {
                iG1 = relSmpGeno->smpGenos1[row][k];
                iG2 = relSmpGeno->smpGenos2[row][k];
                jG1 = relSmpGeno->smpGenos1[col][k];
                jG2 = relSmpGeno->smpGenos2[col][k];

                isOi = ~(iG1 ^ iG2); /// homozygous SNPs of sample i
                isOj = ~(jG1 ^ jG2);
                isXi = ~iG1 & iG2;   /// hetrozygous SNPs of sample i: 01. 10 for missing
                isXj = ~jG1 & jG2;

                homoSameBits = ~(iG1 ^ jG2) & isOi & isOj;
                homoDiffBits = (iG1 ^ jG2) & isOi & isOj;
                hetroSameBits = isXi & isXj;
                ibs1Bits = (isOi & isXj) | (isOj & isXi);

                homoSame  += __builtin_popcountl(homoSameBits);
                homoDiff  += __builtin_popcountl(homoDiffBits);
                hetroSame += __builtin_popcountl(hetroSameBits);
                ibs1      += __builtin_popcountl(ibs1Bits);

                if (chkWordNum && k+1 == chkWordNum && homoDiff > minHomoDiff && !pedRelNo) {
                    /// Assumming PO MR 0.02.  We are 99.99% sure this is not PO.  Quit loop.
                    isRel = 0;
                    break;
                }
            }
            totComps += k+1;

            allSame = homoSame + hetroSame;
            allDiff = homoDiff + ibs1;
            if (smpId1 && smpId2 && isRel &&
                (    homoSame + homoDiff > 400 /// Enough number of homogenous SNPs to compare
                || (allSame > 50 && (allDiff * 4 < allSame || pedRelNo == 1 || pedRelNo == 2) )
                ) /// Enough SNPs for ID agmr < 0.2
                ) {
                /// Don't normalize HGMR. Use the expected HGMR for each relationship to calculate p-value

                double hgmr = homoDiff * 100.0 / (homoSame + homoDiff);
                if (hgmr < maxHgmr || pedRelNo) {
                    agmr = allDiff * 100.0 / (allSame + allDiff);

                    char* line = (char*)malloc(sizeof(char) * 256);
                    sprintf(line, "%d\t%s\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%f\t%d\t%d\t%f\t%d\n",
                    phs1, smpName1, sbjName1, sex1, phs2, smpName2, sbjName2, sex2, allSame, allDiff, agmr, homoSame, homoDiff, hgmr, pedRelNo);
                    threadRelLines[thNo].push_back(line);

                    numRelPairs++;
                }
            }

            pairNo++;

            if (pairNo % 10000000 == 0 && thNo == 0) {
                int mPairNo = pairNo / 1000000;
                printf("\t\tCompared %d M pairs of samples. Found %d pairs of related subjects.  Total %d M pairs to compare.\n", mPairNo, numRelPairs, totMPairs);
                if (debug) {
                    double avgComp = totComps * 1.0 / pairNo;
                    printf("\t\t\tAverage %3.1f numbers per pair compared.\n", avgComp);
                }
            }
        }

        row++;
    }

    if (thNo == 0) {
        if (chkPairs < 5000000) {
            printf("\t\tCompared total %d pairs of samples.", chkPairs);
        }
        else {
            printf("\t\tCompared total %d M pairs of samples.", totMPairs);
        }
        printf(" Found %d pairs of relatives.\n\n", numRelPairs);

        double avgComp = totComps * 1.0 / pairNo;
        if (debug) printf("Compared averaged %3.1f numbers per pair.\n", avgComp);
    }
}

int CombinePlinkSets(char* plinkList, char* outFile, map<int, int> ssToRs, int stDsNo)
{
    int totSmps = 0;

    int i, j, k;

    vector<char*> plinkSets;
    int len = strlen(plinkList);
    char file[512];
    int nameLen = 0;
    for (i = 0; i <= len; i++) {
        char ch = plinkList[i];
        if (ch == ',' || i == len) {
            file[nameLen] = 0;
            if (nameLen > 4 && (
                (file[nameLen-4] == '.' && file[nameLen-3] == 'b' && file[nameLen-2] == 'e' && file[nameLen-1] == 'd') ||
                (file[nameLen-4] == '.' && file[nameLen-3] == 'b' && file[nameLen-2] == 'i' && file[nameLen-1] == 'm') ||
                (file[nameLen-4] == '.' && file[nameLen-3] == 'f' && file[nameLen-2] == 'a' && file[nameLen-1] == 'm') )  ) {
                file[nameLen-4] = 0;
            }

            char* newFile = (char*)malloc(sizeof(char)*(nameLen+1));
            strcpy(newFile, file);
            plinkSets.push_back(newFile);
            nameLen = 0;
        }
        else {
            file[nameLen] = ch;
            nameLen++;
        }
    }
    int numPlinkSets = plinkSets.size();

    if (FileExists(outFile) && stDsNo < 1) {
        cout << "ERROR: Output file " << outFile << " already exists! Please use a different output file or delete the old file.\n";
        return totSmps;
    }

    if (!FileWriteable(outFile)) {
        cout << "ERROR: Couldn't open " << outFile << " for writing!\n";
        return totSmps;
    }

    /// Make sure all the files are there
    vector<string> famFiles;
    vector<string> bimFiles;
    vector<string> bedFiles;
    vector<char*> errLines;
    for (i = 0; i < numPlinkSets; i++) {
        char* baseFile = plinkSets[i];

        char* famFile = (char*)malloc(sizeof(char) * 256);
        char* bimFile = (char*)malloc(sizeof(char) * 256);
        char* bedFile = (char*)malloc(sizeof(char) * 256);

        sprintf(famFile, "%s.fam", baseFile);
        sprintf(bimFile, "%s.bim", baseFile);
        sprintf(bedFile, "%s.bed", baseFile);

        famFiles.push_back(famFile);
        bimFiles.push_back(bimFile);
        bedFiles.push_back(bedFile);

        if (!FileExists(bimFile)) {
            char* errFile = (char*)malloc(sizeof(char) * 256);
            strcpy(errFile, bimFile);
            errLines.push_back(errFile);
        }
        if (!FileExists(bedFile)) {
            char* errFile = (char*)malloc(sizeof(char) * 256);
            strcpy(errFile, bedFile);
            errLines.push_back(errFile);
        }
        if (!FileExists(famFile)) {
            char* errFile = (char*)malloc(sizeof(char) * 256);
            strcpy(errFile, famFile);
            errLines.push_back(errFile);
        }
    }

    int numFamFiles = famFiles.size();
    int numErrLines = errLines.size();

    if (numErrLines > 0) {
        cout << "ERROR: Didn't find the following " << numErrLines << " files:\n";
        for (i = 0; i < numErrLines; i++) {
            cout << "\t" << errLines[i] << "\n";
        }

        return totSmps;
    }

    /// Check all the fam files to see if all the sample IDs are unique
    map<string, int> allSmpIds;
    map<string, int>::iterator smpIt;
    int numDupSmps = 0;
    int numTotSmps = 0;
    int dsSmpCnts[numFamFiles];

    for (i = 0; i < numFamFiles; i++) {
        string famFile = famFiles[i];
        bool hasErr = false;
        vector<string> dsSmps = FindSamplesInFamFile(famFile.c_str(), &hasErr);
        if (hasErr) return totSmps;

        int numSmps = dsSmps.size();
        dsSmpCnts[i] = numSmps;

        for (j = 0; j < numSmps; j++) {
            string smp = dsSmps[j];
            smpIt = allSmpIds.find(smp);
            if (smpIt == allSmpIds.end()) {
                allSmpIds[smp] = 1;
            }
            else {
                numDupSmps++;
            }

            numTotSmps++;
        }
    }

    cout << "Found total " << numTotSmps << " samples in " << numFamFiles << " fam files.\n";
    if (numDupSmps > 0) cout << numDupSmps << " sample IDs are found in more than one file.\n";
    for (i = 0; i < numFamFiles; i++)
    cout << "\t" << dsSmpCnts[i] << "\t" << famFiles[i] << "\n";

    for (i = 0; i < numFamFiles; i++) {
        int dsNo = i + stDsNo + 1;
        char* famFile = (char*)famFiles[i].c_str();
        char* bimFile = (char*)bimFiles[i].c_str();
        char* bedFile = (char*)bedFiles[i].c_str();

        SampleGeno* smpGeno = new SampleGeno(famFile, dsNo);
        smpGeno->SetMaxSnpMiss(minSmpFpSnps);
        int numSmps = smpGeno->numSmps;
        cout << "Found " << numSmps << " samples in " << famFile << "\n\n";

        bool hasErr = smpGeno->ReadGenotypesFromBedFile(bimFile, bedFile, 200, ssToRs);
        if (hasErr) {
            cout << "Skipping dataset No. " << dsNo << ": " << bedFile << "\n";
        }
        else {
            cout << "Read genotypes from " << bedFile << "\n\n";
            smpGeno->SaveHexGenotypes(outFile);
            totSmps += numSmps;
        }
    }

    return totSmps;
}


vector<string> FindSamplesInFamFile(const char* famFile, bool* hasErr)
{
    vector<string> samples;
    map<string, int> dupSmps;
    map<string, int>::iterator smpIt;

    *hasErr = false;
    int nameLen = 80;
    int lineLen = 300;
    char famLine[lineLen];
    int numFamSmps = 0;

    FILE *ifp = fopen(famFile, "r");
    if(ifp) {
        char fam[nameLen];
        char smp[nameLen];
        char dad[nameLen];
        char mom[nameLen];
        int sex;
        char pheno[16];
        while (fgets(famLine, lineLen, ifp) != NULL && *hasErr != true) {
            sscanf(famLine, "%s %s %s %s %d %d", fam, smp, dad, mom, &sex, pheno);
            if (smp[0]) {
                samples.push_back(smp);
                smpIt = dupSmps.find(smp);
                if (smpIt == dupSmps.end()) {
                    dupSmps[smp] = 1;
                }
                else {
                    cout << "ERROR in " << famFile << ": found duplicated sample " << smp << ".\n";
                    *hasErr = true;
                }
                numFamSmps++;
            }
        }
        fclose(ifp);
    }
    else {
        cout << "ERROR: couldn't open " << famFile << "!\n";
        *hasErr = true;
    }

    return samples;
}


string GetExecutablePath()
{
    char rawPathName[PATH_MAX];
    realpath(PROC_SELF_EXE, rawPathName);

    string exePath = string(rawPathName);
    size_t slashPos = exePath.find_last_of("/\\");
    string exeDir = exePath.substr(0, slashPos);

    return exeDir;
}
