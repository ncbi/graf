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

#include <fstream>
#include "Graf.h"

SampleInfo::SampleInfo()
{
    phs = 0;
    smpId = 0;
    sbjId = 0;
    strcpy(smpName, "");
    strcpy(sbjName, "");
    sex = 0;
}

SampleInfo::SampleInfo(char* smp, char* sbj, int gender, int sampleId, int subjectId, int std)
{
    phs = std;
    smpId = sampleId;
    sbjId = subjectId;
    strcpy(smpName, smp);
    strcpy(sbjName, sbj);
    sex = gender;
}

void SampleGeno::SetConstants()
{
    smpNameLen      = 128;
    numFpSnps       = 10000;
    numFpHexChars   = numFpSnps / 2;
    genoHexLen      = numFpHexChars + 2;
    genoFileLineLen = genoHexLen + smpNameLen * 2 + 30;
    numFpGenoInts   = (numFpSnps - 1) / 64 + 1;

    unsigned long base = 1;
    for (int i = 0; i < 64; i++) {
        baseNums[i] = base;
        base = base << 1;
    }
}

void SampleGeno::SetMaxSnpMiss(int minFpSnps)
{
    maxTotSnpMiss = 10000 - minFpSnps; /// Maximum SNPs with missing genotypes to be read from the geno file
}

SampleGeno::SampleGeno(int totSmps)
{
    maxTotSnpMiss = 9500;
    SetConstants();
    numSmps = totSmps;
    smpGenos1 = (unsigned long**)malloc(sizeof(unsigned long*) * numSmps);
    smpGenos2 = (unsigned long**)malloc(sizeof(unsigned long*) * numSmps);
    samples = (SampleInfo*)malloc(sizeof(SampleInfo) * numSmps);

    for (int i = 0; i < numSmps; i++) {
        smpGenos1[i] = NULL;
        smpGenos2[i] = NULL;
        samples[i].smpId = i;
        samples[i].sbjId = i;
    }

    snpInfo = new SnpInfo(numFpSnps);
    fpSnps = new FpSnps(numFpSnps);
    hgmrAdj = 1.0;
}

SampleGeno::SampleGeno(char* famFile, int filePhs)
{
    vector<string> smpNames;
    vector<int> sexes;

    maxTotSnpMiss = 9500;
    SetConstants();

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
        while (fgets(famLine, lineLen, ifp) != NULL) {
            sscanf(famLine, "%s %s %s %s %d %d", fam, smp, dad, mom, &sex, pheno);
            if (smp[0]) {
                smpNames.push_back(smp);
                sexes.push_back(sex);
                numFamSmps++;
            }
        }
        fclose(ifp);

        numSmps = numFamSmps;

        smpGenos1 = (unsigned long**)malloc(sizeof(unsigned long*) * numSmps);
        smpGenos2 = (unsigned long**)malloc(sizeof(unsigned long*) * numSmps);
        samples = (SampleInfo*)malloc(sizeof(SampleInfo) * numSmps);

        for (int i = 0; i < numSmps; i++) {
            smpGenos1[i] = NULL;
            smpGenos2[i] = NULL;

            int smpId = i + 1;
            samples[i].smpId = smpId;
            samples[i].sbjId = smpId;
            samples[i].phs = filePhs;
            samples[i].sex = sexes[i];

            strcpy(samples[i].smpName, smpNames[i].c_str());
            strcpy(samples[i].sbjName, smpNames[i].c_str());
        }

        snpInfo = new SnpInfo(numFpSnps);
        fpSnps = new FpSnps(numFpSnps);
        hgmrAdj = 1.0;
    }
}

SampleGeno::SampleGeno(char* genoFile)
{
    vector<string> smpNames;
    vector<int> stds;

    SetConstants();

    int nameLen = 80;
    int lineLen = 5120;
    char genoLine[lineLen];
    int numGenoSmps = 0;

    FILE *ifp = fopen(genoFile, "r");
    if(ifp) {
        int phs;
        char smp[nameLen];
        while (fgets(genoLine, lineLen, ifp) != NULL) {
            sscanf(genoLine, "%d %s", &phs, smp);
            if (smp[0]) {
                smpNames.push_back(smp);
                stds.push_back(phs);
                numGenoSmps++;
            }
        }
        fclose(ifp);

        numSmps = numGenoSmps;
        cout << "Found " << numSmps << " samples in " << genoFile << "\n";

        smpGenos1 = (unsigned long**)malloc(sizeof(unsigned long*) * numSmps);
        smpGenos2 = (unsigned long**)malloc(sizeof(unsigned long*) * numSmps);
        samples = (SampleInfo*)malloc(sizeof(SampleInfo) * numSmps);

        for (int i = 0; i < numSmps; i++) {
            smpGenos1[i] = NULL;
            smpGenos2[i] = NULL;

            int smpId = i + 1;
            samples[i].smpId = smpId;
            samples[i].sbjId = smpId;
            samples[i].phs = stds[i];
            samples[i].sex = 0;
            strcpy(samples[i].smpName, smpNames[i].c_str());
            strcpy(samples[i].sbjName, smpNames[i].c_str());
        }

        snpInfo = new SnpInfo(numFpSnps);
        fpSnps = new FpSnps(numFpSnps);
        hgmrAdj = 1.0;
    }
}

char SampleGeno::GetCompAllele(char a)
{
    char c = '0';
    if      (a == 'A') c = 'T';
    else if (a == 'T') c = 'A';
    else if (a == 'G') c = 'C';
    else if (a == 'C') c = 'G';

    return c;
}

bool SampleGeno::ReadGenotypesFromBedFile(char* bimFile, char* bedFile, int minFpSnps, map<int, int> ssToRs, bool rsIsInt)
{
    /// Check the bim file to find fingerprinting SNPs, and whether the alleles need to be flipped
    vector<int> fpSnpIdxs;    /// Array of SNP posision at bim file (0 to total num SNPs -1 in bim file)
    vector<int> fpSnpIds;     /// Array of SNP posision at FP database (0 to 9999)
    vector<int> snpFlips;     /// Whether alleles need to be flipped or not (AG => TC)
    vector<int> snpSwaps;     /// Whether alleles need to be swapped or not (TC => CT, GA => TC)
    vector<string> errLines;
    map<int, int> snpIdBimFp; /// FP SNPs appear in bim file not in order. We first read genos from the bed file and save them to
    /// memory using SNP order in bim file.  Then we need to reorder the SNPs using the FP SNP order.
    bool hasErr = false;
    int nameLen = 64;
    int lineLen = 512;
    char bimLine[lineLen];
    int numBimSnps = 0;
    int numBimFps = 0;
    int numFlips = 0;
    int numNoFlips = 0;
    int numMismats = 0;
    int numGenoInts = (numFpSnps - 1 ) / 64 + 1;

    /// In case of same SNPs in multiple rows, include one SNP only once
    int genoFpSnps[numFpSnps];
    for (int i = 0; i < numFpSnps; i++) genoFpSnps[i] = 0;

    map<int, int>::iterator ssRsIt;
    int numSsMaps = ssToRs.size();

    FILE *bimFilePtr = fopen(bimFile, "r");
    if(bimFilePtr) {
        int chr;
        char snp[nameLen];
        int junk;
        int pos;
        char a1Str[nameLen];
        char a2Str[nameLen];
        while (fgets(bimLine, lineLen, bimFilePtr) != NULL) {
            sscanf(bimLine, "%d %s %f %d %s %s", &chr, snp, &junk, &pos, a1Str, a2Str);

            int rs = 0;
            if (rsIsInt) {
                rs = atoi(snp);
            }
            else {
                int snpLen = strlen(snp);
                if (snpLen > 2) {
                    if ((snp[0] == 'r' || snp[0] == 'R') && (snp[1] == 's' || snp[1] == 'S')) {
                        rs = atoi(&snp[2]);
                    }
                    else if (numSsMaps && (snp[0] == 's' || snp[0] == 'S') && (snp[1] == 's' || snp[1] == 'S')) {
                        int ss = atoi(&snp[2]);
                        ssRsIt = ssToRs.find(ss);
                        if (ssRsIt != ssToRs.end()) rs = ssRsIt->second;
                    }
                }
            }

            if (rs) {
                int snpNo = fpSnps->FindRs(rs);
                if (snpNo > -1) {
                    int snpId = fpSnps->snps[snpNo].snpId;
                    char dbA1 = fpSnps->snps[snpNo].a1;
                    char dbA2 = fpSnps->snps[snpNo].a2;

                    char a1 = '0';
                    char a2 = '0';
                    if (a1Str[0] && !a1Str[1]) a1 = a1Str[0];
                    if (a2Str[0] && !a1Str[1]) a2 = a2Str[0];

                    if (a1 != '0' && a1 != '0') {
                        bool flip;
                        bool swap;
                        bool mismatch = false;
                        char c1 = GetCompAllele(a1);
                        char c2 = GetCompAllele(a2);

                        if      ((a1 == dbA1 && a2 == dbA2) || (a1 == dbA2 && a2 == dbA1) ) {
                            flip = false;
                            if (a1 == dbA1) swap = false;
                            else            swap = true;
                            numNoFlips++;
                        }
                        else if ((c1 == dbA1 && c2 == dbA2) || (c1 == dbA2 && c2 == dbA1) ) {
                            flip = true;
                            if (c1 == dbA1) swap = false;
                            else            swap = true;
                            numFlips++;
                        }
                        else {
                            char errLine[64];
                            sprintf(errLine, "rs%d:\t%c%c => %c%c", rs, dbA1, dbA2, a1, a2);
                            errLines.push_back((string)errLine);
                            mismatch = true;
                            numMismats++;
                        }

                        if (!mismatch && !genoFpSnps[snpNo]) {
                            fpSnpIdxs.push_back(numBimSnps);
                            fpSnpIds.push_back(snpId);
                            snpFlips.push_back(flip);
                            snpSwaps.push_back(swap);
                            snpIdBimFp[snpId] = numBimFps;
                            genoFpSnps[snpNo] = 1;
                            numBimFps++;
                        }
                    }
                }
            }

            numBimSnps++;
        }

        cout << "\nReading genotypes from PLINK set.\n";
        printf("\tTotal %d SNPs in the bim file. %d SNPs are FP SNPs.\n", numBimSnps, numBimFps);
        printf("\t%d SNPs need to be flipped.  %d SNPs don't need to be flipped.  %d SNPs have allele conflicts with dbSNP.\n",
        numFlips, numNoFlips, numMismats);

        if (numMismats > 0) {
            cout << "\nERROR: ";
            cout << numMismats << " SNPs have allele conflicts with dbSNP.\n";
            for (int i = 0; i < numMismats; i++) {
                cout << "\t" << errLines[i] << "\n";
                if (i > 10) {
                    cout << "\t...\n";
                    break;
                }
            }
            cout << "\n";
        }

        fclose(bimFilePtr);
    }
    else {
        cout << "ERROR: Couldn't open bim file " << bimFile << "!\n";
    }

    if (numBimFps < minFpSnps) {
        cout << "Not enough fingerprinting SNPs found in bim file!\n";
        return true;
    }

    /// Read genotypes from the bed file
    char header[2];
    char mode[1];

    ifstream bedFilePtr (bedFile, ios::in | ios::binary);

    long snpNumBytes = (numSmps - 1) / 4 + 1;
    long expFileLen = snpNumBytes * numBimSnps + 3;

    bedFilePtr.seekg (0, bedFilePtr.end);
    long fileLen = bedFilePtr.tellg();
    bedFilePtr.seekg(0, bedFilePtr.beg);

    bedFilePtr.read (header, 2);
    bedFilePtr.read (mode, 1);

    if (header[0] != 108 || header[1] != 27) {
        cout << "ERROR: File " << bedFile << " is not a valid PLINK bed file!\n";
        hasErr = true;
    }
    else if (mode[0] != 1) {
        cout << "ERROR: File " << bedFile << " is not in SNP mode!\n";
        hasErr = true;
    }
    else if (fileLen != expFileLen) {
        cout << "ERROR: Number of genotypes in bed file doesn't match fam and bim File!\n";
        cout << "Fam file has " << numSmps << " samples.  Bim file has " << numBimSnps
        << " SNPs. Each SNP should have " << snpNumBytes << " bytes.  Expected total " << expFileLen << " bytes.\n";
        cout << "Bed file has " << fileLen << " bytes.\n";
        hasErr = true;
    }

    if (hasErr) return hasErr;

    /// Read genotypes to memory first, in array of SNPs
    char* snpGenoStrs[numFpSnps];      /// Each time read genos for all samples at one SNP
    /// fpSnpIdxs[snpNo] is the SNP position in the whole bim file.
    char buff[snpNumBytes];            /// Reusable memory to keep the genotypes
    int bimFpSnpNo = 0;
    for (int i = 0; i < numBimSnps; i++) {
        bedFilePtr.read (buff, snpNumBytes);
        if (i == fpSnpIdxs[bimFpSnpNo]) {
            int snpId = fpSnpIds[bimFpSnpNo];
            char* snpGenoStr = (char*)malloc(sizeof(char) * snpNumBytes);
            for (int j = 0; j < snpNumBytes; j++) snpGenoStr[j] = buff[j];
            assert(bimFpSnpNo < numFpSnps);
            snpGenoStrs[bimFpSnpNo] = snpGenoStr;


            bimFpSnpNo++;
        }
    }

    bedFilePtr.close();

    cout << "\n\tRecoding genotypes of " << numSmps << " samples.\n";

    /// Now transpose the matrix, copy the genotypes to array of samples
    for (int smpNo = 0; smpNo < numSmps; smpNo++) {
        unsigned long* genoInts1 = (unsigned long*)malloc(sizeof(unsigned long) * numGenoInts);
        unsigned long* genoInts2 = (unsigned long*)malloc(sizeof(unsigned long) * numGenoInts);

        int smpBedByteNo   = smpNo / 4;          /// One bed byte stores geno for 4 samples: 3 2 1 0   7 6 5 4 ...
        int smpBedBytePos1 = (smpNo % 4) * 2;    /// Counting from the least siginificant digit
        int smpBedBytePos2 = smpBedBytePos1 + 1;

        /// Put genotypes of 64 SNPs (16 bed bytes) into one long integer.
        for (int genoIntNo = 0; genoIntNo < numGenoInts; genoIntNo++) {    /// 0 - 156
            unsigned long genoInt1 = 0;
            unsigned long genoInt2 = 0;

            for (int bitNo = 0; bitNo < 64; bitNo++) {
                int genoIntBitPos = 64 - bitNo - 1;            /// Position at the 64-bit integer, from most sig digit
                int snpId = genoIntNo * 64 + bitNo;
                bimFpSnpNo = -1;
                map<int, int>::iterator iter = snpIdBimFp.find(snpId);
                if (iter != snpIdBimFp.end()) {
                    bimFpSnpNo = iter->second;
                }

                /// Set genotype to all missing (85 = 01010101.  Read from right to left) when an FP SNPs is not found in the array
                char bedByte = 85;
                bool swap = false;
                if (bimFpSnpNo > -1) {
                    swap = snpSwaps[bimFpSnpNo];
                    bedByte = snpGenoStrs[bimFpSnpNo][smpBedByteNo];
                }

                int bit1 = bedByte & baseNums[smpBedBytePos1];
                int bit2 = bedByte & baseNums[smpBedBytePos2];

                if (bit1 && !bit2) {       /// Missing genotype
                    genoInt1 = genoInt1 | baseNums[genoIntBitPos];
                }
                else if (!bit1 && bit2) {  /// Hetrozygous
                    genoInt2 = genoInt2 | baseNums[genoIntBitPos];
                }
                else {
                    if (swap) {
                        if (bit1) {           /// 00 => 11
                            genoInt1 = genoInt1 | baseNums[genoIntBitPos];
                            genoInt2 = genoInt2 | baseNums[genoIntBitPos];
                        }
                    }
                    else {
                        if (!bit1) {            /// Set 11
                            genoInt1 = genoInt1 | baseNums[genoIntBitPos];
                            genoInt2 = genoInt2 | baseNums[genoIntBitPos];
                        }
                    }
                }

                int gapGeno = 3;
                if      ( bit1 &&  bit2) gapGeno = 2;
                else if (!bit1 &&  bit2) gapGeno = 1;
                else if (!bit1 && !bit2) gapGeno = 0;

                if (snpId < numFpSnps) snpInfo->AddSnpGeno(snpId, gapGeno);
            }

            genoInts1[genoIntNo] = genoInt1;
            genoInts2[genoIntNo] = genoInt2;
        }

        smpGenos1[smpNo] = genoInts1;
        smpGenos2[smpNo] = genoInts2;

        int smpCnt = smpNo + 1;
        if (smpCnt % 10000 == 0) cout << "\t\tRecoded genotypes for " << smpCnt << " samples.\n";
    }

    for (int i = 0; i < numBimFps; i++) delete snpGenoStrs[i];
    cout << "\n\tBed file has genotypes of " << numBimSnps << " SNPs. Read genotypes of "
    << numBimFps << " fingerprinting SNPs for " << numSmps << " samples.\n";

    snpInfo->Summarize();


    return hasErr;
}


void SampleGeno::ShowSamples(int numShowSmps)
{
    cout << "Showing " << numShowSmps << " samples total " << numSmps << " samples.\n";
    if (numShowSmps > numSmps) numShowSmps = numSmps;

    for (int i = 0; i < numShowSmps; i++) {
        printf("No. %5d: Sample ID %5d name %s Subject ID %5d name %s sex %d phs %3d\n",
        i, samples[i].smpId, samples[i].smpName, samples[i].sbjId, samples[i].sbjName, samples[i].sex, samples[i].phs);
    }
    cout << "Total " << numSmps << " samples.\n";
}

SampleGeno::~SampleGeno()
{
    DeleteSamples();
}

void SampleGeno::DeleteSamples()
{
    for (int i = 0; i < numSmps; i++) {
        if (smpGenos1[i]) delete smpGenos1[i];
        if (smpGenos2[i]) delete smpGenos2[i];
    }
    delete smpGenos1;
    delete smpGenos2;
    delete samples;
    delete snpInfo;
    delete fpSnps;

    numSmps = 0;
}

int SampleGeno::ReadGenotypesFromHexGenoFile(char* genoFile, bool hasId, int* stPos, int* edPos)
{
    int i, j, k;

    /// Summary of SNP genotypes: counts of AA, AB, BB for each SNP
    printf("Reading genotypes from file %s.\n", genoFile);

    FILE *ifp = fopen(genoFile, "r");
    if(!ifp) {
        cout << "Could not open genotype file " << genoFile << "!\n";
        DeleteSamples();
        return 0;
    }

    char genoLine[genoFileLineLen];
    char genoStr[genoHexLen];
    int lineNo;

    for (i = 0; i < numSmps; i++) {
        char* smpName = (char*)malloc(smpNameLen * sizeof(char));
        char* sbjName = (char*)malloc(smpNameLen * sizeof(char));
        samples[i].phs = 0;
        samples[i].smpId = 0;
        samples[i].sbjId = 0;
    }

    /// Initilize the genotype arrays
    for (i = 0; i < numSmps; i++) {
        unsigned long* genoInts1;
        unsigned long* genoInts2;

        genoInts1 = (unsigned long*)malloc(numFpGenoInts * sizeof(unsigned long));
        genoInts2 = (unsigned long*)malloc(numFpGenoInts * sizeof(unsigned long));
        assert(genoInts1 != NULL);
        assert(genoInts2 != NULL);
        for (j = 0; j < numFpGenoInts; j++) {
            genoInts1[j] = 0;
            genoInts2[j] = 0;
        }

        smpGenos1[i] = genoInts1;
        smpGenos2[i] = genoInts2;
    }

    lineNo = 0;
    int genoSmpNo = 0;
    int iSmpId, iSbjId, iPhs, iSex;
    char iSmpName[smpNameLen], iSbjName[smpNameLen];
    int numMissSmps = 0;

    while (fgets(genoLine, genoFileLineLen, ifp) != NULL) {
        if (hasId) {
            sscanf(genoLine, "%d %d %s %d %s %d %s", &iPhs, &iSmpId, iSmpName, &iSbjId, iSbjName, &iSex, genoStr);
        }
        else {
            sscanf(genoLine, "%d %s %s", &iPhs, iSmpName, genoStr);
        }

        int numMiss = 0;
        int bitSnpPos = 0;
        for (i = 0; i < numFpGenoInts; i++) {
            /// Each pair 64-bit integers stores genotypes of 64 SNPs
            unsigned long gInt1 = 0;
            unsigned long gInt2 = 0;

            /// Each hex letter contains genotypes of two SNPs.  Read 32 hex letters and build a 64-bit long integer
            for (j = 0; j < 32; j++) {
                /// Starting from the most significant digit, with 1 representing 2^63
                int baseNo = 64 - j*2 - 1;
                unsigned long base1 = baseNums[baseNo];
                unsigned long base2 = baseNums[baseNo-1];

                int hexNo = i * 32 + j; /// The index of the hex char in the sequence
                char ch = 'F';          /// Last integer doesn't have 64 SNPs.  Append with NULLs.
                if (hexNo < numFpHexChars) ch = genoStr[hexNo];
                int snpNo = hexNo * 2;  /// Index of the first SNP of the two encoded in the hex letter

                /// Split each hex into two integers, with 0 = BB; 1 = AB; 2 = AA; 3 = missing
                int val = ch - '0';
                if (val < 0 || val > 9) {
                    val = 10 + ch - 'A';
                }

                assert (val > -1 && val < 16);
                short g1 = val / 4;
                short g2 = val % 4;

                if (snpNo < numFpSnps) {
                    snpInfo->AddSnpGeno(snpNo, g1);
                    snpInfo->AddSnpGeno(snpNo+1, g2);
                }

                bitSnpPos++;
                if      (g1 == 1) { /// AB => 01
                    gInt2 += base1;
                }
                else if (g1 == 2) { /// BB => 11
                    gInt1 += base1;
                    gInt2 += base1;
                }
                else if (g1 == 3) { /// XX => 10
                    gInt1 += base1;
                    if (bitSnpPos < 10000) numMiss++;
                }

                bitSnpPos++;
                if      (g2 == 1) { /// AB => 01
                    gInt2 += base2;
                }
                else if (g2 == 2) { /// BB => 11
                    gInt1 += base2;
                    gInt2 += base2;
                }
                else if (g2 == 3) { /// XX => 10
                    gInt1 += base2;
                    if (bitSnpPos < 10000) numMiss++;
                }
            }

            smpGenos1[genoSmpNo][i] = gInt1;
            smpGenos2[genoSmpNo][i] = gInt2;
        }

        lineNo++;

        if (numMiss < maxTotSnpMiss) {
            samples[genoSmpNo].phs = iPhs;
            strcpy(samples[genoSmpNo].smpName, iSmpName);
            if (hasId) {
                samples[genoSmpNo].smpId = iSmpId;
                samples[genoSmpNo].sbjId = iSbjId;
                samples[genoSmpNo].sex = iSex;
                strcpy(samples[genoSmpNo].sbjName, iSbjName);
            }
            else {
                samples[genoSmpNo].smpId = i + 1;
                samples[genoSmpNo].sbjId = i + 1;
                samples[genoSmpNo].sex = 0;
                strcpy(samples[genoSmpNo].sbjName, iSmpName);
            }

            genoSmpNo++;
        }
        else {
            numMissSmps++;
        }

        if (lineNo % 50000 == 0) cout << "\tRead genotypes for " << genoSmpNo << " samples from  " << lineNo << " lines\n";
    }
    fclose(ifp);

    cout << "\tRead genotypes for total " << genoSmpNo << " samples from  " << lineNo << " lines.\n";
    if (numMissSmps > 0) cout << "\tSkipped " << numMissSmps << " samples with low genotyping rates.\n";

    snpInfo->Summarize();


    int actSmps = genoSmpNo;
    numSmps = actSmps;

    return actSmps;
}

// Reorder SNPs by MAF and missing rate so that the more informative SNPs are in the front
void SampleGeno::ReorderSnps()
{
    int i, j, k;

    /// Compute expected HGMRs for each relative group
    double unHms = 0.0, unHmt = 0.0;
    double d3Hms = 0.0, d3Hmt = 0.0;
    double d2Hms = 0.0, d2Hmt = 0.0;
    double fsHms = 0.0, fsHmt = 0.0;

    numValidSnps = 0;

    SnpMaf snpMafs[numFpSnps];

    int mafStSmpNo = 0;
    int mafEdSmpNo = 0;

    double p2q2 = 0;
    double p4q4 = 0;
    double p3q3 = 0;
    double p2aq2 = 0;
    double p2q = 0;
    double pq2 = 0;

    for (i = 0; i < numFpSnps; i++) {
        int numAlleles =  snpInfo->aaCnts[i] + snpInfo->abCnts[i] + snpInfo->bbCnts[i];
        double maf = snpInfo->mafs[i];
        double cr = 0, mafCr = 0;

        if (numAlleles > 0) {
            double p = maf;
            double q = 1 - p;
            cr = snpInfo->crs[i];
            mafCr = snpInfo->mafCrs[i];

            p2q2  += 2 * p*p * q*q;
            p4q4  += p*p*p*p + q*q*q*q;;
            p3q3  += p*p*p   + q*q*q;
            p2aq2 += p*p     + q*q;

            p2q += p*p * q;
            pq2 += q*q * p;

            double hms = 2 * p*p * q*q;
            double p4q4 = p*p*p*p + q*q*q*q;
            double p3q3 = p*p*p   + q*q*q;
            double p2q2 = p*p     + q*q;

            unHms += hms;
            unHmt += p4q4;

            d3Hms += hms*3/4;
            d3Hmt += p4q4*3/4 + p3q3/4;

            d2Hms += hms/2;
            d2Hmt += p4q4/2 + p3q3/2;

            fsHms += hms/4;
            fsHmt += p4q4/4 + p3q3/2 + p2q2/4;

            numValidSnps++;
        }
        snpMafs[i].snpNo = i;
        snpMafs[i].maf = maf;
        snpMafs[i].cr = cr;
        snpMafs[i].mafCr = mafCr;
    }

    /// Sort SNPs by Call rate
    qsort(snpMafs, numFpSnps, sizeof(struct SnpMaf), CompareSnpMaf);

    SnpInfo* tmpSnpInfo = snpInfo;
    snpInfo = new SnpInfo(numFpSnps);
    for (i = 0; i < numFpSnps; i++) {
        int snpNo = snpMafs[i].snpNo;
        snpInfo->aaCnts[i]  = tmpSnpInfo->aaCnts[snpNo];
        snpInfo->abCnts[i]  = tmpSnpInfo->abCnts[snpNo];
        snpInfo->bbCnts[i]  = tmpSnpInfo->bbCnts[snpNo];
        snpInfo->ooCnts[i]  = tmpSnpInfo->ooCnts[snpNo];
        snpInfo->refCnts[i] = tmpSnpInfo->refCnts[snpNo];
        snpInfo->altCnts[i] = tmpSnpInfo->altCnts[snpNo];
        snpInfo->mafs[i]    = tmpSnpInfo->mafs[snpNo];
        snpInfo->crs[i]     = tmpSnpInfo->crs[snpNo];
        snpInfo->mafCrs[i]  = tmpSnpInfo->mafCrs[snpNo];
    }


    delete tmpSnpInfo;
    cout << "\tSorted genotypes by MAF and call rate.\n";

    /// Reorder the genotypes, using the SNPs ordered by MAF
    unsigned long** smpOdGenos1;
    unsigned long** smpOdGenos2;
    smpOdGenos1 = (unsigned long**)malloc(sizeof(unsigned long*) * numSmps);
    smpOdGenos2 = (unsigned long**)malloc(sizeof(unsigned long*) * numSmps);
    for (i = 0; i < numSmps; i++) {
        unsigned long* genoInts1;
        unsigned long* genoInts2;

        genoInts1 = (unsigned long*)malloc(numFpGenoInts * sizeof(unsigned long));
        genoInts2 = (unsigned long*)malloc(numFpGenoInts * sizeof(unsigned long));
        assert(genoInts1 != NULL);
        assert(genoInts2 != NULL);

        for (j = 0; j < numFpGenoInts; j++) {
            genoInts1[j] = 0;
            genoInts2[j] = 0;

            int stSnp = j * 64;
            /// Pick up 64 bits from different long int numbers
            for (k = 0; k < 64; k++) {
                int snpIdx = stSnp + k;
                int bit1 = 0;
                int bit2 = 0;
                if (snpIdx < numFpSnps) {
                    int pickSnp = snpMafs[snpIdx].snpNo;
                    double pickMaf = snpMafs[snpIdx].maf;

                    /// Pick the oldBitNum's bit of the oldIntNum's geno integer
                    int oldIntNum = pickSnp / 64;
                    int oldBitNum = pickSnp % 64;

                    if (smpGenos1[i][oldIntNum] & baseNums[64-oldBitNum-1]) bit1 = 1;
                    if (smpGenos2[i][oldIntNum] & baseNums[64-oldBitNum-1]) bit2 = 1;
                }
                else {
                    /// Missing genotype: 10
                    bit1 = 1;
                    bit2 = 0;
                }

                if (bit1) genoInts1[j] = genoInts1[j] | baseNums[64-k-1];
                if (bit2) genoInts2[j] = genoInts2[j] | baseNums[64-k-1];
            }
        }

        smpOdGenos1[i] = genoInts1;
        smpOdGenos2[i] = genoInts2;
    }

    /// hgmrAdj is out of date.
    double unExpHgmr = p2q2 / (p2q2 + p4q4);
    double avgUnMr = 0.2482;
    hgmrAdj = avgUnMr / unExpHgmr;

    /// Set expected mismatch rates for each relationship
    double pz0 = (2 * p2q2 + p4q4) / numValidSnps;
    double pz1 = (p3q3 + p2q + pq2) / numValidSnps;
    double pz2 = 1;

    poExpAgmr = 100.0 * (1 - pz1);
    fsExpAgmr = 100.0 * (1 - (0.25 * pz0 + 0.50 * pz1 + 0.25 * pz2));
    d2ExpAgmr = 100.0 * (1 - (0.50 * pz0 + 0.50 * pz1));
    d3ExpAgmr = 100.0 * (1 - (0.75 * pz0 + 0.25 * pz1));

    fsExpHgmr = 100.0 * 0.25 * p2q2 / (0.25 * (p2q2 + p4q4) + 0.50 * p3q3 + 0.25 * p2aq2);
    d2ExpHgmr = 100.0 * 0.50 * p2q2 / (0.50 * (p2q2 + p4q4) + 0.50 * p3q3);
    d3ExpHgmr = 100.0 * 0.75 * p2q2 / (0.75 * (p2q2 + p4q4) + 0.25 * p3q3);

    for (i = 0; i < numSmps; i++) {
        delete smpGenos1[i];
        delete smpGenos2[i];
        smpGenos1[i] = smpOdGenos1[i];
        smpGenos2[i] = smpOdGenos2[i];
    }

    delete smpOdGenos1;
    delete smpOdGenos2;
}

int SampleGeno::GetNumberStudies()
{
    map<int, int> allStds;
    for (int i = 0; i < numSmps; i++) {
        int phs = samples[i].phs;
        allStds[phs] = 1;
    }

    int numStds = allStds.size();

    return numStds;
}

bool SampleGeno::SaveHexGenotypes(char* outFile)
{
    bool saved = true;

    FILE *ifp = fopen(outFile, "a");
    if(ifp) {
        int i, j, k;
        for (i = 0; i < numSmps; i++) {
            char* smp = samples[i].smpName;
            int phs = samples[i].phs;
            char line[numFpHexChars+1];

            for (j = 0; j < numFpGenoInts; j++) {
                unsigned long genoInt1 = smpGenos1[i][j];
                unsigned long genoInt2 = smpGenos2[i][j];







                for (k = 0; k < 32; k++) {
                    int charNo = j * 32 + k;

                    int snpNo1 = 64 - k*2 - 1;
                    int snpNo2 = snpNo1 - 1;;
                    int genoIntBitPos = 64 - k - 1;

                    bool snp1b1 = genoInt1 & baseNums[snpNo1];
                    bool snp1b2 = genoInt2 & baseNums[snpNo1];
                    int g1 = 3;
                    if (!snp1b1 || snp1b2) {
                        g1 = snp1b1 + snp1b2;
                    }
                    //	  int g1 = GetSnpGenoInt(snp1b1, snp1b2);


                    bool snp2b1 = genoInt1 & baseNums[snpNo2];
                    bool snp2b2 = genoInt2 & baseNums[snpNo2];
                    int g2 = 3;
                    if (!snp2b1 || snp2b2) {
                        g2 = snp2b1 + snp2b2;
                    }



                    int gNum = g1 * 4 + g2;
                    char gHex = '0' + gNum;
                    if (gNum > 9) gHex = 'A' + gNum - 10;

                    if (charNo < numFpHexChars) line[charNo] = gHex;

                }

            }

            line[numFpHexChars] = 0;

            fprintf(ifp, "%d\t%s\t%s\n", phs, smp, line);
        }
    }
    else {
        cout << "ERROR: Can't open " << outFile << " for writing!\n";
        return false;
    }

    fclose(ifp);

    cout << "Saved genotypes of " << numSmps << " samples to " << outFile << ".\n";
    return true;
}

bool SampleGeno::SavePopulationResults(char* outFile)
{
    bool saved = true;
    int numPops = fpSnps->numChkPops;

    GenoDist ecAll, fcAll, acAll, spTmp;

    fcAll.f = fpSnps->popExpPfvals[0];
    fcAll.e = fpSnps->popExpPevals[0];
    fcAll.a = fpSnps->popExpPavals[0];

    ecAll.f = fpSnps->popExpPfvals[1];
    ecAll.e = fpSnps->popExpPevals[1];
    ecAll.a = fpSnps->popExpPavals[1];

    acAll.f = fpSnps->popExpPfvals[2];
    acAll.e = fpSnps->popExpPevals[2];
    acAll.a = fpSnps->popExpPavals[2];

    /// Dump sample point for TransformAllDists. Only the three vertices need to be transformed
    spTmp.f = fcAll.f;
    spTmp.e = fcAll.e;
    spTmp.a = fcAll.a;

    SampleGenoDist allGd(&ecAll, &fcAll, &acAll, &spTmp); /// The three vertices assuming all SNPs have genotypes
    allGd.TransformAllDists();

    printf("Positions (x, y, z coordinates) of the three vertices:\n");
    allGd.ShowPositions(0);

    FILE *ifp = fopen(outFile, "w");
    if(ifp) {
        int i, j;

        char header[256];
        sprintf(header, "%s\t%s\t%s\tGD1\tGD2\tGD3\tGD4\tF(\%)\tE(\%)\tA(\%)", "DS No.", "Sample", "#SNPs");
        for (j = 0; j < numPops-2; j++) {  /// Don't save Chinese and Japanese values for now
            sprintf(header, "%s\t%s", header, fpSnps->checkPops[j].c_str());
        }
        fprintf(ifp, "%s\n", header);

        for (i = 0; i < numSmps; i++) {
            char* smp = samples[i].smpName;
            int phs = samples[i].phs;
            char line[256];

            sprintf(line, "%d\t%s\t%d", phs, smp, samples[i].popNumSnps);

            /// Use barycentric coordinates to map the first three D values (to AA, E, EA) to the position
            /// expected for the case that all SNPs have genotypes

            GenoDist ecSmp, fcSmp, acSmp, spSmp;

            ecSmp.f = samples[i].popExpPfvals[1];
            ecSmp.e = samples[i].popExpPevals[1];
            ecSmp.a = samples[i].popExpPavals[1];

            fcSmp.f = samples[i].popExpPfvals[0];
            fcSmp.e = samples[i].popExpPevals[0];
            fcSmp.a = samples[i].popExpPavals[0];

            acSmp.f = samples[i].popExpPfvals[2];
            acSmp.e = samples[i].popExpPevals[2];
            acSmp.a = samples[i].popExpPavals[2];

            spSmp.f = samples[i].popMeanPvals[0];
            spSmp.e = samples[i].popMeanPvals[1];
            spSmp.a = samples[i].popMeanPvals[2];

            SampleGenoDist smpGd(&ecSmp, &fcSmp, &acSmp, &spSmp);
            smpGd.TransformAllDists();
            smpGd.CalculateBaryCenters();

            /// Show rotated x, y, z values as GD1, GD2, GD3
            double xVal = smpGd.eWt * allGd.ePt.x + smpGd.fWt * allGd.fPt.x + smpGd.aWt * allGd.aPt.x;
            double yVal = smpGd.eWt * allGd.ePt.y + smpGd.fWt * allGd.fPt.y + smpGd.aWt * allGd.aPt.y;
            double zVal = smpGd.sPt.z;

            /// GD4 = D_mexican - D_india_pakistani
            double gd4 = samples[i].popMeanPvals[3] - samples[i].popMeanPvals[4];

            double fjWt = smpGd.fWt > 0 ? smpGd.fWt : 0;
            double ejWt = smpGd.eWt > 0 ? smpGd.eWt : 0;
            double ajWt = smpGd.aWt > 0 ? smpGd.aWt : 0;
            double totWt = fjWt + ejWt + ajWt;
            double fPct = fjWt * 100 / totWt;
            double ePct = ejWt * 100 / totWt;
            double aPct = ajWt * 100 / totWt;

            sprintf(line, "%s\t%7.6f\t%7.6f\t%7.6f\t%7.6f\t%6.2f\t%6.2f\t%6.2f", line, xVal, yVal, zVal, gd4, fPct, ePct, aPct);

            for (j = 0; j < numPops-2; j++) { /// Don't save Chinese and Japanese values for now
                sprintf(line, "%s\t%7.6f", line, samples[i].popMeanPvals[j]);
            }

            fprintf(ifp, "%s\n", line);
        }
    }
    else {
        cout << "ERROR: Can't open " << outFile << " for writing!\n";
        return false;
    }

    fclose(ifp);

    cout << "Saved population results of " << numSmps << " samples to " << outFile << ".\n";
    return true;
}

void SampleGeno::InitPopPvalues()
{
    int numPops = fpSnps->numChkPops;

    for (int smpNo = 0; smpNo < numSmps; smpNo++) {
        samples[smpNo].popPvalues   = (double*)malloc(sizeof(double)*numPops);
        samples[smpNo].popMeanPvals = (double*)malloc(sizeof(double)*numPops);
        samples[smpNo].popExpPfsums = (double*)malloc(sizeof(double)*3);
        samples[smpNo].popExpPesums = (double*)malloc(sizeof(double)*3);
        samples[smpNo].popExpPasums = (double*)malloc(sizeof(double)*3);
        samples[smpNo].popExpPfvals = (double*)malloc(sizeof(double)*3);
        samples[smpNo].popExpPevals = (double*)malloc(sizeof(double)*3);
        samples[smpNo].popExpPavals = (double*)malloc(sizeof(double)*3);
        samples[smpNo].popNumSnps = 0;
    }
}

void SampleGeno::SetPopulationPvalues(int thNo)
{
    int meanThSmps = int(numSmps / numThreads);
    int rmSmps = numSmps % numThreads;
    int chkThSmps = meanThSmps;
    if (thNo+1 <= rmSmps) chkThSmps = meanThSmps + 1;

    int stSmp = thNo * chkThSmps;
    if (thNo >= rmSmps) stSmp = thNo * chkThSmps + rmSmps;
    int edSmp = stSmp + chkThSmps - 1;


    int j, k;
    int smpCnt = 0;

    for (int smpNo = stSmp; smpNo <= edSmp; smpNo++) {
        int snpNo = 0;

        /// Read sample genotypes (0, 1, 2)
        bool bitGeno1[10000];
        bool bitGeno2[10000];
        if (smpNo >= 0 && smpNo < numSmps) {
            for (j = 0; j < numFpGenoInts; j++) {
                unsigned long genoInt1 = smpGenos1[smpNo][j];
                unsigned long genoInt2 = smpGenos2[smpNo][j];

                for (k = 0; k < 32; k++) {
                    int snpNo1 = 64 - k*2 - 1;
                    int snpNo2 = snpNo1 - 1;;
                    int genoIntBitPos = 64 - k - 1;

                    if (snpNo < 10000) {
                        bitGeno1[snpNo] = genoInt1 & baseNums[snpNo1];
                        bitGeno2[snpNo] = genoInt2 & baseNums[snpNo1];
                        snpNo++;
                    }

                    if (snpNo < 10000) {
                        bitGeno1[snpNo] = genoInt1 & baseNums[snpNo2];
                        bitGeno2[snpNo] = genoInt2 & baseNums[snpNo2];
                        snpNo++;
                    }
                }
            }
        }

        int numPops = fpSnps->numChkPops;

        for (int popId = 0; popId < numPops; popId++) {
            samples[smpNo].popPvalues[popId] = 0;   /// The raw log p-value
            samples[smpNo].popMeanPvals[popId] = 0; /// Genetic distancs from the sample to each ref population
        }

        /// Initialize expected values for each vertex population (E, F, A)
        for (int vetId = 0; vetId < 3; vetId++) {
            samples[smpNo].popExpPfsums[vetId] = 0;
            samples[smpNo].popExpPesums[vetId] = 0;
            samples[smpNo].popExpPasums[vetId] = 0;
            samples[smpNo].popExpPfvals[vetId] = 0;
            samples[smpNo].popExpPevals[vetId] = 0;
            samples[smpNo].popExpPavals[vetId] = 0;
        }

        int refPopSnps[7]; /// Counts of SNPs with freqs for each ref population
        for (int popId = 0; popId < numPops; popId++) {
            refPopSnps[popId] = 0;
        }

        int vetNumSnps[3]; /// Counts of SNPs with freqs for each ref population
        for (int vetId = 0; vetId < 3; vetId++) {
            vetNumSnps[vetId] = 0;
        }

        for (snpNo = 0; snpNo < 10000; snpNo++) {
            bool gb1 = bitGeno1[snpNo];
            bool gb2 = bitGeno2[snpNo];
            bool hasRefPv = false;
            bool hasVetPv = false;

            double v0p = fpSnps->snps[snpNo].vertexAf[0];
            double v1p = fpSnps->snps[snpNo].vertexAf[1];
            double v2p = fpSnps->snps[snpNo].vertexAf[2];

            if (!gb1 || gb2) {
                for (int popId = 0; popId < numPops; popId++) {
                    /// For the first three ref pops, since the barycentric coordinates need to be calculated for each sample,
                    /// only the SNPs on which all of the three vertex populations have genotypes are counted
                    double pv = fpSnps->snps[snpNo].popAf[popId]; /// reference population allele freq p value

                    if (pv > 0 && pv < 1 && (popId > 3 || (v0p > 0 && v0p < 1 && v1p > 0 && v1p < 1 && v2p > 0 && v2p < 1)) ) {
                        double qv = 1 - pv;
                        double bbPv = log(qv) * 2;
                        double abPv = log(pv * qv * 2);
                        double aaPv = log(pv) * 2;

                        if      (!gb1 && !gb2) {
                            samples[smpNo].popPvalues[popId] += bbPv;
                        }
                        else if (!gb1 &&  gb2) {
                            samples[smpNo].popPvalues[popId] += abPv;
                        }
                        else if ( gb1 &&  gb2) {
                            samples[smpNo].popPvalues[popId] += aaPv;
                        }

                        hasRefPv = true;
                        refPopSnps[popId]++;
                    }
                }

                /// Calculate expected P values for each of the vertices E, F, A
                for (int vetId = 0; vetId < 3; vetId++) {
                    double vpv = fpSnps->snps[snpNo].vertexAf[vetId];  /// Vertex population allele freq p value
                    double pfv = fpSnps->snps[snpNo].popAf[0];         /// Ref population 0 (Black) allele freq
                    double pev = fpSnps->snps[snpNo].popAf[1];         /// Ref population 1 (White) allele freq
                    double pav = fpSnps->snps[snpNo].popAf[2];         /// Ref population 2 (Asian) allele freq

                    /// Count only the SNPs on which all vertices have allele freqs, assuming ref pops always have freqs
                    if (v0p > 0 && v0p < 1 && v1p > 0 && v1p < 1 && v2p > 0 && v2p < 1) {
                        double vqv = 1 - vpv;
                        double qfv = 1 - pfv;
                        double qev = 1 - pev;
                        double qav = 1 - pav;

                        double aaPfv = log(pfv) * 2;
                        double bbPfv = log(qfv) * 2;
                        double abPfv = log(pfv * qfv * 2);

                        double aaPev = log(pev) * 2;
                        double bbPev = log(qev) * 2;
                        double abPev = log(pev * qev * 2);

                        double aaPav = log(pav) * 2;
                        double bbPav = log(qav) * 2;
                        double abPav = log(pav * qav * 2);

                        samples[smpNo].popExpPfsums[vetId] += aaPfv * vpv * vpv + bbPfv * vqv * vqv + abPfv * 2 * vpv * vqv;
                        samples[smpNo].popExpPesums[vetId] += aaPev * vpv * vpv + bbPev * vqv * vqv + abPev * 2 * vpv * vqv;
                        samples[smpNo].popExpPasums[vetId] += aaPav * vpv * vpv + bbPav * vqv * vqv + abPav * 2 * vpv * vqv;

                        hasVetPv = true;
                        vetNumSnps[vetId]++;
                    }
                }
            }

            if (hasRefPv) samples[smpNo].popNumSnps++;
        }

        for (int popId = 0; popId < numPops; popId++) {
            if (refPopSnps[popId] > 0) {
                samples[smpNo].popMeanPvals[popId] = -1 * samples[smpNo].popPvalues[popId]/refPopSnps[popId];
            }
        }

        for (int vetId = 0; vetId < 3; vetId++) {
            if (vetNumSnps[vetId] > 0) {
                samples[smpNo].popExpPfvals[vetId]  = -1 * samples[smpNo].popExpPfsums[vetId]/vetNumSnps[vetId];
                samples[smpNo].popExpPevals[vetId]  = -1 * samples[smpNo].popExpPesums[vetId]/vetNumSnps[vetId];
                samples[smpNo].popExpPavals[vetId]  = -1 * samples[smpNo].popExpPasums[vetId]/vetNumSnps[vetId];
            }
        }

        smpCnt++;

    }
}

int SampleGeno::GetSnpGenoInt(bool b1, bool b2)
{
    int g = 3;

    if (!b1 && !b2) {
        g = 0;
    }
    else if (!b1 && b2) {
        g = 1;
    }
    else if (b1 && b2) {
        g = 2;
    }

    return g;
}

void SampleGeno::SetSampleSubjectIds (map<string, string> ssm)
{
    map<string, string>::iterator sbjIt;

    for (int i = 0; i < numSmps; i++) {
        string smp = samples[i].smpName;
        int sex = samples[i].sex;
        sbjIt = ssm.find(smp);
        if (sbjIt != ssm.end()) {
            string sbj = ssm[smp];
            strcpy(samples[i].sbjName, sbj.c_str());

        }
    }
}

int SampleGeno::SetSampleGenders (map<string, int> sbjSexes)
{
    int numUpds = 0;

    map<string, string>::iterator sbjIt;
    map<string, int>::iterator sexIt;

    for (int i = 0; i < numSmps; i++) {
        string smp = samples[i].smpName;
        string sbj = samples[i].sbjName;
        int sex = samples[i].sex;
        sexIt = sbjSexes.find(sbj);
        if (sexIt != sbjSexes.end()) {
            int newSex = sbjSexes[sbj];
            if (newSex && newSex != sex) {
                samples[i].sex = newSex;

                numUpds++;
            }
        }
    }

    return numUpds;
}

int CompareSnpMaf(const void *a, const void *b)
{
    struct SnpMaf *ia = (struct SnpMaf *)a;
    struct SnpMaf *ib = (struct SnpMaf *)b;

    return (int)(10000.0*ib->mafCr - 10000.0*ia->mafCr);
}

