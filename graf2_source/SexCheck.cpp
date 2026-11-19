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
#include "FamSubjects.h"
#include "SexCheck.h"
#include <bitset>

using namespace std;


SexCheck::SexCheck(string file, XchrFpSnps* fpSnps, BimFileSnps* plinkSnps, FamSubjects* sbjs)
{
    ASSERT(sbjs,  "sbjs is NULL!");

    numBimSnps = 0;
    numSubjects = 0;
    numConflicts = 0;
    numMatches = 0;
    numErrSbjs = 0;
    numOtherSbjs = 0;

    minFpSnps = 16;
    maxMaleHetroRate = 0.08;
    minFemaleHetroRate = 0.2;
    maxFemaleHetroRate = 0.6;

    bedFile = file;
    sexFpSnps = fpSnps;
    bimSnps = plinkSnps;
    famSbjs = sbjs;
    numSubjects = sbjs->subjects.size();

    unsigned long base = 1;
    for (int i = 0; i < 64; i++) {
        baseNums[i] = base;
        base = base << 1;
    }
}

SexCheck::~SexCheck()
{
    /// Nothing to delete for now.
}

void SexCheck::SetCutoffValues(int minSnps, float maxMaleHr, float minFemaleHr, float maxFemaleHr)
{
    minFpSnps = minSnps;
    maxMaleHetroRate = maxMaleHr;
    minFemaleHetroRate = minFemaleHr;
    maxFemaleHetroRate = maxFemaleHr;
}


bool SexCheck::CalculateHetroRates()
{
    bool hasErr = false;

    /// Read genotypes from the bed file
    char header[2];
    char mode[1];

    ifstream bedFilePtr (bedFile.c_str(), ios::in | ios::binary);

    int numBimSnps = bimSnps->GetNumSnps();
    int numFpSnps = sexFpSnps->GetNumXchrFpSnps();

    long snpNumBytes = (numSubjects - 1) / 4 + 1;
    long expFileLen = snpNumBytes * numBimSnps + 3;

    bedFilePtr.seekg (0, bedFilePtr.end);
    long fileLen = bedFilePtr.tellg();
    bedFilePtr.seekg(0, bedFilePtr.beg);

    bedFilePtr.read (header, 2);
    bedFilePtr.read (mode, 1);

    if (header[0] != BYTE1_IN_BED_FILE || header[1] != BYTE2_IN_BED_FILE) {
        cout << "ERROR: File " << bedFile << " is not a valid PLINK bed file!\n";
        hasErr = true;
    }
    else if (mode[0] != BYTE_OF_SNP_MODE) {
        cout << "ERROR: File " << bedFile << " is not in SNP mode!\n";
        hasErr = true;
    }
    else if (fileLen != expFileLen) {
        cout << "ERROR: Number of genotypes in bed file don't match fam and bim File!\n";
        cout << "\tFam file has " << numSubjects << " samples.  Bim file has " << numBimSnps
        << " SNPs. Each SNP should have " << snpNumBytes << " bytes.  Expected total " << expFileLen << " bytes.\n";
        cout << "\tBed file has " << fileLen << " bytes.\n\n";
        hasErr = true;
    }

    if (hasErr) return hasErr;

    vector<int> sbjHomoCnts;
    vector<int> sbjHetroCnts;
    vector<int> sbjOtherCnts;

    for (int i = 0; i < numSubjects; i++) {
        sbjHomoCnts.push_back(0);
        sbjHetroCnts.push_back(0);
        sbjOtherCnts.push_back(0);
    }

    char buff[snpNumBytes];            /// Reusable memory to keep the genotypes
    int bimFpSnpNo = 0;

    /// For each sample, read genotypes of all SNPs
    for (int i = 0; i < numBimSnps; i++) {
        bedFilePtr.read (buff, snpNumBytes);
        int fpSnpId = bimSnps->GetFpSnpIdGivenBimSnpPos(i);

        /// Check genotypes if SNP is an X-chr FP SNP
        if (fpSnpId >= 0) {
            int smpNo = 0;

            /// Check each byte, which stores genotypes of 4 samples
            for (int j = 0; j < snpNumBytes; j++) {
                char bedByte = buff[j];


                /// Check genotypes of the for samples stored in this byte, from right to left
                for (int byteSmpNo = 0; byteSmpNo < 4; byteSmpNo++) {
                    if (smpNo < numSubjects) {
                        int bit1 = bedByte & baseNums[byteSmpNo * 2];
                        int bit2 = bedByte & baseNums[byteSmpNo * 2 + 1];

                        if (bit1 && !bit2) {       /// 10 = Missing genotype
                            sbjOtherCnts[smpNo]++;
                        }
                        else if (!bit1 && bit2) {  /// 01 = Hetrozygous
                            sbjHetroCnts[smpNo]++;
                        }
                        else {
                            sbjHomoCnts[smpNo]++;  /// Homozygous
                        }
                    }

                    smpNo++;
                }
            }

            ASSERT(bimFpSnpNo < numFpSnps, "bimFpSnpNo = " << bimFpSnpNo << " numFpSnps = " << numFpSnps);

            bimFpSnpNo++;
        }
    }

    bedFilePtr.close();

    cout << "Read genotypes of " << numSubjects << " samples " << numBimSnps << " SNPs\n";

    for (int i = 0; i < numSubjects; i++) {
        famSbjs->subjects[i].numHetros = sbjHetroCnts[i];
        famSbjs->subjects[i].numFpSnps = sbjHetroCnts[i] + sbjHomoCnts[i];
    }

    Summarize();

    return hasErr;
}

void SexCheck::Summarize()
{
    numMatches = 0;
    numErrSbjs = 0;
    numConflicts = 0;

    int sbjNo = 0;

    for (int i = 0; i < numSubjects; i++) {
        PedigreeSubject sbj = famSbjs->subjects[i];

        PhenotypeSex phenoSex = sbj.phenoSex;
        GenotypeSex genoSex = GenotypeSex::NoInfo;



        if (sbj.numFpSnps >= minFpSnps) {
            float hetroRate = sbj.numHetros * 1.0 / sbj.numFpSnps;
            famSbjs->subjects[i].hetroRate = hetroRate;

            if (hetroRate < maxMaleHetroRate) {
                genoSex = GenotypeSex::Male;
            }
            else if (hetroRate < minFemaleHetroRate) {
                genoSex = GenotypeSex::NoInfo;
            }
            else if (hetroRate < maxFemaleHetroRate) {
                genoSex = GenotypeSex::Female;
            }
            else {
                genoSex = GenotypeSex::Error;
            }
        }

        PhenoGenoSexCheck match = PhenoGenoSexCheck::NoInfo;

        /// Valid means either male or female
        if (SexIsValid(phenoSex) && SexIsValid(genoSex)) {
            if (GenoSexMatchesPhenoSex(genoSex, phenoSex)) {
                match = PhenoGenoSexCheck::Match;
                numMatches++;
            }
            else if (genoSex == GenotypeSex::Error) {
                match = PhenoGenoSexCheck::Error;
                numErrSbjs++;
            }
            else {
                match = PhenoGenoSexCheck::Mismatch;
                numConflicts++;
            }
        }
        else {
            match = PhenoGenoSexCheck::NoInfo;
            numOtherSbjs++;
        }

        /// For debugging. cout doesn't take enum class values
        /// cout << "sex " << static_cast<int>(phenoSex) << " geno " << static_cast<int>(genoSex)
        ///      << " match " << static_cast<int>(match) << "\n";

        famSbjs->subjects[i].genoSex = genoSex;
        famSbjs->subjects[i].sexMatch = match;

        sbjNo++;
    }
}


void SexCheck::ShowSummary()
{
    cout << "\nTotal " << numSubjects << " samples checked\n";
    cout << "\t" << numMatches << " samples have sex matches\n";
    cout << "\t" << numConflicts << " samples have sex discrepancies\n";
    cout << "\t" << numOtherSbjs << " samples have missing phenotype sexes or no sexes determined using genotypes\n";
    cout << "\n";
}

int SexCheck::SaveResults(string outputFile)
{
    FILE *ofp = fopen(outputFile.c_str(), "w");

    if(ofp) {
        fputs("sample_id\tpheno_sex\tgeno_sex\tmatch\theterozygous_snps\tcompared_snps\theterozygosity_rate\n", ofp);
        for (int i = 0; i < numSubjects; i++) {
            PedigreeSubject sbj = famSbjs->subjects[i];

            fprintf(ofp, "%s\t%d\t%d\t%d\t%d\t%d\t%5.4f\n",
            sbj.name.c_str(), sbj.phenoSex, sbj.genoSex, sbj.sexMatch,
            sbj.numHetros, sbj.numFpSnps, sbj.hetroRate);
        }
    }
    else {
        cout << "ERROR: Can't open " << outputFile << " for writing!\n";
        return false;
    }

    fclose(ofp);

    return numSubjects;
}
