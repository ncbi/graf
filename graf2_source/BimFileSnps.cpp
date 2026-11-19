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
#include <sstream>

using namespace std;


BimFileSnps::BimFileSnps(int minSnps = 20)
{
    minSexFpSnps = minSnps;
    filename = "";
    numBimSnps = 0;
    numXchrFpSnps = 0;
    numGoodFpSnps = 0;
}

int BimFileSnps::GetRsNumFromSnpName(string snpName)
{
    int rsNum = 0;
    int chNo = 0;
    bool isValid = true;

    while(isValid && chNo < snpName.length()) {
        if (chNo == 0) {
            if (snpName[chNo] != 'r' && snpName[chNo] != 'R') {
                isValid = false;
            }
        }
        else if (chNo == 1) {
            if (snpName[chNo] != 's' && snpName[chNo] != 'S') {
                isValid = false;
            }
        }
        else {
            int value = snpName[chNo] - '0';
            if (value >= 0 && value < 10) {
                rsNum = rsNum * 10 + value;
            }
            else {
                isValid = false;
            }
        }

        chNo++;
    }

    if (!isValid) rsNum = -1;

    return rsNum;
}

char BimFileSnps::FlipAllele(char allele)
{
    char flipAllele = '0';

    switch(allele) {
        case 'A': flipAllele = 'T'; break;
        case 'T': flipAllele = 'A'; break;
        case 'G': flipAllele = 'C'; break;
        case 'C': flipAllele = 'G'; break;
    }

    return flipAllele;
}

int BimFileSnps::CompareFpSnpAlleles(const string& str1, const string& str2, const char expA1, const char expA2)
{
    int match = 0;

    if (str1.length() == 1 && str2.length() == 1) {
        char a1 = str1[0];
        char a2 = str2[0];

        char fa1 = FlipAllele(a1);
        char fa2 = FlipAllele(a2);

        if (a1 == expA1 && a2 == expA2) {
            match = 1;
        }
        else if (a1 == expA2 && a2 == expA1) {
            match = 2;
        }
        else if (fa1 == expA1 && fa2 == expA2) {
            match = -1;
        }
        else if (fa1 == expA2 && fa2 == expA1) {
            match = -2;
        }
    }

    return match;
}

int BimFileSnps::ReadXchrFpSnpsFromFile(string bimFile, XchrFpSnps* xchrFpSnps)
{
    cout << "Reading SNPs from file " << bimFile << "\n";

    ASSERT(FileExists(bimFile.c_str()), "File " << bimFile << " does not exist.");

    filename = bimFile;

    int lineLen = 300;
    char fpLine[lineLen];

    FILE *ifp = fopen(bimFile.c_str(), "r");
    ASSERT(ifp, "Could not open " << bimFile << "\n");

    bool fileIsValid = true;

    int rs, pos;
    string chrStr, rsStr, cm, a1, a2;

    while (fgets(fpLine, lineLen, ifp) != NULL && fileIsValid == true) {
        // sscanf(fpLine, "%s %s %s %d %s %s", chrStr, rsStr, cm, &pos, a1, a2);

        istringstream iss(fpLine);
        iss >> chrStr >> rsStr >> cm >> pos >> a1 >> a2;

        int rsNum = GetRsNumFromSnpName(rsStr);
        int fpSnpId = xchrFpSnps->FindSnpIdGivenRs(rsNum);

        if (fpSnpId > -1) {
            SexFpSnp fpSnp = xchrFpSnps->GetSexFpSnp(fpSnpId);
            char expA1 = fpSnp.allele1;
            char expA2 = fpSnp.allele2;

            int match = CompareFpSnpAlleles(a1, a2, expA1, expA2);

            numXchrFpSnps++;

            if (match) {
                numGoodFpSnps++;
            }
            else {
                fpSnpId = -1;
            }
        }

        bimSnpFpSnpIds.push_back(fpSnpId);

        numBimSnps++;
    }
    fclose(ifp);

    return numBimSnps;
}

bool BimFileSnps::HasEnoughXchrFpSnps()
{
    return numGoodFpSnps >= minSexFpSnps;
}

void BimFileSnps::ShowSummary()
{
    int numBadFpSnps = numXchrFpSnps - numGoodFpSnps;

    cout << "Total " << numBimSnps << " SNPs in bim file.\n";
    cout << "\t" << numXchrFpSnps << " SNPs are X-chr FP SNPs.\n";
    cout << "\t" << numGoodFpSnps << " SNPs have expected alleles.\n";

    if (numBadFpSnps > 0) {
        cout << numBadFpSnps << " SNPs do not have expected alleles.\n";
    }
    cout << "\n";
    cout << "Minimum FP SNPs needed to do sex check: " << minSexFpSnps << "\n";

    if (HasEnoughXchrFpSnps()) {
        cout << "There are enough FP SNPs to do sex checking.\n";
    }
    else {
        cout << "There are not enough FP SNPs to do sex checking.\n";
    }

    cout << "\n";
}
