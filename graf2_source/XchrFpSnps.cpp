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

using namespace std;


SexFpSnp::SexFpSnp(int id, int rsNum, int g37, int g38, char a1, char a2)
{
    snpId = id;
    rs = rsNum;
    posG37 = g37;
    posG38 = g38;
    allele1 = a1;
    allele2 = a2;
}

XchrFpSnps::XchrFpSnps()
{
    numXchrFpSnps = 0;
}

int XchrFpSnps::ReadXchrFpSnpsFromFile(string sexFpSnpFile)
{
    ASSERT(FileExists(sexFpSnpFile.c_str()), "File " << sexFpSnpFile << " does not exist!\n");

    int lineLen = 300;
    char fpLine[lineLen];

    FILE *ifp = fopen(sexFpSnpFile.c_str(), "r");
    ASSERT(ifp, "Couldn't open file " << sexFpSnpFile << ".\n");

    int lineNo = 0;
    bool fileIsValid = true;

    int numFpSnps = 0;
    int rsNum, chr, g37, g38;
    char a1, a2;

    while (fgets(fpLine, lineLen, ifp) != NULL && fileIsValid == true) {
        if (lineNo == 0) {
            if (fpLine[0] != 'r' || fpLine[1] != 's' || fpLine[2] != '#') {
                fileIsValid = false;
            }
        }
        else {
            sscanf(fpLine, "%d %d %d %d %c %c", &rsNum, &chr, &g37, &g38, &a1, &a2);

            if (chr == 23) {
                SexFpSnp fpSnp(numXchrFpSnps, rsNum, g37, g38, a1, a2);
                fpSnps.push_back(fpSnp);

                rsToFpSnpId[rsNum] = numXchrFpSnps;
                numXchrFpSnps++;
            }

            numFpSnps++;
        }

        lineNo++;
    }
    fclose(ifp);

    cout << "Read " << numXchrFpSnps << " X-chr FP SNPs from file " << sexFpSnpFile << "\n\n";
    ASSERT(numXchrFpSnps == 1000, "numXchrFpSnps = " << numXchrFpSnps << ".\n");
    return numXchrFpSnps;
}

int XchrFpSnps::FindSnpIdGivenRs(int rsNum)
{
    int snpId = -1;

    if (rsToFpSnpId.find(rsNum) != rsToFpSnpId.end()) {
        snpId = rsToFpSnpId[rsNum];
    }

    return snpId;
}

SexFpSnp XchrFpSnps::GetSexFpSnp(int snpId)
{
    return fpSnps[snpId];
}

void XchrFpSnps::SaveXchrFpSnps(string outputFile)
{
    int numFpSnps = fpSnps.size();

    FILE *ofp = fopen(outputFile.c_str(), "a");
    if(ofp) {
        int snpNo = 0;
        for (int i = 0; i < numFpSnps; i++) {
            SexFpSnp fpSnp = fpSnps[i];
            char line[128];
            sprintf(line, "%d\t%d\t%d\t%d\t%c\t%c", i, fpSnp.rs, fpSnp.posG37, fpSnp.posG38, fpSnp.allele1, fpSnp.allele2);
            fprintf(ofp, "%s\n", line);
        }
        fclose(ofp);
    }

    cout << "Saved " << numFpSnps << " FP SNPs to " << outputFile << "\n";
}

void XchrFpSnps::ShowXchrFpSnps()
{
    int numFpSnps = fpSnps.size();

    cout << "Total " << numXchrFpSnps << " FP SNPs on X chromosome. #FP SNPs in array: " << numFpSnps << "\n";
}
