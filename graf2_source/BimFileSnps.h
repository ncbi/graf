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

#ifndef BIM_FILE_SNPS_H
#define BIM_FILE_SNPS_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
#include <map>
#include <vector>

using namespace std;

class BimFileSnps
{
    int minSexFpSnps;

    string filename;
    int numBimSnps;
    int numXchrFpSnps;
    int numGoodFpSnps;

    map<int, int> rsToFpSnpId;

    int GetRsNumFromSnpName(string);
    char FlipAllele(char);

    /// For each bim SNP, if it is a FP SNP, the FP SNP ID is saved here
    /// FP SNP ID is 0-based. If a bim SNP is not an FP SNP, the SNP ID is -1
    vector<int> bimSnpFpSnpIds;

public:
    BimFileSnps(int);
    int ReadXchrFpSnpsFromFile(string, XchrFpSnps*);
    int CompareFpSnpAlleles(const string&, const string&, const char, const char);
    int GetNumSnps() { return numBimSnps; };
    int GetGoodFpSnps() { return numGoodFpSnps; };
    int GetXchrFpSnps() { return numXchrFpSnps; };
    int GetFpSnpIdGivenBimSnpPos(int bimSnpPos) {
        return bimSnpPos >= 0 && bimSnpPos < bimSnpFpSnpIds.size() ? bimSnpFpSnpIds[bimSnpPos] : -1;
    };
    bool HasEnoughXchrFpSnps();
    void ShowSummary();
};


#endif
