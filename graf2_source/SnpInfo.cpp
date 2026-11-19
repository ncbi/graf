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

using namespace std;


SnpInfo::SnpInfo(int snps)
{
    numSnps = snps;
    numSmps = 0;
    aaCnts = (int*) malloc(numSnps * sizeof(int));
    abCnts = (int*) malloc(numSnps * sizeof(int));
    bbCnts = (int*) malloc(numSnps * sizeof(int));
    ooCnts = (int*) malloc(numSnps * sizeof(int));
    refCnts = (int*) malloc(numSnps * sizeof(int));
    altCnts = (int*) malloc(numSnps * sizeof(int));
    mafs = (double*) malloc(numSnps * sizeof(double));
    crs  = (double*) malloc(numSnps * sizeof(double));
    mafCrs = (double*) malloc(numSnps * sizeof(double));
    
    for (int i = 0; i < numSnps; i++) {
        aaCnts[i] = 0;
        abCnts[i] = 0;
        bbCnts[i] = 0;
        ooCnts[i] = 0;
        refCnts[i] = 0;
        altCnts[i] = 0;
        mafs[i] = 0;
        crs[i] = 0;
        mafCrs[i] = 0;
    }
}

SnpInfo::~SnpInfo()
{
    delete aaCnts;
    delete abCnts;
    delete bbCnts;
    delete ooCnts;
    delete refCnts;
    delete altCnts;
    delete mafs;
    delete crs;
    delete mafCrs;
}

void SnpInfo::AddSnpGeno(int snpNo, int geno)
{
    /// geno: 0 = aa; 1 = ab; 2 = bb
    assert(snpNo > -1 && snpNo < numSnps);
    if (geno == 0) {
        aaCnts[snpNo]++;
        refCnts[snpNo] += 2;
    }
    else if (geno == 1) {
        abCnts[snpNo]++;
        refCnts[snpNo] += 1;
        altCnts[snpNo] += 1;
    }
    else if (geno == 2) {
        bbCnts[snpNo]++;
        altCnts[snpNo] += 2;
    }
}

void SnpInfo::Summarize()
{
    for (int i = 0; i < numSnps; i++) {
        double maf = -1, cr = 0, mafCr = 0;
        int callGenos = aaCnts[i] + abCnts[i] + bbCnts[i];
        int totGenos = callGenos + ooCnts[i];
        
        refCnts[i] = aaCnts[i] * 2 + abCnts[i];
        altCnts[i] = bbCnts[i] * 2 + abCnts[i];
        
        if (callGenos > 0) {
            maf = refCnts[i] < altCnts[i] ? refCnts[i] / (callGenos * 2.0) : altCnts[i] / (callGenos * 2.0);
            cr = callGenos / totGenos;
            mafCr = maf * cr;
        }
        mafs[i] = maf;
        crs[i] = cr;
        mafCrs[i] = mafCr;
    }
}

void SnpInfo::Show(int startSnp, int endSnp, int step)
{
    if (endSnp > numSnps) endSnp = numSnps;
    for (int snpNo = startSnp; snpNo <= endSnp; snpNo += step) {
        printf("SNP No. %5d: MAF %7.5f   CR %7.5f   MAFCR %7.5f\n", mafs[snpNo], crs[snpNo], mafCrs[snpNo]);
    }
}
