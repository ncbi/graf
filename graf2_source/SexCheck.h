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

#ifndef SEX_CHECK_H
#define SEX_CHECK_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
#include <map>
#include <vector>

#include "XchrFpSnps.h"
#include "BimFileSnps.h"
#include "FamSubjects.h"
#include "SexCheck.h"

using namespace std;
static const int BYTE1_IN_BED_FILE = 108;
static const int BYTE2_IN_BED_FILE = 27;
static const int BYTE_OF_SNP_MODE  = 1;

class SexCheck
{
    int minFpSnps;              /// Minimum SNPs needed for determining pheno sexes
    unsigned long baseNums[64]; /// bits 1, 10, 100 ... for decoding genos in bed file

    /// Cutoff hetrozygosity rates to determine geno sexes
    float maxMaleHetroRate;
    float minFemaleHetroRate;
    float maxFemaleHetroRate;

    int numBimSnps;
    int numSubjects;
    int numConflicts; /// Subjects with discrepancies between pheno and geno sexes
    int numMatches;   /// Subjects with geno sexes same as pheno sexes
    int numErrSbjs;   /// Subjects with potential genotype errors
    int numOtherSbjs; /// Subjects without enough information to determine pheno or geno sexes

    string bedFile;
    XchrFpSnps *sexFpSnps;
    BimFileSnps *bimSnps;
    FamSubjects *famSbjs;

    bool SexIsValid(PhenotypeSex sex) {
        return sex == PhenotypeSex::Male || sex == PhenotypeSex::Female ? true : false;
    };

    bool SexIsValid(GenotypeSex sex) {
        return sex == GenotypeSex::Male || sex == GenotypeSex::Female ? true : false;
    };

    bool GenoSexMatchesPhenoSex(GenotypeSex genoSex, PhenotypeSex phenoSex) {
        bool match = false;

        if ( (phenoSex == PhenotypeSex::Male   && genoSex == GenotypeSex::Male) ||
            (phenoSex == PhenotypeSex::Female && genoSex == GenotypeSex::Female) )
            match = true;

            return match;
        };

        void Summarize();

public:
        SexCheck(string, XchrFpSnps*, BimFileSnps*, FamSubjects*);
        ~SexCheck();
        void SetCutoffValues(int, float, float, float);
        bool CalculateHetroRates();
        void ShowSummary();
        int SaveResults(string);
    };


    #endif
