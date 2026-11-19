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

#ifndef FAM_SUBJECTS_H
#define FAM_SUBJECTS_H

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

enum class PhenotypeSex
{
    Male = 1,
    Female = 2,
    Unknown = 0
};

enum class GenotypeSex
{
    Male = 1,
    Female = 2,
    NoInfo = 0,
    Error = -1
};

enum class PhenoGenoSexCheck
{
    Mismatch = 0,
    Match = 1,
    NoInfo = -1,
    Error = -2
};

class PedigreeSubject
{
public:
    string father;
    string mother;
    int numFpSnps;
    int numHetros;
    float hetroRate;

    string name;
    PhenotypeSex phenoSex;
    GenotypeSex genoSex;
    PhenoGenoSexCheck sexMatch;

public:
    PedigreeSubject(string, string, string, PhenotypeSex);
};

class FamSubjects
{
    string filename;
    int numFamSbjs;
    int numMales;
    int numFemales;
    int numUnkSexSbjs;

    bool Summarize();

public:
    vector<PedigreeSubject> subjects;

    FamSubjects(string);
    int ReadSubjectsFromFile();
    int GetNumFamSubjects() {return numFamSbjs;};
    void ShowSummary();
    void SaveResults(string);
};


#endif
