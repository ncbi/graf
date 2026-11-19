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
#include "FamSubjects.h"

using namespace std;


PedigreeSubject::PedigreeSubject(string sbj, string dad, string mom, PhenotypeSex gender)
{
    name = sbj;
    father = dad;
    mother = mom;
    phenoSex = gender;

    numFpSnps = 0;
    numHetros = 0;
    hetroRate = -1;
    genoSex = GenotypeSex::NoInfo;
    sexMatch = PhenoGenoSexCheck::NoInfo;
}

FamSubjects::FamSubjects(string file)
{
    filename = file;
    numFamSbjs = 0;
    numMales = 0;
    numFemales = 0;
    numUnkSexSbjs = 0;
}

bool FamSubjects::Summarize()
{
    bool success = false;
    int arraySize = subjects.size();

    if (arraySize == numFamSbjs) {
        for (int i = 0; i < numFamSbjs; i++) {
	    PedigreeSubject sbj = subjects[i];
	    if (sbj.phenoSex == PhenotypeSex::Male) {
	        numMales++;
	    }
	    else if (sbj.phenoSex == PhenotypeSex::Female) {
	        numFemales++;
	    }
	    else {
	        numUnkSexSbjs++;
	    }
	}

	success = true;
    }

    return success;
}

void FamSubjects::ShowSummary()
{
    bool success = false;
    cout << "Total " << numFamSbjs << " samples in fam file " << filename << ".\n";
    cout << "\t" << numMales << " males\n";
    cout << "\t" << numFemales << " females\n";
    if (numUnkSexSbjs > 0) cout << "\t" << numUnkSexSbjs << " subjects don't hase sex information\n";
    cout << "\n";
}

int FamSubjects::ReadSubjectsFromFile()
{
    int numFileSbjs = 0;

    ASSERT(FileExists(filename.c_str()), "File " << filename << " does not exist.");

    int lineLen = 300;
    char fpLine[lineLen];

    FILE *ifp = fopen(filename.c_str(), "r");
    ASSERT(ifp, "Couldn't open file " << filename << ".\n");

    int lineNo = 0;
    bool fileIsValid = true;

    int numSbjs = 0;
    int sex, pheno;
    char famId[80], sbjId[80], dadId[80], momId[80];

    while (fgets(fpLine, lineLen, ifp) != NULL && fileIsValid == true) {
        sscanf(fpLine, "%s %s %s %s %d %d", famId, sbjId, dadId, momId, &sex, &pheno);

	if (!sex) sex = 0;

	if (sbjId) {
	    PedigreeSubject subject(sbjId, dadId, momId, (PhenotypeSex)sex);
	    subjects.push_back(subject);
	    numFileSbjs++;
	}

        lineNo++;
    }
    fclose(ifp);

    numFamSbjs = numFileSbjs;

    Summarize();

    return numFileSbjs;
}

