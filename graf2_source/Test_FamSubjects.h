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


#ifndef __TEST_FAMSUBJECTS_H
#define __TEST_FAMSUBJECTS_H

#include "Graf.h"
#include "FamSubjects.h"

#include <cxxtest/TestSuite.h>
#include <iostream>

using namespace std;

class Test_FamSubjects : public CxxTest::TestSuite {
    FamSubjects * famSbjs;
    int result;

public:
    Test_FamSubjects() : famSbjs(0), result(0) {
        cout << "\n\nTesting FamSubjects function\n";
    }
    virtual ~Test_FamSubjects() { delete famSbjs; }

    static Test_FamSubjects *createSuite() { return new Test_FamSubjects(); }
    static void destroySuite(Test_FamSubjects *suite) { delete suite; }

    void setUp()
    {
        string plinkFile = "data/affy_hapmap.fam";
        famSbjs = new FamSubjects(plinkFile);
        result = famSbjs->ReadSubjectsFromFile();
    }

    void tearDown()
    {
        if (famSbjs != 0) {
            delete famSbjs;
        }

        famSbjs = 0;
    }

    void test_read_subjects()
    {
        int expected = 277;
        int actual = famSbjs->GetNumFamSubjects();
        cout << "Subjects read: " << actual << ", expected: " << expected << "\n";;
        TS_ASSERT_EQUALS(expected, actual);
    }
};

#endif /// __TEST_FAMSUBJECTS_H
