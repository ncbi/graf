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


#ifndef __TEST_BIMFILESNPS_H
#define __TEST_BIMFILESNPS_H

#include <cxxtest/TestSuite.h>
#include <iostream>

#include "Graf.h"
#include "XchrFpSnps.h"
#include "BimFileSnps.h"

using namespace std;

class Test_BimFileSnps : public CxxTest::TestSuite {

    BimFileSnps * bimFile;
    char a1[2], a2[2];
    char e1 = 'A';
    char e2 = 'G';


public:
    void setUp()
    {
        bimFile = new BimFileSnps(1);
        strcpy(a1, "G");
        strcpy(a2, "A");
    }

    void tearDown()
    {
        delete bimFile;
        bimFile = 0;
    }

    Test_BimFileSnps() : bimFile(0) {
        cout << "\n\nTesting CompareFpSnpAlleles function\n";
    }

    virtual ~Test_BimFileSnps() {}

    static Test_BimFileSnps *createSuite() { return new Test_BimFileSnps(); }
    static void destroySuite(Test_BimFileSnps *suite) { delete suite; }

    void test_ag_ga_2() {

        int expected = 2;
        int actual = 0;

        actual = bimFile->CompareFpSnpAlleles(a1, a2, e1, e2);
        cout << "\t" << e1 << e2 << " / " << a1 << a2 << " => " << actual << "\n";

        TS_ASSERT_EQUALS(expected, actual)
    }

    void test_ag_ag_1() {
        int expected = 1;
        int actual = 0;

        a1[0] = 'A';
        a2[0] = 'G';
        actual = bimFile->CompareFpSnpAlleles(a1, a2, e1, e2);
        cout << "\t" << e1 << e2 << " / " << a1 << a2 << " => " << actual << "\n";

        TS_ASSERT_EQUALS(expected, actual)

    }

    void test_ag_tc_negative_1() {
        int expected = -1;
        int actual = 0;

        a1[0] = 'T';
        a2[0] = 'C';

        actual = bimFile->CompareFpSnpAlleles(a1, a2, e1, e2);
        cout << "\t" << e1 << e2 << " / " << a1 << a2 << " => " << actual << "\n";

        TS_ASSERT_EQUALS(expected, actual)
    }


    void test_ag_ct_negative_2() {
        int expected = -2;
        int actual = 0;

        a1[0] = 'C';
        a2[0] = 'T';
        actual = bimFile->CompareFpSnpAlleles(a1, a2, e1, e2);
        cout << "\t" << e1 << e2 << " / " << a1 << a2 << " => " << actual << "\n";

        TS_ASSERT_EQUALS(expected, actual)
    }


    void test_ag_ac_0() {
        int expected = 0;
        int actual = -100;

        a1[0] = 'A';
        a2[0] = 'C';
        actual = bimFile->CompareFpSnpAlleles(a1, a2, e1, e2);
        cout << "\t" << e1 << e2 << " / " << a1 << a2 << " => " << actual << "\n";

        TS_ASSERT_EQUALS(expected, actual)
    }


    void test_ag_aa_0() {
        int expected = 0;
        int actual = -100;

        a1[0] = 'A';
        a2[0] = 'G';

        a1[1] = 'A';
        a1[2] = '\0';

        actual = bimFile->CompareFpSnpAlleles(a1, a2, e1, e2);
        cout << "\t" << e1 << e2 << " / " << a1 << a2 << " => " << actual << "\n";

        TS_ASSERT_EQUALS(expected, actual)
    }
};

#endif // __TEST_BIMFILESNPS_H
