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


#ifndef __TEST_XCHRFPSNPS_H
#define __TEST_XCHRFPSNPS_H

#include "Graf.h"
#include "XchrFpSnps.h"

#include <cxxtest/TestSuite.h>
#include <iostream>

using namespace std;

class Test_XchrFpSnps : public CxxTest::TestSuite {
    XchrFpSnps * xchrFpSnps;
    int result;

public:
    Test_XchrFpSnps() : xchrFpSnps(0), result(0) {
        cout << "\n\nTesting FindSnpIdGivenRs function\n";
    }
    virtual ~Test_XchrFpSnps() { delete xchrFpSnps; }

    static Test_XchrFpSnps *createSuite() { return new Test_XchrFpSnps(); }
    static void destroySuite(Test_XchrFpSnps *suite) { delete suite; }

    void setUp()
    {
        xchrFpSnps = new XchrFpSnps();
        string fpSnpFile = "data/FP_SNPs.txt";
        result = xchrFpSnps->ReadXchrFpSnpsFromFile(fpSnpFile);
    }

    void tearDown()
    {
        if (xchrFpSnps != 0) {
            delete xchrFpSnps;
        }

        xchrFpSnps = 0;
    }

    void test_fp_snp1()
    {
        int expected = 5;
        int rs = 1205545;
        int actual = xchrFpSnps->FindSnpIdGivenRs(rs);
        TS_ASSERT_EQUALS(expected, actual);
        cout << "rs" << rs << " is an FP SNP on X-chromosome.\n";;
    }

    void test_fp_snp2()
    {
        int expected = -1;
        int rs = 804;
        int actual = xchrFpSnps->FindSnpIdGivenRs(rs);
        TS_ASSERT_EQUALS(expected, actual);
        cout << "rs" << rs << " is not an FP SNP on X-chromosome.\n";;
    }
};

#endif /// __TEST_XCHRFPSNPS_H
