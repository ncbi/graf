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

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
#include "SampleGenoDist.h"

int main(int argc, char* argv[])
{
  if (argc < 2) {
    cout << "Usage: transTest <input_file>\n";;
    exit(0);
  }

  char* inFile = argv[1];

  cout << "Reading data from " << inFile << "\n";

  GenoDist ec, fc, ac, sp;

  if (inFile) {
    char line[256];
    int numLines = 0;
    FILE *inFilePtr = fopen(inFile, "r");
    if(inFilePtr) {
      char smp[80];
      double eDist, fDist, aDist;;
      while (fgets(line, 256, inFilePtr) != NULL) {
	sscanf(line, "%s %lf %lf %lf", smp, &eDist, &fDist, &aDist);
	if (smp && aDist) {
	  if      (smp[0] == 'W') {
	    ec.e = eDist;
	    ec.f = fDist;
	    ec.a = aDist;
	  }
	  else if (smp[0] == 'B') {
	    fc.e = eDist;
	    fc.f = fDist;
	    fc.a = aDist;
	  }
	  else if (smp[0] == 'A') {
	    ac.e = eDist;
	    ac.f = fDist;
	    ac.a = aDist;
	  }
	  else if (smp[0] == 'S') {
	    sp.e = eDist;
	    sp.f = fDist;
	    sp.a = aDist;
	  }
	  numLines++;
	}
      }
      fclose(inFilePtr);
      cout << "Read " << numLines << " dist numbers from " << inFile << ".\n";

      if (ec.e && fc.e && ac.e && sp.e) {
	SampleGenoDist d1(&ec, &fc, &ac, &sp);
	d1.TransformAllDists();
	d1.CalculateBaryCenters();
	cout << "Set dist\n";
      }
    }
  }
}
