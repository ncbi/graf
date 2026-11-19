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

bool FileExists (char* filename) {
  FILE* fp = fopen(filename, "r");
  if (fp) {
    fclose(fp);
    return true;
  }
  else {
    return false;
  }
}

int GetGapGenoValue(unsigned long bit1, unsigned long bit2)
{
  int val = 0;
  if      ( bit1 &&  bit2) val = 2;
  else if (!bit1 &&  bit2) val = 1;
  else if (!bit1 && !bit2) val = 0;
  else                     val = 3;

  return val;
}

int main(int argc, char* argv[])
{
  if (argc < 3) {
    printf("Usage: LoadFpGenotypes <PLINK set> <output_file>\n");
    printf("\n");
    exit(0);
  }

  int i, j, k;

  char* baseFile = argv[1];
  char* outFile = argv[2];
  char bimFile[512];
  char bedFile[512];
  char famFile[512];

  sprintf(bimFile, "%s.bim", baseFile);
  sprintf(bedFile, "%s.bed", baseFile);
  sprintf(famFile, "%s.fam", baseFile);

  if (!FileExists(bimFile)) {
    cout << "ERROR: file " << bimFile << " doesn't exist.\n";
    exit(0);
  }
  if (!FileExists(bedFile)) {
    cout << "ERROR: file " << bedFile << " doesn't exist.\n";
    exit(0);
  }
  if (!FileExists(famFile)) {
    cout << "ERROR: file " << famFile << " doesn't exist.\n";
    exit(0);
  }

  SampleGeno* smpGeno = new SampleGeno(famFile, 0);
  bool hasErr = smpGeno->ReadGenotypesFromBedFile(bimFile, bedFile, 200);
  if (hasErr) {
    cout << "ERROR: failed to read genotypes from " << bedFile << "\n\n";
    exit(0);
  }

  FILE *ofp = fopen(outFile, "w");
  if (!ofp) {
    cout << "ERROR: failed to open " << outFile << " for writing!\n\n";
    exit(0);
  }

  int numFpGenoInts = (smpGeno->numFpSnps - 1) / 64 + 1;
  for (i = 0; i < smpGeno->numSmps; i++) {
    int smpId = smpGeno->samples[i].smpId;
    char* smpName = smpGeno->samples[i].smpName;
    //printf("No. %4d: ID %5d name %6s ", i, smpId, smpName);

    char genoStr[smpGeno->numFpHexChars+2+1];
    int hexNo = 0;
    for (j = 0; j < smpGeno->numFpGenoInts; j++) {
      unsigned long genoInt1 = smpGeno->smpGenos1[i][j];
      unsigned long genoInt2 = smpGeno->smpGenos2[i][j];

      // Save genotypes of two SNPs into one hex number. One pair of 64-bit integers can be saved into 32 hex numbers.
      for (k = 0; k < 32; k++) {
	int snpPos1 = 64 - k*2 - 1;  // #63 is the most siginificant digit
	int snpPos2 = snpPos1 - 1;

	unsigned long snp1bit1 = genoInt1 & smpGeno->baseNums[snpPos1];
	unsigned long snp1bit2 = genoInt2 & smpGeno->baseNums[snpPos1];

	unsigned long snp2bit1 = genoInt1 & smpGeno->baseNums[snpPos2];
	unsigned long snp2bit2 = genoInt2 & smpGeno->baseNums[snpPos2];

	int genoVal1 = GetGapGenoValue(snp1bit1, snp1bit2);
	int genoVal2 = GetGapGenoValue(snp2bit1, snp2bit2);

	int hexVal = genoVal1 * 4 + genoVal2;
	char hex = '0' + hexVal;
	if (hexVal > 9) hex = 'A' + hexVal - 10;

	if (hexNo < smpGeno->numFpHexChars+3) genoStr[hexNo] = hex;
	hexNo++;
      }
    }

    if (hexNo < smpGeno->numFpHexChars+3) {
      genoStr[hexNo] = 'F';
      genoStr[hexNo+1] = 'F';
    }
    genoStr[smpGeno->numFpHexChars+2] = 0;

    fprintf(ofp, "%s\t%s\n", smpName, genoStr);
  }

  fclose(ofp);
  printf("Saved genotypes of %d samples to %s.\n\n", smpGeno->numSmps, outFile);

  return 1;
}
