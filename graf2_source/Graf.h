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

using namespace std;

#ifndef NDEBUG
#   define ASSERT(condition, message) \
    if (!(condition)) {	\
        cerr << "Assertion `" #condition "` failed in " << __FILE__  << " line " << __LINE__ << ": " << message << "\n"; \
        terminate(); \
    } 
#else
#   define ASSERT(condition, message) {}
#endif

#ifndef GRAF_H
#define GRAF_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <map>
#include <vector>
#include <unistd.h>
#include <thread>
#include <mutex>

struct GenoDist
{
    double e;
    double f;
    double a;
};

struct Point
{
    double x;
    double y;
    double z;
};

struct FpSnp
{
    int snpId;
    int rs;
    int chr;
    int pos;
    char a1;
    char a2;
    float maf;
    float popAf[7];
    float vertexAf[3];
};

class SampleGenoDist
{
private:
    GenoDist eDist;
    GenoDist fDist;
    GenoDist aDist;
    GenoDist sDist;

public:
    Point ePt;
    Point fPt;
    Point aPt;
    Point sPt;

    double eWt;
    double fWt;
    double aWt;

    SampleGenoDist(GenoDist*, GenoDist*, GenoDist*, GenoDist*);

    void SetDistWithDist(GenoDist*, GenoDist*);
    void SetPointWithDist(Point*, GenoDist*);
    void SetPointWithPoint(Point*, Point*);
    void MovePoint(Point*, double, double, double);
    void RotatePointOnx(Point*, double);
    void RotatePointOny(Point*, double);
    void RotatePointOnz(Point*, double);
    Point TransformGenoDist(GenoDist);
    void TransformAllDists();
    void CalculateBaryCenters();
    void ShowPositions(int);
};

class FpSnps
{
public:
    int numFpSnps;
    int numChkPops;
    string checkPops[7];
    double popExpPfvals[7];
    double popExpPevals[7];
    double popExpPavals[7];
    int popNumSnps[7];

    map<int, int> rsToId;
    FpSnp* snps;

    FpSnps(int);
    ~FpSnps();
    int FindRs(int);
    void SetPopulationAFs(void);
    void SetVertexAFs(void);
    void SetPopulationExpectedPvalues(void);
    void SaveSnps(char*);
};

// For sorting genotypes by MAF and call rate
struct SnpMaf
{
    int snpNo;
    double maf, cr, mafCr;
};

// For sorting relative pairs by mismatch rate
struct RelMissRate
{
    int pairNo;
    double mr;
};

class SnpInfo
{
public:
    int numSnps;
    int numSmps;

    int* aaCnts;
    int* abCnts;
    int* bbCnts;
    int* ooCnts;
    int* refCnts;
    int* altCnts;
    double* mafs;
    double* crs;
    double* mafCrs;

public:
    SnpInfo(int);
    ~SnpInfo();
    void AddSnpGeno(int, int);
    void Summarize();
    void Show(int, int, int);
};


class SampleInfo
{
public:
    int  smpId;
    int  sbjId;
    char smpName[80];
    char sbjName[80];
    int  phs;
    int  sex;

    double* popPvalues;       /// P-values calculate using genotypes against each reference population
    int popNumSnps;           /// Number of valid SNPs checked for above p-values
    double* popMeanPvals;     /// P-values averaged over valid SNPs
    double* popExpPfsums;     /// Expected Pf-values for a given vertex population (E, F, A)
    double* popExpPesums;     /// Expected Pe-values for a given vertex population (E, F, A)
    double* popExpPasums;     /// Expected Pa-values for a given vertex population (E, F, A)
    double* popExpPfvals;     /// Expected Pf-values averaged over valid SNPs, for a given vertex population
    double* popExpPevals;     /// Expected Pe-values averaged over valid SNPs, for a given vertex population
    double* popExpPavals;     /// Expected Pa-values averaged over valid SNPs, for a given vertex population

    SampleInfo();
    SampleInfo(char*, char*, int, int=0, int=0, int=0);
};


class SampleGeno
{
public:
    int numSmps;
    int numFpSnps;             /// Number of SNPs used for fingerprinting
    int numFpHexChars;         /// Number of hex chars of the FP SNPs.  One char for two SNPs.
    int numFpGenoInts;         /// Number of long intergers to code the FP genotypes.  One integer for 64 SNPs
    int numThreads;            /// Number of threads for parallel computing

    unsigned long** smpGenos1; /// First  bit of the genotype for each SNP
    unsigned long** smpGenos2; /// Second bit of the genotype for each SNP
    unsigned long baseNums[64];
    SampleInfo* samples;
    SnpInfo* snpInfo;          /// MAF and missing rate at each SNP position
    double hgmrAdj;            /// Coefficient for normalizing HGMR
    FpSnps* fpSnps;

public:
    double poExpAgmr;          /// Expected genotype mismatch rates for each relationship
    double fsExpAgmr;
    double d2ExpAgmr;
    double d3ExpAgmr;
    double fsExpHgmr;
    double d2ExpHgmr;
    double d3ExpHgmr;
    int numValidSnps;

    SampleGeno();
    SampleGeno(int);
    SampleGeno(char*);         /// Construct with .geno file
    SampleGeno(char*, int);    /// Construct with .fam file
    ~SampleGeno();
    void DeleteSamples();
    bool ReadGenotypesFromBedFile(char*, char*, int, map<int, int>, bool=false);
    int ReadGenotypesFromHexGenoFile(char*, bool, int* = NULL, int* = NULL);
    void ReorderSnps();
    void ShowSamples(int);
    void SetSampleSubjectIds (map<string, string>);
    void SetMaxSnpMiss(int);
    int  GetNumberStudies();
    int  SetSampleGenders (map<string, int>);
    bool SaveHexGenotypes(char*);
    bool SavePopulationResults(char*);
    void SetPopulationPvalues(int);
    void InitPopPvalues();

private:
    int genoFileLineLen; /// Max length of one sample geno file (with sample info)
    int smpNameLen;
    int genoHexLen;
    int maxTotSnpMiss;

    char GetCompAllele(char);
    void SetConstants();
    int  GetSnpGenoInt(bool, bool);
};

class RelativePair
{
public:
    SampleInfo sample1, sample2;
    int homoSame, homoDiff;
    int allSame, allDiff;
    double hgmr, agmr;
    int pedRelation;
    int genoRelation;
    double pvalue;

    RelativePair();
};

class Relatives
{
public:
    int numPairs;
    vector<RelativePair> relatives;
    string genoRelTypes[6];
    string pedRelTypes[13];
    double idCutMr, poCutMr, fsCutMr, d2CutMr, d3CutMr;
    int numGapGenoSnps;  /// Average number of FP SNPs with genotypes in dbGaP

    double poMux, poSdx, poMuy, poSdy, poCor;
    double fsMux, fsSdx, fsMuy, fsSdy, fsCor;
    double d2Mux, d2Sdx, d2Muy, d2Sdy, d2Cor;
    double d3Mux, d3Sdx, d3Muy, d3Sdy, d3Cor;
    double poHgmrK;

    Relatives(double);
    double GetBivarProp(double, double, int, int);
    double GetPoBivarProp(double, double, int);
    void SetExpectedMismatchRates(double, double, double, double, double, double, double);
    void AddOnePair(SampleInfo, SampleInfo, int, int, double, int, int, double, int);
    void AddOnePair(RelativePair);
    bool SaveRelativesToFile(char*, int);
    void ShowSummary();
};

class Pedigree
{
public:
    int numGenoSmps;
    int numSsmSmps;
    int numPedSbjs;
    vector<string> genoSmps;
    vector<string> pedSbjList;

    bool ssmIsGap; // SSM file is in dbGaP format:
                   // the two colums are subject ID and sample ID, no matter what the headers are

    map<string, int> sbjSexes;
    map<string, string> ssm;
    map<string, string> sbjDads;
    map<string, string> sbjMoms;
    map<string, int> pedSbjs;
    map<string, string> sbjTwins;
    map<string, int> smpRels;
    map<string, int> sbjRels;
    string pedRelTypes[13];

    Pedigree();
    Pedigree(SampleInfo*, int);
    bool ReadSmpSbjMapping(char*);
    bool ReadPedigree(char*, bool);
    void ShowSummary(int);
    void StringToUpper(char*);
    string FindAncestor(map<string, string>*, string);
    int GetPedigreeRelatives();
    int GetPedLineValue(char*, int, char*);

    void SetSsmIsGapFormat(bool isGap) {ssmIsGap = isGap;};
};

struct RelParams
{
    char* plinkBase;
    char* genoFile;
    char* outFile;
    char* outPopFile;
    char* outSexFile;
    char* ssmFile;
    char* fpSnpFile;
    char* pedigreeFile;
    int checkRelType;
    int chkWordNum;
    int minHomoDiff;
    double maxHgmr;
    bool isUserHgmr;
    int expMrCalc;
    int ssmIsGap;
};

int CompareSnpMaf(const void*, const void*);
int ComparePairMissRate(const void*, const void*);

void ShowNumberInBits (unsigned long, int);
void ShowTimeDiff(struct timeval, struct timeval);
bool FileExists (const char*);
bool FileWriteable (const char*);

int CombinePlinkSets(char*, char*, map<int, int>, int);
vector<string> FindSamplesInFamFile(const char*, bool*);
Relatives* FindRelatedSubjects(struct RelParams, bool, map<int, int>, int);
void FindRelatedSubjectsForThread(int);
void SaveLog(vector<string>);

int CheckSubjectPopulations(struct RelParams, bool, map<int, int>, int);
int CheckSubjectSexes(struct RelParams, int);
string GetExecutablePath();

const double pi = 3.1415926;

#endif
