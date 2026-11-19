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

int ComparePairMissRate(const void *a, const void *b)
{
    struct RelMissRate *ia = (struct RelMissRate *)a;
    struct RelMissRate *ib = (struct RelMissRate *)b;
    
    return (int)(10000.0*ia->mr - 10000.0*ib->mr);
}


RelativePair::RelativePair()
{
    homoSame = 0;
    homoDiff = 0;
    allSame = 0;
    allDiff = 0;
    hgmr = -1;
    agmr = -1;
    pedRelation = 0;
    pvalue = -1;
}

Relatives::Relatives(double hgmrAdj)
{
    numPairs = 0;
    string genoTypes[] = {"ID", "PO", "FS", "D2", "D3", "UN"};
    string pedTypes[]  = {"", "DP", "MT", "PO", "FS", "SB", "HS", "GP", "AV", "FC", "HF", "HA", "HC"};
    
    for (int i = 0; i < 6;  i++) genoRelTypes[i] = genoTypes[i];
    for (int i = 0; i < 13; i++) pedRelTypes[i]  = pedTypes[i];
    
    idCutMr = 14;
    poCutMr = 1.4;
    fsCutMr = 8.2;
    d2CutMr = 14.0;
    d3CutMr = 20.0;
    numGapGenoSnps = 9000;
    
    poMux =  0.037292;
    poSdx =  0.073727;
    poMuy = 39.477589;
    poSdy =  1.468133;
    poCor =  0.104462;
    
    fsMux =  4.863680;
    fsSdx =  1.015553;
    fsMuy = 33.592671;
    fsSdy =  2.322399;
    fsCor =  0.784207;
    
    d2Mux = 11.162160;
    d2Sdx =  1.353599;
    d2Muy = 47.304737;
    d2Sdy =  1.096585;
    d2Cor =  0.802548;
    
    d3Mux = 17.320116;
    d3Sdx =  1.354330;
    d3Muy = 51.069507;
    d3Sdy =  0.958527;
    d3Cor =  0.830226;
    
    poHgmrK = 10.0;
    
    relatives.clear();
}

void Relatives::SetExpectedMismatchRates(double poAg, double fsAg, double d2Ag, double d3Ag,
double fsHg, double d2Hg, double d3Hg)
{
    poMuy = poAg;
    fsMuy = fsAg;
    d2Muy = d2Ag;
    d3Muy = d3Ag;
    
    fsMux = fsHg;
    d2Mux = d2Hg;
    d3Mux = d3Hg;
}

double Relatives::GetBivarProp(double x, double y, int numGenoSnps, int type)
{
    double dx = 0.0, dy = 0.0, rho = 0.0, sigx = 0.0, sigy = 0.0;
    double p = 0;
    
    if (type == 1) {          /// FS
        dx = x - fsMux;
        dy = y - fsMuy;
        sigx = fsSdx;
        sigy = fsSdy;
        rho = fsCor;
    }
    else if (type == 2) {    /// D2
        dx = x - d2Mux;
        dy = y - d2Muy;
        sigx = d2Sdx;
        sigy = d2Sdy;
        rho = d2Cor;
    }
    else if (type == 3) {    /// D3
        dx = x - d3Mux;
        dy = y - d3Muy;
        sigx = d3Sdx;
        sigy = d3Sdy;
        rho = d3Cor;
    }
    else {
        return p;
    }
    
    double sdAdj = sqrt(numGapGenoSnps * 1.0 / numGenoSnps);  /// Adjust SD based on number of FP SNPs with genos
    sigx = sigx * sdAdj;
    sigy = sigy * sdAdj;
    
    double z = dx*dx/(sigx*sigx) - 2 * rho * dx * dy / (sigx * sigy) + dy*dy/(sigy*sigy);
    p = 1 / (2 * pi * sigx * sigy * sqrt(1 - rho*rho) * exp(z / (2 * (1 - rho*rho))));
    
    return p;
}

double Relatives::GetPoBivarProp(double x, double y, int numGenoSnps)
{
    double p = 0.0;
    
    double sigy = poSdy * sqrt(numGapGenoSnps * 1.0 / numGenoSnps);
    p = poHgmrK * exp(-1*poHgmrK * x - (y - poMuy)*(y - poMuy)/(2 * sigy*sigy)) / (sigy * sqrt(2 * pi));
    
    return p;
}


void Relatives::AddOnePair(SampleInfo smp1, SampleInfo smp2, int hs, int hd, double hgmr, int aSame, int aDiff, double agmr, int pr)
{
    RelativePair pair;
    pair.sample1 = smp1;
    pair.sample2 = smp2;
    pair.homoSame = hs;
    pair.homoDiff = hd;
    pair.allSame = aSame;
    pair.allDiff = aDiff;
    pair.hgmr = hgmr;
    pair.agmr = agmr;
    pair.pedRelation = pr;
    
    int numGenoSnps = aSame + aDiff;
    double poProb = GetPoBivarProp(hgmr, agmr, numGenoSnps);
    double fsProb = GetBivarProp(hgmr, agmr, numGenoSnps, 1);
    double d2Prob = GetBivarProp(hgmr, agmr, numGenoSnps, 2);
    double d3Prob = GetBivarProp(hgmr, agmr, numGenoSnps, 3);
    double allProb = poProb + fsProb + d2Prob + d3Prob;
    
    int genoRel = 5;
    if (agmr < idCutMr) {
        genoRel = 0;
    }
    else {
        genoRel = 1;
        double maxProb = poProb;
        
        if (fsProb > maxProb) {
            genoRel = 2;
            maxProb = fsProb;
        }
        
        if (d2Prob > maxProb) {
            genoRel = 3;
            maxProb = d2Prob;
        }
        
        if (d3Prob > maxProb) {
            genoRel = 4;
            maxProb = d3Prob;
        }
    }
    
    double pval = -1;
    if (genoRel == 1) {
        pval = (fsProb + d2Prob + d3Prob) / allProb;
    }
    else if (genoRel == 2) {
        pval = (poProb + d2Prob + d3Prob) / allProb;
    }
    else if (genoRel == 3) {
        pval = (poProb + fsProb + d3Prob) / allProb;
    }
    
    pair.genoRelation = genoRel;
    pair.pvalue = pval;
    
    AddOnePair(pair);
}

void Relatives::AddOnePair(RelativePair pair)
{
    relatives.push_back(pair);
    numPairs++;
}

void Relatives::ShowSummary()
{
    assert(numPairs == relatives.size());
    
    int totGenoRelTypes = 6;
    int numIdMiss = 0, numPoMiss = 0, numMissDp = 0, numMissPo = 0;
    
    if (numPairs > 0) {
        int numRelPairs = 0;
        int relTypePairs[totGenoRelTypes];
        
        for (int i = 0; i < totGenoRelTypes; i++) relTypePairs[i] = 0;
        
        for (int i = 0; i < numPairs; i++) {
            RelativePair pair = relatives[i];
            int genoRelCode = pair.genoRelation;
            int pedRelCode = pair.pedRelation;
            if (genoRelCode < 5) {
                relTypePairs[genoRelCode]++;
                numRelPairs++;
            }
            
            if      (genoRelCode == 0 && pedRelCode > 2)   numIdMiss++;                    /// G = ID; P <> DP, MT
            else if (genoRelCode == 1 && pedRelCode > 0 && pedRelCode != 3)  numPoMiss++;  /// G = PO; P <> PO
            else if (pedRelCode  == 1 && genoRelCode != 0)  numMissDp++;                   /// P = DP; G <> ID
            else if (pedRelCode  == 3 && genoRelCode != 1)  numMissPo++;                   /// P = PO; G <> PO
        }
        
        cout << "Total " << numRelPairs << " pairs of related samples or subjects found.\n";;
        for (int i = 0; i < totGenoRelTypes-1; i++) {
            printf("    %s: %6d\n", genoRelTypes[i].c_str(), relTypePairs[i]);
        }
        
        if (numIdMiss + numPoMiss + numMissDp + numMissPo > 0) {
            cout << "\nERROR:\n";
            if (numIdMiss > 0)
                cout << "    " << numIdMiss << " pairs of samples with identical genotypes are reported as other relationships in the pedigree.\n";
                if (numPoMiss > 0)
                cout << "    " << numPoMiss << " parent/offspring pairs are reported as other relationships in the pedigree.\n";
                if (numMissDp > 0)
                cout << "    " << numMissDp << " pairs of reported duplicated samples don't have same genotypes.\n";
                if (numMissPo > 0)
                cout << "    " << numMissPo << " parent/offspring pairs reported in the pedigree cannot be confirmed by genotypes.\n";
            }
            cout << "\n";
        }
        else {
        cout << "No relatives found.\n";
    }
}

bool Relatives::SaveRelativesToFile(char* filename, int numPhs)
{
    assert(numPairs == relatives.size());
    
    vector<RelMissRate> idPairs;
    vector<RelMissRate> otherPairs;
    
    if (numPairs > 0) {
        int numRelPairs = 0, numIdPairs = 0, numOthPairs = 0;
        
        for (int i = 0; i < numPairs; i++) {
            RelativePair pair = relatives[i];
            
            int genoRelCode = pair.genoRelation;
            if (genoRelCode < 5) {
                if (genoRelCode == 0) {
                    RelMissRate idMr;
                    idMr.pairNo = i;
                    idMr.mr = pair.agmr;
                    idPairs.push_back(idMr);
                    numIdPairs++;
                }
                else {
                    RelMissRate othMr;
                    othMr.pairNo = i;
                    othMr.mr = pair.hgmr;
                    
                    otherPairs.push_back(othMr);
                    numOthPairs++;
                }
                
                numRelPairs++;
            }
        }
        
        if (numRelPairs > 0) {
            FILE *ofp = fopen(filename, "w");
            if (!ofp) {
                cout << "ERROR: failed to open " << filename << " for writing!\n\n";
                return false;
            }
            
            fprintf(ofp, "# HG match: number of SNPs with matched genotypes when only homozygous SNPs are counted\n");
            fprintf(ofp, "# HG miss: number of SNPs with mismatched genotypes when only homozygous SNPs are counted\n");
            fprintf(ofp, "# HGMR: Homozygous Genotype Mismatch Rate (%)\n");
            fprintf(ofp, "# AG match: number of SNPs with matched genotypes when all SNPs are counted\n");
            fprintf(ofp, "# AG miss: number of SNPs with mismatched genotypes when all SNPs are counted\n");
            fprintf(ofp, "# AGMR: All Genotype Mismatch Rate (%)\n");
            fprintf(ofp, "# Geno rel: relationship determined by sample genotypes. ID = identical subjects or twins; PO = parent-offspring; FS = full sibling; D2 = second degree relative; D3 = third degree relative; UN = unrelated or remotely related subjects\n");
            fprintf(ofp, "# Ped rel: Relationship derived from pedigree file. DP = duplicated samples; MT = monozygotic twins; PO = parent-offspring; FS = full sibling; SB = sibling with one parent missing in pedigree; HS = half sibling; GP = grandparent-grandchild; AV = avuncular; HA = half avuncular; FC = first cousin; HC = half first cousin; HF = half sibling + first cousin (half siblings whose fathers or mothers are full siblings)");
            fprintf(ofp, "# p_value: Probability that the genetic relationship is NOT the predicted type\n");
            fprintf(ofp, "#\n");
            fprintf(ofp, "# Expected HGMR and AGMR values (%) for different types of relationships:\n");
            fprintf(ofp, "# FS HGMR: %6.4f\n", fsMux);
            fprintf(ofp, "# D2 HGMR: %6.4f\n", d2Mux);
            fprintf(ofp, "# D3 HGMR: %6.4f\n", d3Mux);
            fprintf(ofp, "# PO AGMR: %6.4f\n", poMuy);
            fprintf(ofp, "# FS AGMR: %6.4f\n", fsMuy);
            fprintf(ofp, "# D2 AGMR: %6.4f\n", d2Muy);
            fprintf(ofp, "# D3 AGMR: %6.4f\n", d3Muy);
            
            if (numPhs > 1)
                fprintf(ofp, "DS1\tsample1\tDS2\tsample2\tsubject1\tsubject2\tsex1\tsex2\tHG match\tHG miss\tHGMR\tAG match\tAG miss\tAGMR\tgeno relation\tped relation\tp_value\n");
                else
                fprintf(ofp, "sample1\tsample2\tsubject1\tsubject2\tsex1\tsex2\tHG match\tHG miss\tHGMR\tAG match\tAG miss\tAGMR\tgeno relation\tped relation\tp_value\n");
                
                cout << "Total " << numRelPairs << " pairs of related samples or subjects found.\n";
                if (numIdPairs > 0) {
                RelMissRate* idArray = (RelMissRate*)malloc(sizeof(struct RelMissRate) * numIdPairs);
                for (int i = 0; i < numIdPairs; i++) {
                    idArray[i].pairNo = idPairs[i].pairNo;
                    idArray[i].mr     = idPairs[i].mr;
                }
                
                qsort(idArray, numIdPairs, sizeof(struct RelMissRate), ComparePairMissRate);
                
                for (int i = 0; i < numIdPairs; i++) {
                    int pairNo = idArray[i].pairNo;
                    RelativePair pair = relatives[pairNo];
                    string genoRel = genoRelTypes[pair.genoRelation];
                    string pedRel  = pedRelTypes[pair.pedRelation];
                    
                    char pvalStr[8];
                    if (pair.pvalue > -0.5) {
                        if (pair.pvalue < 0.001) {
                            sprintf(pvalStr, "%2.1e", pair.pvalue);
                        }
                        else {
                            sprintf(pvalStr, "%5.4f", pair.pvalue);
                        }
                    }
                    else {
                        pvalStr[0] = 0;
                    }
                    
                    if (numPhs > 1)
                        fprintf(ofp, "%d\t%s\t%d\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%4.2f\t%d\t%d\t%4.2f\t%s\t%s\t%s\n",
                        pair.sample1.phs, pair.sample1.smpName, pair.sample2.phs, pair.sample2.smpName,
                        pair.sample1.sbjName, pair.sample2.sbjName,
                        pair.sample1.sex,     pair.sample2.sex,
                        pair.homoSame, pair.homoDiff, pair.hgmr,
                        pair.allSame, pair.allDiff, pair.agmr, genoRel.c_str(), pedRel.c_str(), pvalStr);
                        else
                        fprintf(ofp, "%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%4.2f\t%d\t%d\t%4.2f\t%s\t%s\t%s\n",
                        pair.sample1.smpName, pair.sample2.smpName,
                        pair.sample1.sbjName, pair.sample2.sbjName,
                        pair.sample1.sex,     pair.sample2.sex,
                        pair.homoSame, pair.homoDiff, pair.hgmr,
                        pair.allSame, pair.allDiff, pair.agmr, genoRel.c_str(), pedRel.c_str(), pvalStr);
                        
                    }
                    delete idArray;
                }
                
                if (numOthPairs > 0) {
                RelMissRate* otherArray = (RelMissRate*)malloc(sizeof(struct RelMissRate) * numOthPairs);
                for (int i = 0; i < numOthPairs; i++) {
                    otherArray[i].pairNo = otherPairs[i].pairNo;
                    otherArray[i].mr = otherPairs[i].mr;
                }
                
                qsort(otherArray, numOthPairs, sizeof(struct RelMissRate), ComparePairMissRate);
                
                for (int i = 0; i < numOthPairs; i++) {
                    int pairNo = otherArray[i].pairNo;
                    RelativePair pair = relatives[pairNo];
                    string genoRel = genoRelTypes[pair.genoRelation];
                    string pedRel  = pedRelTypes[pair.pedRelation];
                    
                    char pvalStr[8];
                    if (pair.pvalue > -0.5) {
                        if (pair.pvalue < 0.001) {
                            sprintf(pvalStr, "%2.1e", pair.pvalue);
                        }
                        else {
                            sprintf(pvalStr, "%5.4f", pair.pvalue);
                        }
                    }
                    else {
                        pvalStr[0] = 0;
                    }
                    
                    if (numPhs > 1)
                        fprintf(ofp, "%d\t%s\t%d\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%4.2f\t%d\t%d\t%4.2f\t%s\t%s\t%s\n",
                        pair.sample1.phs, pair.sample1.smpName, pair.sample2.phs, pair.sample2.smpName,
                        pair.sample1.sbjName, pair.sample2.sbjName,
                        pair.sample1.sex,     pair.sample2.sex,
                        pair.homoSame, pair.homoDiff, pair.hgmr,
                        pair.allSame, pair.allDiff, pair.agmr, genoRel.c_str(), pedRel.c_str(), pvalStr);
                        else
                        fprintf(ofp, "%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%4.2f\t%d\t%d\t%4.2f\t%s\t%s\t%s\n",
                        pair.sample1.smpName, pair.sample2.smpName,
                        pair.sample1.sbjName, pair.sample2.sbjName,
                        pair.sample1.sex,     pair.sample2.sex,
                        pair.homoSame, pair.homoDiff, pair.hgmr,
                        pair.allSame, pair.allDiff, pair.agmr, genoRel.c_str(), pedRel.c_str(), pvalStr);
                    }
                    delete otherArray;
                }
                
                fclose(ofp);
                cout << "Saved " << numRelPairs << " pairs of relatives to " << filename << "\n";
            }
            else {
            cout << "No relative pairs found.\n";
        }
    }
    
    return true;
}
