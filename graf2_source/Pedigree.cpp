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
#include <Graf.h>

Pedigree::Pedigree()
{
    numGenoSmps = 0;
    numSsmSmps = 0;
    numPedSbjs = 0;
    
    ssmIsGap = false;

    string pedTypes[]  = {"", "DP", "MT", "PO", "FS", "SB", "HS", "GP", "AV", "FC", "HF", "HA", "HC"};
    
    for (int i = 0; i < 13;  i++) pedRelTypes[i] = pedTypes[i];
    
    genoSmps.clear();
    ssm.clear();
    pedSbjs.clear();
    sbjSexes.clear();
    sbjDads.clear();
    sbjMoms.clear();
    sbjTwins.clear();
    smpRels.clear();
    pedSbjList.clear();
}

Pedigree::Pedigree(SampleInfo* smps, int numSmps)
{
    numGenoSmps = numSmps;
    numSsmSmps = 0;
    numPedSbjs = 0;
    
    string pedTypes[]  = {"", "DP", "MT", "PO", "FS", "SB", "HS", "GP", "AV", "FC", "HF", "HA", "HC"};
    for (int i = 0; i < 13;  i++) pedRelTypes[i] = pedTypes[i];
    
    ssm.clear();
    genoSmps.clear();
    ssm.clear();
    pedSbjs.clear();
    sbjSexes.clear();
    sbjDads.clear();
    sbjMoms.clear();
    sbjTwins.clear();
    smpRels.clear();
    pedSbjList.clear();
    
    for (int i = 0; i < numSmps; i++) {
        string smp = smps[i].smpName;
        string sbj = smp;
        ssm[smp] = sbj;
        genoSmps.push_back(smp);
    }
}

void Pedigree::ShowSummary(int numShowPairs)
{
    int maxShowSmps = 10;
    
    cout << "SSM and pedigree summary:\n";
    cout << "\tTotal " << numSsmSmps << " samples from SSM file have genotypes.\n";
    map<string, string>::iterator iter;
    if (numShowPairs > 0) {
        int smpNo = 0;
        for (iter = ssm.begin(); iter != ssm.end(); iter++) {
            if (smpNo < numShowPairs) {
                cout << "\t\t" << iter->first << " => " << iter->second << "\n";
            }
            else {
                cout << "\t\t...\n";
                break;
            }
            smpNo++;
        }
    }
    
    int numSbjDads = sbjDads.size();
    int numSbjMoms = sbjDads.size();
    int males = 0, females = 0, unSexes = 0;
    int numTwinSbjs = sbjTwins.size();
    
    map<string, int>::iterator pedIter;
    for (pedIter = sbjSexes.begin(); pedIter != sbjSexes.end(); pedIter++) {
        if      (pedIter->second == 1) males++;
        else if (pedIter->second == 2) females++;
        else                           unSexes++;
    }
    
    cout << "\tTotal " << numPedSbjs << " subjects from pedigree file.\n";
    printf("\t\t%5d subjects have father IDs in pedigree.\n", numSbjDads);
    printf("\t\t%5d subjects have mother IDs in pedigree.\n", numSbjMoms);
    printf("\t\t%5d subjects have twin IDs in pedigree.\n", numTwinSbjs);
    printf("\t\t%5d subjects are males.\n", males);
    printf("\t\t%5d subjects are females.\n", females);
    printf("\t\t%5d subjects have unknown sex values.\n", unSexes);
    
    int numRelTypes = 13;
    int relTypeCnts[numRelTypes];
    for (int i = 0; i < numRelTypes; i++) relTypeCnts[i] = 0;
    
    for (pedIter = smpRels.begin(); pedIter != smpRels.end(); pedIter++) {
        int type = pedIter->second;
        assert(type < numRelTypes);
        if (type > 0) relTypeCnts[type]++;
    }
    
    int totRels = smpRels.size();
    cout << "\tTotal " << totRels << " pairs of relatives.\n";
    for (int i = 0; i < numRelTypes; i++) {
        string relation = pedRelTypes[i];
        int cnt = relTypeCnts[i];
        if (cnt > 0) printf("\t\t%2d (%s): %5d\n", i, relation.c_str(), cnt);
        
        if (cnt > 0 && numShowPairs > 0) {
            int pairNo = 0;
            for (pedIter = smpRels.begin(); pedIter != smpRels.end(); pedIter++) {
                int type = pedIter->second;
                if (type == i) {
                    if (pairNo < numShowPairs) {
                        cout << "\t\t\t" << pedIter->first << "\n";
                    }
                    else {
                        cout << "\t\t\t...\n";
                        break;
                    }
                    pairNo++;
                }
            }
        }
    }
    cout << "\n";
}

void Pedigree::StringToUpper(char* str)
{
    int len = strlen(str);
    for (int i = 0; i < len; i++) {
        str[i] = toupper(str[i]);
    }
}

bool Pedigree::ReadSmpSbjMapping(char* ssmFile)
{
    bool flag = true;
    
    int nameLen = 80;
    int lineLen = 300;
    char fileLine[lineLen];
    int fileType = 0;
    
    int totRels = 0;
    map<string, string>::iterator iter;
    
    /// Store whether a geno sample is found in SSM or not
    map<string, bool> smpInSsm;
    for (int i = 0; i < numGenoSmps; i++) {
        string smp = genoSmps[i];
        smpInSsm[smp] = false;
    }
    map<string, bool>::iterator genoIt;
    vector<string> ssmOnlySmps;         /// Samples found in SSM file but not in genotype list
    ssm.clear();
    
    cout << "Reading subject sample mapping from " << ssmFile << "\n";

    FILE *ifp = fopen(ssmFile, "r");
    if(ifp) {
        fgets(fileLine, lineLen, ifp);
        char var1[nameLen];
        char var2[nameLen];
        sscanf(fileLine, "%s %s", var1, var2);
        StringToUpper(var1);
        StringToUpper(var2);
        
	if (ssmIsGap) {
            // This is default order specified in submission guide:
            // https://www.ncbi.nlm.nih.gov/gap/docs/submissionguide/#10-how-do-i-create-subject-sampl
	    fileType = 2;
	}
	else {
            if (strstr(var1, "SAMP") && strstr(var2, "SUBJ")) {
                fileType = 1;
            }
            else if (strstr(var2, "SAMP") && strstr(var1, "SUBJ")) {
                fileType = 2;
            }
	}
        
        if (fileType < 1) {
	    cout << "ERROR: Invalid SSM file. The first two columns should be SAMPLE_ID and SUBJECT_ID.\n";
            flag = false;
            return flag;
        }
        
        char val1[nameLen];
        char val2[nameLen];
        char* smp;
        char* sbj;
        while (fgets(fileLine, lineLen, ifp) != NULL) {
            sscanf(fileLine, "%s %s", val1, val2);
            if (fileType == 1) {
                smp = val1;
                sbj = val2;
            }
            else {
                smp = val2;
                sbj = val1;
            }
            
            genoIt = smpInSsm.find(smp);     /// Check if the sample is in genotyped sample list
            if (genoIt == smpInSsm.end()) {
                ssmOnlySmps.push_back(smp);
            }
            else {
                iter = ssm.find(smp);
                if (iter == ssm.end()) {
                    ssm[smp] = sbj;
                    smpInSsm[smp] = true;
                    numSsmSmps++;
                }
                else {
                    if (strcmp(sbj, iter->second.c_str()) != 0) {
                        cout << "ERROR in " << ssmFile << ": Sample " << smp << " maps to two different subjects: "
                        << iter->second << " and " << sbj << "\n";
                        numSsmSmps = 0;
                        flag = false;
                        
                        break;
                    }
                }
            }
        }
        
        fclose(ifp);
    }
    else {
        cout << "ERROR: Couldn't open file " << ssmFile << "!\n";
    }
    
    vector<string> genoOnlySmps;
    for (genoIt = smpInSsm.begin(); genoIt != smpInSsm.end(); genoIt++) {
        if (genoIt->second == false) genoOnlySmps.push_back(genoIt->first);
    }
    
    int numSsmOnlySmps = ssmOnlySmps.size();
    int numGenoOnlySmps = genoOnlySmps.size();
    
    if (numSsmOnlySmps > 0) cout  << "\t******* NOTE ****** " << numSsmOnlySmps
        << " samples in SSM file are not found in the genotype dataset.\n";
        if (numGenoOnlySmps > 0) cout << "\t***** WARNING ***** " << numGenoOnlySmps
        << " genotyped samples are not found in SSM file. "
        << "Subject IDs are set to the same as sample IDs for these samples.\n";
        for (int i = 0; i < numGenoOnlySmps; i++) {
        string smp = genoOnlySmps[i];
        ssm[smp] = smp;
    }
    
    cout << "\tFound " << numSsmSmps << " genotyped samples in SSM file.\n\n";
    
    return flag;
}

int Pedigree::GetPedLineValue(char* line, int linePos, char* val)
{
    val[0] = 0;
    
    int valPos = 0;
    while (line[linePos] && line[linePos] != '\t') {
        /// Ignore spaces
        if (line[linePos] != ' ' && line[linePos] != '\n' && line[linePos] != '"' && line[linePos] != 39) {
            val[valPos] = line[linePos];
            valPos++;
        }
        linePos++;
    }
    
    val[valPos] = 0;
    
    if (line[linePos]) linePos++;
    
    return linePos;
}

bool Pedigree::ReadPedigree(char* pedFile, bool isPlink)
{
    int nameLen = 80;
    int lineLen = 300;
    char famLine[lineLen];
    int numFamSmps = 0;
    bool flag = false;
    
    FILE *ifp = fopen(pedFile, "r");
    if(ifp) {
        char fam[nameLen];
        char sbj[nameLen];
        char dad[nameLen];
        char mom[nameLen];
        char sexStr[nameLen];
        char twin[nameLen];
        int pheno;
        int sex;
        bool hasTwin = 0;
        bool swapDadMom = false;
        
        if (!isPlink) {
            fgets(famLine, lineLen, ifp);
            sscanf(famLine, "%s %s %s %s %s %s", fam, sbj, dad, mom, sexStr, twin);
            int linePos = 0;
            linePos = GetPedLineValue(famLine, linePos, fam);
            linePos = GetPedLineValue(famLine, linePos, sbj);
            linePos = GetPedLineValue(famLine, linePos, dad);
            linePos = GetPedLineValue(famLine, linePos, mom);
            linePos = GetPedLineValue(famLine, linePos, sexStr);
            linePos = GetPedLineValue(famLine, linePos, twin);
            StringToUpper(fam);
            StringToUpper(sbj);
            StringToUpper(dad);
            StringToUpper(mom);
            StringToUpper(sexStr);
            StringToUpper(twin);
            
            if (!(strstr(fam, "FAM") &&
                (strstr(dad, "FATHER") || strstr(mom, "FATHER")) &&
                (strstr(mom, "MOTHER") || strstr(dad, "MOTHER")) &&
                strstr(sexStr, "SEX") ) ) {
                cout << "ERROR: Invalid pedigree file " << pedFile << ".  Expected the following columns:\n"
                << "\tFamily_ID\n"
                << "\tSubject_ID\n"
                << "\tFather\n"
                << "\tMother\n"
                << "\tSex\n"
                << "\tMZ_Twin_ID (optional)\n";
                
                cout << "Fam = " << fam << "\n";
                cout << "Sbj = " << sbj << "\n";
                cout << "Dad = " << dad << "\n";
                cout << "Mom = " << mom << "\n";
                cout << "Sex = " << sexStr << "\n";
                cout << "Twin = " << twin << "\n";
                
                return flag;
            }
            
            if (strstr(twin, "TWIN")) hasTwin = 1;
            if (strstr(dad, "MOTHER")) swapDadMom = 1;
        }
        
        cout << "Reading pedigree from " << pedFile << "\n";
       
        
        while (fgets(famLine, lineLen, ifp) != NULL) {
            fam[0] = 0;
            sbj[0] = 0;
            dad[0] = 0;
            mom[0] = 0;
            sexStr[0] = 0;
            twin[0] = 0;
            
            /// Trime the line breakers
            int len = strlen(famLine);
            if (len > 1 && famLine[len-1] == 10) {
                if (famLine[len-2] == 13) { /// CRLF
                    famLine[len-2] = 0;
                }
                else {                      /// New line
                    famLine[len-1] = 0;
                }
            }
            
            if (isPlink) {
                sscanf(famLine, "%s %s %s %s %d %d", fam, sbj, dad, mom, &sex, &pheno);
            }
            else {
                int linePos = 0;
                linePos = GetPedLineValue(famLine, linePos, fam);
                linePos = GetPedLineValue(famLine, linePos, sbj);
                if (swapDadMom) {
                    linePos = GetPedLineValue(famLine, linePos, mom);
                    linePos = GetPedLineValue(famLine, linePos, dad);
                }
                else {
                    linePos = GetPedLineValue(famLine, linePos, dad);
                    linePos = GetPedLineValue(famLine, linePos, mom);
                }
                
                linePos = GetPedLineValue(famLine, linePos, sexStr);
                linePos = GetPedLineValue(famLine, linePos, twin);
                
                sex = 0;
                char sexChar = toupper(sexStr[0]);
                if      (sexChar == '1' || sexChar == 'M') sex = 1;
                else if (sexChar == '2' || sexChar == 'F') sex = 2;
            }
            
            pedSbjs[sbj] = 1;
            
            if (dad[0] && dad[0] != '0' && dad[0] != ' ') {
                sbjDads[sbj] = dad;
                pedSbjs[dad] = 1;
            }
            
            if (mom[0] && mom[0] != '0' && mom[0] != ' ') {
                sbjMoms[sbj] = mom;
                pedSbjs[mom] = 1;
            }
            
            sbjSexes[sbj] = sex;
            if (hasTwin && twin[0] && twin[0] != '0' && twin[0] != ' ') {
                sbjTwins[sbj] = twin;
            }
        }
        fclose(ifp);
        
        map<string, int>::iterator sbjIt;
        for (sbjIt = pedSbjs.begin(); sbjIt != pedSbjs.end(); sbjIt++) {
            string sbj = sbjIt->first;
            pedSbjList.push_back(sbj);
        }
        
        numPedSbjs = pedSbjs.size();
        
        flag = 1;
    }
    else {
        cout << "ERROR: couldn't open pedigree file " << pedFile << "!\n";
        flag = 0;
    }
    
    return flag;
}

string Pedigree::FindAncestor(map<string, string>* sbjAns, string sbj)
{
    string ans = "";
    map<string, string>::iterator it = sbjAns->find(sbj);
    if (it != sbjAns->end()) ans = (*sbjAns)[sbj];
    
    return ans;
}

int Pedigree::GetPedigreeRelatives()
{
    int i, j;
    int totRels = 0;
    
    map<string, string> sbjYayas;  /// Father's father, or yaya in Chinese
    map<string, string> sbjNanas;  /// Father's mother, or nana in Chinese
    map<string, string> sbjWagos;  /// Mother's father, or waigon in Chinese
    map<string, string> sbjWapos;  /// Mother's mother, or waipo in Chinese
    
    vector<string> pedChildren;    /// Only need to check the subjects with parents
    
    map<string, int>::iterator sbjIt;
    map<string, int>::iterator pairIt;
    for (sbjIt = pedSbjs.begin(); sbjIt != pedSbjs.end(); sbjIt++) {
        string sbj = sbjIt->first;
        string dad = FindAncestor(&sbjDads, sbj);
        string mom = FindAncestor(&sbjMoms, sbj);
        
        if (dad != "" || mom != "") pedChildren.push_back(sbj);
        
        if (dad != "") {
            string  po = sbj + "\t" + dad;
            sbjRels[po] = 3;  /// 3 = PO
            
            string yaya = FindAncestor(&sbjDads, dad);
            string nana = FindAncestor(&sbjMoms, dad);
            
            if (yaya != "") {
                string  gp = sbj + "\t" + yaya;
                sbjRels[gp] = 7;  /// 7 = GP
                sbjYayas[sbj] = yaya;
            }
            
            if (nana != "") {
                string  gp = sbj + "\t" + nana;
                sbjRels[gp] = 7;  /// 7 = GP
                sbjNanas[sbj] = nana;
            }
        }
        
        if (mom != "") {
            string  po = sbj + "\t" + mom;
            sbjRels[po] = 3;  /// 3 = PO
            
            string wago = FindAncestor(&sbjDads, mom);
            string wapo = FindAncestor(&sbjMoms, mom);
            
            if (wago != "") {
                string  gp = sbj + "\t" + wago;
                sbjRels[gp] = 7;  /// 7 = GP
                sbjWagos[sbj] = wago;
            }
            
            if (wapo != "") {
                string  gp = sbj + "\t" + wapo;
                sbjRels[gp] = 7;  /// 7 = GP
                sbjWapos[sbj] = wapo;
            }
        }
    }
    
    int numRelPairs = 0;
    int numChildren = pedChildren.size();
    int totSbjPairs = numChildren * (numChildren - 1) / 2;
    cout << "\tTotal " << numPedSbjs << " subjects in pedigree file. " << numChildren << " subjects have parent IDs.\n\n";
    cout << "Checking " << totSbjPairs << " pairs of subjects in pedigree for reported relationships.\n";
    int pairNo = 0;
    for (i = 0; i < numChildren-1; i++) {
        string sbj1  = pedChildren[i];
        
        string  dad1 = FindAncestor(&sbjDads,  sbj1);
        string  mom1 = FindAncestor(&sbjMoms,  sbj1);
        string yaya1 = FindAncestor(&sbjYayas, sbj1);
        string nana1 = FindAncestor(&sbjNanas, sbj1);
        string wago1 = FindAncestor(&sbjWagos, sbj1);
        string wapo1 = FindAncestor(&sbjWapos, sbj1);
        
        for (j = i+1; j < numChildren; j++) {
            string sbj2  = pedChildren[j];
            
            string  dad2 = FindAncestor(&sbjDads,  sbj2);
            string  mom2 = FindAncestor(&sbjMoms,  sbj2);
            string yaya2 = FindAncestor(&sbjYayas, sbj2);
            string nana2 = FindAncestor(&sbjNanas, sbj2);
            string wago2 = FindAncestor(&sbjWagos, sbj2);
            string wapo2 = FindAncestor(&sbjWapos, sbj2);
            
            string pair1 = sbj1 + "\t" + sbj2;
            string pair2 = sbj2 + "\t" + sbj1;
            
            int relation = 0;
            
            if (dad1 != "" && dad2 != "" && mom1 != "" && mom2 != "") {
                if (dad1 == dad2 && mom1 == mom2) {
                    relation = 4;    /// FS
                }
                else if (dad1 == dad2 && mom1 != mom2) { /// Check if moms are siblings
                    if ((wago1 != "" && wago1 == wago2) || (wapo1 != "" && wapo1 == wapo2)) {
                        relation = 10; /// HF
                    }
                    else {
                        relation = 6;  /// HS
                    }
                }
                else if (dad1 != dad2 && mom1 == mom2) { /// Check if dads are siblings
                    if ((yaya1 != "" && yaya1 == yaya2) || (nana1 != "" && nana1 == nana2)) {
                        relation = 10; /// HF
                    }
                    else {
                        relation = 6;  /// HS
                    }
                }
            }
            else if ( (dad1 != "" && dad2 != "" && dad1 == dad2) || (mom1 != "" && mom2 != "" && mom1 == mom2) ) {
                relation = 5;      /// SB
            }
            
            if (relation < 1) {
                if      (dad1 != "" && yaya2 != "" && dad1 == yaya2) {
                    if    (mom1 != "" && nana2 != "" && mom1 != nana2) {
                        relation = 11; /// HA
                    }
                    else {
                        relation = 8;  /// AV
                    }
                }
                else if (dad2 != "" && yaya1 != "" && dad2 == yaya1) {
                    if    (mom2 != "" && nana1 != "" && mom2 != nana1) {
                        relation = 11;
                    }
                    else {
                        relation = 8;
                    }
                }
                else if (mom1 != "" && nana2 != "" && mom1 == nana2) {
                    if    (dad1 != "" && yaya2 != "" && dad1 != yaya2) {
                        relation = 11;
                    }
                    else {
                        relation = 8;
                    }
                }
                else if (mom2 != "" && nana1 != "" && mom2 == nana1) {
                    if    (dad2 != "" && yaya1 != "" && dad2 != yaya1) {
                        relation = 11;
                    }
                    else {
                        relation = 8;
                    }
                }
                else if (dad1 != "" && wago2 != "" && dad1 == wago2) {
                    if    (mom1 != "" && wapo2 != "" && mom1 != wapo2) {
                        relation = 11;
                    }
                    else {
                        relation = 8;
                    }
                }
                else if (dad2 != "" && wago1 != "" && dad2 == wago1) {
                    if    (mom2 != "" && wapo1 != "" && mom2 != wapo1) {
                        relation = 11;
                    }
                    else {
                        relation = 8;
                    }
                }
                else if (mom1 != "" && wapo2 != "" && mom1 == wapo2) {
                    if    (dad1 != "" && wago2 != "" && dad1 != wago2) {
                        relation = 11;
                    }
                    else {
                        relation = 8;
                    }
                }
                else if (mom2 != "" && wapo1 != "" && mom2 == wapo1) {
                    if    (dad2 != "" && wago1 != "" && dad2 != wago1) {
                        relation = 11;
                    }
                    else {
                        relation = 8;
                    }
                }
                else if ( (yaya1 != "" && yaya2 != "" && yaya1 == yaya2) || (nana1 != "" && nana2 != "" && nana1 == nana2) ) {
                    relation = 9; /// FC
                }
                else if ( (wago1 != "" && wago2 != "" && wago1 == wago2) || (wapo1 != "" && wapo2 != "" && wapo1 == wapo2) ) {
                    relation = 9;
                }
                else if ( (wago1 != "" && yaya2 != "" && wago1 == yaya2) || (wapo1 != "" && nana2 != "" && wapo1 == nana2) ) {
                    relation = 9;
                }
                else if ( (yaya1 != "" && wago2 != "" && yaya1 == wago2) || (nana1 != "" && wapo2 != "" && nana1 == wapo2) ) {
                    relation = 9;
                }
            }
            
            if (relation == 9) {
                if (  (yaya1 != "" && yaya1 == yaya2) && (nana1 != "" && nana2 != "" && nana1 != nana2) ||
                    (yaya1 != "" && yaya1 == wago2) && (nana1 != "" && wapo2 != "" && nana1 != wapo2) ||
                    (wago1 != "" && wago1 == yaya2) && (wapo1 != "" && nana2 != "" && wapo1 != nana2) ||
                    (wago1 != "" && wago1 == wago2) && (wapo1 != "" && wapo2 != "" && wapo1 != wapo2) ||
                    (nana1 != "" && nana1 == nana2) && (yaya1 != "" && yaya2 != "" && yaya1 != yaya2) ||
                    (nana1 != "" && nana1 == wapo2) && (yaya1 != "" && wago2 != "" && yaya1 != wago2) ||
                    (wapo1 != "" && wapo1 == nana2) && (wago1 != "" && yaya2 != "" && wago1 != yaya2) ||
                    (wapo1 != "" && wapo1 == wapo2) && (wago1 != "" && wago2 != "" && wago1 != wago2)  ) {
                    relation = 12; /// HC
                }
            }
            
            if (relation) {
                map<string, int>::iterator pairIt1 = sbjRels.find(pair1);
                map<string, int>::iterator pairIt2 = sbjRels.find(pair2);
                if (pairIt1 == sbjRels.end() && pairIt2 == sbjRels.end()) sbjRels[pair1] = relation;
            }
            
            pairNo++;
            if (pairNo % 1000000 == 0) {
                numRelPairs = sbjRels.size();
                cout << "\tChecked " << pairNo << " pairs.  Found " << numRelPairs << " pairs of related subjects.\n";
            }
        }
    }
    cout << "\tChecked " << pairNo << " pairs of subjects.  Found " << numRelPairs << " pairs of related subjects.\n";
    
    /// Find twins
    vector<string> twinSbjs;
    map<string, string>::iterator twinIt;
    for (twinIt = sbjTwins.begin(); twinIt != sbjTwins.end(); twinIt++) {
        string sbj = twinIt->first;
        twinSbjs.push_back(sbj);
    }
    int numTwinSbjs = twinSbjs.size();
    
    int numMzTwins = 0, numAllTwins = 0;
    int twinSexMis = 0, twinParMis = 0;
    vector<string> twinParConf;
    
    for (i = 0; i < numTwinSbjs-1; i++) {
        string sbj1 = twinSbjs[i];
        for (j = i+1; j < numTwinSbjs; j++) {
            string sbj2 = twinSbjs[j];
            
            if (sbjTwins[sbj1] == sbjTwins[sbj2]) {
                bool isMZ = true;
                
                if (sbjSexes[sbj1] && sbjSexes[sbj2] && sbjSexes[sbj1] != sbjSexes[sbj2]) {
                    twinSexMis++;
                    isMZ = false;
                }
                if (sbjDads[sbj1] != sbjDads[sbj2] || sbjMoms[sbj1] != sbjMoms[sbj2]) {
                    twinParMis++;
                    isMZ = false;
                    string err = "Twin " + sbj1 + "/" + sbj2 + ": parents " + sbjDads[sbj1] + "/" + sbjMoms[sbj1]
                    + " vs. " +    sbjDads[sbj2] + "/" + sbjMoms[sbj2];
                    twinParConf.push_back(err);
                }
                
                if (isMZ) {
                    string pair1 = sbj1 + "\t" + sbj2;
                    sbjRels[pair1] = 2; /// 2 = MT
                    numMzTwins++;
                }
                numAllTwins++;
            }
        }
    }
    
    if (twinSexMis > 0) cout << "\n********* NOTE ******** " << twinSexMis
        << " pairs of twins have different genders, indicating that they are not MZ twins\n\n";
        if (twinParMis > 0) cout << "ERROR: " << twinParMis << " pairs of twins have different parents!\n";
        for (i = 0; i < twinParMis; i++) {
        if (i < 10) {
            cout << "\t" << twinParConf[i] << "\n";
        }
        else {
            cout << "\t...\n";
            break;
        }
    }
    
    totRels = sbjRels.size();
    cout << "\n\tFound " << numMzTwins << " pairs of MZ twins. Total " << numTwinSbjs << " subjects have twin IDs.\n";
    if (numAllTwins > numMzTwins) cout << "Total " << numAllTwins << " pairs of twins in the pedigree file.\n";
    cout << "\nChecked " << numPedSbjs << " subjects in SSM and pedigree file. Found " << totRels << " relative pairs.\n";
    cout << "\n";
    
    /// Get sample pair relations
    vector<string> smpList;
    map<string, string>::iterator iter;
    int smpNo = 0;
    for (iter = ssm.begin(); iter != ssm.end(); iter++) {
        smpList.push_back(iter->first);
        smpNo++;
    }
    
    for (i = 0; i < smpNo-1; i++) {
        string smp1 = smpList[i];
        string sbj1 = ssm[smp1];
        for (j = i+1; j < smpNo; j++) {
            string smp2 = smpList[j];
            string sbj2 = ssm[smp2];
            
            string smpPair  = smp1 + "\t" + smp2;
            string sbjPair1 = sbj1 + "\t" + sbj2;
            string sbjPair2 = sbj2 + "\t" + sbj1;
            
            if (sbj1 == sbj2) {
                smpRels[smpPair] = 1;
            }
            else {
                int rel = 0;
                pairIt = sbjRels.find(sbjPair1);
                if (pairIt != sbjRels.end()) rel = pairIt->second;
                if (!rel) {
                    pairIt = sbjRels.find(sbjPair2);
                    if (pairIt != sbjRels.end()) rel = pairIt->second;
                }
                
                if (rel) smpRels[smpPair] = rel;
            }
        }
    }
    
    return totRels;
}

