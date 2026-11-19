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

bool FileExists (const char* filename)
{
    FILE* fp = fopen(filename, "r");
    if (fp) {
        fclose(fp);
        return true;
    }
    else {
        return false;
    }
}

bool FileWriteable (const char* filename)
{
    bool pre_exists = FileExists(filename);

    FILE* fp = fopen(filename, "a");

    if (fp) {
        fclose(fp);
	// If a file doesn't pre-exits, fopen creates new file.
	// it is removed here in order not to add extra files to user's directory.
	if (!pre_exists) remove (filename);
        return true;
    }
    else {
        return false;
    }
}

void ShowTimeDiff(struct timeval t1, struct timeval t2)
{
    int usec = t2.tv_usec - t1.tv_usec;
    int sec = t2.tv_sec - t1.tv_sec;
    if (usec < 0) {
        usec += 1000000;
        sec += 1;
    }
    int min = 0;
    if (sec > 60) {
        min = sec / 60;
        sec = sec % 60;
    }
    int hour = 0;
    if (min > 60) {
        hour = min / 60;
        min = min % 60;
    }

    printf("Time used: ");
    if (hour > 0) { printf("%d hours ", hour); }
    if (min > 0) { printf("%d minutes ", min); }
    if (sec > 0) { printf("%d seconds ", sec); }
    printf("%d microseconds\n\n", usec);
}

void ShowNumberInBits (unsigned long n, int len)
{
    char bitStr[len+1];
    int i;
    for (i = 0; i < len; i++) {
        bitStr[i] = '0';
    }
    bitStr[len] = 0;

    int pos = len-1;
    while(n > 0) {
        int r = n % 2;
        n = n >> 1;
        if (r) bitStr[pos] = '1';
        pos--;
    }

    printf("%s\n", bitStr);
}
