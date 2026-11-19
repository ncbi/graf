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

#include <fstream>
#include "Graf.h"

void SampleGenoDist::SetDistWithDist(GenoDist *p, GenoDist *d) {
    p->e = d->e;
    p->f = d->f;
    p->a = d->a;
}

void SampleGenoDist::SetPointWithDist(Point *p, GenoDist *d) {
    p->x = d->a;
    p->y = d->e;
    p->z = d->f;
}

void SampleGenoDist::SetPointWithPoint(Point *p, Point *s) {
    p->x = s->x;
    p->y = s->y;
    p->z = s->z;
}

SampleGenoDist::SampleGenoDist(GenoDist *ec, GenoDist *fc, GenoDist *ac, GenoDist *sp) {
    SetDistWithDist(&eDist, ec);
    SetDistWithDist(&fDist, fc);
    SetDistWithDist(&aDist, ac);
    SetDistWithDist(&sDist, sp);
}

void SampleGenoDist::MovePoint(Point *p, double dx, double dy, double dz)
{
    p->x += dx;
    p->y += dy;
    p->z += dz;
}

void SampleGenoDist::RotatePointOnx(Point *p, double theta)
{
    theta = theta * pi / 180;
    double y = p->y;
    double z = p->z;
    p->y = y * cos(theta) - z * sin(theta);
    p->z = y * sin(theta) + z * cos(theta);
}

void SampleGenoDist::RotatePointOny(Point *p, double theta)
{
    theta = theta * pi / 180;
    double x = p->x;
    double z = p->z;
    p->x = x * cos(theta) + z * sin(theta);
    p->z = x * -1 * sin(theta) + z * cos(theta);
}

void SampleGenoDist::RotatePointOnz(Point *p, double theta)
{
    theta = theta * pi / 180;
    double x = p->x;
    double y = p->y;
    p->x = x * cos(theta) - y * sin(theta);
    p->y = x * sin(theta) + y * cos(theta);
}

void SampleGenoDist::TransformAllDists()
{
    SetPointWithDist(&ePt, &eDist);
    SetPointWithDist(&fPt, &fDist);
    SetPointWithDist(&aPt, &aDist);
    SetPointWithDist(&sPt, &sDist);
    
    /// Move points so that the f point is at the origin
    double dx = fPt.x;
    double dy = fPt.y;
    double dz = fPt.z;
    
    MovePoint(&ePt, -dx, -dy, -dz);
    MovePoint(&fPt, -dx, -dy, -dz);
    MovePoint(&aPt, -dx, -dy, -dz);
    MovePoint(&sPt, -dx, -dy, -dz);
    
    /// Rotate on z-axis
    double theta1 = atan2(aPt.y, aPt.x) * -180 / pi;
    RotatePointOnz(&ePt, theta1);
    RotatePointOnz(&fPt, theta1);
    RotatePointOnz(&aPt, theta1);
    RotatePointOnz(&sPt, theta1);
    
    /// Rotate on y-axis
    double theta2 = atan2(aPt.z, aPt.x) * 180 / pi;
    RotatePointOny(&ePt, theta2);
    RotatePointOny(&fPt, theta2);
    RotatePointOny(&aPt, theta2);
    RotatePointOny(&sPt, theta2);
    
    /// Rotate on y-axis
    double theta3 = atan2(ePt.y, ePt.z) * 180 / pi;
    theta3 -= 90;
    RotatePointOnx(&ePt, theta3);
    RotatePointOnx(&fPt, theta3);
    RotatePointOnx(&aPt, theta3);
    RotatePointOnx(&sPt, theta3);
    
    /// Move points back
    MovePoint(&ePt, dx, dy, dz);
    MovePoint(&fPt, dx, dy, dz);
    MovePoint(&aPt, dx, dy, dz);
    MovePoint(&sPt, dx, dy, dz);
    
    int debug = 0;
    if (debug) {
        printf("E: x: %6.4f  y: %6.4f  z: %6.4f\n", ePt.x, ePt.y, ePt.z);
        printf("F: x: %6.4f  y: %6.4f  z: %6.4f\n", fPt.x, fPt.y, fPt.z);
        printf("A: x: %6.4f  y: %6.4f  z: %6.4f\n", aPt.x, aPt.y, aPt.z);
        printf("J: x: %6.4f  y: %6.4f  z: %6.4f\n", sPt.x, sPt.y, sPt.z);
        
        assert(0);
    }
}

void SampleGenoDist::CalculateBaryCenters()
{
    double x1 = ePt.x;
    double y1 = ePt.y;
    double x2 = fPt.x;
    double y2 = fPt.y;
    double x3 = aPt.x;
    double y3 = aPt.y;
    double xp = sPt.x;
    double yp = sPt.y;
    
    double det = (y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3);
    
    eWt = ((y2 - y3)*(xp - x3) + (x3 - x2)*(yp - y3))/det;
    fWt = ((y3 - y1)*(xp - x3) + (x1 - x3)*(yp - y3))/det;
    aWt = 1 - eWt - fWt;
    
    if (0) {
        cout << "Weights\n";
        cout << "\tWhite: " << eWt << "\n";
        cout << "\tBlack: " << fWt << "\n";
        cout << "\tAsian: " << aWt << "\n";
    }
}

void SampleGenoDist::ShowPositions(int showSmp)
{
    printf("\tE: %6.4f  %6.4f  %6.4f\n", ePt.x, ePt.y, ePt.z);
    printf("\tF: %6.4f  %6.4f  %6.4f\n", fPt.x, fPt.y, fPt.z);
    printf("\tA: %6.4f  %6.4f  %6.4f\n", aPt.x, aPt.y, aPt.z);
    
    if (showSmp) {
        printf("\tSMP: %6.4f  %6.4f  %6.4f\n", sPt.x, sPt.y, sPt.z);
    }
    printf("\n");
}
