#ifndef CAPTURE_LIGHTING_H
#define CAPTURE_LIGHTING_H

#include <bits/stdc++.h>
#include "1605055_objects.h"
#include "1605055_glut_display.h"
#include "1605055_bitmap_image.hpp"

using namespace std;

/*
 * CAPTURE VECTORS
 */

double deg2Rad(double deg) {
    return pi*deg/180;
}

vector<vector<point>>capturePts;

void updateCapturePts() {
    if(capturePts.empty()) {
        capturePts = vector<vector<point>>(pxCount, vector<point>(pxCount));
    }

    double scrDist = (screenHt/2.0) / tan(deg2Rad(fovY)/2.0);
    point scrMid = posCam + lVec*scrDist;
    double dx = (double)screenWd/pxCount, dy = (double)screenHt/pxCount;
    point scrCorner = scrMid - rVec*screenWd*0.5 - uVec*screenHt*0.5 + rVec*dx*0.5 + uVec*dy*0.5;

    for(int xPx = 0; xPx < pxCount; ++xPx) {
        for(int yPx = 0; yPx < pxCount; ++yPx) {
            capturePts[xPx][yPx] = scrCorner + (rVec * dx * xPx) + (uVec * dy * yPx);
        }
    }
}

/*
 * REFLECTION
 */

reflection reflect(const ray& rayNow) {  //warning: rayNow MUST be normalized
    if(!rayNow.isNormalized())
        cout << "WARNING: rayNow in reflect() not normalized" << endl;

    double c = inf;
    ray reflRay;
    object* reflObj;

    for(object *obj:things) {
        //if(obj->type == "curve")  //curve not shown in glut
            //continue;
        double cTmp = obj->intersect(rayNow);
        if(cTmp > 0 && cTmp < c) {
            reflRay = obj->reflect(rayNow);
            c = cTmp;
            reflObj = obj;
        }
    }

    return {reflRay, c, reflObj};
}

/*
 * LIGHTING
 */

bool allAmbient = true;

void rayTraceColor(const reflection& reflNow, color& colorNow, int level = 0) {
    //send first reflection here with "0" color

    if(reflNow.isInValid())
        return;

    if(level < recurLevel)
        rayTraceColor(reflect(reflNow.reflRay), colorNow, level+1);

    const object* reflObj = reflNow.reflObj;
    const vector<double>& coeff = reflObj->coeff;  //coEff: 0->ambient, 1->diffuse, 2->specular, 3->recursive

    const point& cutPt = reflNow.reflRay.origin;
    const color& cutColor = reflObj->getColor(cutPt);

    //recursive component
    colorNow *= coeff[3];

    //ambient component
    if(level == 0 || allAmbient)
        colorNow += cutColor*coeff[0];

    for(const light& lt: lights) {
        ray ltRay(cutPt, lt.pt, false, false);
        if(!reflect(ltRay).isInValid())  //shading
            continue;

        //diffusion component
        double cosD = abs(ltRay.dir.dot( reflObj->normal(cutPt) ));
        colorNow += lt.ltColor * coeff[1] * cosD * cutColor;
        //specular component
        double cosS = pow(max(ltRay.dir.dot(reflNow.reflRay.dir),0.0), reflObj->shine);
        colorNow += lt.ltColor * coeff[2] * cosS * cutColor;
    }
}

/*
 * BMP generation
 */

int bmpFileCount = 0;
string imgPath = "IO\\output";

void genImage() {
    updateCapturePts();
    bitmap_image image(pxCount,pxCount);

    for(int xPx = 0; xPx < pxCount; ++xPx) {
        for(int yPx = 0; yPx < pxCount; ++yPx) {
            color pxColor;
            reflection refl1st = reflect({posCam, capturePts[xPx][yPx],false,false});

            rayTraceColor(refl1st, pxColor,0);

            vector<int> pxRGB = pxColor.intColor();
            image.set_pixel(xPx, pxCount-1 - yPx, pxRGB[0], pxRGB[1], pxRGB[2]);
            //image.set_pixel(xPx, pxCount-1 - yPx, pxColor.r*255, pxColor.g*255, pxColor.b*255);
        }
    }

    image.save_image(imgPath + to_string(++bmpFileCount) + ".bmp");

    cout << "Image generated" << endl;
}

#endif //CAPTURE_LIGHTING_H
