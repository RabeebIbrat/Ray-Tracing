#ifndef GLUT_DISPLAY_H
#define GLUT_DISPLAY_H

#ifdef freeGlut
    #include <GL/glut.h>
#else
    #include <windows.h>
    #include <glut.h>
#endif //freeGlut

#include <bits/stdc++.h>
#include "1605055_objects.h"

using namespace std;

/*
 * LOAD DATA
 */

string inputPath = "IO\\input.txt";

//inputs to be loaded
int recurLevel;
int pxCount;
vector<object*> things;
vector<light> lights;

//display params
int maxRecurLevel = 6;
extern bool allAmbient;

//axes
object* axisLines;
bool showAxes = false;

void loadStaticData() {
    axisLines = new axes;
    axisLines->coeff = axisCoeff;
    axisLines->shine = axisShine;

    if(showBoard) {
        object *obj = new board;
        obj->coeff = boardCoeff;
        obj->shine = boardShine;
        things.push_back(obj);
    }
}

void loadData(const string& file = inputPath) {
    loadStaticData();

    /*
     * OBJECTS
     */
    ifstream read(inputPath);

    read >> recurLevel;
    read >> pxCount;

    if(recurLevel > maxRecurLevel)
        maxRecurLevel = recurLevel;

    int count;
    read >> count;
    while(count--) {  //each object
        string type;
        read >> type;
        object* obj;

        if(type == "sphere") {
            double x,y,z, R;
            read >> x >> y >> z >> R;
            obj = new sphere(point(x,y,z),R);
        }

        else if(type == "triangle") {
            vector<point> pts;
            double x,y,z;
            for(int i = 0; i < 3; ++i) {
                read >> x >> y >> z;
                pts.emplace_back(x,y,z);
            }
            obj = new triangle(pts);
        }

        else if(type == "general") {
            vector<double> params;
            double db;
            for(int i = 0; i < 10; ++i) {
                read >> db;
                params.push_back(db);
            }

            double x,y,z;
            read >> x >> y >> z;

            vector<double> len;
            for(int i = 0; i < 3; ++i) {
                read >> db;
                len.push_back(db);
            }

            obj = new curve(params,point(x,y,z),len);
        }

        else {
            cout << "object type invalid in loadData(). Exiting...";
            exit(1);
        }

        double r,g,b;
        read >> r >> g >> b;
        obj->color1 = color(r,g,b);

        double coeff;
        for(int i = 0; i < 4; ++i) {
            read >> coeff;
            obj->coeff.push_back(coeff);
        }
        read >> obj->shine;

        things.push_back(obj);
    }

    /*
     * LIGHTS
     */

    read >> count;
    while(count--) {  //each light
        double x,y,z,r,g,b;
        read >> x >> y >> z >> r >> g >> b;
        lights.emplace_back(point(x,y,z),color(r,g,b));
    }
}

void printData() {
    cout << "recursion level: " << recurLevel << "\n";
    cout << "pixel count: " << pxCount << "\n";

    cout << "\nTHINGS: \n\n";
    axisLines->print();
    cout << "\n";
    for(object* obj: things) {
        obj->print();
        cout << "\n";
    }

    cout << "LIGHTS: " << "\n";
    for(const light& lt: lights)
        lt.print();
    cout << endl;
}

/*
 ***** GLUT DISPLAY *****
 */

/*
 * CAMERA PARAMS
 */

#define moveSen 1
#define camRotSen 1


point posCam(111.022,8.39314,55.5365);
point lVec(-0.916695,0.0679805,-0.393763);
point rVec(0.110341,0.990172,-0.0859341);
point uVec(-0.384051,0.122223,0.915186);

#define fovY 80
#define aspect 1
#define zNear 1
#define zFar 1000

#define screenWd 500
#define screenHt 500

/*
 * LASER
 */

color laserColor(1,0,0);
double laserLength = 1000;
int laserReflect = 4;
#define showHitRad 3

point laserCam, laserL;
bool showLaser = false;
double hitRad = 0;

struct reflection {
    ray reflRay;
    double reflDist;
    object* reflObj;

    reflection(const ray &reflRay, double reflDist, object *reflObj) :
                reflRay(reflRay), reflDist(reflDist), reflObj(reflObj) {}

    bool isInValid() const {
        return reflRay.isInvalid();
    }
};

extern reflection reflect(const ray& rayNow);

void drawLaser(ray laser = ray{laserCam, laserL}, int level = laserReflect) {
    if(level == 0)
        return;

    reflection refl = reflect(laser);
    ray reflected = refl.reflRay;
    double c = refl.reflDist;

    glColor3f(laserColor);
    laser.draw(min(c,laserLength));

    sphere hit((laser.origin), hitRad);
    if(level == laserReflect-1) {
        hit.color1 = color(1,0,1);
        hit.draw();
    }
    else if(level == laserReflect-2) {
        hit.color1 = color(0,1,0);
        hit.draw();
    }

    if(!reflected.isInvalid())
        drawLaser(reflected, level-1);
}

/*
 * GRID
 */

color gridColor(0,0,1);
#define gridPxCount 9

point gridCam;
bool showGrid = false;
vector<vector<point>> gridPts;

extern vector<vector<point>> capturePts;
extern void updateCapturePts();

void drawGrid() {
    glColor3f(gridColor);
    int offset = pxCount / gridPxCount;
    for(int x = 0; x < gridPxCount; ++x) {
        for(int y = 0; y < gridPxCount; ++y) {
            ray(gridCam, gridCam.vct(gridPts[x*offset][y*offset])).draw();
        }
    }
}

/*
 * GLUT CORE CODE
 */

void rotate(point& ref, point& rt, point& up, double angle) {  //angle in degree
    angle = angle*pi/180;
    point temp_rt =  rt*cos(angle) + up*sin(angle);
    up = up*cos(angle) - rt*sin(angle);
    rt = temp_rt;
}

//image generation function from "1605055_capture_lighting.h"
extern void genImage();

void keyboardListener(unsigned char key, int x,int y) {

    switch(key) {
        case '1':
            rotate(uVec, rVec, lVec, camRotSen);
            break;

        case '2':
            rotate(uVec, rVec, lVec, -camRotSen);
            break;

        case '3':
            rotate(rVec, lVec, uVec, camRotSen);
            break;

        case '4':
            rotate(rVec, lVec, uVec, -camRotSen);
            break;

        case '5':
            rotate(lVec, uVec, rVec, camRotSen);
            break;

        case '6':
            rotate(lVec, uVec, rVec, -camRotSen);
            break;

        /*
         * DEBUG
         */

        case 'v':
            posCam.print("posCam");
            cout << "\n";
            lVec.print("lVec");
            cout << "\n";
            rVec.print("rVec");
            cout << "\n";
            uVec.print("uVec");
            cout << "\n" << endl;
            break;

        case 'l':
            showLaser = !showLaser;
            if(showLaser) {
                laserCam = posCam;
                laserL = lVec;
                hitRad = showHitRad;
            }
            break;

        case 'g':
            showGrid = !showGrid;
            if(showGrid) {
                gridCam = posCam;
                updateCapturePts();
                gridPts = move(capturePts);
            }
            break;

        /*
         * IMAGE GEN
         */

        case 'c':
            genImage();
            break;

        case '+':
            if(recurLevel < maxRecurLevel)
                recurLevel++;
            cout << "recur level: " << recurLevel << endl;
            break;

        case '-':
            if(recurLevel > 0)
                recurLevel--;
            cout << "recur level: " << recurLevel << endl;
            break;

        case 'a':
            allAmbient = !allAmbient;
            cout << "allAmbient: " << allAmbient << endl;
            break;

        default:
            break;
    }
}

void specialKeyListener(int key, int x,int y) {

    switch(key) {
        case GLUT_KEY_DOWN:
            posCam -= lVec * moveSen;
            break;

        case GLUT_KEY_UP:
            posCam += lVec * moveSen;
            break;

        case GLUT_KEY_RIGHT:
            posCam += rVec * moveSen;
            break;

        case GLUT_KEY_LEFT:
            posCam -= rVec * moveSen;
            break;

        case GLUT_KEY_PAGE_UP:
            posCam += uVec * moveSen;
            break;

        case GLUT_KEY_PAGE_DOWN:
            posCam -= uVec * moveSen;
            break;

        default:
            break;
    }
}

void mouseListener(int button, int state, int x, int y){	//x,y is the x-y of the screen (2D)
    switch(button){
        case GLUT_RIGHT_BUTTON:
            if(state == GLUT_DOWN)
                showAxes=!showAxes;
            break;

        default:
            break;
    }
}

void display(){

    //clear the display
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0,0,0,0);	//color black
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    //load the correct matrix -- MODEL-VIEW matrix
    glMatrixMode(GL_MODELVIEW);
    //initialize the matrix
    glLoadIdentity();

    /*
     * CAMERA UPDATE
     */

    point eye = posCam, look = eye + lVec;
    gluLookAt(posCam.x, posCam.y, posCam.z, look.x, look.y, look.z, uVec.x, uVec.y, uVec.z);

    //again select MODEL-VIEW
    glMatrixMode(GL_MODELVIEW);

    /*
     * MY OBJECTS
     */
    if(showAxes)
        axisLines->draw();
    for(object* obj: things)
        obj->draw();
    for(light lt: lights)
        lt.draw();

    if(showLaser)
        drawLaser();
    if(showGrid)
        drawGrid();

    //ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
    glutSwapBuffers();
}

void animate(){
    glutPostRedisplay();
}

void init(){
    //clear the screen
    glClearColor(0,0,0,0);

    //load the PROJECTION matrix
    glMatrixMode(GL_PROJECTION);

    //initialize the matrix
    glLoadIdentity();

    //give PERSPECTIVE parameters
    gluPerspective(fovY, aspect, zNear, zFar);
}

void glutStart(int argc, char** argv) {
    glutInit(&argc,argv);
    glutInitWindowSize(screenWd, screenHt);
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB | GLUT_MULTISAMPLE);	//Depth, Double buffer, RGB color

    glutCreateWindow("Ray Tracing");

    init();

    glEnable(GL_DEPTH_TEST);	//enable Depth Testing

    glutDisplayFunc(display);	//display callback function
    glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    glutMouseFunc(mouseListener);

    glutMainLoop();		//The main loop of OpenGL
}

#endif //GLUT_DISPLAY_H
