#ifndef OBJECTS_H
#define OBJECTS_H

#ifdef freeGlut
    #include <GL/glut.h>
#else
    #include <windows.h>
    #include <glut.h>
#endif //freeGlut

#include <bits/stdc++.h>
using namespace std;

/*
 * POINT,COLOR,RAY
 */

struct point {
    double x,y,z;

    point(): x(0), y(0), z(0) {};

    point(double x, double y, double z) : x(x), y(y), z(z) {}

    point operator+(const point& rhs) const {
        return {x+rhs.x, y+rhs.y, z+rhs.z};
    }

    point& operator+=(const point&rhs) {
        *this = *this+rhs;
        return *this;
    }

    point operator*(const double scalar) const {
        return {x*scalar, y*scalar, z*scalar};
    }

    point& operator*=(const double scalar) {
        *this = *this*scalar;
        return *this;
    }

    point operator-(const point&rhs) const {
        return *this + rhs*-1;
    }

    point& operator-=(const point& rhs) {
        *this = *this-rhs;
        return *this;
    }

    bool operator==(const point &rt) const {
        return x == rt.x && y == rt.y && z == rt.z;
    }

    point vct(const point& to) const {
        return to - *this;
    }

    double distSq(const point& rt) const {
        return (*this-rt).lengthSq();
    }

    double dist(const point& rt) const {
        return sqrt(distSq(rt));
    }

    /*
     * VCT methods
     */

    double lengthSq() const {
        return this->dot(*this);
    }

    double length() const {
        return sqrt(lengthSq());
    }

    bool normalize() {  //return false if already normalized
        double lenSq = lengthSq();
        if(abs(lenSq - 1) < eps)
            return false;

        if(abs(lenSq) < eps) {
            cout << "warning: trying to normalize 0 vector in point::normalize()" << endl;
            return false;
        }

        (*this) *= (1/sqrt(lenSq));
        return true;
    }

    bool isNormalized() const {
        double lenSq = lengthSq();
        if(abs(lenSq - 1) < eps)
            return true;
        return false;
    }

    double dot(const point& rt) const {
        return x*rt.x + y*rt.y + z*rt.z;
    }

    point cross(const point& rt) const {
        return {y*rt.z-z*rt.y, z*rt.x-x*rt.z, x*rt.y-y*rt.x};
    }

    point reflect(point normal) const {  //returns reflected ray (with unchanged length)
        double c = this->dot(normal);
        return *this + normal*c*(-2);
    }

    /*
     * <VCT END>
     */

    vector<double> toVector() {  //not a VCT function!
        return vector<double>({x,y,z});
    }

    void print(const string& name="") const {
        cout << name << "(" << x << "," << y << "," << z << ")" << flush;
    }
};

void glVertex3f(const point& p) {
    glVertex3f(p.x,p.y,p.z);
}

struct color {
    double r,g,b;

    color(): r(0), g(0), b(0) {};

    color(double r, double g, double b) : r(r), g(g), b(b) {}

    color operator+(const color& rhs) const {
        return {r+rhs.r, g+rhs.g, b+rhs.b};
    }

    color& operator+=(const color&rhs) {
        *this = *this+rhs;
        return *this;
    }

    color operator*(const double scalar) const {
        return {r*scalar, g*scalar, b*scalar};
    }

    color& operator*=(const double scalar) {
        *this = *this*scalar;
        return *this;
    }

    color operator*(const color rt) const {
        return {r*rt.r, g*rt.g, b*rt.b};
    }

    color& operator*=(const color rt) {
        *this = *this * rt;
        return *this;
    }

    vector<int> intColor() const {
        vector<int> rgb;
        rgb.push_back(r>=1 ? 255 : max(r*256,0.0));
        rgb.push_back(g>=1 ? 255 : max(g*256,0.0));
        rgb.push_back(b>=1 ? 255 : max(b*256,0.0));

        return rgb;
    }

    void print(const string& name="") const {
        cout << name << "->RGB(" << r << "," << g << "," << b << ")" << flush;
    }
};

struct ray {
    point origin, dir;

    ray() = default;

    ray(const point &origin, const point &dir, bool dirIsVct = true, bool normalized = true) : origin(origin){
        if(dirIsVct)
            this->dir = dir;
        else
            this->dir = origin.vct(dir);

        if(!normalized)
            normalize();
    }

    bool normalize() {  //return false if already normalized
        return dir.normalize();
    }

    bool isNormalized() const {
        return dir.isNormalized();
    }

    point fwd(double c) const {
        return origin + dir*c;
    }

    bool isInvalid() const {
        return dir == point();
    }

    void draw(double c=1) const {  //must set color beforehand
        glBegin(GL_LINES); {
            glVertex3f(origin);
            glVertex3f(fwd(c));
        } glEnd();
    }

    void print() const {
        origin.print("origin");
        cout << ", ";
        dir.print("dir");
    }
};

void glColor3f(color rgb) {
    glColor3f(rgb.r,rgb.g,rgb.b);
}

/*
 * ALGEBRAIC CALC
 */

vector<double> quadRealSolve(double b, double c) {  //smaller solution inserted first
    double disc = b*b-4*c;
    vector<double> soln;
    if(disc < 0)
        return soln;
    disc = sqrt(disc);
    soln.push_back( (-b-disc)/2.0 );
    soln.push_back( (-b+disc)/2.0 );
    return soln;
}

vector<double> quadRealSolve(double a, double b, double c) {  //smaller solution inserted first
    return quadRealSolve(b/a,c/a);
}

vector<double> linearEqn2(const vector<double>& a, const vector<double>& b, const vector<double>& c) {
    //eqn. form: ax+by=c -> result: (x,y)
    double denom = a[0]*b[1] - a[1]*b[0];
    return vector<double>({ (c[0]*b[1] - c[1]*b[0]) / denom, (a[0]*c[1] - a[1]*c[0]) / denom });
}

/*
 * OBJECTS
 */

//sphere drawing params
#define radPrec 150
#define htPrec 50
const double pi = 2*acos(0.0);

//axis params
#define axisLength 400
color xAxisColor(1,0,0);
color yAxisColor(0,0,1);
vector<double> axisCoeff({0,0,0,0});
double axisShine = 0;

//board params
#define showBoard true
#define boardRowCol 30
#define tileLength 20
color tileColor1(1,1,1), tileColor2(0,0,0);
vector<double> boardCoeff({0.5,0.25,0.4,0.5});
double boardShine = 30;

//light sources
struct light;
extern vector<light> lights;

#define lightRad 2

struct object {
    string type;
    vector<point> pt;
    vector<double> len;
    //color & lighting -> set yourself
    color color1;
    color color2;
    vector<double> coeff;
    double shine;

    virtual void draw() const {};

    void glutColor(int index = 1) const {
        if(index == 1)
            glColor3f(color1);
        else
            glColor3f(color2);
    };

    virtual double intersect(const ray& in) const { return -1; }

    virtual point normal(const point& p) const { return {}; }  //unit normal vct; doesn't check if point is on body

    ray reflect(const ray& in) const {  //-1 in case of not reflect
        double c = intersect(in);
        if(c <= 0)
            return {};
        point intPt = in.fwd(c-reflectEps);
        return { intPt, in.dir.reflect(normal(intPt)) };
    }

    virtual color getColor(const point& p) const {
        return color1;
    }

    virtual void print() const {
        cout << "[DEFAULT]" << "\n";
        cout << type << "\n";
        printPt();
        printLen();
        printColor();
        printCoeff();
    };

    void printPt(const string& ptName = "pt") const {
        if(pt.size() == 1) {
            pt[0].print(ptName);
            cout << endl;
            return;
        }

        for(int i = 0; i < pt.size(); ++i) {
            pt[i].print(ptName + to_string(i));
            cout << ", ";
        }
        cout << "\b\b " << endl;
    }

    void printLen(const string& lenName = "len") const {
        cout << lenName << "(";
        for(double db: len) {
            cout << db << ",";
        }
        cout << (len.empty()? "" : "\b") << ")" << endl;
    }

    void printColor(int colorN = 1) const {
        color1.print("color");
        cout << "\n";

        if(colorN > 1) {
            color2.print("color2");
            cout << "\n";
        }
        cout << flush;
    }

    void printCoeff() const {
        cout << "coeff(";
        for(double db: coeff) {
            cout << db << ",";
        }
        cout << (coeff.empty()? "" : "\b") << ")";
        cout << ", shine=" << shine << endl;
    }
};

struct sphere: object {
    sphere(point center, double radius) {
        type = "sphere";
        pt.push_back(center);
        len.push_back(radius);
    }

    //override draw
    void draw() const override {
        int stacks = htPrec, slices = radPrec;

        double h,r, R = len[0];
        point sphere[stacks+1][slices+1];
        for(int i = 0; i <= stacks; i++) {
            h = R * cos( ((double)i/stacks)*(pi/2) );
            r = R * sin( (double)i/stacks*(pi/2) );
            for(int j = 0; j <= slices; j++) {
                sphere[i][j].x = r * cos( (double)j/slices*2*pi );
                sphere[i][j].y = r * sin( (double)j/slices*2*pi );
                sphere[i][j].z = h;
            }
        }

        glutColor();
        glPushMatrix();
        glTranslatef(pt[0].x, pt[0].y, pt[0].z);
        for(int sign = -1; sign <= 1; sign+= 2)
            for(int i = 0; i < stacks; i++) {
                for(int j = 0; j < slices; j++) {
                    glBegin(GL_QUADS);
                    {
                        glVertex3f(sphere[i][j].x, sphere[i][j].y, sign*sphere[i][j].z);
                        glVertex3f(sphere[i][j + 1].x, sphere[i][j + 1].y, sign*sphere[i][j + 1].z);
                        glVertex3f(sphere[i + 1][j + 1].x, sphere[i + 1][j + 1].y, sign*sphere[i + 1][j + 1].z);
                        glVertex3f(sphere[i + 1][j].x, sphere[i + 1][j].y, sign*sphere[i + 1][j].z);
                    }
                    glEnd();
                }
            }
        glPopMatrix();
    }

    double intersect(const ray &in) const override {
        point R0 =  pt[0].vct(in.origin), Rd = in.dir;
        double b = R0.dot(Rd)*2;
        double c = R0.lengthSq() - len[0]*len[0];
        vector<double> soln = quadRealSolve(b,c);
        for(double db: soln) {
            if(db > 0)
                return db;
        }
        return -1;
    }

    point normal(const point& p) const override {
        point norm = pt[0].vct(p);
        norm.normalize();
        return norm;
    }

    void print() const override {
        cout << type << "\n";
        printPt("center");
        printLen("radius");
        printColor();
        printCoeff();
    }
};

struct triangle: object {
    point normalVct;

    explicit triangle(const vector<point>& vertex) {
        type = "triangle";
        pt = vector<point>(vertex);

        normalVct = (pt[1]-pt[0]).cross(pt[2]-pt[0]);
        normalVct.normalize();
    }

    void draw() const override {
        glutColor();
        glBegin(GL_TRIANGLES); {
            for(point p: pt)
                glVertex3f(p);
        } glEnd();
    }

    double intersect(const ray &in) const override {
        double c = (pt[0]-in.origin).dot(normalVct) / in.dir.dot(normalVct);
        if(c < 0)
            return -1;

        point intPt = in.fwd(c);
        vector<double> soln = linearEqn2( (pt[1]-pt[0]).toVector(), (pt[2]-pt[0]).toVector(),
                                          (intPt-pt[0]).toVector() );

        double sum = 0;
        for(double db:soln) {
            if(db < 0)
                return -1;
            sum += db;
        }
        if(sum > 1)
            return 0;

        return c;
    }

    point normal(const point &p) const override {
        return normalVct;
    }

    void print() const override {
        cout << type << "\n";
        printPt("v");
        printColor();
        printCoeff();
    }
};

struct curve: object {
    vector<double> params;

    // len: 1st 3 -> cubeDim, then params
    curve(const vector<double>& params, point cubeRef, const vector<double>& cubeDim) {
        type = "curve";
        this->params = vector<double>(params);
        len = vector<double>(cubeDim);
        pt.push_back(cubeRef);
    }

    double intersect(const ray &in) const override {
        point R0 =  pt[0].vct(in.origin), Rd = in.dir;
        vector<double> quad({0,0,0});  //a,b,c

        for(int i = 0; i < 2; i++) {  //0->a, 1->c
            point& R = ( (i == 0)? Rd:R0 );
            double temp = params[0] * R.x * R.x + params[1] * R.y * R.y + params[2] * R.z * R.z
                          + params[3]*R.x*R.y + params[4]*R.x*R.z + params[5]*R.y*R.z;
            if(i == 1)
                temp += params[6]*R.x + params[7]*R.y + params[8]*R.z + params[9];
            quad[(i == 0)? 0:2] = temp;
        }

        quad[1] = Rd.dot({params[0]*R0.x, params[1]*R0.y, params[2]*R0.z})*2
                    + point(R0.x*Rd.y + R0.y*Rd.x, R0.x*Rd.z + R0.z*Rd.x, R0.y*Rd.z + R0.z*Rd.y)
                            .dot({params[3], params[4], params[5]}) +
                            Rd.dot({params[6], params[7], params[8]});


        vector<double> soln = quadRealSolve(quad[0],quad[1],quad[2]);
        for(double db: soln) {
            if(db > 0) {
                point refVct = in.fwd(db) - pt[0];
                if( (len[0] == 0 || refVct.x >= 0 && refVct.x <= len[0])
                    && (len[1] == 0 || refVct.y >= 0 && refVct.y <= len[1])
                    && (len[2] == 0 || refVct.z >= 0 && refVct.z <= len[2]) )
                return db;
            }
        }
        return -1;
    }

    point normal(const point &p) const override {
        vector<double> pV({p.x,p.y,p.z});
        vector<double> normV({0,0,0});
        for(int i = 0; i < 3; i++) {
            normV[i] += 2*params[i]*pV[i];
            normV[i] += params[6+i];
        }
        point norm(normV[0], normV[1], normV[2]);
        norm += point(params[3]*pV[1] + params[4]*pV[2], params[3]*pV[0] + params[5]*pV[2],
                      params[4]*pV[0] + params[5]*pV[1]);
        norm.normalize();

        return norm;
    }

    void print() const override {
        cout << type << "\n";
        printPt("cubePt");
        printLen("cubeDim");

        cout << "params(";
        for(double db: params)
            cout << db << ",";
        cout << "\b)" << "\n";

        printColor();
        printCoeff();
    }
};

struct axes: object {

    explicit axes(point origin = point(0,0,0), double axisLen = axisLength) {
        type = "axes";
        pt.push_back(origin);
        len.push_back(axisLen);

        color1 = xAxisColor;
        color2 = yAxisColor;
    }

    void draw() const override {
        double axisLen = len[0];
        glBegin(GL_LINES);{
            glutColor(1);
            glVertex3f( axisLen,0,0);
            glVertex3f(-axisLen,0,0);

            glutColor(2);
            glVertex3f(0,0,axisLen);
            glVertex3f(0,0,-axisLen);
        }glEnd();
    }

    void print() const override {
        cout << type << "\n";
        printPt("origin");
        printLen("axisLength");
        printColor(2);
        printCoeff();
    }
};

struct board: object {
    int rows, cols;
    point normalVct{0,0,1};

    explicit board(double tileLen = tileLength, int rows = boardRowCol, int cols = boardRowCol) {
        type = "board";
        pt.emplace_back(-cols*tileLen/2.0, -rows*tileLen/2.0, 0);
        len.push_back(tileLen);
        this->rows = rows;
        this->cols = cols;

        color1 = tileColor1;
        color2 = tileColor2;
    }

    void draw() const override {
        glPushMatrix();
        glTranslatef(pt[0].x, pt[0].y, pt[0].z);
        double tileLen = len[0];

        for(int col = 0; col < cols; col++) {
            for(int row = 0; row < rows; row++) {
                glPushMatrix();
                glTranslatef(col*tileLen,row*tileLen,0);
                glutColor(2 - (row+col)%2);

                glBegin(GL_QUADS); {
                    glVertex3f(0,0,0);
                    glVertex3f(tileLen,0,0);
                    glVertex3f(tileLen,tileLen,0);
                    glVertex3f(0,tileLen,0);
                }glEnd();
                glPopMatrix();
            }
        }
        glPopMatrix();
    }

    double intersect(const ray &in) const override {
        double c = -in.origin.z/in.dir.z;
        point fwdPt = in.fwd(c);
        double dx = fwdPt.x-pt[0].x, dy = fwdPt.y-pt[0].y;

        if(dx >= 0 && dy >= 0 && dx <= len[0]*cols && dy <= len[0]*rows)
            return c;
        return -1;
    }

    point normal(const point& p) const override {  //warning -> there are two sides to a board
        return normalVct;
    }

    color getColor(const point& p) const override {  //doesn't check if point is on the board
        int xSlot = (p.x-pt[0].x)/len[0], ySlot = (p.y-pt[0].y)/len[0];
        if((xSlot+ySlot) % 2)
            return color1;
        return color2;
    }

    void print() const override {
        cout << type << "\n";
        printPt("corner");
        printLen("tileLength");
        cout << "rows(" << rows << "), cols(" << cols << ")" << "\n";
        printColor(2);
        printCoeff();
    }
};

/*
 * LIGHT
 */

struct light {
    point pt;
    color ltColor;
    object* ltObj;

    light(const point &pt, const color &ltColor) : pt(pt), ltColor(ltColor) {
        ltObj = new sphere(pt,lightRad);
        ltObj->color1 = ltColor;
    }

    void draw() const {
        ltObj->draw();
    }

    void print() const {
        pt.print("src");
        ltColor.print();
        cout << endl;
    }
};

#endif //OBJECTS_H
