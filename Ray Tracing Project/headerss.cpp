#include <iostream>
#include <fstream>
#include <cassert>
#include <iomanip>
#include <cmath>
#include <stack>
#include <queue>
#include <glut.h>
#include<math.h>
#include "bitmap_image.hpp"
using namespace std;

#define pi (2*acos(0.0))
#define epsilon (1.0e-6)

class homogeneous_point
{
public:
    double x, y, z, w;

    // set the three coordinates, set w to 1
    homogeneous_point(double x, double y, double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
        this->w = 1;
    }

    /*
    default constructor. does nothing. allows declarations like below:
        matrix m;
    therefore, usage is dangerous
    */
    homogeneous_point() {
    }

    // constructs a homogeneous point with given coordinates. forces w to be 1.0
    // if w is zero, raises error
    homogeneous_point(double x, double y, double z, double w)
    {

        this->x = x;
        this->y = y;
        this->z = z;
        this->w = 1;
    }

    // adds two points. returns a point forcing w to be 1.0
    homogeneous_point operator+ (const homogeneous_point& point)
    {
        double x = this->x + point.x;
        double y = this->y + point.y;
        double z = this->z + point.z;
        double w = this->w + point.w;
        homogeneous_point p(x, y, z, w);
        return p;
    }

    // subtracts one point from another. returns a point forcing w to be 1.0
    homogeneous_point operator- (const homogeneous_point& point)
    {
        double x = this->x - point.x;
        double y = this->y - point.y;
        double z = this->z - point.z;
        double w = this->w - point.w;
        homogeneous_point p(x, y, z, w);
        return p;
    }

    // Print the coordinates of a point. exists for testing purpose.
    void print()
    {
        cout << "Point: " << endl;
        cout << x << " " << y << " " << z << " " << w << endl;
    }

};


class Vector
{
public:
    double x, y, z;

    // constructs a vector with given components
    Vector()
    {

    }
    Vector(double x, double y, double z) {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    // keeps the direction same. recalculates the vector to be unit.
    void normalize()
    {
        double r = sqrt(x*x + y*y + z*z);
        x = x / r;
        y = y / r;
        z = z / r;
    }
    double mod()
    {
        return sqrt(x*x + y*y + z*z);
    }
    // add two vectors
    Vector operator+(const Vector& v)
    {
        Vector v1(x+v.x, y+v.y, z+v.z);
        return v1;
    }

    // subtract one vector from another
    Vector operator-(const Vector& v)
    {
        Vector v1(x-v.x, y-v.y, z-v.z);
        return v1;
    }

    homogeneous_point operator+(const homogeneous_point& v)
    {
        homogeneous_point v1(x+v.x, y+v.y, z+v.z,1);
        return v1;
    }

    // subtract one vector from another
    homogeneous_point operator-(const homogeneous_point& v)
    {
        homogeneous_point v1(x-v.x, y-v.y, z-v.z,1);
        return v1;
    }
    // scale a vector with a given coefficient
    Vector operator* (double m)
    {
        Vector v(x*m, y*m, z*m);
        return v;
    }

    // get the dot product of two vectors
    static double dot(Vector a, Vector b)
    {
        return a.x*b.x + a.y*b.y + a.z*b.z;
    }

    // get the cross product of two vectors
    static Vector cross(Vector a, Vector b)
    {
        Vector v(a.y*b.z - a.z*b.y, b.x*a.z - b.z*a.x, a.x*b.y - a.y*b.x);
        return v;
    }

    // print a vector. only for testing purposes.
    void print ()
    {
        cout << "Vector" << endl;
        cout << x << " " << y << " " << z << endl;
    }
};
class Ray{

public:
    homogeneous_point p;
    Vector v;

    Ray(){

    }

};

/*
The matrices are forced to be 4x4. This is because in this assignment, we will deal with points in triangles.
Maximum # of points that we will deal with at once is 3. And all the standard matrices are 4x4 (i.e. scale, translation, rotation etc.)
*/
class matrix
{
public:
    double values[4][4];
    int num_rows, num_cols;

    // only set the number of rows and cols
    matrix(int rows, int cols)
    {
        assert (rows <= 4 && cols <= 4);
        num_rows = rows;
        num_cols = cols;
    }

    // prepare an nxn square matrix
    matrix(int n)
    {
        assert (n <= 4);
        num_rows = num_cols = n;
    }

    // prepare and return an identity matrix of size nxn
    static matrix make_identity(int n)
    {
        assert (n <= 4);
        matrix m(n);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (i == j)
                    m.values[i][j] = 1;
                else
                    m.values[i][j] = 0;
            }
        }
        return m;
    }

    // print the matrix. exists for testing purposes
    void print()
    {
        cout << "Matrix:" << endl;
        for (int i = 0; i < num_rows; i++)
        {
            for (int j = 0; j < num_cols; j++)
            {
                cout << values[i][j] << "\t";
            }
            cout << endl;
        }
    }

    // add the two matrices. Raise error if dimension mismatches
    matrix operator+ (const matrix& m)
    {
        assert (this->num_rows == m.num_rows);
        assert (this->num_cols == m.num_cols);

        matrix m1(num_rows, num_cols);
        for (int i = 0; i < num_rows; i++)
        {
            for (int j = 0; j < num_cols; j++)
            {
                m1.values[i][j] = values[i][j] + m.values[i][j];
            }
        }
        return m1;
    }

    // subtract a matrix from another. raise error if dimension mismatches
    matrix operator- (const matrix& m)
    {
        assert (this->num_rows == m.num_rows);
        assert (this->num_cols == m.num_cols);

        matrix m1(num_rows, num_cols);
        for (int i = 0; i < num_rows; i++)
        {
            for (int j = 0; j < num_cols; j++)
            {
                m1.values[i][j] = values[i][j] - m.values[i][j];
            }
        }
        return m1;
    }

    // multiply two matrices. allows statements like m1 = m2 * m3; raises error is dimension mismatches
    matrix operator* (const matrix& m)
    {
        assert (this->num_cols == m.num_rows);
        matrix m1(this->num_rows, m.num_cols);

        for (int i = 0; i < m1.num_rows; i++) {
            for (int j = 0; j < m1.num_cols; j++) {
                double val = 0;
                for (int k = 0; k < this->num_cols; k++) {
                    val += this->values[i][k] * m.values[k][j];
                }
                m1.values[i][j] = val;
            }
        }
        return m1;
    }

    // multiply a matrix with a constant
    matrix operator* (double m)
    {
        matrix m1(this->num_rows, this->num_cols);
        for (int i = 0; i < num_rows; i++) {
            for (int j = 0; j < num_cols; j++) {
                m1.values[i][j] = m * this->values[i][j];
            }
        }
        return m1;
    }

    // multiply a 4x4 matrix with a homogeneous point and return the resulting point.
    // usage: homogeneous_point p = m * p1;
    // here, m is a 4x4 matrix, intended to be the transformation matrix
    // p1 is the point on which the transformation is being made
    // p is the resulting homogeneous point
    homogeneous_point operator* (const homogeneous_point& p)
    {
        assert (this->num_rows == this->num_cols && this->num_rows == 4);

        matrix m(4, 1);
        m.values[0][0] = p.x;
        m.values[1][0] = p.y;
        m.values[2][0] = p.z;
        m.values[3][0] = p.w;

        matrix m1 = (*this)*m;
        homogeneous_point p1(m1.values[0][0], m1.values[1][0], m1.values[2][0], m1.values[3][0]);
        return p1;
    }

    // return the transpose of a matrix
    matrix transpose()
    {
        matrix m(num_cols, num_rows);
        for (int i = 0; i < num_rows; i++) {
            for (int j = 0; j < num_cols; j++) {
                m.values[j][i] = values[i][j];
            }
        }
        return m;
    }

};

/*
A simple class to hold the color components, r, g, b of a certain shade.
*/
class color {
public:
    double r, g, b;
    color(double r, double g, double b) {
        this->r = r;
        this->g = g;
        this->b = b;
    }
    color() {
        this->r = 0;
        this->g = 0;
        this->b = 0;
    }
};




class drawable{

public:
    string name;
    homogeneous_point center;
    double radius;
    homogeneous_point top,a,b,c,d;
    double ambient, diffuse, specular, reflection;
    int shininess;
    bool texer;
    color **textureBuffer;
    int Theight, Twidth;

    double width;
    color col;
    void load()
    {

        bitmap_image b_img ("tiger.bmp");
        Theight = b_img.height();
        Twidth = b_img.width();
        textureBuffer = new color* [Twidth];
        for (int i = 0; i < Twidth; i++) {
        textureBuffer[i] = new color [Theight];
            for (int j = 0; j < Theight; j++) {
            unsigned char r, g, b;
            b_img.get_pixel(i, j, r, g, b);
            color c(r/255.0, g/255.0, b/255.0);
            textureBuffer[i][j] = c;
            }
        }
    }
    drawable()
    {
        ambient= diffuse= specular= reflection=shininess=0;
        texer = false;

    }

    void print()
    {
        cout<<this->name<<endl;
       if(this->name=="sphere"){
            cout<<"S "<< radius <<endl;
            center.print();


        }
        else if ( this->name == "pyramid")
        {
            cout<<"P"<<endl;

            a.print();b.print();c.print();d.print();top.print();

        }
        else{
            cout<<"BOARD"<<this->width<<endl;
        }
        cout<<ambient<<" "<<diffuse<<" "<<specular<<" "<<reflection<<" "<<shininess<<endl;
    }
    void draw()
    {
        if(this->name=="sphere"){
            glColor3f(this->col.r, this->col.g, this->col.b);
            glPushMatrix();
            glTranslatef(this->center.x, this->center.y, this->center.z);
            glutSolidSphere(this->radius, 150, 150);
            glPopMatrix();

        }
        else if ( this->name == "pyramid")
        {
            glColor3f(this->col.r, this->col.g, this->col.b);
            glBegin(GL_QUADS);{
                //upper hemisphere
                glVertex3f(a.x,a.y,a.z);
                glVertex3f(b.x,b.y,b.z);
                glVertex3f(c.x,c.y,c.z);
                glVertex3f(d.x,d.y,d.z);

            }glEnd();

            glBegin(GL_TRIANGLES);
            {
                glVertex3f(a.x,a.y,a.z);
                glVertex3f(b.x,b.y,b.z);
                glVertex3f(top.x, top.y, top.z);
            }glEnd();

            glBegin(GL_TRIANGLES);
            {
                glVertex3f(b.x,b.y,b.z);
                glVertex3f(c.x,c.y,c.z);
                glVertex3f(top.x, top.y, top.z);
            }glEnd();

            glBegin(GL_TRIANGLES);
            {
                glVertex3f(c.x,c.y,c.z);
                glVertex3f(d.x,d.y,d.z);
                glVertex3f(top.x, top.y, top.z);
            }glEnd();

            glBegin(GL_TRIANGLES);
            {
                glVertex3f(a.x,a.y,a.z);
                glVertex3f(d.x,d.y,d.z);
                glVertex3f(top.x, top.y, top.z);
            }glEnd();

        }
        else if( this->name == "cb"){

            for(int i =-100 ; i<=100 ;i++)
            {

                for(int j = -100;j<=100;j++){
                    double col = (abs(i+j)%2)*1.0;
                    //if(i*j <0 )col = 1-col;
                    glColor3f(col, col, col);
                    glBegin(GL_QUADS);{
                //upper hemisphere
                        glVertex3f(i*this->width,j*this->width,0);
                        glVertex3f(i*this->width+this->width,j*this->width,0);
                        glVertex3f(i*this->width+this->width,j*this->width+this->width,0);
                        glVertex3f(i*this->width,j*this->width+this->width,0);

                    }glEnd();



                }
            }

        }

    }

    double getint(Ray R,color & col,homogeneous_point & intp, Vector & normal)
    {
        if(this->name=="cb")
        {
            double t = -1*R.p.z/R.v.z;
            homogeneous_point p = R.v*t+R.p;
            int c = abs(int(p.x/this->width)+int(p.y/this->width))%2;
            c=1-c;
            //cout<<"colour"<<c<<endl<<endl;
            if(p.x*p.y<0 ){
                c=1-c;
            }

            col.b = c;
            col.r = c;
            col.g = c;
            homogeneous_point point = R.v*t+R.p;
            intp.x = point.x;   intp.y = point.y;   intp.z = point.z;
            normal.x = 0 ; normal.y = 0; normal.z = 1;
            if(texer)
            {
                double i,j;
                i = abs(p.x)/this->width;
                j = abs(p.y)/this->width;
                i = i - int(i);
                j = j - int(j);
                int x , y ;
                x = int(i*Twidth);
                y = int(j*Theight);
                if(p.x<0 ) x = Twidth - x- 1;
                if(p.y<0) y = Theight - y - 1;
                //cout<<x<<"    "<<y<<" "<<Twidth<<" "<<Theight<<endl;
                col.b = textureBuffer[x][y].b;
                col.r = textureBuffer[x][y].r;
                col.g = textureBuffer[x][y].g;

            }
            return t;

        }
        else if(this->name == "sphere")
        {
           homogeneous_point r0 = R.p-this->center;
           double a =1;
           double b = 2*(r0.x*R.v.x+ r0.y*R.v.y + r0.z*R.v.z);
           double c = r0.x*r0.x+r0.y*r0.y+r0.z*r0.z- this->radius*this->radius;
           double d = b*b-4*a*c;

           if(d<0){ return -1;}

           d = sqrt(d);
           col.r = this->col.r;     col.g = this->col.g;     col.b = this->col.b;
           double t = min( (-1*b+d)/(2*a), (-1*b-d)/(2*a));
           double t2 = max( (-1*b+d)/(2*a), (-1*b-d)/(2*a)  );
           if (t*t2<0){
            t = t2;
           }
           homogeneous_point point = R.v*t+R.p;
           intp.x = point.x;   intp.y = point.y;   intp.z = point.z;
           normal.x = intp.x - this->center.x ; normal.y = intp.y - this->center.y  ; normal.z = intp.z - this->center.z;
           return t;


        }
        else if (this->name == "pyramid"){
            double p[5];
            Vector n[5];
            p[0]=this->raytoTriangle(R,this->a,this->b,this->top,n[0]);
            p[1]=this->raytoTriangle(R,this->c,this->b,this->top,n[1]);
            p[2]=this->raytoTriangle(R,this->c,this->d,this->top,n[2]);
            p[3]=this->raytoTriangle(R,this->a,this->d,this->top,n[3]);
            n[4].x=0;n[4].y=0;n[4].z=-1;
            p[4] = max(this->raytoTriangle(R,this->a,this->b,this->c,n[4]),this->raytoTriangle(R,this->b,this->c,this->d,n[4]));
            double res=1000000;
            int flag = 1;
            for(int i = 0;i<5;i++)
            {
                if(p[i]>0 && p[i]<res)
                {
                    res=p[i];
                    normal.x = n[i].x ; normal.y = n[i].y; normal.z = n[i].z;
                }
            }
            int sign = 1;
            if(normal.z<0 && res!=p[4]) sign = -1;
            normal.x =normal.x*sign ; normal.y = normal.y*sign; normal.z = normal.z*sign;
            col.r = this->col.r;     col.g = this->col.g;     col.b = this->col.b;
            homogeneous_point point = R.v*res+R.p;
            intp.x = point.x;   intp.y = point.y;   intp.z = point.z;
            return res;


        }
        else -1;
    }

    double raytoTriangle(Ray R,homogeneous_point p0,homogeneous_point p1,homogeneous_point p2,Vector & normal)
    {
        Vector e1,e2,n,q,s,r;
        e1.x = p1.x - p0.x; e1.y = p1.y - p0.y; e1.z = p1.z - p0.z;
        e2.x = p2.x - p0.x; e2.y = p2.y - p0.y; e2.z = p2.z - p0.z;

        n = e1.cross(e1,e2);
        n.normalize();
        normal.x = n.x; normal.y = n.y; normal.z = n.z;
        q = R.v.cross(R.v,e2);
        double a = e1.dot(e1,q);
        if(abs(a)<=.00001) return -1;
        s.x = R.p.x - p0.x; s.y = R.p.y - p0.y; s.z = R.p.z - p0.z;
        s = s*(1/a);
        r = s.cross(s,e1);
        double alpa,beta,theta;
        alpa = s.dot(s,q);
        beta = r.dot(r,R.v);
        theta = 1.0 - alpa - beta;
        if ( alpa<0 || beta < 0 || theta<0) return -1;

        return e2.dot(e2,r);
    }

};




class light{

public:
    string name;
    homogeneous_point pos, dir_pon;
    double falloff;
    double angle;

    light(){
        angle = 400;
    }
    void draw()
    {
       glColor3f(1,1,1);
        glPushMatrix();
        glTranslatef(this->pos.x, this->pos.y, this->pos.z);
        glutSolidSphere(10, 150, 150);
        glPopMatrix();

    }
    bool checkspot(homogeneous_point point)
    {
        if(this->angle==400) return true;
        Vector sp, dir ;
        sp.x = point.x -this->pos.x; sp.y = point.y - this->pos.y; sp.z = point.z - this->pos.z;
        sp.normalize();
        dir.x = this->dir_pon.x - this->pos.x; dir.y = this->dir_pon.y - this->pos.y; dir.z = this->dir_pon.z - this->pos.z;
        dir.normalize();
        double angle = acos(sp.dot(sp,dir))*180/3.141592653;
        //if(angle <= this->angle)cout<<angle<<this->name;
        if (angle > this->angle) return false;
        else return true;
    }

};

class buffer{
public:
    color col;
    homogeneous_point p;
    buffer(){
    }

};







/**double eye_x, eye_y, eye_z;
double look_x, look_y, look_z;
double up_x, up_y, up_z;
double fov_x, fov_y, aspectRatio, near, far;
color backgroud;
int screen_x, screen_y;**/



/*void scan_convert() {
    ifstream stage3;
    stage3.open("stage3.txt");

    color** pixels = new color*[screen_x];
    double** zs = new double*[screen_x];
    for (int i = 0; i < screen_x; i++) {
        pixels[i] = new color [screen_y];
        for (int j = 0; j < screen_y; j++) {
            pixels[i][j] = backgroud;
        }
        zs[i] = new double [screen_y];
        for (int j = 0; j < screen_y; j++) {
            zs[i][j] = +20; // a very large value intended as +INFINITY
        }
    }

    // perform scan conversion, populate the 2D array pixels
    // the array zs is the z-buffer.


    // the following code generates a bmp image. do not change this.
    bitmap_image image(screen_x, screen_y);
    for (int x = 0; x < screen_x; x++) {
        for (int y = 0; y < screen_y; y++) {
            image.set_pixel(x, y, pixels[x][y].r, pixels[x][y].g, pixels[x][y].b);
        }
    }
    image.save_image("out.bmp");

    // free the dynamically allocated memory

}**/



