hon#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <iostream>
#include <fstream>
#include <windows.h>
#include <glut.h>
//#include"headers.cpp"
#define pi (2*acos(0.0))


//all library functions are here
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

        bitmap_image b_img ("texture.bmp");
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
            //c=1-c;
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

// ends here
double cameraHeight;
double cameraAngle, n, f;
double fov_y,ration;
int rec_lev,pic_dim;
int drawaxes;

vector<drawable> drawings;
vector<light> lights;
homogeneous_point pos;
Vector l,u,r;
void drawAxes()
{
	if(drawaxes==1)
	{
		glColor3f(6.0, 8.0, 1.0);
		glBegin(GL_LINES);{
			glVertex3f( 1000,0,10);
			glVertex3f(-1000,0,10);

			glVertex3f(0,-1000,10);
			glVertex3f(0, 1000,10);

			glVertex3f(0,0, 1000);
			glVertex3f(0,0,-1000);
		}glEnd();
	}
}
buffer ** getbuffer()
{
    buffer ** pixels = new buffer*[pic_dim];
    for( int i = 0;i<pic_dim;i++)
    {
        pixels[i] = new buffer[pic_dim];
    }
    homogeneous_point m = l+pos;
    double y =n*tan(pi*fov_y/(2*180));
    double x = ration*y;
    homogeneous_point top_left =  r*-1*x + u*y+m;
    double xr = 2*x/pic_dim;
    double yr = 2*y/pic_dim;
    cout<<pic_dim<<endl;
    for(int i =0;i<pic_dim;i++)
    {
        for(int j = 0 ; j<=pic_dim;j++)
        {
            homogeneous_point p =  r*i*xr + u*j*-1*yr+top_left ;

            pixels[i][j].p = p;

        }
    }
    return pixels;

}

color getcolor(Ray R, int depth, int & r_id)
{
    color tempc,final_col;
    homogeneous_point intp,tempp;
    Vector normal,tempn;
    Vector need = R.v;
    double temp=1000000.0,t;
    int d_id=-1;
    for(int k =0;k<drawings.size();k++)
    {
        t = drawings[k].getint(R,tempc,tempp,tempn);
        if(t<temp && t>0){
            temp=t;
            final_col.r = tempc.r; final_col.g = tempc.g; final_col.b = tempc.b;
            intp.x = tempp.x; intp.y = tempp.y; intp.z = tempp.z;
            normal.x = tempn.x ; normal.y = tempn.y; normal.z = tempn.z;
            d_id = k;
        }
    }
    normal.normalize();
    if(d_id==-1) return final_col;
    double lambert=0, fhong = 0;
    r_id = d_id;
    for( int j =0;j<=lights.size();j++)
    {
        Ray R ;
        R.v.x = lights[j].pos.x-intp.x; R.v.y = lights[j].pos.y-intp.y; R.v.z = lights[j].pos.z-intp.z;
        R.p.x = intp.x ; R.p.y = intp.y ; R.p.z = intp.z ;
        R.p = R.v*.0001+R.p;
        double limit = R.v.mod();
        R.v.normalize();
        int flag = 0;
        for(int k =0;k<drawings.size();k++)
        {
            t = drawings[k].getint(R,tempc,tempp,tempn);
            if(t<limit && t>0) flag = 1;
        }
        if(flag==1) continue;
        if(!lights[j].checkspot(intp)){
                continue;
        }

        double scaling_factor = exp(-1*limit*limit*lights[j].falloff);
        lambert += max(R.v.dot(R.v,normal),0.0)*scaling_factor;
        Vector r =  ((normal * R.v.dot(R.v*-1, normal) * 2.0) - R.v*-1);
        r.normalize();

        fhong += pow( max(r.dot(r,need),0.0), drawings[d_id].shininess)*scaling_factor;
        /*Vector r =  ( need - (normal * need.dot(need, normal) * 2.0) );
        r.normalize();

        fhong += pow( max(r.dot(r,need),0.0), drawings[d_id].shininess)*scaling_factor;*/

    }
    Vector refl = ( R.v - (normal * R.v.dot(R.v, normal) * 2.0) );

    //drawings[d_id].diffuse = 1- drawings[d_id].ambient-drawings[d_id].specular;
    final_col.r = (drawings[d_id].ambient+lambert*drawings[d_id].diffuse+fhong*drawings[d_id].specular) *final_col.r ;
    final_col.g = (drawings[d_id].ambient+lambert*drawings[d_id].diffuse+fhong*drawings[d_id].specular) *final_col.g ;
    final_col.b = (drawings[d_id].ambient+lambert*drawings[d_id].diffuse+fhong*drawings[d_id].specular) *final_col.b ;

    if(depth!=0){
        Ray ray ;
        refl.normalize();
        ray.v = refl;
        int mY_id;
        ray.p.x = intp.x;   ray.p.y = intp.y;   ray.p.z = intp.z;
        ray.p = ray.v*.0001+ray.p ;
        color reflected_col = getcolor(ray,depth -1,mY_id );
        final_col.r +=  reflected_col.r*drawings[d_id].reflection;
        final_col.g +=  reflected_col.g*drawings[d_id].reflection;
        final_col.b +=  reflected_col.b*drawings[d_id].reflection;
    }

    return final_col;
}
void mainloop()
{
    buffer ** pixels = getbuffer();
    for(int i =0;i<pic_dim;i++)
    {
        for(int j = 0 ; j<pic_dim;j++)
        {
            Ray R;
            R.p = pixels[i][j].p;
            R.v.x = R.p.x - pos.x; R.v.y = R.p.y - pos.y; R.v.z = R.p.z - pos.z;
            R.v.normalize();
            int fake;
            color final_col = getcolor(R,rec_lev,fake);
            /*homogeneous_point intp,tempp;
            Vector normal,tempn;
            double temp=10000,t;
            for(int k =0;k<drawings.size();k++)
            {
                t = drawings[k].getint(R,tempc,tempp,tempn);
                if(t<temp && t>0){
                    temp=t;
                    final_col.r = tempc.r; final_col.g = tempc.g; final_col.b = tempc.b;
                    intp.x = tempp.x; intp.y = tempp.y; intp.z = tempp.z;
                    normal.x = tempn.x ; normal.y = tempn.y; normal.z = tempn.z;
                }
            }*/

            pixels[i][j].col.r = final_col.r; pixels[i][j].col.g = final_col.g; pixels[i][j].col.b = final_col.b;


        }
    }
    bitmap_image image(pic_dim, pic_dim);
    for (int x = 0; x < pic_dim; x++) {
        for (int y = 0; y < pic_dim; y++) {
            if(pixels[x][y].col.r>1 || pixels[x][y].col.g > 1 || pixels[x][y].col.b > 1 ){
                double fac = max( max(pixels[x][y].col.g,pixels[x][y].col.r),pixels[x][y].col.b);
                pixels[x][y].col.g = pixels[x][y].col.g/fac;
                pixels[x][y].col.b = pixels[x][y].col.b/fac;
                pixels[x][y].col.r = pixels[x][y].col.r/fac;
            }
            image.set_pixel(pic_dim-x-1, y, 255*pixels[x][y].col.r, 255*pixels[x][y].col.g,255* pixels[x][y].col.b);
        }
    }
    image.save_image("out.bmp");
    cout<<"done"<<endl;


}

void drawGrid()
{
	int i;

		glColor3f(0.6, 0.6, 0.6);	//grey
		glBegin(GL_LINES);{
			for(i=-8;i<=8;i++){

				if(i==0)
					continue;	//SKIP the MAIN axes

				//lines parallel to Y-axis
				glVertex3f(i*10, -90, 0);
				glVertex3f(i*10,  90, 0);

				//lines parallel to X-axis
				glVertex3f(-90, i*10, 0);
				glVertex3f( 90, i*10, 0);
			}
		}glEnd();

}

void drawSquare(double a)
{
    //glColor3f(1.0,0.0,0.0);
	glBegin(GL_QUADS);{
		glVertex3f( a, a,2);
		glVertex3f( a,-a,2);
		glVertex3f(-a,-a,2);
		glVertex3f(-a, a,2);
	}glEnd();
}










void keyboardListener(unsigned char key, int xx,int yy){
    double x,y,z;
    double rate = 0.01;
	switch(key){

		case '1':

			{
            x=l.x;y=l.y;z=l.z;
			l.x = l.x*cos(rate)+1*r.x*sin(rate);
			l.y = l.y*cos(rate)+r.y*sin(rate);
			l.z = l.z*cos(rate)+r.z*sin(rate);

			r.x = r.x*cos(rate)-x*sin(rate);
			r.y = r.y*cos(rate)-y*sin(rate);
			r.z = r.z*cos(rate)-z*sin(rate);}
			break;
        case '2':

			{
            x=l.x;y=l.y;z=l.z;
			l.x = l.x*cos(-rate)+r.x*sin(-rate);
			l.y = l.y*cos(-rate)+r.y*sin(-rate);
			l.z = l.z*cos(-rate)+r.z*sin(-rate);

			r.x = r.x*cos(-rate)-x*sin(-rate);
			r.y = r.y*cos(-rate)-y*sin(-rate);
			r.z = r.z*cos(-rate)-z*sin(-rate);
			}
			break;
        case '3':
            x=l.x;y=l.y;z=l.z;
			l.x = l.x*cos(rate)+u.x*sin(rate);
			l.y = l.y*cos(rate)+u.y*sin(rate);
			l.z = l.z*cos(rate)+u.z*sin(rate);

			u.x = u.x*cos(rate)-x*sin(rate);
			u.y = u.y*cos(rate)-y*sin(rate);
			u.z = u.z*cos(rate)-z*sin(rate);
			break;
        case '4':
            x=l.x;y=l.y;z=l.z;
			l.x = l.x*cos(-rate)+1*u.x*sin(-rate);
			l.y = l.y*cos(-rate)+u.y*sin(-rate);
			l.z = l.z*cos(-rate)+u.z*sin(-rate);

			u.x = u.x*cos(-rate)-x*sin(-rate);
			u.y = u.y*cos(-rate)-y*sin(-rate);
			u.z = u.z*cos(-rate)-z*sin(-rate);
			break;
        case '6':
            x=r.x;y=r.y;z=r.z;
			r.x = r.x*cos(rate)+u.x*sin(rate);
			r.y = r.y*cos(rate)+u.y*sin(rate);
			r.z = r.z*cos(rate)+u.z*sin(rate);

			u.x = u.x*cos(rate)-x*sin(rate);
			u.y = u.y*cos(rate)-y*sin(rate);
			u.z = u.z*cos(rate)-z*sin(rate);
			break;
        case '5':
            x=r.x;y=r.y;z=r.z;
			r.x = r.x*cos(-rate)+u.x*sin(-rate);
			r.y = r.y*cos(-rate)+u.y*sin(-rate);
			r.z = r.z*cos(-rate)+u.z*sin(-rate);

			u.x = u.x*cos(-rate)-x*sin(-rate);
			u.y = u.y*cos(-rate)-y*sin(-rate);
			u.z = u.z*cos(-rate)-z*sin(-rate);
			break;
        case '0':
            mainloop();
            break;
        case ' ':
            for(int i = 0 ; i<drawings.size();i++)
            {
               if(drawings[i].name == "cb")
               {
                  drawings[i].texer = true;
               }

            }
            break;


		default:
			break;
	}

}


void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_UP:		//down arrow key
			pos.x+=l.x;
			pos.y+=l.y;
			pos.z+=l.z;
			break;
		case GLUT_KEY_DOWN:		// up arrow key
			pos.x-=l.x;
			pos.y-=l.y;
			pos.z-=l.z;
			break;

		case GLUT_KEY_LEFT :
			pos.x+=r.x;
			pos.y+=r.y;
			pos.z+=r.z;
			break;
		case GLUT_KEY_RIGHT :
			pos.x-=r.x;
			pos.y-=r.y;
			pos.z-=r.z;
			break;

		case GLUT_KEY_PAGE_UP:
		    pos.x+=u.x;
			pos.y+=u.y;
			pos.z+=u.z;
			break;
		case GLUT_KEY_PAGE_DOWN:
            pos.x-=u.x;
			pos.y-=u.y;
			pos.z-=u.z;
			break;

		case GLUT_KEY_INSERT:
			break;

		case GLUT_KEY_HOME:
			break;
		case GLUT_KEY_END:
			break;

		default:
			break;
	}
}


void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_LEFT_BUTTON:
		    if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
			}
			break;


		case GLUT_RIGHT_BUTTON:
		    if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
				drawaxes=1-drawaxes;

			}
			break;
			//........

		case GLUT_MIDDLE_BUTTON:
			//........
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

	/********************
	/ set-up camera here
	********************/
	//load the correct matrix -- MODEL-VIEW matrix
	glMatrixMode(GL_MODELVIEW);

	//initialize the matrix
	glLoadIdentity();

	//now give three info
	//1. where is the camera (viewer)?
	//2. where is the camera looking?
	//3. Which direction is the camera's UP direction?

	//gluLookAt(100,100,100,	0,0,0,	0,0,1);
	//gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
	gluLookAt(pos.x,pos.y,pos.z,	pos.x+l.x,pos.y+l.y,pos.z+l.z,	u.x,u.y,u.z);


	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects



    //draw_cylinder(2,500,segment);
    //100*sin(angle_c_x)+(500-100*cos(angle_c_x))*sin(angle_t_x);


    drawAxes();

    for(int i = 0 ; i<drawings.size();i++)
    {
        drawings[i].draw();
    }
    for(int i = 0 ; i<lights.size();i++)
    {
        lights[i].draw();
    }


	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}


void animate(){
	glutPostRedisplay();
}

void init(){
	//codes for initialization

	drawaxes=1;


	//clear the screen
	glClearColor(0,0,0,0);

	/************************
	/ set-up projection here
	************************/
	//load the PROJECTION matrix
	glMatrixMode(GL_PROJECTION);

	//initialize the matrix
	glLoadIdentity();

	//give PERSPECTIVE parameters
	gluPerspective(fov_y,	ration,	n,	f);
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
}

int main(int argc, char **argv){

    ifstream scene;
    scene.open ("description.txt");
    string command;

    scene>>n>>f;
    scene>>fov_y>>ration;

    scene>>rec_lev>>pic_dim;
    //pic_dim = 50;
    int i ;

    drawable b ;
    b.name = "cb";
    scene>>b.width;
    scene>>b.ambient>>b.diffuse>>b.reflection;
    b.load();
    drawings.push_back(b);

    scene>>i;
    while(i!=0)
    {
        i--;
        scene>>command;
        cout<<command<<endl;
        if(command=="sphere")
        {
           homogeneous_point center;double r;color c;
           scene>>center.x>>center.y>>center.z;
           scene>>r;
           scene>>c.r>>c.g>>c.b;
           drawable d;
           d.name = "sphere";
           d.center = center;
           d.radius = r;
           d.col = c;
           scene>>d.ambient>>d.diffuse>>d.specular>>d.reflection;
           scene>>d.shininess;
           drawings.push_back(d);

        }
        else if(command=="pyramid")
        {
            homogeneous_point top,a,b,c,d;
            double h,w;
            color col;
            scene>>a.x>>a.y>>a.z;
            scene>>w>>h;
            b = a+homogeneous_point(w,0,0);
            d = a+homogeneous_point(0,w,0);
            c = a+homogeneous_point(w,w,0);
            top = a+homogeneous_point(w/2,w/2,0)+homogeneous_point(0,0,h);
            scene>>col.r>>col.g>>col.b;
            drawable D;
            D.name = "pyramid";
            D.top = top ; D.a = a; D.b = b; D.c = c;D.d= d;
            D.col = col;
            scene>>D.ambient>>D.diffuse>>D.specular>>D.reflection;
            scene>>D.shininess;
            drawings.push_back(D);
        }
    }
    for(int i = 0 ; i<drawings.size();i++)
    {
        drawings[i].print();

    }
    scene>>i;
    while(i!=0)
    {
        i--;
        light l ;
        scene>>l.pos.x>>l.pos.y>>l.pos.z>>l.falloff;
        l.name = "normal";
        lights.push_back(l);

    }

    scene>>i;
    while(i!=0)
    {
        i--;
        light l ;
        scene>>l.pos.x>>l.pos.y>>l.pos.z>>l.falloff;
        scene>>l.dir_pon.x>>l.dir_pon.y>>l.dir_pon.z>>l.angle;
        l.name = "spot";
        lights.push_back(l);

    }
    for( int i = 0 ;i<lights.size();i++)
    {
        cout<<lights[i].name<<" "<<lights[i].angle<<endl;
    }
    pos.x=0;
    pos.y=0;
    pos.z=50;
    l.x=0;u.x=1;r.x=0;
    l.y=0;u.y=0;r.y=1;
    l.z=-1;u.z=0;r.z=0;
	glutInit(&argc,argv);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("My OpenGL Program");


	init();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();		//The main loop of OpenGL
    scene.close();
	return 0;
}
