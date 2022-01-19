#include <bits/stdc++.h>
#include <stdlib.h>
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include <sstream>
using namespace std;
const double PI = 3.14159;
struct Pixel{
    double z;
    int rgb[3];
    Pixel(){
        z=0;
        set_color(0,0,0);
    }
    Pixel(double value){
        z=value;
        set_color(0,0,0);
    }
   void set_color(double r,double g,double b){
        rgb[0] = r;
        rgb[1] = g;
        rgb[2] = b;
    }
    void set_z(double value){
        z = value;
    }
    void print(){
        cout<<"z value: "<<z<<" r: "<<rgb[0]<<" g: "<<rgb[1]<<" b: "<<rgb[2]<<"\t";
    }
    friend ostream &operator<<( ostream &output, const Pixel &p ) {
         output <<p.z<<"\t";
         return output;
      }
};
struct Z_buffer{
    int row,col;
    Pixel * buffer;

Z_buffer(int a,int b,double z){
    row = a;
    col = b;
    cout<<row*col<<endl;
    buffer = new Pixel[row*col];
        for (int i = 0; i < row; i++)
        {
            for (int j = 0; j < col; j++)
                {
                    buffer[i*row + j].set_z(z);
                }
        }
    }
void print(){
    for (int i = 0; i < row; i++)
        {
            for (int j = 0; j < col; j++)
                {
                    buffer[i*row+j].print();
                }
            cout << endl;
    }
}
};

struct point
{
	double x,y,z,w;
	point(){
        x=y=z=0;
        w=1;
        }
	point(double a,double b,double c){
	    x=a;
	    y=b;
	    z=c;
	    w=1;
	}
	point(double a,double b,double c,double d){
	    x=a;
	    y=b;
	    z=c;
	    w=d;
	}

	point operator+ (point p) {
	    point res;
	    res.x = x + p.x;
	    res.y = y + p.y;
	    res.z = z + p.z;
	    return res; }
    void print(){
        cout<<x<<"\t"<<y<<"\t"<<z<<"\t"<<w<<endl;
    }
    void normalize(){
        if(w!=1){
            x /= w;
            y /= w;
            z /= w;
            w /= w;
        }
    }
    void precision(){
        x =roundf(x * 10000) / 10000.0 ;
        y =roundf(y * 10000) / 10000.0;
       z = roundf(z * 10000) / 10000.0;
    }

};

struct triangle
{
	point a,b,c;
	double boundary_values[4];
	int color[3];
	triange(){
        color[0]=rand()%256;
        color[1]=rand()%256;
        color[2]=rand()%256;
        }
	triangle(point x,point y,point z){
	    a=x;
	    b=y;
	    c=z;
	    color[0]=rand()%256;
        color[1]=rand()%256;
        color[2]=rand()%256;
        set_max_values();
	}
void	set_max_values(){
        boundary_values[0]=a.x;
        boundary_values[1]=a.x;
        boundary_values[2]=a.y;
        boundary_values[3]=a.y;

        if(boundary_values[0] > b.x && b.x < c.x )
            boundary_values[0]=b.x;
        else if(boundary_values[0] > c.x)
            boundary_values[0]=c.x;

        if(boundary_values[1] < b.x && b.x > c.x )
            boundary_values[1]=b.x;
        else if(boundary_values[1] < c.x)
            boundary_values[1]=c.x;

        if(boundary_values[3] < b.y && b.y > c.y )
            boundary_values[3]=b.y;
        else if(boundary_values[3] < c.y)
            boundary_values[3]=c.y;

        if(boundary_values[2] > b.y && b.y < c.y )
            boundary_values[2]=b.y;
        else if(boundary_values[2] > c.y)
            boundary_values[2]=c.y;
    }
friend ostream &operator<<( ostream &output, const triangle &t ) {
         output << fixed<<setprecision(6)<< t.a.x << " " << t.a.y << " "<<t.a.z<<"\n"<< t.b.x << " " << t.b.y << " "<<t.b.z<<"\n"<<  t.c.x << " " << t.c.y << " "<<t.c.z<<"\n";
         return output;
      }
    void print(){
        cout<<"triangle printing"<<endl;
        a.print();
        b.print();
        c.print();
        cout<<endl;
    }
    void precision(){
        a.precision();
        b.precision();
        c.precision();
        }
};

struct matrix{
    int row,col;
    double * mat;
    matrix(int a,int b){
    row = a;
    col = b;
    mat = new double[row*col];
        for (int i = 0; i < row; i++)
        {
            for (int j = 0; j < col; j++)
                {
                    mat[i*col+j] =0;
                }
    }
    }
matrix(int a,int b,double z,string str){
    row = a;
    col = b;
    mat = new double[row*col];
        for (int i = 0; i < row; i++)
        {
            for (int j = 0; j < col; j++)
                {
                    mat[i*col+j] =z;
                }
    }
    }
    matrix(int a,int b,string str){
    row = a;
    col = b;
    mat = new double[row*col];
        for (int i = 0; i < row; i++)
        {
            for (int j = 0; j < col; j++)
                {
                    if(i == j)
                        mat[i*col+j] =1;
                    else
                    mat[i*col+j] =0;
                }
    }
    }
    void print(){
    for (int i = 0; i < row; i++)
        {
            for (int j = 0; j < col; j++)
                {
                    cout << mat[i*col+j] << " ";
                }
            cout << endl;
    }
}

matrix operator*(matrix m2){
    matrix res(row,col);
    int i, j, k;
    for (i = 0; i < row; i++) {
        for (j = 0; j < col; j++) {
            for (k = 0; k < col; k++)
                res.mat[i*col+j] += mat[i*col+k] * m2.mat[k*col+j];
        }
    }
    return res;
}

point operator*(point p){
    matrix res(4,1);
    double arr[4];
    arr[0]=p.x;
    arr[1]=p.y;
    arr[2]=p.z;
    arr[3]=p.w;
    int i, j;
    for (i = 0; i < row; i++) {
        for (j = 0; j < col; j++) {
                res.mat[i] += mat[i*col+j] * arr[j];
        }
    }
    point r(res.mat[0],res.mat[1],res.mat[2],res.mat[3]);
    r.normalize();
    return r;
}
triangle operator*(triangle t){
    t.a = *(this) * t.a;
    t.b = *(this) * t.b;
    t.c = *(this) * t.c;
    return t;
}
bool operator==(matrix m2){
    int i,j;
        for (i = 0; i < row; i++) {
            for (j = 0; j < col; j++) {
                if(mat[i*col+j] != m2.mat[i*col+j])
                    return false;
        }
    }
    return true;
}
bool is_identity_matrix(){
    int i,j;
        for (i = 0; i < row; i++) {
            for (j = 0; j < col; j++) {
                if((i==j && mat[i*col+j] !=1) || (i != j && mat[i*col+j] !=0))
                    return false;
        }
    }
    return true;
}
};
