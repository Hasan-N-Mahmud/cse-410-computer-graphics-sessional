#include "custom.hpp"
#include <bits/stdc++.h>
#include <stdlib.h>
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include <sstream>
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
        color[0]=rand();
        color[1]=rand();
        color[2]=rand();
        }
	triangle(point x,point y,point z){
	    a=x;
	    b=y;
	    c=z;
	    color[0]=rand();
        color[1]=rand();
        color[2]=rand();
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
int main(){
            double temp,temp2,dz,scan_lx,scan_rx,scan_ty,scan_zu,scan_zd;
            triangle t(point(0,0,1),point(-2,-2,2),point(3,-3,3));
            scan_lx = 0;
            scan_rx = 0;
            scan_zu = 0;
            scan_zd = 0;
            scan_ty = -0.5;
        if((t.b.y - t.a.y) !=0){
            temp = t.a.x + (scan_ty - t.a.y) * (t.b.x -t.a.x) / (t.b.y - t.a.y);
            temp2 = t.a.z + (scan_ty - t.a.y) * (t.b.z -t.a.z) / (t.b.y - t.a.y);
            if(temp >= t.boundary_values[0] && temp <= t.boundary_values[1] ){
                if(temp < scan_lx)
                {
                    scan_lx = temp;
                    scan_zd = temp2;
                }
                else{
                    scan_rx = temp;
                    scan_zu = temp2;
                }
            }
        }
        if((t.b.y - t.c.y) !=0){
            temp = t.c.x + (scan_ty - t.c.y) * (t.b.x -t.c.x) / (t.b.y - t.c.y);
            temp2 = t.c.z + (scan_ty - t.c.y) * (t.b.z -t.c.z) / (t.b.y - t.c.y);
            if(temp >= t.boundary_values[0] && temp <= t.boundary_values[1] ){
                if(temp < scan_lx)
                {
                    scan_lx = temp;
                    scan_zd = temp2;
                }
                else{
                    scan_rx = temp;
                    scan_zu = temp2;
                }
        }
    }
        if((t.c.y - t.a.y) !=0){
            temp = t.a.x + (scan_ty - t.a.y) * (t.c.x -t.a.x) / (t.c.y - t.a.y);
            temp2 = t.a.z + (scan_ty - t.a.y) * (t.c.z -t.a.z) / (t.c.y - t.a.y);
            if(temp >= t.boundary_values[0] && temp <= t.boundary_values[1] ){
                if(temp < scan_lx)
                {
                    scan_lx = temp;
                    scan_zd = temp2;
                }
                else{
                    scan_rx = temp;
                    scan_zu = temp2;
                }
        }
}
        cout<<scan_lx<<" "<<scan_rx<<" "<<scan_zd<<" "<<scan_zu<<endl;
        }
