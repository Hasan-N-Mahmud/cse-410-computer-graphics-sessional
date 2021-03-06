#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include <windows.h>
#include <GL/glut.h>

#define pi (2*acos(0.0))

double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
double angle;

struct point
{
	double x,y,z;
	point operator+ (point p) {
	    point res;
	    res.x = x + p.x;
	    res.y = y + p.y;
	    res.z = z + p.z;
	    return res; }
	    point operator*(double k) {
	    point res;
	    res.x = x * k;
	    res.y = y * k;
	    res.z = z * k;
	    return res; }
	    point operator-(point p) {
	    point res;
	    res.x = x - p.x;
	    res.y = y - p.y;
	    res.z = z - p.z;
	    return res; }
	    bool operator==(point p) {
	   if((x == p.x)&&(y == p.y)&&(z == p.z))
        return true;
        else{
            return false;
        }}
};
point pos,u,r,l;
double k =2.0 ;
double a = 4;
double b,c,d,e;
double dotProduct(point vect_A, point vect_B)
{
    double product = 0;
    product += vect_A.x * vect_B.x;
    product += vect_A.y * vect_B.y;
    product += vect_A.z * vect_B.z;
    product+=1;
    product -=1;
    return product;
}

point crossProduct(point vect_A, point vect_B)
  {
     point cross_P;
     cross_P.x = vect_A.y * vect_B.z - vect_A.z * vect_B.y;
     cross_P.y = vect_A.z * vect_B.x - vect_A.x * vect_B.z;
     cross_P.z = vect_A.x * vect_B.y - vect_A.y * vect_B.x;
     return cross_P;
}
void drawAxes()
{
	if(drawaxes==1)
	{
		//glColor3f(1.0, 1.0, 1.0);
		glBegin(GL_LINES);{
		    glColor3f(1,0,0);
			glVertex3f( 100,0,0);
			glVertex3f(-100,0,0);

            glColor3f(0,1,0);
			glVertex3f(0,-100,0);
			glVertex3f(0, 100,0);

            glColor3f(0,0,1);
			glVertex3f(0,0, 100);
			glVertex3f(0,0,-100);
		}glEnd();
	}
}


void drawGrid()
{
	int i;
	if(drawgrid==1)
	{
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


void drawCircle(double radius,int segments)
{
    int i;
    struct point points[100];
    glColor3f(0.7,0.7,0.7);
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)*2*pi);
        points[i].y=radius*sin(((double)i/(double)segments)*2*pi);
    }
    //draw segments using generated points
    for(i=0;i<segments;i++)
    {
        glBegin(GL_LINES);
        {
			glVertex3f(points[i].x,points[i].y,0);
			glVertex3f(points[i+1].x,points[i+1].y,0);
        }
        glEnd();
    }
}

void drawCone(double radius,double height,int segments)
{
    int i;
    double shade;
    struct point points[100];
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)*2*pi);
        points[i].y=radius*sin(((double)i/(double)segments)*2*pi);
    }
    //draw triangles using generated points
    for(i=0;i<segments;i++)
    {
        //create shading effect
        if(i<segments/2)shade=2*(double)i/(double)segments;
        else shade=2*(1.0-(double)i/(double)segments);
        glColor3f(shade,shade,shade);

        glBegin(GL_TRIANGLES);
        {
            glVertex3f(0,0,height);
			glVertex3f(points[i].x,points[i].y,0);
			glVertex3f(points[i+1].x,points[i+1].y,0);
        }
        glEnd();
    }
}


void drawSphere(double radius,int slices,int stacks)
{
	struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=radius*cos(((double)i/(double)stacks)*(pi/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
			points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
			points[i][j].z=h;
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
        glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
		for(j=0;j<slices;j++)
		{
			glBegin(GL_QUADS);{
			    //upper hemisphere
				glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                //lower hemisphere
                glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
			}glEnd();
		}
	}
}

void drawCylinder(double radius,int slices,int stacks)
{
	struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=radius;
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
			points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
			points[i][j].z=h;
		}
	}
	//draw quads using generated point
	bool stripe = true;
	for(i=0;i<stacks;i++)
	{

		for(j=0;j<slices;j++)
		{
			glBegin(GL_QUADS);
			{
                //lower hemisphere
                if(stripe){
//        glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
                glColor3f(1,1,1);
                stripe = false;
                }
        else{
            glColor3f(0,0,0);
            stripe = true;
        }
                glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);

				glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
			}glEnd();
		}
	}
}

void drawUpperSphere(double radius,int slices,int stacks)
{
	struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=radius*cos(((double)i/(double)stacks)*(pi/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
			points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
			points[i][j].z=h;
		}
	}
	//draw quads using generated points
	bool stripe = true;
	for(i=0;i<stacks;i++)
	{

		for(j=0;j<slices;j++)
		{
			glBegin(GL_QUADS);{
			    if(stripe){
//        glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
                glColor3f(1,1,1);
                stripe = false;
                }
        else{
            glColor3f(0,0,0);
            stripe = true;
        }
			    //upper hemisphere
				glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);

			}glEnd();
		}
	}
}

void drawLowerSphere(double radius,int slices,int stacks)
{
	struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=radius*cos(((double)i/(double)stacks)*(pi/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
			points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
			points[i][j].z=h;
		}
	}
	//draw quads using generated points
	bool stripe = true;
	for(i=0;i<stacks;i++)
	{

		for(j=0;j<slices;j++)
		{
			glBegin(GL_QUADS);
			{
                //lower hemisphere
                if(stripe){
//        glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
                glColor3f(1,1,1);
                stripe = false;
                }
        else{
            glColor3f(0,0,0);
            stripe = true;
        }
                glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
			}glEnd();
		}
	}
}

void drawFunel(double radius,int slices,int stacks)
{
	struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=2*radius-radius*cos(((double)i/(double)stacks)*(pi/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
			points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
			points[i][j].z=h;
		}
	}
	//draw quads using generated points
	bool stripe = true;
	for(i=0;i<stacks;i++)
	{

		for(j=0;j<slices;j++)
		{
			glBegin(GL_QUADS);{
			    if(stripe){
//        glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
                glColor3f(1,1,1);
                stripe = false;
                }
        else{
                glColor3f(0,0,0);
            stripe = true;
        }
			     glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
			}glEnd();
		}
	}
}
void drawSS()
{
    glColor3f(1,0,0);
    drawSquare(20);

    glRotatef(angle,0,0,1);
    glTranslatef(110,0,0);
    glRotatef(2*angle,0,0,1);
    glColor3f(0,1,0);
    drawSquare(15);

    glPushMatrix();
    {
        glRotatef(angle,0,0,1);
        glTranslatef(60,0,0);
        glRotatef(2*angle,0,0,1);
        glColor3f(0,0,1);
        drawSquare(10);
    }
    glPopMatrix();

    glRotatef(3*angle,0,0,1);
    glTranslatef(40,0,0);
    glRotatef(4*angle,0,0,1);
    glColor3f(1,1,0);
    drawSquare(5);
}

void drawBullet(){

}
point rotateAxis(point k,point v,double a){
    point res;
    res = crossProduct(k,v) * sin(pi/180*a);
    res = res + v * cos(pi/180*a);
    return res;
}
void keyboardListener(unsigned char key, int x,int y){
	switch(key){

		case '1':
		    r = rotateAxis(u,r,a);
		    l = rotateAxis(u,l,a);
		    //printf("%f %f %f",crossProduct(r,l).x,crossProduct(r,l).y,crossProduct(r,l).z);
			break;
        case '2':
            r = rotateAxis(u,r,-a);
		    l = rotateAxis(u,l,-a);
			break;
        case '3':
            l = rotateAxis(r,l,a);
		    u = rotateAxis(r,u,a);
			break;
        case '4':
             l = rotateAxis(r,l,-a);
		    u = rotateAxis(r,u,-a);
			break;
        case '5':
             r = rotateAxis(l,r,-a);
		    u = rotateAxis(l,u,-a);
			break;
        case '6':
		     r = rotateAxis(l,r,a);
		    u = rotateAxis(l,u,a);
			break;
        case 'q':
            if(b < 45)
		     b += 1;
			break;
        case 'w':
            if(b > -45)
            b -= 1;
			break;
        case 'e':
            if(c < 45)
            c += 1;
			break;
        case 'r':
            if(c > -45)
            c -= 1;
			break;
        case 'a':
            if(d < 30)
            d += 1;
            break;
        case 's':
            if(d > -30)
            d -= 1;
            break;
        case 'd':
            e += 1;
            break;
        case 'f':
            e -= 1;
            break;
		default:
			break;
	}
}


void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
			//cameraHeight -= 3.0;
//			pos.x -= k * l.x ;
//			pos.y -= k*l.y;
//			pos.z -= k*l.z;
                pos = pos - l * k;
			break;
		case GLUT_KEY_UP:		// up arrow key
			//cameraHeight += 3.0;
//			pos.x += k * l.x ;
//			pos.y += k*l.y;
//			pos.z += k*l.z;
                pos = pos + l * k;
			break;

		case GLUT_KEY_RIGHT:
			//cameraAngle += 0.03;
//			pos.x += k * r.x ;
//			pos.y += k*r.y;
//			pos.z += k*r.z;
			pos = pos + r * k;
			break;
		case GLUT_KEY_LEFT:
			//cameraAngle -= 0.03;
//			pos.x -= k * r.x ;
//			pos.y -= k*r.y;
//			pos.z -= k*r.z;
                pos = pos - r * k;
			break;

		case GLUT_KEY_PAGE_UP:
//		    pos.x += k * u.x ;
//			pos.y += k*u.y;
//			pos.z += k*u.z;
                pos = pos + u * k;
			break;
		case GLUT_KEY_PAGE_DOWN:
//		    pos.x -= k * u.x ;
//			pos.y -= k*u.y;
//			pos.z -= k*u.z;
			pos = pos - u*k;
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
                drawBullet();
			}
			break;

		case GLUT_RIGHT_BUTTON:
			if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
				drawaxes=1-drawaxes;
				drawgrid =1-drawgrid;
			}
			break;

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
	//gluLookAt(0,0,200,	0,0,0,	0,1,0);
    gluLookAt(pos.x, pos.y, pos.z, pos.x + l.x, pos.y + l.y, pos.z + l.z, u.x, u.y, u.z);

	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects

	drawAxes();
	drawGrid();

    //glColor3f(1,0,0);
    //drawSquare(10);

   // drawSS();

    //drawCircle(30,24);

    //drawCone(20,50,24);

	//drawSphere(30,24,20);

	glPushMatrix();
{
    glRotatef(b,0,1,0);
    drawUpperSphere(30,24,20);
    glRotatef(c,1,0,0);
    drawLowerSphere(30,24,20);
    glTranslatef(0,0,-45);
    glRotatef(d,1,0,0);
    glRotatef(e,0,0,1);
    drawUpperSphere(15,24,20);
    glTranslatef(0,0,-15);
    drawCylinder(15,24,10);
    glTranslatef(0,0,-30);
    drawCylinder(15,24,10);
    glTranslatef(0,0,-15);
    drawFunel(15,24,20);
}
	glPopMatrix();

    glTranslatef(0,0,-200);
    glColor3f(0.698,0.745,0.7);
    drawSquare(80);
	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}


void animate(){
	angle+=0.05;
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void init(){
	//codes for initialization
	drawgrid=1;
	drawaxes=1;
	cameraHeight=150.0;
	cameraAngle=1.0;
	angle=0;

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
	gluPerspective(80,	1,	1,	1000.0);
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
	pos.x = 100;
	pos.y = 100;
	pos.z = 0;
    u.x = 0;
    u.y = 0;
    u.z = 1;
    r.x = -(1.0/sqrt(2));
    r.y = (1.0 / sqrt(2));
    r.z = 0;
    l.x = -(1.0/sqrt(2));
    l.y = -(1.0/sqrt(2));
    l.z = 0;
}

int main(int argc, char **argv){
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

	return 0;
}

