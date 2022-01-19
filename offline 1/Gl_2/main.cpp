#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<dos.h>
#include<iostream>
#include <windows.h>
#include <GL/glut.h>
#include<chrono>
#define pi (2*acos(0.0))
using namespace std;
using namespace std::chrono;

double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
bool created;
bool paused;
bool insideCircle[5];
bool overlap[5][5];
double smallRadius = 10;
double largeRadius = 60;
int numberofBubbles = 5;
time_t startTime,currentTime;
double speedConstant = 0.5;
struct point
{
	double x,y,z;
	point(){
	x=0;
	y=0;
	z=0;
	}
	point(double a,double b,double c){
	    x = a;
	    y = b;
	    z = c;
	}
    void print(){
        cout<<"X: "<<x<<" Y : "<<y<<" Z: "<<z<<endl;
    }
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
point pos;
double length ;
double a = 4;
point position[5];
point direction[5];
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
	glBegin(GL_LINES);{
		glVertex3f( a, a,2);
		glVertex3f( a,-a,2);

		glVertex3f(-a,-a,2);
		glVertex3f(-a, a,2);

		glVertex3f( -a, -a,2);
		glVertex3f( a,-a,2);

		glVertex3f( a, a,2);
		glVertex3f( -a,a,2);
	}glEnd();
}


void drawCircle(double radius,int segments)
{
    int i;
    struct point points[100];
   // glColor3f(0.7,0.7,0.7);
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


void keyboardListener(unsigned char key, int x,int y){
	switch(key){
    case 'p':
        paused = !paused;
        break;
		default:
			break;
	}
}


void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
			//cameraHeight -= 3.0;
			if(speedConstant > 0.2)
                speedConstant -= 0.1;

			break;
		case GLUT_KEY_UP:		// up arrow key
			//cameraHeight += 3.0;
            if(speedConstant < 3)
                speedConstant += 0.1;
			break;

		case GLUT_KEY_RIGHT:
			//cameraAngle += 0.03;

			break;
		case GLUT_KEY_LEFT:
			//cameraAngle -= 0.03;

			break;

		case GLUT_KEY_PAGE_UP:

			break;
		case GLUT_KEY_PAGE_DOWN:

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
    gluLookAt(pos.x, pos.y, pos.z,   0, 0, 0,   0, 1, 0);

	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects

	drawAxes();
	drawGrid();

    glColor3f(0,1,0);
    drawSquare(length);

    glColor3f(1,0,0);
    drawCircle(60,50);

    glColor3f(0.7,0.7,0.7);
    time(&currentTime);
    float duration = difftime (currentTime,startTime);
    //cout<<duration<<endl;
    if(duration < 5) {
        numberofBubbles =1+duration;
    }
    for(int i=0;i<numberofBubbles;i++){
     //glColor3f(0.1*(i+1),0.1*(i+1),0.1*(i+1));
    glPushMatrix();
    glTranslatef(position[i].x,position[i].y,position[i].z);
    drawCircle(smallRadius,20);
    glPopMatrix();
    }
    created = true;
   // drawSS();

    //drawCircle(30,24);

    //drawCone(20,50,24);

	//drawSphere(30,24,20);


	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}
double distance(point p1,point p2){
    return sqrt(((p1.x - p2.x) * (p1.x - p2.x)) + ((p1.y -p2.y) * (p1.y - p2.y)));
}

point normalize(point p){
            double value = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
            point res;
            res.x = p.x  / value;
            res.y = p.y / value;
            res.z = p.z / value;
            return res;
}

point getReflected(point d,point n){
    point r;
    //cout<<n.x<<endl<<n.y<<endl<<n.z<<endl;
    double c = (2 * dotProduct(d,n) );
    r = d - n * c ;
    return r;
}

void animate(){
    if(!paused){
        for(int i=0;i<numberofBubbles;i++){
            if(sqrt((position[i].x * position[i].x) + (position[i].y * position[i].y) ) <= (60 - smallRadius)){
                insideCircle[i] = true;}
            if(!insideCircle[i]){
                if(position[i].x >= (length-smallRadius) && abs(int(position[i].x)) > (length-smallRadius)){
                    direction[i] = getReflected(direction[i], point(-1,0,0));}
                else if(position[i].x <= -(length-smallRadius) && abs(int(position[i].x)) > (length-smallRadius)){
                    direction[i] = getReflected(direction[i], point(1,0,0));}
                if(position[i].y >= (length-smallRadius) && abs(int(position[i].y)) > (length-smallRadius)){
                    direction[i] = getReflected(direction[i], point(0,-1,0));}
                else if(position[i].y <= -(length-smallRadius) && abs(int(position[i].y)) > (length-smallRadius)){
                    direction[i] = getReflected(direction[i], point(0,1,0));}
    }
            if(insideCircle[i]){
//                if(overlap[i][i] && (sqrt((position[i].x * position[i].x) + (position[i].y * position[i].y) ) < (60 - smallRadius)))
//                   overlap[i][i] = false;
                if(sqrt((position[i].x * position[i].x) + (position[i].y * position[i].y) ) >= (60 - smallRadius) && !overlap[i][i]){
//                    overlap[i][i] = true;
                    point temp;
                    temp.x =  position[i].x;
                    temp.y =  position[i].y;
                    temp.z = 0;
                    direction[i] = getReflected(direction[i],normalize(temp));
                    while(distance(position[i],point(0,0,0)) >=  (60 - smallRadius)){
                        position[i].x +=  direction[i].x;
                        position[i].y +=  direction[i].y;
                    }
                    }
//                     if((sqrt((position[i].x * position[i].x) + (position[i].y * position[i].y) ) > (60 - smallRadius)) && !overlap[i][i])
//                        overlap[i][i] = true;

                for(int j=0;j<numberofBubbles;j++){
                    if(j == i)
                        continue;
                    if((distance(position[i] , position[j]) - (2.0 * smallRadius) < -1 ) && !overlap[i][j]){
                        overlap[i][j] = true;
                        overlap[j][i] =true;
                       }
                       else if(overlap[i][j]  && (distance(position[i] , position[j]) - (2.0 * smallRadius) > 0 )){
                       overlap[i][j] =false;
                       overlap[j][i] =false;
                       }
                    if((distance(position[i] , position[j]) - (2.0 * smallRadius) <= 0  )  && insideCircle[j]  && !overlap[i][j]  ){
                        //cout<<i<<" "<<j<<endl;
                        point normal;
                        normal = position[i] - position[j];
                        normal = normalize(normal);
//                        cout<<"before: ";
//                        position[i].print();
//                        position[j].print();
//                        cout<<"Direction: ";
//                        direction[i].print();
//                        cout<<"Normal: ";
//                        normal.print();
                        direction[i] = getReflected(direction[i],normal);
//                        cout<<"After: ";
//                        direction[i].print();
//                        cout<<endl;
                    }
                }
            }

    }
    //cout<< "End of loop" <<endl;
    for(int i=0;i<numberofBubbles;i++){
        position[i].x += speedConstant * direction[i].x;
        position[i].y += speedConstant * direction[i].y;
    }
    }
   //printf("%f %f\n",direction[i].x,direction[i].y);

	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void init(){
	//codes for initialization
	drawgrid=0;
	drawaxes=0;
	cameraHeight=150.0;
	cameraAngle=1.0;


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
	pos.x = 0;
	pos.y = 0;
	pos.z = 150;
	length=100;
    srand(time(0));
	for(int i=0;i<numberofBubbles;i++){
    position[i].x = -(length-smallRadius);
    position[i].y = -(length-smallRadius);
    direction[i].x = rand()*0.1/RAND_MAX;
    direction[i].y = rand()*0.1/RAND_MAX;
    printf("%f %f\n",direction[i].x,direction[i].y);
    insideCircle[i] = false;
	}
    created = false;
    paused = false;
    time(&startTime);
//    direction[0].x = 0.02;
//    direction[0].y = 0.01;

}

int main(int argc, char **argv){
	glutInit(&argc,argv);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("Task 3");

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

