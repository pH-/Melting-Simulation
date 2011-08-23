//============================================================================
// Name        : MeltingIce.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "includes.h"
using namespace std;
object *iceCube;
kdtree *kdTree;
heatsrc *heatSource;
int global=0;
int iteratorCounter=0;
double camMotionStep=0.5;
double eyeX=0, eyeY=15, eyeZ=50;
double lookatx=0, lookaty=15, lookatz=0;
//ofstream outFileVariable;

void init(void)
{
	glClearColor (0.0,0.0,0.0, 0.0);
	glLoadIdentity();
}

void display(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();
	gluLookAt(eyeX,eyeY,eyeZ, lookatx,lookaty,lookatz, 0,1,0);
	deque<particle>::iterator iter;
	deque<particle> particleList = iceCube->getParticleList();
	stringstream data;
	stringstream fileName;
	if(!iteratorCounter)
	{
		fileName<<"IceMeltFrames/MeltFrame"<<iteratorCounter<<".txt";
		ofstream outFile(fileName.str().c_str());
		//outFile<<"\"Frame\",0,";
		outFile.flush();
		//outFileVariable.seekp(outFile.tellp());
		//outFile<<"          ,"<<endl;
		//outFile.flush();
		for(unsigned int i=0; i<particleList.size(); i++)
		{
			data<<particleList[i].getId()<<","<<particleList[i].getPosition()[0]<<","<<particleList[i].getPosition()[1]<<","
			       <<particleList[i].getPosition()[2]<<","<<particleList[i].getTemp()<<","<<particleList[i].getStatus()<<","<<endl;
		}
		//double numBytes = outFile.tellp()-outFileVariable.tellp();
		//cout<<"Number of bytes written"<<numBytes<<endl;
//		outFile<<data.str().length();
//		outFile<<","<<endl;
		outFile<<data.str()<<endl;
		outFile.close();
		//outFileVariable<<outFile.tellp()-outFileVariable.tellp();
		//outFileVariable.flush();
	}
	glBegin(GL_LINES);
		glColor3f(1,0,0);
		glVertex3f(0,0,0);
		glVertex3f(10,0,0);
		glColor3f(0,1,0);
		glVertex3f(0,0,0);
		glVertex3f(0,10,0);
		glColor3f(0,0,1);
		glVertex3f(0,0,0);
		glVertex3f(0,0,10);
	glEnd();
	iceCube->calcDensities();
	heatSource->clearPhotonList();
	for(int pc=0; pc<PHOTONSPERTIMESTEP; pc++)
		heatSource->shootPhotons();
	deque<photon> photons;
	photons = heatSource->getPhotonList();
	for(unsigned int i=0;i<photons.size(); i++)
	{
		//heatSource->tracePhoton(photons[i]);
		heatSource->tracephotonByRay(photons[i]);
	}
	iceCube->calcTemp();
//	iceCube->calcDensities();
	for(iter = particleList.begin();iter!=particleList.end();++iter)
	{
//		if(iter->getParticleId()==-1)
//			continue;
		double colr[3];
		colr[0]=iter->getColor()[0];
		colr[1]=iter->getColor()[1];
		colr[2]=iter->getColor()[2];
		glPushMatrix();
			glColor3dv(iter->getColor());
			glTranslated(iter->getPosition()[0]/DISPLAYSCALE,iter->getPosition()[1]/DISPLAYSCALE, iter->getPosition()[2]/DISPLAYSCALE);
//			cout<<iter->getId()<<":"<<iter->getTemp()<<endl;
			glutSolidSphere (0.2, 16, 16);
		glPopMatrix();
		//cout<<"colors"<<colr[0]<<":"<<colr[1]<<":"<<colr[2]<<endl;
	}
	glPushMatrix();
		glColor3d(1.0,1.0,0.0);
		glTranslated(HEATSRCPOSX/DISPLAYSCALE,HEATSRCPOSY/DISPLAYSCALE,HEATSRCPOSZ/DISPLAYSCALE);
		glutSolidSphere(.4,16,16);
	glPopMatrix();
	fileName.str("");
	iceCube->calcForces();
	iceCube->calcPosition();
	iteratorCounter++;
	fileName<<"IceMeltFrames/MeltFrame"<<iteratorCounter<<".txt";
	ofstream outFile(fileName.str().c_str());
	//outFile<<"\"Frame\","<<iteratorCounter<<",";
	//outFileVariable.seekp(outFile.tellp());
	//outFile<<"          ,"<<endl;
	//outFile.flush();
	data.str("");
	for(unsigned int i=0; i<particleList.size(); i++)
	{
		//if(particleList[i].getStatus()==LIQUID || !i)
			data<<particleList[i].getId()<<","<<particleList[i].getPosition()[0]<<","<<particleList[i].getPosition()[1]<<","
			       <<particleList[i].getPosition()[2]<<","<<particleList[i].getTemp()<<","<<particleList[i].getStatus()<<","<<endl;
	}
	//outFile.flush();
//	outFileVariable<<outFile.tellp()-outFileVariable.tellp();
//	outFileVariable.flush();
	cout<<"frame count"<<iteratorCounter<<endl;
	//cout<<"bytes read"<<data.str().length();
	outFile<<data.str()<<endl;
	data.str("");
	outFile.close();
	glutSwapBuffers();
	if(iteratorCounter==53)
	{
		cout<<"browse for particle 105"<<endl;
	}
	if(iteratorCounter>700)
	{
	//	exit(0);
	}
	if(global)
		glutPostRedisplay();
}

void keyboard(unsigned char key, int x, int y)
{
	switch (key) {
		case 'c':
		case 'C':
			global=1;
			glutPostRedisplay();
			break;
		case 's':
			global = 0;
			glutPostRedisplay();
			break;
		case 27:
			exit(0);
			break;
		case 'x':
			eyeX-=camMotionStep;
			glutPostRedisplay();
			break;
		case 'X':
			eyeX+=camMotionStep;
			glutPostRedisplay();
			break;
		case 'u':
			lookatx-=camMotionStep;
			glutPostRedisplay();
			break;
		case 'U':
			lookatx+=camMotionStep;
			glutPostRedisplay();
			break;
		case 'Y':
			eyeY+=camMotionStep;
			glutPostRedisplay();
			break;
		case 'y':
			eyeY-=camMotionStep;
			glutPostRedisplay();
			break;
		case 'v':
			lookaty-=camMotionStep;
			glutPostRedisplay();
			break;
		case 'V':
			lookaty+=camMotionStep;
			glutPostRedisplay();
			break;
		case 'z':
			eyeZ-=camMotionStep;
			glutPostRedisplay();
			break;
		case 'Z':
			eyeZ+=camMotionStep;
			glutPostRedisplay();
			break;
		case 'w':
			lookatz-=camMotionStep;
			glutPostRedisplay();
			break;
		case 'W':
			lookatz+=camMotionStep;
			glutPostRedisplay();
			break;
		default:
			break;
	}
}


void reshape(int w, int h)
{
	/*
	glViewport(0, 0, (GLsizei) w, (GLsizei) h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective (50.0, (GLdouble)w/(GLdouble)h, 3.0, 50.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();*/
	glViewport (0, 0, (GLsizei) w, (GLsizei) h);
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity();
	gluPerspective (50.0, (GLdouble)w/(GLdouble)h, 0.0, 100.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

int main(int argc, char** argv)
{
	iceCube = new object();
	kdTree = new kdtree();
	double heatSourcePos[3]={HEATSRCPOSX,HEATSRCPOSY,HEATSRCPOSZ};
	double photonDir[3];
	double sizx=10,sizy=10,sizz=10;
	cout<<"hello World"<<endl;
	srand ( time(NULL) );
	//iceCube->setDims(sizx,sizy,sizz);
	iceCube->insertParticles(sizx,sizy,sizz);
	MINUS(photonDir,iceCube->getBboxDims()[1],heatSourcePos);
	heatSource = new heatsrc();
	heatSource->setHeatSrcPos(heatSourcePos);
	heatSource->setCubeDimensions(iceCube->getBboxDims());
	heatSource->setPhotonDirNorm(photonDir);
	//iceCube->insertParticles();
	kdTree->buildTree(iceCube->getParticleList(),0);
	cout<<kdTree->getRoot()->getKeyValue()<<endl;
	//outFileVariable.open("IceMeltAnim.txt");
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition (100, 100);
	glutCreateWindow("Variational Tetrahedralization");
	init();
	glutReshapeFunc(reshape);
	glutDisplayFunc(display);
	glutKeyboardFunc (keyboard);
	glutMainLoop();
	return 0;
}
