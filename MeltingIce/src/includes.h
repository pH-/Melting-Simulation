/*
 * includes.h
 *
 *  Created on: May 15, 2011
 *      Author: phanindrabhagavatula
 */

#ifndef INCLUDES_H_
#define INCLUDES_H_



#include <iostream>
#include <deque>
#include <algorithm>
#include <list>
#include <cmath>
#include <limits>
#include <fstream>
#include <sstream>

#ifdef __APPLE__
	#include <GL/glut.h>
	#include <GL/gl.h>
	#include <GL/glu.h>
#endif
#ifdef __linux__
	#include<GL/glut.h>
	#include<GL/gl.h>
	#include<GL/glu.h>
#endif

#ifdef _WIN32 || _WIN64
	#include<GL/glut.h>
	#include<GL/gl.h>
	#include<GL/glu.h>
#endif

enum trippleBool{falseVal, indeterminateVal, trueVal};
using namespace std;

#include "particle.h"
#include "object.h"
#include "kdtree.h"
#include "heatSource.h"


#define CROSS(O,A,B)		 {(O)[0] = (A)[1]*(B)[2]-(A)[2]*(B)[1]; \
							  (O)[1] = (A)[2]*(B)[0]-(A)[0]*(B)[2]; \
							  (O)[2] = (A)[0]*(B)[1]-(A)[1]*(B)[0];}

#define DOT(A,B) ((A)[0]*(B)[0] + (A)[1]*(B)[1]+(A)[2]*(B)[2])

#define PLUS(O,A,B)			 {(O)[0] = (A)[0]+(B)[0]; \
							  (O)[1] = (A)[1]+(B)[1]; \
							  (O)[2] = (A)[2]+(B)[2];}

#define MINUS(O,A,B)		 {(O)[0] = (A)[0]-(B)[0]; \
							  (O)[1] = (A)[1]-(B)[1]; \
							  (O)[2] = (A)[2]-(B)[2];}

#define ASSIGN(O,A)			 {(O)[0] = (A)[0]; \
							  (O)[1] = (A)[1]; \
							  (O)[2] = (A)[2];}

#define MULT(O,A,c)			 {(O)[0] = (A)[0]*c; \
							  (O)[1] = (A)[1]*c; \
							  (O)[2] = (A)[2]*c;}

#define ZERO(O)				 {(O)[0]=0.0;\
							  (O)[1]=0.0;\
							  (O)[2]=0.0;}


#define XAXIS 0
#define YAXIS 1
#define ZAXIS 2
#define X 0
#define Y 1
#define Z 2
#define GLOBALYMIN 0
#define CUBESIZE .001
#define PARTICLEDIST .001        //meters
#define KERNELRANGE 0.0011      //PARTICLEDIST+0.0001
#define MINSEPARATION 0.001
#define PARTICLEMASS 0.0000009162  // in kgs
#define HEATSRCPOSX .011
#define HEATSRCPOSY .013
#define HEATSRCPOSZ .011
#define SOLID    1
#define LIQUID   0
#define AIRTEMP 272.8
#define ICETHERMCONDUCTIVITY 2.22    // W/mK at 273 Kelvin
#define ICESPHEAT 2050                  // units: J/KgK
#define TIMESTEP .003				// units: s
#define PI 3.14159268
#define ICEDIFFCONST 0.117132803     // units: m^2/s
#define WATERDIFFCONST 0.000000137963844   // units: m^2/s
#define ICEDENSITY 916.2
#define INITTEMP 272.6
#define ICELATENTHEAT 334000            // units: J/KG
#define INTERFACIALWATER .0000005
#define INTERFACIALICE   1
#define AIRGASCONST 286.9
#define DISPLAYSCALE 0.001
#define USEHEATSOURCE 1
#define PHOTONSPERTIMESTEP 100
#define STEPHANBOLTZCONST .00000005670373      //units: W/m^2*K^4
#define TEMPHEATSOURCE 1000.0                    //units: K
#define HEATSRCAREA .05                    //units: m^2
#define EMMISSIVITYOFHEATSRC 1
inline double euclidDistSq(double point1[3],double point2[3])
{
	//cout<<"points"<<point1[0]<<":"<<point1[1]<<":"<<point1[2]<<"&&&&"<<point2[0]<<":"<<point2[1]<<":"<<point2[2]<<endl;
	return (pow(point1[0]-point2[0],2))+(pow(point1[1]-point2[1],2))+(pow(point1[2]-point2[2],2));
}
inline double densityKernel(double distanceSq, double effRad)
{
	if(distanceSq > effRad*effRad)
		return 0;
	return 315*(pow((effRad*effRad - distanceSq),3)/(64*PI*pow(effRad,9)));
}
inline double pressureKernel(double r, double h)
{
	return 15*(pow((h-r),3))/(PI*pow(h,6));
}
inline double norm3d(double *vect)//normalizes a 3d vector
{
	double length;
	length = sqrt(vect[0]*vect[0] + vect[1]*vect[1] + vect[2]*vect[2]);

	if(length)
    {//normalize the vector if its non zero
  		vect[0]=vect[0] / length;
  		vect[1]=vect[1] / length;
		vect[2]=vect[2] / length;
    }
	return length;
}
inline void matrixMult(double rotationMatrix[3][3], double coordMatrix[3])
{
	double tempVar = 0.0;
	double newCoords[3];
	for(int i=0; i<3; i++)
	{
		for(int j=0; j<3; j++)
		{
			tempVar+= rotationMatrix[i][j]*coordMatrix[j];
		}
		newCoords[i]=tempVar;
		tempVar=0.0;
	}
	ASSIGN(coordMatrix,newCoords);
}
inline void rotateCoords(double& x, double& y, double& z, double angle, int axis)
{
	angle = angle*PI/180;
	double coordMatrix[3] = {x,y,z};
	double matrixZ[3][3] = {{cos(angle),-sin(angle), 0}, {sin(angle), cos(angle), 0}, {0,0,1}};
	if(axis == ZAXIS)
	{
		matrixMult(matrixZ,coordMatrix);
		x=coordMatrix[0];
		y=coordMatrix[1];
		z=coordMatrix[2];
	}
}
inline void translateCoords(double& x,double& y,double& z,double deltaX, double deltaY, double deltaZ)
{
	x=x+deltaX;
	y=y+deltaY;
	z=z+deltaZ;
}
inline bool chkParticle(double center[3], double radius, double particleCenter[3])
{
	return(euclidDistSq(particleCenter,center)<radius*radius);
}

inline double randnum(void)           // random number in range [0,1]
{
	return ((double)rand())/((double)RAND_MAX+double(1.0));
}
#endif /* INCLUDES_H_ */
