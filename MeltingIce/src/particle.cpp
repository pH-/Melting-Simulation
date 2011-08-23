/*
 * particle.cpp
 *
 *  Created on: May 15, 2011
 *      Author: phanindrabhagavatula
 */

#include "includes.h"
int idParticle = 0;

particle::particle(double posx,double posy, double posz, double radius, double surfaceA, double mass, double temp, double heatValue, int status, double density)
{
	pos[0]=posx; pos[1]=posy; pos[2]=posz;
	this->radius=radius;
	this->surfaceArea=surfaceA;
	this->mass=mass;
	this->temp=temp;
	this->heatVal=heatValue;
	this->status=status;
	this->particleId=idParticle++;
	this->color[0]=0.0;
	this->color[1]=0.0;
	this->color[2]=1.0;
	this->diffusionConst=ICEDIFFCONST;
	this->density=density;
	ZERO(this->velocity);
	this->obstaclePos=pos[1];
	this->dropletCount=1;
	this->photonCount=0;
}
double particle::getHeatVal()
{
	return heatVal;
}

double particle::getMass()
{
	return mass;
}

int particle::getParticleId()
{
	return particleId;
}

double* particle::getPosition()
{
	return pos;
}

double particle::getRadius()
{
	return radius;
}

int particle::getStatus()
{
	return status;
}

int particle::getOpenFaces()
{
	return openFaces;
}

double particle::getTemp()
{
	return temp;
}
double particle::getDensity()
{
	return density;
}

int particle::getId()
{
	return particleId;
}
double* particle::getColor()
{
	return color;
}
double particle::getDiffConst()
{
	return diffusionConst;
}
double* particle::getForces()
{
	return forceVector;
}
bool particle::liesInRegion(double region[2][3])
{
	return pos[0]>=region[0][0] && pos[0]<=region[1][0]
	    && pos[1]>=region[0][1] && pos[1]<=region[1][1]
	    && pos[2]>=region[0][2] && pos[2]<=region[1][2];

}
double* particle::getVelocity()
{
	return velocity;
}
double particle::getobstaclePos()
{
	return obstaclePos;
}
int particle::getDropletCnt()
{
	return dropletCount;
}

int particle::getPhotonCount()
{
	return photonCount;
}
//// SETTERS  ////
void particle::setParticleId(int id)
{
	particleId = id;
}
void particle::setOpenFaces(int openFaces)
{
	this->openFaces = openFaces;
}
void particle::setTemp(double tempVal)
{
	this->temp=tempVal;
}
void particle::setColor()
{
	double r=0,g=0,b=0;
	if(temp>273)
	{
		g = (AIRTEMP-temp)/(AIRTEMP-273);
		r = (temp-273.0)/(AIRTEMP-273.0);
	}
	else
	{
		b = (273-temp)/(273-INITTEMP);
		g = (temp-INITTEMP)/(273-INITTEMP);
	}
	this->color[0]=r*255;this->color[1]=g*255;this->color[2]=b*255;
}
void particle::setHeat(double value)
{
	this->heatVal = value;
}
void particle::setStatus(int value)
{
	this->status = value;
}
void particle::setDiffConst(double value)
{
	this->diffusionConst=value;
}
void particle::setForces(double forceVector[3])
{
	ASSIGN(this->forceVector,forceVector);
}
void particle::setVelocity(double vel[3])
{
	ASSIGN(velocity,vel);
}
void particle::setDensity(double densityVal)
{
	this->density=densityVal;
}
void particle::setPosition(double xpos, double ypos, double zpos)
{
	pos[0]=xpos;
	pos[1]=ypos;
	pos[2]=zpos;
}
void particle::setObstaclePos(double y)
{
	obstaclePos=y;
}
void particle::clearVelocity()
{
	ZERO(velocity);
}
void particle::clearVelocity(int axis)
{
	velocity[axis]=0;
}
void particle::incrDropLetCount()
{
	dropletCount++;
	mass+=PARTICLEMASS;
}
void particle::incrPhotonCount()
{
	photonCount++;
}

void particle::clearPhotonCount()
{
	photonCount=0;
}
