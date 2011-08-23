/*
 * particle.h
 *
 *  Created on: May 14, 2011
 *      Author: phanindrabhagavatula
 */

#ifndef PARTICLE_H_
#define PARTICLE_H_
class particle {

public:
	particle(double posx,double posy, double posz, double radius, double surfaceA, double mass, double temp, double heatValue, int status, double density);
	void setParticleId(int);
	void setPosition(double xpos, double ypos, double zpos);
	void setTemp(double);
	void setHeat(double);
	void setStatus(int );
	void setOpenFaces(int);
	void setColor();
	void setDiffConst(double);
	void setForces(double[3]);
	void setVelocity(double[3]);
	void setDensity(double);
	void setObstaclePos(double);
	void clearVelocity();
	void clearVelocity(int);
	void incrDropLetCount();
	void incrPhotonCount();
	void clearPhotonCount();

	double getHeatVal();
	double getMass();
	double* getPosition();
	double getRadius();
	double getTemp();
	double getDensity();
	double getDiffConst();
	int    getParticleId();
	int    getOpenFaces();
	int    getStatus();
	int    getId();
	double* getColor();
	double* getForces();
	double* getVelocity();
	double getobstaclePos();
	int    getDropletCnt();
	int    getPhotonCount();
	bool   liesInRegion(double[2][3]);
private:
	double pos[3];
	double radius;
	double surfaceArea;
	double mass;
	double temp;
	double heatVal;
	double color[3];
	double density;
	double diffusionConst;
	double forceVector[3];
	double velocity[3];
	int    openFaces;
	int   status;
	double obstaclePos;
	int   particleId;
	int	  dropletCount;
	int   photonCount;

};

#endif /* PARTICLE_H_ */
