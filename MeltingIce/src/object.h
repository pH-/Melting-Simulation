/*
 * object.h
 *
 *  Created on: May 14, 2011
 *      Author: phanindrabhagavatula
 */

#ifndef OBJECT_H_
#define OBJECT_H_

using namespace std;

class object {
public:
	void calcTemp();
	void calcForces();
	void insertParticles(double,double,double);
	void insertParticles();
	void calcPosition();
	particle* findParticleatPos(double, double, double);
	void setDims(double, double, double);
	deque<particle>& getParticleList();
	//void particlesAtDist(deque<particle>*,particle*,double);
	void waterParticlesAtDist(deque<deque<particle>::pointer>&,particle&,double);
	void getPressureForce(double[3]);
	void getViscocityForce(double[3]);
	double (&getBboxDims())[2][3];
	double weightValue(double,double);
	double calcWaterDensity(particle&);
	void calcDensities();

private:
	double dims[3];
	double bboxDims[2][3];
	deque<particle> particlesList;

};


#endif /* OBJECT_H_ */
