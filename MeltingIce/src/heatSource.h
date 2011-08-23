/*
 * heatSource.h
 *
 *  Created on: Aug 18, 2011
 *      Author: phanindrabhagavatula
 */

#ifndef HEATSOURCE_H_
#define HEATSOURCE_H_


class photon
{
public:
	photon(double[3],double[3]);
	void setOrigin(double[3]);
	void setDirection(double[3]);
	double* getOrigin();
	double* getDirection();
private:
	double origin[3];
	double direction[3];
};
class heatsrc
{
public:
	void setHeatSrcPos(double[3]);
	void shootPhotons();
	void setCubeDimensions(double[2][3]);
	void tracePhoton(photon&);
	void tracephotonByRay(photon&);
	void setPhotonDirNorm(double[3]);
	double* getPhotonDirNorm();
	deque<photon>& getPhotonList();
	void clearPhotonList();
private:
	bool bboxChk(double[3],double[3]);
	bool photonIntersectsTri(double[3],double&,double*,double*,double[3][3]);
	double position[3];
	double bboxDims[2][3];
	double photonDirNorm[3];
	deque<photon> photonList;

};


#endif /* HEATSOURCE_H_ */
