/*
 * object.cpp
 *
 *  Created on: May 15, 2011
 *      Author: phanindrabhagavatula
 */

#include "includes.h"
using namespace std;
extern kdtree *kdTree;
//// SETTERS ////

void object::setDims(double x, double y, double z)
{
	dims[0]=x;
	dims[1]=y;
	dims[2]=z;
}
void object::insertParticles(double dimx,double dimy,double dimz)
{
	double i,j,k;
	int count=0,openFaceCount=0;
	double x,y,z;
	//double cubeLimi = dims[0], cubeLimj = dims[1], cubeLimk = dims[2];
	for(i=0;i<dimx;i++)
		for(j=0;j<dimy;j++)
			for(k=0;k<dimz;k++)
			{
				x = (double)i*PARTICLEDIST+0.0005;
				y = (double)j*PARTICLEDIST+0.0005;
				z = (double)k*PARTICLEDIST+0.0005;
				//rotateCoords(x,y,z,45,ZAXIS);
				translateCoords(x,y,z,0,2*PARTICLEDIST,0);
				particle *newParticle = new particle(x,y,z,KERNELRANGE,6.0,PARTICLEMASS,INITTEMP,0,SOLID,ICEDENSITY);
				particlesList.push_back(*newParticle);
				count++;
				openFaceCount=0;
				if(!i && !j && !k)
				{
					bboxDims[0][0]=x;
					bboxDims[0][1]=y;
					bboxDims[0][2]=z;
				}
			}
	bboxDims[1][0]=x;
	bboxDims[1][1]=y;
	bboxDims[1][2]=z;
	//cout<<"list ready"<<endl;
	//cout<<"particle list size:"<<particlesList.size()<<endl;
	cout.flush();
}

void object::insertParticles()
{
	double icicleLen = .01;
	double icicleMaxRad = .003;
	double icicleMinRad = .001;
	double icicleLayerRadDiff = (icicleMaxRad-icicleMinRad)/((icicleLen)/0.001);
	double icicleAxis[] = {0.0,0.0};
	double icicleTopY = 0.01;
	int	   spheresAtLevel=0;
	double particlePos[3];
	double levelCenter[3];
	double levelRadius;
	int twiceNumRows;
	for(int i=0;i<10;i++)
	{
		spheresAtLevel = pow(((icicleMaxRad-icicleLayerRadDiff*i)/0.0002),2);
		levelRadius = icicleMaxRad-i*icicleLayerRadDiff;
		twiceNumRows = floor(levelRadius/0.001);
		levelCenter[0]=icicleAxis[0]+(pow(-1.0,rand()%20))*(rand()%64/1000000);
		levelCenter[1]=icicleTopY-i*0.001;
		levelCenter[2]=icicleAxis[1]+(pow(-1.0,rand()%20))*(rand()%64/1000000);
		for(int j=-twiceNumRows; j<twiceNumRows; j++)
		{
			for(int k=-twiceNumRows; k<twiceNumRows; k++)
			{
				particlePos[0]=levelCenter[0]+j*PARTICLEDIST;
				particlePos[1]=levelCenter[1];
				particlePos[2]=levelCenter[2]+k*PARTICLEDIST;
				if(chkParticle(levelCenter,levelRadius,particlePos))
				{
					particle *newParticle = new particle(particlePos[0],particlePos[1],particlePos[2],KERNELRANGE,6.0,PARTICLEMASS,INITTEMP,0,SOLID,ICEDENSITY);
					particlesList.push_back(*newParticle);
				}
			}
		}
	}
	cout<<"list ready"<<endl;
	cout<<"particle list size:"<<particlesList.size()<<endl;
	cout.flush();
}

void object::calcTemp()
{
	deque<particle>::iterator it;
	deque<deque<particle>::pointer > neighbours;
	double dTps,dT;
	double tempSum=0;
	double heatRcvd=0;
	double rangeVal;
	double openFaces;
	int loopcount =0;
	for(it=particlesList.begin(); it!=particlesList.end(); it++)
	{

		if(it->getId()==-1)
			continue;
		tempSum=0;
		heatRcvd=0;
		if(USEHEATSOURCE)
			heatRcvd+=it->getPhotonCount()*(EMMISSIVITYOFHEATSRC*STEPHANBOLTZCONST*(pow(TEMPHEATSOURCE,4))*HEATSRCAREA*TIMESTEP)/PHOTONSPERTIMESTEP;
		loopcount++;
		double currTemp = it->getTemp();
		if(currTemp!=currTemp)
		{
			cout<<"found Nan"<<endl;
		}
		rangeVal = it->getRadius();
		neighbours.clear();
		kdTree->searchNodes(neighbours,kdTree->getRoot(),it->getRadius(),it->getPosition());
		//cout<<"neighbour size tempcalc:"<<neighbours.size()<<endl;
		cout.flush();
		if(neighbours.size()<8)
			openFaces = 6.0 - (neighbours.size()-1);
		else
			openFaces = 0;
		it->setOpenFaces(openFaces);
		if(openFaces>0)
		{
			//if(rand()%10>6)
				//continue;
			heatRcvd += ICETHERMCONDUCTIVITY*(AIRTEMP-currTemp)*(pow(CUBESIZE,2)*openFaces);
			//cout<<"Count:"<<loopcount<<"currTemp"<<currTemp<<endl;
			if(currTemp >= 273)
			{
				if(it->getStatus()==SOLID)
				{
					if(it->getHeatVal()>=ICELATENTHEAT*it->getMass())

					{
						it->setStatus(LIQUID);
						it->setHeat(0.0);
						it->setDiffConst(WATERDIFFCONST);
					}
					else
					{
						it->setHeat(it->getHeatVal()+heatRcvd);
					}
				}
				else
				{
					dTps = heatRcvd/(ICESPHEAT*it->getMass());
					dT = dTps*TIMESTEP;
					if(dT!=dT)
						cout<<"problem"<<endl;
					it->setTemp(currTemp+dT);
					it->setColor();
				}
			}
			else
			{
				dTps = heatRcvd/(ICESPHEAT*it->getMass());
				dT = dTps*TIMESTEP;
				it->setTemp(currTemp+dT);
				it->setColor();
			}
			//cout<<"dT is:"<<dT;
		}
		else
		{
			deque<deque<particle>::pointer>::iterator iter;
			for(iter=neighbours.begin(); iter!=neighbours.end(); ++iter)
			{
				if((*iter)->getId()==it->getId())
					continue;
				double partDistance,partDensity;
				partDistance = euclidDistSq(it->getPosition(),(*iter)->getPosition());
				partDensity = (*iter)->getDensity();
				tempSum += (*iter)->getDiffConst()*(*iter)->getMass()*((*iter)->getTemp()-it->getTemp())*weightValue(sqrt(partDistance),sqrt(rangeVal))*(1/partDensity);
			}
			dTps = tempSum;
			dT = dTps*TIMESTEP;
			if(it->getTemp()>=273)
			{
				if(it->getStatus()==SOLID)
				{
					if(it->getHeatVal()>=ICELATENTHEAT)
					{
						it->setStatus(LIQUID);
						it->setHeat(0);
						it->setDiffConst(WATERDIFFCONST);
					}
					else
					{
						heatRcvd=dT*ICESPHEAT*it->getMass();
						it->setHeat(it->getHeatVal()+heatRcvd);
					}
				}
				else
				{
					currTemp+=dT;
					if(currTemp!=currTemp)
						cout<<"problem"<<endl;
					it->setTemp(currTemp);
					it->setColor();
				}
			}
			else
			{
				currTemp+=dT;
				if(currTemp!=currTemp)
					cout<<"problem"<<endl;
				it->setTemp(currTemp);
				it->setColor();
			}
			//cout<<"dT is:"<<dT<<endl;
		}
	}
}

void object::calcForces()
{
	deque<particle>::iterator iter;
	double surfaceForces[3];
	double dragForce[3];
	//double pressureForce[3];
	deque<int> deleteList;
	//double viscosityForce[3];
	double forceRange;
	double totalForce[3];
	double direction[3];
	double normalDir[3];
	double flowDir[3];
	double flowDirUnit[3];
	double gravity[3];
//	double gravityDir[3];
	double tempDotPr=0.0;
	double tempVector[3];
	//bool obstacle=false;
	unsigned int closedFaces=0;
	for(iter=particlesList.begin(); iter!=particlesList.end(); ++iter)
	{
		if(iter->getId()==-1)
			continue;
		ZERO(surfaceForces);
		ZERO(totalForce);
		ZERO(gravity);
		ZERO(dragForce);
//		ZERO(gravityDir);
		ZERO(normalDir);
		ZERO(flowDir);
		flowDir[1]=-1;
		gravity[1]=-8.8;
		closedFaces=0;
//		gravityDir[1]=-1;
		if(iter->getStatus()==LIQUID)
		{
			deque<deque<particle>::pointer> neighbours;
			neighbours.clear();
//			if((iter->getId()%100)/10==9 && iter->getId()/100<8 &&iter->getId()/100>3)
//				cout<<"stop to chk"<<endl;
			forceRange = iter->getRadius();
			kdTree->searchNodes(neighbours,kdTree->getRoot(),forceRange*5,iter->getPosition());
			//				cout<<"neighbours of the liquid prticle:"<<neighbours.size();
			//				cout.flush();
			for(deque<deque<particle>::pointer>::iterator niter = neighbours.begin(); niter!=neighbours.end(); ++niter)
			{
				if((*niter)->getId()==iter->getId() || (*niter)->getId()==-1)
					continue;
				if(euclidDistSq((*niter)->getPosition(),iter->getPosition())<=0.3*forceRange*forceRange)
				{
					closedFaces++;
					if((*niter)->getPosition()[1]<iter->getPosition()[1])
					{
						if((*niter)->getStatus()==SOLID)
						{
							iter->setObstaclePos(iter->getPosition()[1]);
							iter->clearVelocity(YAXIS);
						}
						else
							iter->setObstaclePos(GLOBALYMIN);
					}
				}
				else
					iter->setObstaclePos(GLOBALYMIN);
				MINUS(direction, (*niter)->getPosition(), iter->getPosition());
				double dirLenSq= direction[0]*direction[0]+direction[1]*direction[1]+direction[2]*direction[2];
				//MULT(direction,direction,(1/dirLenSq));
				if((*niter)->getStatus()==LIQUID)
				{
					if(dirLenSq<MINSEPARATION*MINSEPARATION)
					{
						if((*niter)->getPosition()[1]<iter->getPosition()[1])
						{
							(*niter)->incrDropLetCount();
							//iter->setParticleId(-1);

						}
						else
						{
							iter->incrDropLetCount();
							//(*niter)->setParticleId(-1);
						}

					}
					else
					{
						//MULT(direction,direction,(1/dirLenSq));
						MULT(direction,direction,INTERFACIALWATER);
						PLUS(surfaceForces,surfaceForces,direction);
					}
				}
				else
				{
					if(dirLenSq<forceRange*forceRange)
						PLUS(normalDir,normalDir,direction);
				}

			}
			if(closedFaces>5)
				iter->setOpenFaces(0);
			MULT(gravity,gravity,iter->getMass());
			if(norm3d(normalDir))
			{
				MULT(normalDir,normalDir,0.71);
				tempDotPr = DOT(normalDir,gravity);
				MULT(tempVector,normalDir,tempDotPr);
				MINUS(flowDir,gravity,tempVector);
				MULT(dragForce,flowDir,-1*INTERFACIALICE);
			}
			PLUS(surfaceForces,surfaceForces,dragForce);
			PLUS(totalForce,surfaceForces,gravity);
			ASSIGN(flowDirUnit,flowDir);
			norm3d(flowDirUnit);
			tempDotPr = DOT(totalForce,flowDirUnit);
			MULT(totalForce,flowDirUnit,tempDotPr);
			//getPressureForce(pressureForce);
			//PLUS(totalForce,totalForce,pressureForce)
			//MULT(totalForce,totalForce,(1/iter->getMass()));
			iter->setForces(totalForce);
		}

	}
//	for(deleteIndex=0; deleteIndex<particlesList.size(); deleteIndex++)
//	{
//		if(particlesList[deleteIndex].getId()==-1)
//		{
//			deleteList.push_back(deleteIndex);
//		}
//	}
//	for(unsigned int i=0; i<deleteList.size(); i++)
//	{
//		particlesList.erase(particlesList.begin()+(deleteList[i]-i));
//	}
}

void object::calcPosition()
{
	deque<particle>::iterator piter;
	double newPos[3];
	double oldPos[3];
	double acceleration[3];
	double *forces;
	double *velocity;
	forces = (double*)malloc(sizeof(double)*3);
	velocity = (double*)malloc(sizeof(double)*3);
	for(piter=particlesList.begin(); piter!=particlesList.end();piter++)
	{
		piter->clearPhotonCount();
		if(piter->getStatus()==LIQUID && piter->getOpenFaces())
		{
			oldPos[0] = piter->getPosition()[0];
			oldPos[1] = piter->getPosition()[1];
			oldPos[2] = piter->getPosition()[2];
			forces = piter->getForces();
			MULT(acceleration,forces,1/piter->getMass());
			velocity = piter->getVelocity();
			velocity[0] += acceleration[0]*TIMESTEP;
			velocity[1] += acceleration[1]*TIMESTEP;
			velocity[2] += acceleration[2]*TIMESTEP;
			newPos[0] = velocity[0]*TIMESTEP;
			newPos[1] = velocity[1]*TIMESTEP;
			newPos[2] = velocity[2]*TIMESTEP;
			//PLUS(newPos,newPos,oldPos);
			//		cout<<"\tPosition:"<<newPos[0]<<":"<<newPos[1]<<":"<<newPos[2]<<endl;
			if(newPos[1]+oldPos[1]>GLOBALYMIN)
			{
				PLUS(newPos,newPos,oldPos);
			}
			else
			{
				if(oldPos[1]>GLOBALYMIN)
				{
					ZERO(velocity);
				}
				newPos[1]=0;
				PLUS(newPos,newPos,oldPos);
				newPos[1]=GLOBALYMIN;
			}
			piter->setPosition(newPos[0],newPos[1],newPos[2]);
			piter->setVelocity(velocity);
		}
	}
}


//// GETTERS ////

deque<particle>& object::getParticleList()
{
	return particlesList;
}

particle* object::findParticleatPos(double x, double y, double z)
{
	deque<particle>::iterator parIt;
	for(parIt=particlesList.begin();parIt != particlesList.end(); ++parIt)
	{
		if(parIt->getPosition()[0]==x && parIt->getPosition()[1]==y && parIt->getPosition()[2]==z)
		{
			return &(*parIt);
		}
	}
	return NULL;
}

//void object::particlesAtDist(deque<particle>*neighbours, particle *centerParticle, double distSquare)
//{
//	deque<particle>::iterator parIt;double distance;
//	for(parIt=particlesList.begin();parIt != particlesList.end(); ++parIt)
//		{
//			distance = euclidDistSq(centerParticle->getPosition(),parIt->getPosition());
//			if(distance && distance<=distSquare)
//			{
//				neighbours->push_back(*parIt);
//			}
//		}
//}

void object::waterParticlesAtDist(deque<deque<particle>::pointer>& neighbours, particle& centerParticle, double distSquare)
{
	deque<deque<particle>::pointer>::iterator parIt;
//	double distance=0;
	kdTree->searchNodes(neighbours,kdTree->getRoot(),distSquare,centerParticle.getPosition());
	for(parIt=neighbours.begin(); parIt != neighbours.end(); ++parIt)
		{
			if((*parIt)->getStatus()!=LIQUID)
			{
				neighbours.erase(parIt++);
			}
		}
}

double object::calcWaterDensity(particle& centerParticle)
{
	deque<deque<particle>::pointer> neighbours;
	deque<deque<particle>::pointer>::iterator iter;
	double totalDensity=0;
	double distance;
	double rangeVal = centerParticle.getRadius();
	waterParticlesAtDist(neighbours,centerParticle,rangeVal);
	//cout<<"neighbors for density:"<<neighbours.size();
	for(iter=neighbours.begin(); iter!=neighbours.end(); ++iter)
	{
		distance=euclidDistSq(centerParticle.getPosition(),(*iter)->getPosition());
		totalDensity+= (*iter)->getMass()*densityKernel(sqrt(distance),sqrt(rangeVal));
	}
	return totalDensity;
}

void object::calcDensities()
{
	deque<particle>::iterator piter;
	for(piter=particlesList.begin();piter!=particlesList.end();piter++)
	{
		deque<deque<particle>::pointer> neighbours;
		deque<particle>::iterator iter;
		double totalDensity=0;
		int effectiveNeighs=0;
		double distanceSq;
		double rangeVal = piter->getRadius();
		//particlesAtDist(&neighbours,(&(*piter)),rangeVal);
		//cout<<"neighbors for density:"<<neighbours.size();
		kdTree->searchNodes(neighbours,kdTree->getRoot(),rangeVal,piter->getPosition());
//		cout<<"position:"<<piter->getPosition()[0]<<":"<<piter->getPosition()[1]<<":"<<piter->getPosition()[2]<<endl;
//		cout<<"neighs found for density:"<<neighbours.size()<<endl;
//		cout.flush();

		for(unsigned int i=0; i<neighbours.size(); ++i)
		{
			distanceSq=euclidDistSq(piter->getPosition(),neighbours[i]->getPosition());
			if(distanceSq<rangeVal*rangeVal)
			{
				effectiveNeighs+=1;
				totalDensity+= neighbours[i]->getMass()*densityKernel(distanceSq,rangeVal);
			}
		}
		//cout<<"effective neighs for density:"<<effectiveNeighs<<endl;
		piter->setDensity(totalDensity);
	}
}
void object::getPressureForce(double pressureForce[3])
{
	deque<particle>::iterator piter;
	deque<deque<particle>::pointer> neighbours;
	deque<deque<particle>::pointer>::iterator niter;
	double particlePress,nparticlePress;
	double scalarForce=0;
	double dirVector[3];
	double rangeVal = 3*(CUBESIZE)*(CUBESIZE);
	double ndist;
	for(piter=particlesList.begin(); piter!=particlesList.end(); ++piter)
	{
		if(piter->getId()==-1)
			continue;
		neighbours.clear();
		particlePress = AIRGASCONST*calcWaterDensity(*piter);
		waterParticlesAtDist(neighbours,*piter,rangeVal);
		for(niter=neighbours.begin();niter!=neighbours.end();++niter)
		{
			double nDensity = calcWaterDensity(*(*niter));
			MINUS(dirVector,piter->getPosition(),(*niter)->getPosition());
			norm3d(dirVector);
			ndist = euclidDistSq(piter->getPosition(),(*niter)->getPosition());
			nparticlePress = AIRGASCONST*nDensity;
			scalarForce += ((*niter)->getMass()*(particlePress+nparticlePress)*pressureKernel(sqrt(ndist),sqrt(rangeVal)))/2*nDensity;
		}
		MULT(pressureForce,dirVector,-1*scalarForce);
	}
}

double (&object::getBboxDims())[2][3]
{
	return bboxDims;
}

double object::weightValue(double r, double h)
{
	double val;
	val = 45*(h-r)/(PI*pow(h,6));
	return val;
}

