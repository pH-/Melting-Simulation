/*
 * heatSource.cpp
 *
 *  Created on: Aug 18, 2011
 *      Author: phanindrabhagavatula
 */

#include "includes.h"
extern kdtree *kdTree;

photon::photon(double orig[3], double dir[3])
{
	ASSIGN(origin,orig);
	ASSIGN(direction,dir);
}
double* photon::getOrigin()
{
	return origin;
}

double* photon::getDirection()
{
	return direction;
}

void photon::setOrigin(double orig[3])
{
	ASSIGN(origin,orig);
}

void photon::setDirection(double dir[3])
{
	ASSIGN(direction,dir);
}


///class heatsource

void heatsrc::setHeatSrcPos(double pos[3])
{
	ASSIGN(position,pos);
}

void heatsrc::setCubeDimensions(double cubeDims[2][3])
{
	ASSIGN(bboxDims[0],cubeDims[0]);
	ASSIGN(bboxDims[1],cubeDims[1]);
}
void heatsrc::setPhotonDirNorm(double normDir[3])
{
	ASSIGN(photonDirNorm,normDir);
}

double* heatsrc::getPhotonDirNorm()
{
	return photonDirNorm;
}
void heatsrc::shootPhotons()
{
	norm3d(photonDirNorm);
	double vec[3];
	//do{
		do{
			vec[0] = (randnum()-.5)*2.0;
			vec[1] = (randnum()-.5)*2.0;
			vec[2] = (randnum()-.5)*2.0;
		}while((DOT(vec,vec)>1.0));

		norm3d(vec);
		if(DOT(vec,photonDirNorm)<0.0)
		{
			vec[0]=-vec[0];
			vec[1]=-vec[1];
			vec[2]=-vec[2];
		}
	//}while(!bboxChk(position,vec));
	photon *newPhoton = new photon(position,vec);
	photonList.push_back(*newPhoton);
}

void heatsrc::tracephotonByRay(photon& photonToTrace)
{
	double searchPoint[3];
	double tempParticleDist;
	double bestParticleDist=numeric_limits<double>::max();
	int winnerIndex=-1;

	deque<deque<particle>::pointer > neighbours;
	for(int i=0;i<1000;i++)
	{
		MULT(searchPoint,photonToTrace.getDirection(),i*PARTICLEDIST);
		PLUS(searchPoint,searchPoint,photonToTrace.getOrigin());
		kdTree->searchNodes(neighbours,kdTree->getRoot(),KERNELRANGE,searchPoint);
		for(unsigned int iter=0; iter<neighbours.size(); ++iter)
		{
			tempParticleDist=euclidDistSq(searchPoint,neighbours[iter]->getPosition());
			if(tempParticleDist<bestParticleDist)
			{
				bestParticleDist=tempParticleDist;
				winnerIndex=iter;
			}
		}
		if(winnerIndex>-1)
		{
			neighbours[winnerIndex]->incrPhotonCount();
			break;
		}

	}
}
void heatsrc::tracePhoton(photon& photonToTrace)
{
	double triangles[12][3][3];
	double vertices[8][3];
	deque<deque<particle>::pointer > neighbours;
	for(int i=0;i<8;i++)
	{
		vertices[i][0]=bboxDims[(i&4)>>2][0];
		vertices[i][1]=bboxDims[(i&2)>>1][1];
		vertices[i][2]=bboxDims[(i&1)][2];
	}
	ASSIGN(triangles[0][0],vertices[0]);
	ASSIGN(triangles[0][1],vertices[1]);
	ASSIGN(triangles[0][2],vertices[3]);

	ASSIGN(triangles[1][0],vertices[0]);
	ASSIGN(triangles[1][1],vertices[3]);
	ASSIGN(triangles[1][2],vertices[2]);

	ASSIGN(triangles[2][0],vertices[0]);
	ASSIGN(triangles[2][1],vertices[4]);
	ASSIGN(triangles[2][2],vertices[5]);

	ASSIGN(triangles[3][0],vertices[0]);
	ASSIGN(triangles[3][1],vertices[5]);
	ASSIGN(triangles[3][2],vertices[1]);

	ASSIGN(triangles[4][0],vertices[0]);
	ASSIGN(triangles[4][1],vertices[2]);
	ASSIGN(triangles[4][2],vertices[6]);

	ASSIGN(triangles[5][0],vertices[0]);
	ASSIGN(triangles[5][1],vertices[6]);
	ASSIGN(triangles[5][2],vertices[4]);

	ASSIGN(triangles[6][0],vertices[7]);
	ASSIGN(triangles[6][1],vertices[6]);
	ASSIGN(triangles[6][2],vertices[2]);

	ASSIGN(triangles[7][0],vertices[7]);
	ASSIGN(triangles[7][1],vertices[2]);
	ASSIGN(triangles[7][2],vertices[3]);

	ASSIGN(triangles[8][0],vertices[7]);
	ASSIGN(triangles[8][1],vertices[5]);
	ASSIGN(triangles[8][2],vertices[4]);

	ASSIGN(triangles[9][0],vertices[7]);
	ASSIGN(triangles[9][1],vertices[4]);
	ASSIGN(triangles[9][2],vertices[6]);

	ASSIGN(triangles[10][0],vertices[7]);
	ASSIGN(triangles[10][1],vertices[3]);
	ASSIGN(triangles[10][2],vertices[1]);

	ASSIGN(triangles[11][0],vertices[7]);
	ASSIGN(triangles[11][1],vertices[1]);
	ASSIGN(triangles[11][2],vertices[5]);

	double currdist,maxdist;        // distance on the path of photon where it hits an object
	double hitpos[3],hitNorm[3];
	int winnerTriangle=-1;
	int winnerIndex=-1;
	double pos[3],norm[3];    // to retrieve hit position and hit normal ofa  triangle
	ZERO(pos);
	ZERO(norm);
	maxdist=numeric_limits<double>::max();
	double bestParticleDist=maxdist, tempParticleDist;
	for (int j=0; j<12; j++)
	{
		if(photonIntersectsTri(photonToTrace.getDirection(),currdist,pos,norm,triangles[j]))
		{
			if (currdist<maxdist)
			{
				maxdist=currdist;
				winnerTriangle=j;
				ASSIGN(hitpos,pos);
				//interpolate normals to get the normal at point of intersection
				ASSIGN(hitNorm,norm);
			}
		}
	}

	if (maxdist!=numeric_limits<double>::max())
	{
		//deque<deque<particle>::pointer>::iterator iter;
		kdTree->searchNodes(neighbours,kdTree->getRoot(),KERNELRANGE,hitpos);
		for(unsigned int iter=0; iter<neighbours.size(); ++iter)
		{
			tempParticleDist=euclidDistSq(hitpos,neighbours[iter]->getPosition());
			if(tempParticleDist<bestParticleDist)
			{
				bestParticleDist=tempParticleDist;
				winnerIndex=iter;
			}
		}
		if(winnerIndex>-1)
			neighbours[winnerIndex]->incrPhotonCount();

	}

}



bool heatsrc::bboxChk(double origpos[3], double sdirection[3])
{
	double r,x,y,z;
	r = (bboxDims[0][X]-origpos[X])/sdirection[X];
	y = origpos[Y]+(sdirection[Y]*r);
	z = origpos[Z]+(sdirection[Z]*r);

	if((r>=0.0) && (y>=bboxDims[0][Y]) && (y<=bboxDims[1][Y]) && (z>=bboxDims[0][Z]) && (z<=bboxDims[1][Z]))
		return 1;

	r = (bboxDims[1][X]-origpos[X])/sdirection[X];
	y = origpos[Y]+(sdirection[Y]*r);
	z = origpos[Z]+(sdirection[Z]*r);

	if((r>=0.0) && (y>=bboxDims[0][Y]) && (y<=bboxDims[1][Y]) && (z>=bboxDims[0][Z]) && (z<=bboxDims[1][Z]))
		return 1;

	// go along the y direction
	r = (bboxDims[0][Y]-origpos[Y])/sdirection[Y];
	x = origpos[X]+(sdirection[X]*r);
	z = origpos[Z]+(sdirection[Z]*r);

	if((r>=0.0) && (x>=bboxDims[0][X]) && (x<=bboxDims[1][X]) && (z>=bboxDims[0][Z]) && (z<=bboxDims[1][Z]))
		return 1;

	r = (bboxDims[1][Y]-origpos[Y])/sdirection[Y];
	x = origpos[X]+(sdirection[X]*r);
	z = origpos[Z]+(sdirection[Z]*r);

	if((r>=0.0) && (x>=bboxDims[0][X]) && (x<=bboxDims[1][X]) && (z>=bboxDims[0][Z]) && (z<=bboxDims[1][Z]))
		return 1;

	// go along the z direction
	r = (bboxDims[0][Z]-origpos[Z])/sdirection[Z];
	y = origpos[Y]+(sdirection[Y]*r);
	x = origpos[X]+(sdirection[X]*r);

	if((r>=0.0) && (y>=bboxDims[0][Y]) && (y<=bboxDims[1][Y]) && (x>=bboxDims[0][X]) && (x<=bboxDims[1][X]))
		return 1;

	r = (bboxDims[1][Z]-origpos[Z])/sdirection[Z];
	y = origpos[Y]+(sdirection[Y]*r);
	x = origpos[X]+(sdirection[X]*r);

	if((r>=0.0) && (y>=bboxDims[0][Y]) && (y<=bboxDims[1][Y]) && (x>=bboxDims[0][X]) && (x<=bboxDims[1][X]))
		return 1;

	return 0;

}

//photonIntersectsTri(photonToTrace.getDirection(),currdist,pos,norm,triangles[j]))
bool heatsrc::photonIntersectsTri(double *sdirection, double& dist,double *intersectionPos,double *intersectionNorm, double triangle[3][3])
{
	double triNorm[3],triSideA[3],triSideB[3],tempA[3],temp1,temp2,vertxPtVector1[3],vertxPtVector2[3];

	MINUS(triSideA,triangle[1],triangle[0]);
	MINUS(triSideB,triangle[2],triangle[0]);
	CROSS(triNorm,triSideA,triSideB);

	MINUS(tempA,triangle[0],position);        // vector from origin of ray to vertex of triangle
	temp1=DOT(tempA,triNorm);
	temp2=DOT(sdirection,triNorm);
	if (!temp2)                               // if temp2=0 , ray is parallel
	{
		if (!temp1)
			return 0;   // The ray is in the plane.
		else
			return 0;   // The ray is parallel to the plane
	}
	else
	{
		dist = temp1/temp2;
		if(dist<0)
			return 0;
		intersectionPos[X]=position[X]+dist*sdirection[X];
		intersectionPos[Y]=position[Y]+dist*sdirection[Y];
		intersectionPos[Z]=position[Z]+dist*sdirection[Z];

		norm3d(triNorm);
		ASSIGN(intersectionNorm,triNorm);

		MINUS(vertxPtVector1,triangle[0],intersectionPos);
		MINUS(vertxPtVector2,triangle[1],intersectionPos);
		CROSS(tempA,vertxPtVector1,vertxPtVector2);
		temp1=DOT(tempA,triNorm);
		if(temp1<0)
			return 0;         // ray doesnt pass through the triangle
		MINUS(vertxPtVector1,triangle[1],intersectionPos);
		MINUS(vertxPtVector2,triangle[2],intersectionPos);
		CROSS(tempA,vertxPtVector1,vertxPtVector2);
		temp1=DOT(tempA,triNorm);
		if(temp1<0)
			return 0;

		MINUS(vertxPtVector1,triangle[2],intersectionPos);
		MINUS(vertxPtVector2,triangle[0],intersectionPos);
		CROSS(tempA,vertxPtVector1,vertxPtVector2);
		temp1=DOT(tempA,triNorm);
		if(temp1<0)
			return 0;

		return 1;          // point is inside the triangle.
	}

}

deque<photon>& heatsrc::getPhotonList()
{
	return photonList;
}

void heatsrc::clearPhotonList()
{
	photonList.clear();
}
