/*
 * kdtree.cpp
 *
 *  Created on: Jul 19, 2011
 *      Author: phanindrabhagavatula
 */

#include"includes.h"
using namespace std;

bool sortParticlesX(deque<particle>::pointer p1, deque<particle>::pointer p2) { return (p1->getPosition()[0]<p2->getPosition()[0]);}
bool sortParticlesY(deque<particle>::pointer p1, deque<particle>::pointer p2) { return (p1->getPosition()[1]<p2->getPosition()[1]);}
bool sortParticlesZ(deque<particle>::pointer p1, deque<particle>::pointer p2) { return (p1->getPosition()[2]<p2->getPosition()[2]);}

//treeNode Class...
treeNode* treeNode::getLeftChild()
{
	return leftChild;
}

treeNode* treeNode::getRightChild()
{
	return rightChild;
}

double treeNode::getKeyValue()
{
	return keyVal;
}

deque<particle>::pointer treeNode::getPtObject()
{
	return pointObj;
}

int treeNode::getaxis()
{
	return axis;
}

void treeNode::setaxis(int axisToChk)
{
	axis = axisToChk;
}

void treeNode::setLeftChild(treeNode* ilc)
{
	leftChild = ilc;
}

void treeNode::setRightChild(treeNode* irc)
{
	rightChild = irc;
}

void treeNode::setKeyValue(double keyValue)
{
	keyVal=keyValue;
}

void treeNode::setRegion(double boundingRegion[2][3])
{
	for(int i=0; i<3; i++)
		regionLow[i]=boundingRegion[0][i];
	for(int i=0; i<3; i++)
		regionHigh[i]=boundingRegion[1][i];
}

void treeNode::attachPtObject(deque<particle>::pointer partRef)
{
	pointObj=partRef;
}

bool treeNode::isLeaf()
{
	return(!leftChild && !rightChild);
}

trippleBool treeNode::testRegion(double queryRegion[2][3])
{
//	trippleBool returnVal=falseVal;
	char returnBits=0x00;
	for(int j=0; j<3; j++)
	{
		if(regionLow[j]>queryRegion[1][j] || regionHigh[j]<queryRegion[0][j])
			return falseVal;
		if(regionLow[j]>=queryRegion[0][j])
		{
			if(j==0)
				returnBits|=0x20;
			else if(j==1)
				returnBits|=0x10;
			else if(j==2)
				returnBits|=0x08;
		}
		if(regionHigh[j]<=queryRegion[1][j])
		{
			if(j==0)
				returnBits|=0x04;
			else if(j==1)
				returnBits|=0x02;
			else if(j==2)
				returnBits|=0x01;
		}
	}
	if((0x24==(0x24 & returnBits)) && (0x12==(0x12 & returnBits)) && (0x09==(0x09 & returnBits)))
		return trueVal;
	else
		return indeterminateVal;
}

//kdtree class...

treeNode* kdtree::getRoot()
{
	return root;
}

void kdtree::buildTree(deque<particle>& particleList, int level)
{
	double max = numeric_limits<double>::max();
	double regionBox[][3]={{-max, -max, -max}, {max, max, max}};
	for(unsigned int i=0; i<particleList.size(); i++)
		listOfParticles.push_back(&particleList[i]);
	sortList(listOfParticles,0);
//	for(unsigned int i=0; i<listOfParticles.size(); i++)
//	{
//		cout<<i<<"\tparticle is at:"<<listOfParticles[i]->getPosition()[0]<<":"<<listOfParticles[i]->getPosition()[1]<<":"<<listOfParticles[i]->getPosition()[2]<<endl;
//	}
	root = buildTreeInternal(listOfParticles,0,regionBox);
}
void kdtree::searchNodes(deque<deque<particle>::pointer >& nodeList, treeNode *root, double distance, double position[3])
{
	double queryRegion[2][3];
	for(int i=0; i<2; i++)
		for(int j=0; j<3; j++)
			queryRegion[i][j] = position[j]+ pow(-1.0, i+1)*distance;
	searchNodes(nodeList,root,queryRegion);
//	deque<deque<particle>::pointer >::iterator nodeIter;
	deque<int> erasePos;
	for(unsigned int i=0; i<nodeList.size();i++)
	{
		if(euclidDistSq(position, nodeList[i]->getPosition())>distance*distance)
		{
			erasePos.push_back(i);
		}
	}
	for(unsigned int i=0; i<erasePos.size(); i++)
	{
		nodeList.erase(nodeList.begin()+(erasePos[i])-i);
	}
}

void kdtree::searchNodes(deque<deque<particle>::pointer >& nodeList, treeNode *root, double queryRegion[2][3])
{
	trippleBool testResult;
	if(root->isLeaf())
	{
		if(root->getPtObject()->liesInRegion(queryRegion))
			nodeList.push_back(root->getPtObject());
	}
	else
	{
		testResult = root->getLeftChild()->testRegion(queryRegion);
		if(testResult == trueVal)
			reportSubTree(nodeList,root->getLeftChild());
		else if(testResult==indeterminateVal)
			searchNodes(nodeList, root->getLeftChild(), queryRegion);
		testResult = root->getRightChild()->testRegion(queryRegion);
		if(testResult == trueVal)
			reportSubTree(nodeList,root->getRightChild());
		else if(testResult==indeterminateVal)
			searchNodes(nodeList,root->getRightChild(),queryRegion);
	}
}

void kdtree::reportSubTree(deque<deque<particle>::pointer >& nodeList, treeNode *subRoot)
{
	if(subRoot->isLeaf())
		nodeList.push_back(subRoot->getPtObject());
	else
	{
		reportSubTree(nodeList,subRoot->getLeftChild());
		reportSubTree(nodeList,subRoot->getRightChild());
	}
}
deque<particle>::pointer kdtree::nearestNeighbour(deque<particle>::pointer queryPoint)
{
	double max = numeric_limits<double>::max();
	double tempNearest = max;
	treeNode* nnVertex = nearestNeighbourInternal(root->getRightChild(), queryPoint, tempNearest);
	return nnVertex->getPtObject();
}

deque<particle>::pointer kdtree::nearestNeighbour(deque<particle>::pointer queryPoint, double tempNearest)
{
	treeNode* nnVertex = nearestNeighbourInternal(root->getRightChild(), queryPoint, tempNearest);
	return nnVertex->getPtObject();
}

treeNode* kdtree::buildTreeInternal(deque<deque<particle>::pointer >& particleList, int level, double regionBox[2][3])
{
	treeNode *newNode = new treeNode();
	int axisToUse = level%3;
	newNode->setaxis(axisToUse);
	deque<deque<particle>::pointer > leftSubTreeList;
	deque<deque<particle>::pointer > rightSubTreeList;
	if(particleList.size()==1)
	{
		newNode->attachPtObject(particleList[0]);
		newNode->setRegion(regionBox);
	}
	else
	{
		treeNode *leftChild, *rightChild;
		double lcRegionBox[2][3], rcRegionBox[2][3];
		memcpy(lcRegionBox, regionBox, 2*sizeof(*regionBox));
		memcpy(rcRegionBox, regionBox, 2*sizeof(*regionBox));
		unsigned int medianIndex = particleList.size()/2;
		for(unsigned int i=0; i<medianIndex; i++)
			leftSubTreeList.push_back(particleList[i]);

		for(unsigned int i=medianIndex; i<particleList.size(); i++)
			rightSubTreeList.push_back(particleList[i]);

		//double keyVal = vertexList[medianIndex-1]->getCoord((axisToSort)axisToUse) + (vertexList[medianIndex]->getCoord((axisToSort)axisToUse) - vertexList[medianIndex-1]->getCoord((axisToSort)axisToUse))/2;
		double keyVal = particleList[medianIndex-1]->getPosition()[axisToUse];
		lcRegionBox[1][axisToUse]=particleList[medianIndex-1]->getPosition()[axisToUse];
		rcRegionBox[0][axisToUse]=particleList[medianIndex]->getPosition()[axisToUse];
		newNode->attachPtObject(particleList[medianIndex-1]);
		newNode->setKeyValue(keyVal);
		newNode->setRegion(regionBox);
		sortList(leftSubTreeList,(axisToUse+1)%3);
		sortList(rightSubTreeList,(axisToUse+1)%3);
		leftChild = buildTreeInternal(leftSubTreeList,level+1,lcRegionBox);
		rightChild = buildTreeInternal(rightSubTreeList, level+1,rcRegionBox);
		newNode->setLeftChild(leftChild);
		newNode->setRightChild(rightChild);
	}
	return newNode;
}

void kdtree::sortList(deque<deque<particle>::pointer >& list, int axis)
{
	if(axis==0)
		sort(list.begin(), list.end(),sortParticlesX);
	else if(axis == 1)
		sort(list.begin(), list.end(),sortParticlesY);
	else
		sort(list.begin(), list.end(),sortParticlesZ);;
}

treeNode* kdtree::nearestNeighbourInternal(treeNode* root, deque<particle>::pointer queryPoint, double& bestDist)
{
	treeNode* candidateNearest=NULL;
	treeNode* tempNearestVertex;
	double tempNearestDist;
	if(root->isLeaf())
	{
		tempNearestVertex = root;
		if(tempNearestVertex)
		{
			tempNearestDist = euclidDistSq(queryPoint->getPosition(),tempNearestVertex->getPtObject()->getPosition());

			if(tempNearestDist<bestDist)
			{
				bestDist=tempNearestDist;
				candidateNearest=tempNearestVertex;
			}
		}
		return candidateNearest;
	}
	else if(queryPoint->getPosition()[root->getaxis()] < root->getKeyValue())
	{
		tempNearestVertex = nearestNeighbourInternal(root->getLeftChild(), queryPoint, bestDist);

		if(tempNearestVertex)
		{
			tempNearestDist = euclidDistSq(queryPoint->getPosition(),tempNearestVertex->getPtObject()->getPosition());

			if(tempNearestDist<=bestDist)
			{
				bestDist=tempNearestDist;
				candidateNearest=tempNearestVertex;
			}
		}
		if(bestDist > pow((queryPoint->getPosition()[root->getaxis()] - root->getKeyValue()),2))
			tempNearestVertex = nearestNeighbourInternal(root->getRightChild(), queryPoint, bestDist);
		else
			return candidateNearest;

		if(tempNearestVertex)
			return tempNearestVertex;
		else
			return candidateNearest;

//		else
//			return candidateNearest;
	}
	else
	{
		tempNearestVertex = nearestNeighbourInternal(root->getRightChild(), queryPoint, bestDist);

		if(tempNearestVertex)
		{
			tempNearestDist = euclidDistSq(queryPoint->getPosition(),tempNearestVertex->getPtObject()->getPosition());

			if(tempNearestDist<=bestDist)
			{
				bestDist=tempNearestDist;
				candidateNearest=tempNearestVertex;
			}
		}
		if(bestDist > pow((queryPoint->getPosition()[root->getaxis()] - root->getKeyValue()),2))
			tempNearestVertex = nearestNeighbourInternal(root->getLeftChild(), queryPoint, bestDist);
		else
			return candidateNearest;

		if(tempNearestVertex)
			return tempNearestVertex;
		else
			return candidateNearest;

//		else
//			return candidateNearest;
	}
}
