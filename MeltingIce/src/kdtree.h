/*
 * kdtree.h
 *
 *  Created on: Jul 19, 2011
 *      Author: phanindrabhagavatula
 */

#ifndef KDTREE_H_
#define KDTREE_H_

class leafNode;

class treeNode
{
public:
	treeNode():leftChild(NULL),rightChild(NULL),keyVal(0.0){};
	treeNode* getLeftChild();
	treeNode* getRightChild();
	double 	  getKeyValue();
	deque<particle>::pointer	getPtObject();
	int		getaxis();
	void	setaxis(int);
	void 	setLeftChild(treeNode*);
	void    setRightChild(treeNode*);
	void    setKeyValue(double);
	void    setRegion(double[2][3]);
	void    attachPtObject(deque<particle>::pointer);
	bool    isLeaf();
	trippleBool    testRegion(double[2][3]);
private:
	treeNode *leftChild;
	treeNode *rightChild;
	int       axis;
	double    keyVal;
	double    regionLow[3], regionHigh[3];
	deque<particle>::pointer pointObj;
};

class kdtree
{
public:
	kdtree():root(NULL){};
	treeNode* getRoot();
	treeNode* setRoot();
	void      buildTree(deque<particle>&, int);
//	void 	  insertNode(treeNode*);
	void 	  searchNodes(deque<deque<particle>::pointer >&, treeNode*, double, double[3]);
	void 	  searchNodes(deque<deque<particle>::pointer >&, treeNode*, double[2][3]);
	void	  reportSubTree(deque<deque<particle>::pointer >&, treeNode*);
	treeNode* getTreeMin();
	treeNode* getTreeMax();
	deque<particle>::pointer nearestNeighbour(deque<particle>::pointer);
	deque<particle>::pointer nearestNeighbour(deque<particle>::pointer, double);
private:
	treeNode*	buildTreeInternal(deque<deque<particle>::pointer >& , int, double[2][3]);
	void 		sortList(deque<deque<particle>::pointer >& , int);
	treeNode* nearestNeighbourInternal(treeNode*, deque<particle>::pointer, double&);
	treeNode *root;
	deque<deque<particle>::pointer > listOfParticles;
};

//class leafNode
//{
//public:
//	void attachPtObject(void *);
//	void* getPtObject();
//private:
//	void *pointObject;
//};

#endif /* KDTREE_H_ */
