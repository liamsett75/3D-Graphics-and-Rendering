#pragma once
#include <vector>
#include "Gz.h"

class KDNode
{
public:
	bool isLeaf; // Instead of introducing extra inheritance

	/* For leaves */
	std::vector<GzTriangle> contents;

	/* For nodes */
	int dim; // should be X, Y, or Z
	float separator; // value of separator plane
	KDNode *leftChild;
	KDNode *rightChild;
};
