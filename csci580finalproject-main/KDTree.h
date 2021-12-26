#pragma once
#include "AContainer.h"
#include "KDNode.h"
#include <map>
#include <vector>
#include <assert.h>
#include <algorithm>

class KDTree : public AContainer
{
public:
	enum class side { L, R };
	std::vector<GzTriangle> triangles;
	KDNode* root;

	int StoreTriangle(GzTriangle triangle);
	void Build();
	bool GetNearestIntersectingSurface(GzRay ray, GzTriangle& closest_triangle, GzCoord intersection, GzCoord interpNorms);

private:
	bool SearchLeaf(std::vector<GzTriangle> tris, GzRay ray, GzTriangle& closest_triangle, GzCoord intersection, GzCoord interpNorms);
	KDNode* BuildHelper(std::vector<GzTriangle> tris, int depth);
	bool TraverseHelper(KDNode* rt, GzRay ray, float tStart, float tEnd, GzTriangle& closest_triangle, GzCoord intersection, GzCoord interpNorms);
};

