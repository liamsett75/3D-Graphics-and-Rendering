#pragma once
#ifndef NAIVE_CONTAINER
#define NAIVE_CONTAINER

#include "AContainer.h"
#include <vector>

class NaiveContainer : public AContainer {
public:
	std::vector<GzTriangle> triangles; 

	int StoreTriangle(GzTriangle triangle);
	void Build();
	bool GetNearestIntersectingSurface(GzRay ray, GzTriangle& closest_triangle, GzCoord intersection, GzCoord interpNorms);
};

#endif