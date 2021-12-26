#include "Gz.h"

#pragma once
class AContainer
{
public:
	// returns true if ray intersects triangle, false if not
	// t returns distance from origin to intersection
	bool intersects(GzRay ray, GzTriangle triangle, float& t, GzCoord intersection, GzCoord interpNorms);

	virtual int StoreTriangle(GzTriangle triangle) = 0;

	virtual void Build() = 0;

	// returns true if ray intersects a triangle, false if it does not
	// if true, closest_triangle will have the value of the first intersecting triangle
	virtual bool GetNearestIntersectingSurface(GzRay ray, GzTriangle& closest_triangle, GzCoord intersection, GzCoord interpNorms) = 0;
};
