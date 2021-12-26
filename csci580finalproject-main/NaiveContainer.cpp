#include "stdafx.h"
#include "NaiveContainer.h"

int NaiveContainer::StoreTriangle(GzTriangle triangle) {
	triangles.push_back(triangle);
	return 0;
}

void NaiveContainer::Build() {}

bool NaiveContainer::GetNearestIntersectingSurface(GzRay ray, GzTriangle& closest_triangle, GzCoord intersection, GzCoord interpNorms) {
	float min_distance = INT_MAX;
	bool triangle_hit = false;
	for (int i = 0; i < triangles.size(); i++) {
		GzCoord norm, intersect;
		float distance;
		if (intersects(ray, triangles[i], distance, intersect, norm) && distance < min_distance) {
			min_distance = distance;
			closest_triangle = triangles[i];
			memcpy(intersection, intersect, sizeof(GzCoord));
			memcpy(interpNorms, norm, sizeof(GzCoord));
			triangle_hit = true;
		}
	}
	return triangle_hit;
}