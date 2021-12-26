#include "stdafx.h"
#include "KDTree.h"

int KDTree::StoreTriangle(GzTriangle triangle)
{
	triangles.push_back(triangle);
	return 0;
}

float TriangleKey(GzTriangle tri, int dim) {
	return tri.vertices[0][dim];
}

//x = 3;
//(-2, 1, 3) ( 3, 3, 3) (0, 1, 2)
bool IntersectsPlane(GzTriangle tri, int dim, float separator) {
	float diffs[3];
	for (int i = 0; i < 3; i++) {
		diffs[i] = tri.vertices[i][dim] - separator;
	}

	bool allLeqThan = true;
	bool allGeqThan = true;

	for (int i = 0; i < 3; i++) {
		allLeqThan = allLeqThan && diffs[i] <= -0.001;
		allGeqThan = allGeqThan && diffs[i] >= 0.001;
	}

	return !(allLeqThan || allGeqThan);
}

KDNode* KDTree::BuildHelper(std::vector<GzTriangle> tris, int depth)
{
	KDNode* newNode = new KDNode();
	if (depth == 0) {
		newNode->isLeaf = true;
		newNode->contents = tris;
	} else {
		newNode->isLeaf = false;
		newNode->dim = depth % 3;
		int dim;
		dim = newNode->dim;

		newNode->contents = tris; //TODO: REMOVE THIS

		std::vector<GzTriangle> *leftSide = new std::vector<GzTriangle>();
		std::vector<GzTriangle> *rightSide = new std::vector<GzTriangle>();
		
		std::sort(tris.begin(), tris.end(), [dim](GzTriangle t1, GzTriangle t2) -> bool {
			return TriangleKey(t1, dim) > TriangleKey(t2, dim);
			});

		// for each triangle, if it doesn't intersect(), then put it in the vector specified by 

		int half = tris.size() / 2;
		float sep = TriangleKey(tris[half], dim);
		newNode->separator = sep;

		for (int i = 0; i < tris.size(); i++) {
			if (IntersectsPlane(tris[i], dim, sep)) {
				leftSide->push_back(tris[i]);
				rightSide->push_back(tris[i]);
			}
			else if (TriangleKey(tris[i], dim) < sep) {
				leftSide->push_back(tris[i]);
			}
			else {
				rightSide->push_back(tris[i]);
			}
		}

		newNode->leftChild = BuildHelper(*leftSide, depth - 1);
		newNode->rightChild = BuildHelper(*rightSide, depth - 1);
	}
	return newNode;
}

void KDTree::Build()
{

	/* Simplification: Since we don't need to find *the best* solution, we can
	 * just proceed for a fixed number of iterations. 2^9 = 512 cells should give
	 * measurable improvement, but we could go higher (this step can be slow)
	 */
	this->root = BuildHelper(triangles, 9);
	int x = 3;
}

bool KDTree::SearchLeaf(std::vector<GzTriangle> tris, GzRay ray, GzTriangle& closest_triangle, GzCoord intersection, GzCoord interpNorms) {
	float min_distance = INT_MAX;
	bool triangle_hit = false;
	for (int i = 0; i < tris.size(); i++) {
		GzCoord norm, intersect;
		float distance;
		if (intersects(ray, tris[i], distance, intersect, norm) && distance < min_distance) {
			min_distance = distance;
			closest_triangle = tris[i];
			memcpy(intersection, intersect, sizeof(GzCoord));
			memcpy(interpNorms, norm, sizeof(GzCoord));
			triangle_hit = true;
		}
	}
	return triangle_hit;
}

bool KDTree::TraverseHelper(KDNode* rt, GzRay ray, float tStart, float tEnd, GzTriangle& closest_triangle, GzCoord intersection, GzCoord interpNorms)
{
	if (rt->isLeaf) {
		SearchLeaf(rt->contents, ray, closest_triangle, intersection, interpNorms);
	}
	else {
		// Check beginning and end of ray
		// if ray crosses our plane, need to find intersection
		GzCoord start, end;
		for (int i = 0; i < 3; i++) {
			start[i] = ray.origin[i] + ray.direction[i] * tStart;
			end[i] = ray.origin[i] + ray.direction[i] * tEnd;
		}
		
		side startSide = start[rt->dim] < rt->separator ? side::L : side::R;
		side endSide = end[rt->dim] < rt->separator ? side::L : side::R;
		float tMid;

		if (startSide != endSide) {
			tMid = (rt->separator - start[rt->dim]) / ray.direction[rt->dim];
		}

		if (startSide == side::L && endSide == side::R) {
			if (TraverseHelper(rt->leftChild, ray, tStart, tMid, closest_triangle, intersection, interpNorms)) {
				return true;
			}
			else {
				return TraverseHelper(rt->rightChild, ray, tMid, tEnd, closest_triangle, intersection, interpNorms);
			}
		}
		else if (startSide == side::R && endSide == side::L) {
			if (TraverseHelper(rt->rightChild, ray, tStart, tMid, closest_triangle, intersection, interpNorms)) {
				return true;
			}
			else {
				return TraverseHelper(rt->leftChild, ray, tMid, tEnd, closest_triangle, intersection, interpNorms);
			}
		}
		else if (startSide == side::L && endSide == side::L) {
			return TraverseHelper(rt->leftChild, ray, tStart, tEnd, closest_triangle, intersection, interpNorms);
		}
		else if (startSide == side::R && endSide == side::R) {
			return TraverseHelper(rt->rightChild, ray, tStart, tEnd, closest_triangle, intersection, interpNorms);
		}
	}
	assert(1 == 0);
}

bool KDTree::GetNearestIntersectingSurface(GzRay ray, GzTriangle& closest_triangle, GzCoord intersection, GzCoord interpNorms)
{
	return TraverseHelper(this->root, ray, 0.0f, 50.0f, closest_triangle, intersection, interpNorms);
}

//bool AContainer::intersects(GzRay ray, GzTriangle triangle, float &t, GzCoord intersection, GzCoord interpNorms)
