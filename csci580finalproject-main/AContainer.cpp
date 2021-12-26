#include "stdafx.h"
#include "AContainer.h"

void crossProduct(GzCoord a, GzCoord b, GzCoord result) {
	result[0] = a[1] * b[2] - a[2] * b[1];
	result[1] = a[2] * b[0] - a[0] * b[2];
	result[2] = a[0] * b[1] - a[1] * b[0];
}

float dotProduct2(GzCoord a, GzCoord b) {
	return a[X] * b[X] + a[Y] * b[Y] + a[Z] * b[Z];
}

void componentwiseAdd(GzCoord a, GzCoord b, GzCoord output) {
	output[X] = a[X] + b[X];
	output[Y] = a[Y] + b[Y];
	output[Z] = a[Z] + b[Z];
}

void componentwiseSubtract(GzCoord a, GzCoord b, GzCoord output) {
	output[X] = a[X] - b[X];
	output[Y] = a[Y] - b[Y];
	output[Z] = a[Z] - b[Z];
}

void componentwiseMultiply(GzCoord a, GzCoord b, GzCoord output) {
	output[X] = a[X] * b[X];
	output[Y] = a[Y] * b[Y];
	output[Z] = a[Z] * b[Z];
}

void scalarMultiply(GzCoord a, float b, GzCoord output) {
	output[X] = a[X] * b;
	output[Y] = a[Y] * b;
	output[Z] = a[Z] * b;
}

void createEdgeVector2(GzCoord edgeVector, GzCoord vert1, GzCoord vert2) {
	edgeVector[0] = vert2[0] - vert1[0];
	edgeVector[1] = vert2[1] - vert1[1];
	edgeVector[2] = vert2[2] - vert1[2];
}

void threeToTwoTransform(GzCoord v1, GzCoord v2, GzCoord in, GzCoord out) {
	out[0] = 0;
	out[1] = 0;

	for (int i = 0; i < 3; i++) {
		out[0] += in[i] * v1[i];
		out[1] += in[i] * v2[i];
	}
}

void normalizeCoord2(GzCoord coord) {
	float l2_norm = sqrt(((double)coord[X] * coord[X]) + ((double)coord[Y] * coord[Y]) + ((double)coord[Z] * coord[Z]));
	for (int i = 0; i < 3; i++) {
		coord[i] /= l2_norm;
	}
}

bool AContainer::intersects(GzRay ray, GzTriangle triangle, float &t, GzCoord intersection, GzCoord interpNorms)
{
	/* Redefine our plane (without changing the plane itself) so that the plane normal is
     * oriented towards the camera direction. We need to do this so that our cross products
	 * are meaningful when we check whether the point is inside the triangle. 
	 */
	GzCoord _w; // triangle plane w, but perhaps flipped by -1
	float _d; // triangle plane d, but perhaps flipped by -1
	memcpy(_w, triangle.trianglePlane->w, sizeof(GzCoord));
	_d = triangle.trianglePlane->d;
	if (dotProduct2(_w, ray.direction) > 0) {
		 scalarMultiply(_w, -1, _w);
		 _d *= -1;
	}

    //origin + t * direction = intersection 
    t = (-_d - dotProduct2(_w, ray.origin))
        / dotProduct2(_w, ray.direction);

	// plane is behind ray
	if (t < 0) {
		return false;
	}

    scalarMultiply(ray.direction, t, intersection);
    componentwiseAdd(ray.origin, intersection, intersection);

	// based on notation from https://courses.cs.washington.edu/courses/csep557/10au/lectures/triangle_intersection.pdf

	GzCoord ab, aq, ac, bc, bq, ca, cq;
	componentwiseSubtract(triangle.vertices[1], triangle.vertices[0], ab);
	componentwiseSubtract(triangle.vertices[2], triangle.vertices[1], bc);
	componentwiseSubtract(triangle.vertices[0], triangle.vertices[2], ca);

	componentwiseSubtract(intersection, triangle.vertices[0], aq);
	componentwiseSubtract(intersection, triangle.vertices[1], bq);
	componentwiseSubtract(intersection, triangle.vertices[2], cq);

	// used only for baycentric coordinates
	componentwiseSubtract(triangle.vertices[2], triangle.vertices[0], ac);

	GzCoord abXaq, bcXbq, caXcq, abXac;
	crossProduct(ab, aq, abXaq);
	crossProduct(bc, bq, bcXbq);
	crossProduct(ca, cq, caXcq);

	crossProduct(ab, ac, abXac);

	GzCoord baryCoords;
	float denom = dotProduct2(abXac, _w);
	baryCoords[0] = dotProduct2(bcXbq, _w) / denom; // alpha
	baryCoords[1] = dotProduct2(caXcq, _w) / denom; // beta
	baryCoords[2] = dotProduct2(abXaq, _w) / denom; // gamma

	for (int i = 0; i < 3; i++) {
		interpNorms[i] = baryCoords[0] * triangle.normals[0][i] + baryCoords[1] * triangle.normals[1][i] + baryCoords[2] * triangle.normals[2][i];
	}
	normalizeCoord2(interpNorms);

	bool test1, test2, test3;
	test1 = dotProduct2(abXaq, _w) >= 0;
	test2 = dotProduct2(bcXbq, _w) >= 0;
	test3 = dotProduct2(caXcq, _w) >= 0;

	bool is_in_triangle = test1 && test2 && test3;

	if (is_in_triangle) {
		int x;
		x = 3;
	}
    return is_in_triangle;
}
