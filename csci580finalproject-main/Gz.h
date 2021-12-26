/*
 * Gz.h - include file for the cs580 rendering library
 */

/*
 * universal constants
 */
#define GZ_SUCCESS      0
#define GZ_FAILURE      1

#define MAX_INTENSITY 4095

/*
 * name list tokens
 */
#define GZ_NULL_TOKEN			0    /* triangle vert attributes */
#define GZ_POSITION             1
#define GZ_NORMAL               2
#define GZ_TEXTURE_INDEX        3

/* renderer-state default pixel color */
#define GZ_RGB_COLOR            99	

#define GZ_INTERPOLATE			95		/* interpolation mode */

#define GZ_DIRECTIONAL_LIGHT	79	/* directional light */
#define GZ_AMBIENT_LIGHT		78	/* ambient light type */

#define GZ_AMBIENT_COEFFICIENT		1001	/* Ka material property */
#define GZ_DIFFUSE_COEFFICIENT		1002	/* Kd material property */
#define GZ_SPECULAR_COEFFICIENT		1003	/* Ks material property */
#define GZ_DISTRIBUTION_COEFFICIENT	1004	/* specular power of material */

#define	GZ_TEXTURE_MAP	1010		/* texture function ptr */

#define GZ_CONTAINER_TYPE 1100

/*
 * value-list attributes
 */

/* select interpolation mode of the shader */
#define GZ_FLAT			0	/* do flat shading with GZ_RBG_COLOR */
#define	GZ_COLOR		1	/* interpolate vertex color */
#define	GZ_NORMALS		2	/* interpolate normals */

typedef int     GzToken;
typedef void    *GzPointer;
typedef float   GzColor[3];
typedef short   GzIntensity;	/* 0-4095 in lower 12-bits for RGBA */
typedef float	GzCoord[3];
typedef float	GzTextureIndex[2];
typedef float	GzMatrix[4][4];
typedef int	GzDepth;		/* signed z for clipping */

typedef	int	(*GzTexture)(float u, float v, GzColor color);	/* pointer to texture lookup method */
/* u,v parameters [0,1] are defined tex_fun(float u, float v, GzColor color) */

/*
 * Gz camera definition
 */
#ifndef GZCAMERA
#define GZCAMERA
typedef struct  GzCamera
{
  GzMatrix			Xiw;  	/* xform from world to image space */
  GzMatrix			Xpi;     /* perspective projection xform */
  GzCoord			position;  /* position of image plane origin */
  GzCoord			lookat;         /* position of look-at-point */
  GzCoord			worldup;   /* world up-vector (almost screen up) */
  float				FOV;            /* horizontal field of view */
} GzCamera;
#endif

#ifndef DATASTRUCTS
#define DATASTRUCTS
struct PlaneEq {
	GzCoord w;
	long double a;
	long double b;
	long double c;
	long double d; // a * x + b * y + c * z + d = 0

	PlaneEq(GzCoord vert1, GzCoord vert2, GzCoord vert3) {
		GzCoord edge1, edge2;
		createEdgeVector(edge1, vert1, vert2);
		createEdgeVector(edge2, vert1, vert3);

		a = ((long double)edge1[1] * edge2[2]) - ((long double)edge1[2] * edge2[1]);
		b = ((long double)edge1[2] * edge2[0]) - ((long double)edge1[0] * edge2[2]);
		c = ((long double)edge1[0] * edge2[1]) - ((long double)edge1[1] * edge2[0]);
		GzCoord temp = { a, b, c };
		memcpy(w, temp, sizeof(GzCoord));
		this->normalizeCoord(w);
		a = w[0];
		b = w[1];
		c = w[2];

		d = -1 * ((a * vert1[0]) + (b * vert1[1]) + (c * vert1[2]));
		
	}

	void createEdgeVector(GzCoord edgeVector, GzCoord vert1, GzCoord vert2) {
		edgeVector[0] = vert2[0] - vert1[0];
		edgeVector[1] = vert2[1] - vert1[1];
		edgeVector[2] = vert2[2] - vert1[2];
	}

	long double calculateZ(int x, int y) {
		long double numerator = -1 * ((a * x) + (b * y) + d);
		return numerator / c;
	}

private:
	void normalizeCoord(GzCoord coord) {
		float l2_norm = sqrt(((double)coord[0] * coord[0]) + ((double)coord[1] * coord[1]) + ((double)coord[2] * coord[2]));
		for (int i = 0; i < 3; i++) {
			coord[i] /= l2_norm;
		}
	}
};

struct LineEq {
	double a;
	double b;
	double c;

	LineEq(GzCoord vert1, GzCoord vert2) {
		double dX = (double)vert2[0] - (double)vert1[0];
		double dY = (double)vert2[1] - (double)vert1[1];
		a = dY;
		b = -1 * dX;
		c = (dX * vert1[1]) - (dY * vert1[0]);
	}

	float calculateX(float y) {
		double numerator = (-1 * b * y) - c;
		return numerator / a;
	}

	bool isBelowLine(int x, int y) {
		return (a * x) + (b * y) + c < 0;
	}

	bool isAboveLine(int x, int y) {
		return (a * x) + (b * y) + c > 0;
	}

	bool isOnLine(int x, int y) {
		return (a * x) + (b * y) + c == 0;
	}
};

/*
 * Information for a triangle, given in ***world*** space.
 */
typedef struct	GzTriangle
{
	GzCoord			vertices[3];
	GzCoord			normals[3];
	GzTextureIndex	textureCoords[3];
	PlaneEq			*trianglePlane;
	bool			textured;

	GzTriangle() {};

	GzTriangle(GzCoord vertexList[3], GzCoord normalList[3], GzTextureIndex uvList[3]) {
		memcpy(vertices, vertexList, sizeof(vertices));
		memcpy(normals, normalList, sizeof(normals));
		memcpy(textureCoords, uvList, sizeof(textureCoords));

		trianglePlane = new PlaneEq(vertices[0], vertices[1], vertices[2]);
	}
} GzTriangle;

/*
 * Also given in ***world*** space.
 */
typedef struct	GzRay
{
	GzCoord			origin;
	GzCoord			direction; //please make sure this is normalized
} GzRay;
#endif

#ifndef GZLIGHT
#define GZLIGHT
typedef struct  GzLight
{
  GzCoord        direction;    /* vector from surface to light */
  GzColor        color;		/* light color intensity */
} GzLight;
#endif

#ifndef GZINPUT
#define GZINPUT
typedef struct  GzInput
{
	GzCoord         rotation;       /* object rotation */
	GzCoord			translation;	/* object translation */
	GzCoord			scale;			/* object scaling */
	GzCamera		camera;			/* camera */
} GzInput;
#endif

#define RED     0         /* array indicies for color vector */
#define GREEN   1
#define BLUE    2

#define X       0      /* array indicies for position vector */
#define Y       1
#define Z       2

#define U       0       /* array indicies for texture coords */
#define V       1


#ifndef GZ_PIXEL
typedef	struct {
  GzIntensity    red;	
  GzIntensity    green;
  GzIntensity    blue;
  GzIntensity    alpha;
  GzDepth	 z;
} GzPixel;
#define GZ_PIXEL
#endif;

#define	MAXXRES	1024	/* put some bounds on size in case of error */
#define	MAXYRES	1024
