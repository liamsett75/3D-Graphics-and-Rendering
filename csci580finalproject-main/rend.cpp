/* CS580 Homework 3 */

#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"

#define PI (float) 3.14159265358979323846

#define EMPTY_STACK (short) -1

GzColor DEF_BG_COLOR = { 0.5, 0.4375, 0.375 };

float degreeToRad(float degree) {
	return (degree * PI) / 180;
}

float dotProduct(GzCoord a, GzCoord b) {
	return a[X] * b[X] + a[Y] * b[Y] + a[Z] * b[Z];
}

void setMatrixToZero(GzMatrix mat) {
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			mat[i][j] = 0;
		}
	}
}

void copyMatrix(GzMatrix dest, GzMatrix src) {
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			dest[i][j] = src[i][j];
		}
	}
}

void matrix_dot_product(GzMatrix dest, GzMatrix mat_a, GzMatrix mat_b) {
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			float dot_product = 0;
			for (int k = 0; k < 4; k++) {
				dot_product += mat_a[i][k] * mat_b[k][j];
			}
			dest[i][j] = dot_product;
		}
	}
}

void create_identity_matrix(GzMatrix matrix) {
	setMatrixToZero(matrix);
	matrix[0][0] = 1;
	matrix[1][1] = 1;
	matrix[2][2] = 1;
	matrix[3][3] = 1;
}

void normalize_matrix(GzMatrix matrix) {
	double a = matrix[0][0];
	double b = matrix[0][1];
	double c = matrix[0][2];

	double l2_norm = sqrt((a * a) + (b * b) + (c * c));

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			matrix[i][j] /= l2_norm;
		}
	}
}

void createNormTransformMatrix(GzMatrix dest, GzMatrix orig) {
	copyMatrix(dest, orig);
	dest[0][3] = 0;
	dest[1][3] = 0;
	dest[2][3] = 0;
	
	normalize_matrix(dest);
}

int GzRender::GzRotXMat(float degree, GzMatrix mat)
{
/* HW 3.1
// Create rotate matrix : rotate along x axis
// Pass back the matrix using mat value
*/

	float rad = degreeToRad(degree);
	setMatrixToZero(mat);


	mat[0][0] = 1;
	mat[1][1] = cos(rad);
	mat[1][2] = -1 * sin(rad);
	mat[2][1] = sin(rad);
	mat[2][2] = cos(rad);
	mat[3][3] = 1;

	return GZ_SUCCESS;
}

int GzRender::GzRotYMat(float degree, GzMatrix mat)
{
/* HW 3.2
// Create rotate matrix : rotate along y axis
// Pass back the matrix using mat value
*/

	float rad = degreeToRad(degree);
	setMatrixToZero(mat);

	mat[0][0] = cos(rad);
	mat[0][2] = sin(rad);
	mat[1][1] = 1;
	mat[2][0] = -1 * sin(rad);
	mat[2][2] = cos(rad);
	mat[3][3] = 1;

	return GZ_SUCCESS;
}

int GzRender::GzRotZMat(float degree, GzMatrix mat)
{
/* HW 3.3
// Create rotate matrix : rotate along z axis
// Pass back the matrix using mat value
*/

	float rad = degreeToRad(degree);
	setMatrixToZero(mat);

	mat[0][0] = cos(rad);
	mat[0][1] = -1 * sin(rad);
	mat[1][0] = sin(rad);
	mat[1][1] = cos(rad);
	mat[2][2] = 1;
	mat[3][3] = 1;

	return GZ_SUCCESS;
}

int GzRender::GzTrxMat(GzCoord translate, GzMatrix mat)
{
/* HW 3.4
// Create translation matrix
// Pass back the matrix using mat value
*/
	setMatrixToZero(mat);
	
	mat[0][0] = 1;
	mat[1][1] = 1;
	mat[2][2] = 1;
	mat[3][3] = 1;

	mat[0][3] = translate[0];
	mat[1][3] = translate[1];
	mat[2][3] = translate[2];

	return GZ_SUCCESS;
}


int GzRender::GzScaleMat(GzCoord scale, GzMatrix mat)
{
/* HW 3.5
// Create scaling matrix
// Pass back the matrix using mat value
*/

	setMatrixToZero(mat);

	mat[0][0] = scale[0];
	mat[1][1] = scale[1];
	mat[2][2] = scale[2];
	mat[3][3] = 1;

	return GZ_SUCCESS;
}


/***********************************************/
/* HW1 methods: copy here the methods from HW1 */

// helper function clamps Intensity to 0 <= intensity <= 4095
GzIntensity clampIntensity(GzIntensity intensity) {
	if (intensity < 0) {
		return 0;
	}

	if (intensity > 4095) {
		return 4095;
	}

	return intensity;
}

unsigned char intensityToChar(GzIntensity intensity) {
	intensity = intensity >> 4;
	return (unsigned char)intensity;
}

float clampColor(float color) {
	if (color < 0.0) {
		return 0.0;
	}

	if (color > 1.0) {
		return 1.0;
	}

	return color;
}

GzRender::GzRender(int xRes, int yRes)
{
	/* HW1.1 create a framebuffer for MS Windows display:
	 -- set display resolution
	 -- allocate memory for framebuffer : 3 bytes(b, g, r) x width x height
	 -- allocate memory for pixel buffer
	 */

	this->xres = xRes;
	this->yres = yRes;
	this->framebuffer = new char[3 * xRes * yRes];
	this->pixelbuffer = new GzPixel[xRes * yRes];
	this->timebuffer = new float[xRes * yRes];

	/* HW 3.6
	- setup Xsp and anything only done once
	- init default camera
	*/

	setMatrixToZero(this->Xsp);

	this->Xsp[0][0] = xRes / 2.0;
	this->Xsp[0][3] = xRes / 2.0;
	this->Xsp[1][1] = (-1.0 * yRes) / 2;
	this->Xsp[1][3] = yRes / 2.0;
	this->Xsp[2][2] = MAXINT;
	this->Xsp[3][3] = 1;

	this->m_camera.FOV = DEFAULT_FOV;
	this->m_camera.lookat[0] = 0;
	this->m_camera.lookat[1] = 0;
	this->m_camera.lookat[2] = 0;
	this->m_camera.worldup[0] = 0;
	this->m_camera.worldup[1] = 1;
	this->m_camera.worldup[2] = 0;
	this->m_camera.position[0] = DEFAULT_IM_X;
	this->m_camera.position[1] = DEFAULT_IM_Y;
	this->m_camera.position[2] = DEFAULT_IM_Z;

	this->matlevel = EMPTY_STACK;

}

GzRender::~GzRender()
{
	/* HW1.2 clean up, free buffer memory */
	delete this->framebuffer;
	delete this->pixelbuffer;
	delete this->timebuffer;

}

int GzRender::GzDefault()
{
	/* HW1.3 set pixel buffer to some default values - start a new frame */
	/*
	int pxl_buffer_size = this->xres * this->yres;
	for (int i = 0; i < pxl_buffer_size; i++) {
		this->pixelbuffer[i].red = 2048;
		this->pixelbuffer[i].green = 1792;
		this->pixelbuffer[i].blue = 1536;
		this->pixelbuffer[i].z = INT32_MAX;
	}
	*/
	return GZ_SUCCESS;
}

float norm(GzCoord coord) {
	return sqrt(pow(coord[X], 2) + pow(coord[Y], 2) + pow(coord[Z], 2));
}

void calculateXiw(GzMatrix Xiw, GzCamera camera) {
	GzCoord cl, vec_x, vec_y, vec_z, up_prime;

	for (int i = 0; i < 3; i++) {
		cl[i] = camera.lookat[i] - camera.position[i];
	}

	float cl_norm = norm(cl);
	for (int i = 0; i < 3; i++) {
		vec_z[i] = cl[i] / cl_norm;
	}

	float up_dotproduct_z = dotProduct(camera.worldup, vec_z);
	for (int i = 0; i < 3; i++) {
		up_prime[i] = camera.worldup[i] - (up_dotproduct_z * vec_z[i]);
	}

	float up_prime_norm = norm(up_prime);
	for (int i = 0; i < 3; i++) {
		vec_y[i] = up_prime[i] / up_prime_norm;
	}

	// Y x Z
	vec_x[X] = vec_y[Y] * vec_z[Z] - vec_y[Z] * vec_z[Y];
	vec_x[Y] = vec_y[Z] * vec_z[X] - vec_y[X] * vec_z[Z];
	vec_x[Z] = vec_y[X] * vec_z[Y] - vec_y[Y] * vec_z[X];

	float x_dotproduct_c = dotProduct(vec_x, camera.position);
	float y_dotproduct_c = dotProduct(vec_y, camera.position);
	float z_dotproduct_c = dotProduct(vec_z, camera.position);

	setMatrixToZero(Xiw);
	Xiw[0][0] = vec_x[X];
	Xiw[0][1] = vec_x[Y];
	Xiw[0][2] = vec_x[Z];
	Xiw[0][3] = -1 * x_dotproduct_c;

	Xiw[1][0] = vec_y[X];
	Xiw[1][1] = vec_y[Y];
	Xiw[1][2] = vec_y[Z];
	Xiw[1][3] = -1 * y_dotproduct_c;

	Xiw[2][0] = vec_z[X];
	Xiw[2][1] = vec_z[Y];
	Xiw[2][2] = vec_z[Z];
	Xiw[2][3] = -1 * z_dotproduct_c;

	Xiw[3][3] = 1;

}

void calculateXpi(GzMatrix Xpi, GzCamera camera) {
	setMatrixToZero(Xpi);

	float rad = degreeToRad(camera.FOV);

	float d_inv = tan(rad / 2);

	Xpi[0][0] = 1;
	Xpi[1][1] = 1;
	Xpi[2][2] = d_inv;
	Xpi[3][2] = d_inv;
	Xpi[3][3] = 1;

}

int GzRender::GzBeginRender()
{
/* HW 3.7 
- setup for start of each frame - init frame buffer color,alpha,z
- compute Xiw and projection xform Xpi from camera definition 
- init Ximage - put Xsp at base of stack, push on Xpi and Xiw 
- now stack contains Xsw and app can push model Xforms when needed 
*/ 

	calculateXpi(this->m_camera.Xpi, this->m_camera);
	calculateXiw(this->m_camera.Xiw, this->m_camera);

	if (this->GzPushMatrix(this->Xsp) || this->GzPushMatrix(this->m_camera.Xpi) || this->GzPushMatrix(this->m_camera.Xiw)) {
		return GZ_FAILURE;
	}

	return GZ_SUCCESS;
}

int GzRender::GzPutCamera(GzCamera camera)
{
/* HW 3.8 
/*- overwrite renderer camera structure with new camera definition
*/
	memcpy(&this->m_camera, &camera, sizeof(camera));
	return GZ_SUCCESS;	
}

int GzRender::GzPushMatrix(GzMatrix	matrix)
{
/* HW 3.9 
- push a matrix onto the Ximage stack
- check for stack overflow
*/

	if (this->matlevel == MATLEVELS) {
		return GZ_FAILURE;
	}

	GzMatrix identity_matrix;
	create_identity_matrix(identity_matrix);

	if (this->matlevel == EMPTY_STACK) {
		this->matlevel++;
		copyMatrix(this->Ximage[matlevel], matrix);
		copyMatrix(this->Xnorm[matlevel], identity_matrix);
		return GZ_SUCCESS;
	}

	this->matlevel++;

	matrix_dot_product(this->Ximage[this->matlevel], this->Ximage[this->matlevel - 1], matrix);

	if (this->matlevel == 1) {
		copyMatrix(this->Xnorm[matlevel], identity_matrix);
	}
	else {
		GzMatrix norm_transform_matrix;
		createNormTransformMatrix(norm_transform_matrix, matrix);
		matrix_dot_product(this->Xnorm[this->matlevel], this->Xnorm[this->matlevel - 1], norm_transform_matrix);
	}
	
	return GZ_SUCCESS;
}

int GzRender::GzPopMatrix()
{
/* HW 3.10
- pop a matrix off the Ximage stack
- check for stack underflow
*/

	if (this->matlevel == EMPTY_STACK) {
		return GZ_FAILURE;
	}

	this->matlevel--;

	return GZ_SUCCESS;
}

int GzRender::GzPut(int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z)
{
	/* HW1.4 write pixel values into the buffer */
		// bounds check
	if (i < 0 || j < 0 || i >= this->xres || j >= this->yres) {
		return GZ_SUCCESS;
	}

	// if current z is less than new z || less than 0, dont store new value
	if (z < 0 || this->pixelbuffer[ARRAY(i, j)].z < z) {
		return GZ_SUCCESS;
	}

	this->pixelbuffer[ARRAY(i, j)].red = clampIntensity(r);
	this->pixelbuffer[ARRAY(i, j)].green = clampIntensity(g);
	this->pixelbuffer[ARRAY(i, j)].blue = clampIntensity(b);
	this->pixelbuffer[ARRAY(i, j)].alpha = a;
	this->pixelbuffer[ARRAY(i, j)].z = z;
	return GZ_SUCCESS;
}

int GzRender::GzPutTime(int i, int j, std::chrono::duration<double, std::nano> d)
{
	this->timebuffer[ARRAY(i, j)] = d.count();
	return GZ_SUCCESS;
}


int GzRender::GzGet(int i, int j, GzIntensity* r, GzIntensity* g, GzIntensity* b, GzIntensity* a, GzDepth* z)
{
	/* HW1.5 retrieve a pixel information from the pixel buffer */
	if (i < 0 || j < 0 || i >= this->xres || j >= this->yres) {
		return GZ_FAILURE;
	}

	*r = this->pixelbuffer[ARRAY(i, j)].red;
	*g = this->pixelbuffer[ARRAY(i, j)].green;
	*b = this->pixelbuffer[ARRAY(i, j)].blue;
	*a = this->pixelbuffer[ARRAY(i, j)].alpha;
	*z = this->pixelbuffer[ARRAY(i, j)].z;
	return GZ_SUCCESS;
}


int GzRender::GzFlushDisplay2File(FILE* outfile)
{
	/* HW1.6 write image to ppm file -- "P6 %d %d 255\r" */
	
	char* header_buf = new char[2000];
	sprintf(header_buf, "P6 %d %d 255\r", this->xres, this->yres);
	fputs(header_buf, outfile);
	delete[] header_buf;

	int buffersize = (3 * this->xres * this->yres);
	char* output = new char[buffersize + 1];
	for (int i = 0; i < buffersize; i += 3) {
		int pixel_num = i / 3;
		if (i + 2 > buffersize) {
			printf("Error: overflowing flush to ppm file");
			return GZ_FAILURE;
		}
		output[i] = intensityToChar(this->pixelbuffer[pixel_num].red);
		output[i + 1] = intensityToChar(this->pixelbuffer[pixel_num].green);
		output[i + 2] = intensityToChar(this->pixelbuffer[pixel_num].blue);
	}
	output[buffersize] = '\0';
	fwrite(output, sizeof(unsigned char), buffersize, outfile);
	return GZ_SUCCESS;
}

int GzRender::GzFlushDisplay2FrameBuffer()
{
	/* HW1.7 write pixels to framebuffer:
		- put the pixels into the frame buffer
		- CAUTION: when storing the pixels into the frame buffer, the order is blue, green, and red
		- NOT red, green, and blue !!!
	*/
	int framebuffer_size = 3 * this->xres * this->yres;
	for (int i = 0; i < framebuffer_size; i += 3) {
		int pixel_num = i / 3;
		if (i + 2 > framebuffer_size) {
			printf("Error: overflowing flush to frame buffer");
			return GZ_FAILURE;
		}
		this->framebuffer[i] = intensityToChar(this->pixelbuffer[pixel_num].blue);
		this->framebuffer[i + 1] = intensityToChar(this->pixelbuffer[pixel_num].green);
		this->framebuffer[i + 2] = intensityToChar(this->pixelbuffer[pixel_num].red);
	}
	return GZ_SUCCESS;
}


/***********************************************/
/* HW2 methods: implement from here */

void swapVertices(GzCoord vertexList[3], int indx1, int indx2) {
	GzCoord swap;
	memcpy(swap, vertexList[indx1], sizeof(swap));
	memcpy(vertexList[indx1], vertexList[indx2], sizeof(vertexList[indx1]));
	memcpy(vertexList[indx2], swap, sizeof(vertexList[indx2]));
}

void swapTextures(GzTextureIndex uvList[3], int indx1, int indx2) {
	GzTextureIndex swap;
	memcpy(swap, uvList[indx1], sizeof(swap));
	memcpy(uvList[indx1], uvList[indx2], sizeof(uvList[indx1]));
	memcpy(uvList[indx2], swap, sizeof(uvList[indx2]));
}

void sortVerticesByY(GzCoord vertexList[3], GzCoord normalList[3], GzTextureIndex uvList[3]) {
	for (int i = 0; i < 3; i++) {
		int max_index = i;
		for (int j = i; j < 3; j++) {
			bool is_y_lower = vertexList[j][1] < vertexList[max_index][1];
			if (is_y_lower) {
				max_index = j;
			}
		}
		swapVertices(vertexList, i, max_index);
		swapVertices(normalList, i, max_index);
		swapTextures(uvList, i, max_index);
	}
}

void sortVerticesByX(GzCoord vertexList[3], GzCoord normalList[3], GzTextureIndex uvList[3]) {
	bool is_flat_top = (int)round(vertexList[0][1]) == (int)round(vertexList[1][1]);
	bool is_flat_bottom = (int)round(vertexList[1][1]) == (int)round(vertexList[2][1]);

	if (is_flat_top) {
		if (vertexList[0][0] > vertexList[1][0]) {
			swapVertices(vertexList, 1, 2);
			swapVertices(normalList, 1, 2);
			swapTextures(uvList, 1, 2);
		}
		return;
	}

	if (is_flat_bottom) {
		if (vertexList[1][0] < vertexList[2][0]) {
			swapVertices(vertexList, 1, 2);
			swapVertices(normalList, 1, 2);
			swapTextures(uvList, 1, 2);
		}
		return;
	}

	LineEq line_v1_v3 = LineEq(vertexList[0], vertexList[2]);
	float x1 = line_v1_v3.calculateX(vertexList[1][1]);
	if (x1 > vertexList[1][0]) {
		swapVertices(vertexList, 1, 2);
		swapVertices(normalList, 1, 2);
		swapTextures(uvList, 1, 2);
	}
}

void sortVertices(GzCoord vertexList[3], GzCoord normalList[3], GzTextureIndex uvList[3]) {
	sortVerticesByY(vertexList, normalList, uvList);
	sortVerticesByX(vertexList, normalList, uvList);
}

void getBoundsBox(GzCoord vertexList[3], int* minX, int* maxX, int* minY, int* maxY) {
	*minX = INT32_MAX;
	*maxX = INT32_MIN;
	*minY = INT32_MAX;
	*maxY = INT32_MIN;
	for (int i = 0; i < 3; i++) {
		*minX = min(*minX, round(vertexList[i][0]));
		*maxX = max(*maxX, round(vertexList[i][0]));
		*minY = min(*minY, round(vertexList[i][1]));
		*maxY = max(*maxY, round(vertexList[i][1]));
	}
}

void setColorToZero(GzColor color) {
	color[RED] = 0;
	color[GREEN] = 0;
	color[BLUE] = 0;
}

float coordCoordDotProduct(GzCoord coord_a, GzCoord coord_b) {
	float dot_product = 0;
	for (int i = 0; i < 3; i++) {
		dot_product += coord_a[i] * coord_b[i];
	}
	return dot_product;
}

void calculateR(GzCoord normal, GzCoord R, GzLight light, float n_dot_l) {
	R[X] = (2 * n_dot_l * normal[X]) - light.direction[X];
	R[Y] = (2 * n_dot_l * normal[Y]) - light.direction[Y];
	R[Z] = (2 * n_dot_l * normal[Z]) - light.direction[Z];
}

void flipCoord(GzCoord src, GzCoord dest) {
	dest[X] = -1 * src[X];
	dest[Y] = -1 * src[Y];
	dest[Z] = -1 * src[Z];
}

void normalizeCoord(GzCoord coord) {
	float l2_norm = sqrt(((double)coord[X] * coord[X]) + ((double)coord[Y] * coord[Y]) + ((double)coord[Z] * coord[Z]));
	for (int i = 0; i < 3; i++) {
		coord[i] /= l2_norm;
	}
}

float zsToPersp(float zs) {
	return zs / (MAXINT - zs);
}

void uvToPersp(GzTextureIndex src, GzTextureIndex dst, float zs) {
	float zp = zsToPersp(zs);
	dst[0] = src[0] / (zp + 1);
	dst[1] = src[1] / (zp + 1);
}

void uvToScreen(GzTextureIndex src, GzTextureIndex dst, float zs) {
	float zp = zsToPersp(zs);
	dst[0] = src[0] * (zp + 1);
	dst[1] = src[1] * (zp + 1);
}

void GzRender::shading_equation(GzCoord normals, GzRender* render, GzColor C, GzTextureIndex uv, bool textured, GzColor gourad_term) {
	GzColor diffuse, specular, texture;
	GzCoord E;
	float n_dot_e;

	if (textured) {
		tex_fun(uv[0], uv[1], C);
		return;
	}

	setColorToZero(diffuse);
	setColorToZero(specular);

	E[X] = 0;
	E[Y] = 0;
	E[Z] = -1;

	n_dot_e = coordCoordDotProduct(normals, E);

	for (int i = 0; i < render->numlights; i++) {
		GzCoord R;
		float n_dot_l;

		n_dot_l = coordCoordDotProduct(normals, render->lights[i].direction);
		if (n_dot_l > 0 && n_dot_e > 0) {
			calculateR(normals, R, render->lights[i], n_dot_l);
		}
		else if (n_dot_l < 0 && n_dot_e < 0) {
			GzCoord flipped_normals;
			flipCoord(normals, flipped_normals);

			n_dot_l = coordCoordDotProduct(flipped_normals, render->lights[i].direction);
			calculateR(flipped_normals, R, render->lights[i], n_dot_l);
		}
		else {
			continue;
		}

		normalizeCoord(R);
		float r_dot_e = coordCoordDotProduct(R, E);
		if (r_dot_e < 0) {
			r_dot_e = 0;
		}
		else if (r_dot_e > 1) {
			r_dot_e = 1;
		}

		specular[RED] += render->lights[i].color[RED] * pow(r_dot_e, render->spec);
		specular[GREEN] += render->lights[i].color[GREEN] * pow(r_dot_e, render->spec);
		specular[BLUE] += render->lights[i].color[BLUE] * pow(r_dot_e, render->spec);

		diffuse[RED] += render->lights[i].color[RED] * n_dot_l;
		diffuse[GREEN] += render->lights[i].color[GREEN] * n_dot_l;
		diffuse[BLUE] += render->lights[i].color[BLUE] * n_dot_l;
	}

	if (gourad_term != NULL) {
		for (int i = 0; i < 3; i++) {
			gourad_term[i] = specular[i] + diffuse[i] + render->ambientlight.color[i];
		}
	}

	for (int i = 0; i < 3; i++) {
		C[i] = clampColor((render->Ks[i] * specular[i]) + (render->Kd[i] * diffuse[i]) + (render->Ka[i] * render->ambientlight.color[i]));
	}
}

int GzRender::GzPutAttribute(int numAttributes, GzToken* nameList, GzPointer* valueList)
{
	/* HW 2.1
	-- Set renderer attribute states (e.g.: GZ_RGB_COLOR default color)
	-- In later homeworks set shaders, interpolaters, texture maps, and lights
	*/

	for (int i = 0; i < numAttributes; i++) {
		switch (nameList[i]) {
		case GZ_RGB_COLOR:
			// Potentially valueList should be incremented by sizeof type(?) - info in hw2 ppt
			memcpy(this->flatcolor, (GzColor*)(valueList[i]), sizeof(this->flatcolor));
			this->flatcolor[0] = clampColor(this->flatcolor[0]);
			this->flatcolor[1] = clampColor(this->flatcolor[1]);
			this->flatcolor[2] = clampColor(this->flatcolor[2]);
			break;

		case GZ_INTERPOLATE:
			this->interp_mode = *((int*)(valueList[i]));
			break;

		case GZ_DIRECTIONAL_LIGHT:
			this->lights[numlights] = *(GzLight*)(valueList[i]);
			this->numlights++;
			break;

		case GZ_AMBIENT_LIGHT:
			this->ambientlight = *(GzLight*)(valueList[i]);
			break;

		case GZ_AMBIENT_COEFFICIENT:
			memcpy(this->Ka, (GzColor*)(valueList[i]), sizeof(this->Ka));
			break;

		case GZ_DIFFUSE_COEFFICIENT:
			memcpy(this->Kd, (GzColor*)(valueList[i]), sizeof(this->Kd));
			break;

		case GZ_SPECULAR_COEFFICIENT:
			memcpy(this->Ks, (GzColor*)(valueList[i]), sizeof(this->Ks));
			break;

		case GZ_DISTRIBUTION_COEFFICIENT:
			this->spec = *((float*)(valueList[i]));
			break;

		case GZ_TEXTURE_MAP:
			this->tex_fun = (GzTexture)(valueList[i]);
			break;

		case GZ_CONTAINER_TYPE:
			this->container = (AContainer*)(valueList[i]);
			break;
		}

	}

	return GZ_SUCCESS;
}

int GzRender::GzPutTriangle(int	numParts, GzToken* nameList, GzPointer* valueList)
/* numParts - how many names and values */
{
	/* HW 2.2
	-- Pass in a triangle description with tokens and values corresponding to
		  GZ_NULL_TOKEN:		do nothing - no values
		  GZ_POSITION:		3 vert positions in model space
	-- Invoke the rastrizer/scanline framework
	-- Return error code
	*/
	GzCoord vertexList[3], normalList[3];
	GzTextureIndex uvList[3];

	for (int i = 0; i < numParts; i++) {
		switch (nameList[i]) {
		case GZ_NULL_TOKEN:
			continue;

		case GZ_POSITION:
			memcpy(vertexList, (GzCoord*)valueList[i], sizeof(vertexList));

		case GZ_NORMAL:
			memcpy(normalList, (GzCoord*)valueList[i], sizeof(normalList));

		case GZ_TEXTURE_INDEX:
			memcpy(uvList, (GzTextureIndex*)valueList[i], sizeof(uvList));
		}
	}

	for (int i = 0; i < 3; i++) {
		normalizeCoord(normalList[i]);
	}
	// sortVertices(vertexList, normalList, uvList);
	GzTriangle triangle(vertexList, normalList, uvList);
	triangle.textured = texturing;
	this->container->StoreTriangle(triangle);

	/*
		GzMatrix vertex_list_4d, normal_list_4d, rot_normal_list_4d;
		setMatrixToZero(vertex_list_4d);
		setMatrixToZero(normal_list_4d);
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				vertex_list_4d[j][i] = vertexList[i][j];
				normal_list_4d[j][i] = normalList[i][j];
			}
			vertex_list_4d[3][i] = 1.0;
			normal_list_4d[3][i] = 1.0;
		}

		this->GzPushMatrix(vertex_list_4d);
		for (int i = 0; i < 3; i++) {
			float w = this->Ximage[this->matlevel][3][i];
			for (int j = 0; j < 3; j++) {
				vertexList[i][j] = this->Ximage[this->matlevel][j][i] / w;
			}
		}
		this->GzPopMatrix();

		matrix_dot_product(rot_normal_list_4d, this->Xnorm[this->matlevel], normal_list_4d);

		for (int i = 0; i < 3; i++) {
			float norm_w = rot_normal_list_4d[3][i];
			for (int j = 0; j < 3; j++) {
				normalList[i][j] = rot_normal_list_4d[j][i] / norm_w;
			}
		}

		for (int i = 0; i < 3; i++) {
			normalizeCoord(normalList[i]);
		}

		rasterizeTriangle(vertexList, normalList, uvList, this);
	*/


	return GZ_SUCCESS;
}

void vectorAddition(GzCoord a, GzCoord b, GzCoord res) {
	for (int i = 0; i < 3; i++) {
		res[i] = a[i] + b[i];
	}
}

void vectorSubtraction(GzCoord a, GzCoord b, GzCoord res) {
	for (int i = 0; i < 3; i++) {
		res[i] = a[i] - b[i];
	}
}

void scalarMultiplication(GzCoord a, float s, GzCoord res) {
	for (int i = 0; i < 3; i++) {
		res[i] = a[i] * s;
	}
}

void zeroColor(GzColor color) {
	color[RED] = 0;
	color[GREEN] = 0;
	color[BLUE] = 0;
}

void reflect(GzRay& ray, GzCoord origin, GzCoord norm) {
	float dDOTn = dotProduct(ray.direction, norm);
	GzCoord dDOTnTIMESn;
	scalarMultiplication(norm, 2 * dDOTn, dDOTnTIMESn);
	vectorSubtraction(ray.direction, dDOTnTIMESn, ray.direction);

	memcpy(ray.origin, origin, sizeof(GzCoord));
}


void GzRender::Build() {
	this->container->Build();
}

void GetBaryCoordsXZ(float Px, float Pz, GzCoord* vertices, float* out) {
	out[0] = ((vertices[1][Z] - vertices[2][Z]) * (Px - vertices[2][X]) + (vertices[2][X] - vertices[1][X]) * (Pz - vertices[2][Z]))
		/ ((vertices[1][Z] - vertices[2][Z]) * (vertices[0][X] - vertices[2][X]) + (vertices[2][X] - vertices[1][X]) * (vertices[0][Z] - vertices[2][Z]));
	out[1] = ((vertices[2][Z] - vertices[0][Z]) * (Px - vertices[2][X]) + (vertices[0][X] - vertices[2][X]) * (Pz - vertices[2][Z]))
		/ ((vertices[1][Z] - vertices[2][Z]) * (vertices[0][X] - vertices[2][X]) + (vertices[2][X] - vertices[1][X]) * (vertices[0][Z] - vertices[2][Z]));
	out[2] = 1 - out[0] - out[1];
}

bool GzRender::rayTrace() {
	float rad = degreeToRad(this->m_camera.FOV);
	float d = (1 / tan(rad / 2));

	GzCoord cameraDirection, center;
	vectorSubtraction(this->m_camera.lookat, this->m_camera.position, cameraDirection);
	normalizeCoord(cameraDirection);
	scalarMultiplication(cameraDirection, d, cameraDirection);
	vectorAddition(cameraDirection, this->m_camera.position, center);

	GzCoord screen_right, screen_up;

	screen_right[X] = this->m_camera.worldup[Y] * cameraDirection[Z] - this->m_camera.worldup[Z] * cameraDirection[Y];
	screen_right[Y] = this->m_camera.worldup[Z] * cameraDirection[X] - this->m_camera.worldup[X] * cameraDirection[Z];
	screen_right[Z] = this->m_camera.worldup[X] * cameraDirection[Y] - this->m_camera.worldup[Y] * cameraDirection[X];
	normalizeCoord(screen_right);

	screen_up[X] = cameraDirection[Y] * screen_right[Z] - cameraDirection[Z] * screen_right[Y];
	screen_up[Y] = cameraDirection[Z] * screen_right[X] - cameraDirection[X] * screen_right[Z];
	screen_up[Z] = cameraDirection[X] * screen_right[Y] - cameraDirection[Y] * screen_right[X];
	normalizeCoord(screen_up);

	for (int i = 0; i < xres; i++) {
		for (int j = 0; j < yres; j++) {
			std::chrono::time_point<std::chrono::steady_clock> start, end;
			GzRay ray;
			GzTriangle closestTriangle;
			GzCoord intersection, interpNorm;
			GzColor pixelColor;
			zeroColor(pixelColor);
			float remainingIntensity = 1.0;
			// ray.origin = c + (2 * i / xres - 1) * screen_right - (2 * j / yres - 1) * screen_up
			for (int k = 0; k < 3; k++) {
				ray.origin[k] = center[k] + ((((float)2 * i / xres) - 1) * screen_right[k]) - ((((float)2 * j / yres) - 1) * screen_up[k]);
			}
			vectorSubtraction(ray.origin, this->m_camera.position, ray.direction);
			normalizeCoord(ray.direction);

  		start = std::chrono::steady_clock::now();

			bool hitTriangle = this->container->GetNearestIntersectingSurface(ray, closestTriangle, intersection, interpNorm);

			while (hitTriangle && remainingIntensity > 0.001) {
				GzColor color;
				GzCoord baryCoords;
				GzTextureIndex uv;
				memset(uv, 0, sizeof(GzTextureIndex));

				// shading_equation(closestTriangle.normals[0], this, color); // Flat Color
				GzCoord colorNorm;
				memcpy(colorNorm, interpNorm, sizeof(GzCoord));
				if (colorNorm[2] < 0) {
					colorNorm[2] *= -1;
				}

				GetBaryCoordsXZ(intersection[X], intersection[Z], closestTriangle.vertices, baryCoords);
				for (int vertInd = 0; vertInd < 3; vertInd++) {
					uv[0] += baryCoords[vertInd] * closestTriangle.textureCoords[vertInd][0];
					uv[1] += baryCoords[vertInd] * closestTriangle.textureCoords[vertInd][1];
				}
				shading_equation(colorNorm, this, color, uv, closestTriangle.textured); // phong shading

				if (closestTriangle.textured) {
					pixelColor[RED] += color[RED] * remainingIntensity;
					pixelColor[GREEN] += color[GREEN] * remainingIntensity;
					pixelColor[BLUE] += color[BLUE] * remainingIntensity;

					remainingIntensity = 0;
					break;
				}

				pixelColor[RED] += color[RED] * (0.7 * remainingIntensity);
				pixelColor[GREEN] += color[GREEN] * (0.7 * remainingIntensity);
				pixelColor[BLUE] += color[BLUE] * (0.7 * remainingIntensity);

				remainingIntensity *= 0.3;

				reflect(ray, intersection, interpNorm); // same normal as used in Phong? pull functionality out?
				hitTriangle = this->container->GetNearestIntersectingSurface(ray, closestTriangle, intersection, interpNorm);
			}

			pixelColor[RED] += DEF_BG_COLOR[RED] * remainingIntensity;
			pixelColor[GREEN] += DEF_BG_COLOR[GREEN] * remainingIntensity;
			pixelColor[BLUE] += DEF_BG_COLOR[BLUE] * remainingIntensity;
			
			end = std::chrono::steady_clock::now();
			this->GzPutTime(i, j, end - start);

			GzIntensity r, g, b;
			r = this->ctoi(pixelColor[RED]);
			g = this->ctoi(pixelColor[GREEN]);
			b = this->ctoi(pixelColor[BLUE]);
			this->GzPut(i, j, r, g, b, 1, 0);
		}
	}
	return true;
}


