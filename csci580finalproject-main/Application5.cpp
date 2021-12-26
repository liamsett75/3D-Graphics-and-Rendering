// Application5.cpp: implementation of the Application5 class.
//
//////////////////////////////////////////////////////////////////////

/*
 * application test code for homework assignment #5
*/

#include "stdafx.h"
#include "CS580HW.h"
#include "Application5.h"
#include "Gz.h"
#include "rend.h"
#include "NaiveContainer.h"
#include "KDTree.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

#define INFILE  "ppot.asc"
#define INFLOOR "floor.asc"
#define OUTFILE "output.ppm"
#define HEATFILE "heatmap.ppm"
#define SCALEFILE "scale.ppm"
#define QUANTFILE "quant"

extern int tex_fun(float u, float v, GzColor color); /* image texture function */
extern int ptex_fun(float u, float v, GzColor color); /* procedural texture function */
extern int GzFreeTexture();

void shade(GzCoord norm, GzCoord color);

int time2color(double t, GzColor c) {

	double wavelength = t / 2000 + 380; // change divisor of t to scale colors [0, 

	// Linear spline borrowed from https://stackoverflow.com/questions/1472514/convert-light-frequency-to-rgb

	if (wavelength < 440) {
		c[RED] = -(wavelength - 440) / (440 - 380) * MAX_INTENSITY;
		c[GREEN] = 0;
		c[BLUE] = MAX_INTENSITY;
	}
	else if ((wavelength >= 440) && (wavelength < 490)) {
		c[RED] = 0;
		c[GREEN] = (wavelength - 440) / (490 - 440) * MAX_INTENSITY;
		c[BLUE] = MAX_INTENSITY;
	}
	else if ((wavelength >= 490) && (wavelength < 510)) {
		c[RED] = 0;
		c[GREEN] = MAX_INTENSITY;
		c[BLUE] = -(wavelength - 510) / (510 - 490);
	}
	else if ((wavelength >= 510) && (wavelength < 580)) {
		c[RED] = (wavelength - 510) / (580 - 510) * MAX_INTENSITY;
		c[GREEN] = MAX_INTENSITY;
		c[BLUE] = 0;
	}
	else if ((wavelength >= 580) && (wavelength < 645)) {
		c[RED] = MAX_INTENSITY;
		c[GREEN] = -(wavelength - 645) / (645 - 580) * MAX_INTENSITY;
		c[BLUE] = 0;
	}
	else {
		c[RED] = MAX_INTENSITY;
		c[GREEN] = 0;
		c[BLUE] = 0;
	}

	return GZ_SUCCESS;
}

unsigned char intensityToChar2(GzIntensity intensity) {
	intensity = intensity >> 4;
	return (unsigned char)intensity;
}

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Application5::Application5()
{

}

Application5::~Application5()
{
	Clean();
}

int Application5::Initialize()
{
	GzCamera	camera;  
	int		    xRes, yRes;	/* display parameters */ 

	GzToken		nameListShader[9]; 	    /* shader attribute names */
	GzPointer   valueListShader[9];		/* shader attribute pointers */
	GzToken     nameListLights[10];		/* light info */
	GzPointer   valueListLights[10];
	GzToken		nameListContainer[1];
	GzPointer   valueListContainer[1];
	int			shaderType, interpStyle;
	float		specpower;
	int		status; 
 
	status = 0; 

	/* 
	 * Allocate memory for user input
	 */
	m_pUserInput = new GzInput;

	/* 
	 * initialize the display and the renderer 
	 */ 
 	m_nWidth = 256;		// frame buffer and display width
	m_nHeight = 256;    // frame buffer and display height

	m_pRender = new GzRender(m_nWidth, m_nHeight);
	m_pRender->GzDefault();

	m_pFrameBuffer = m_pRender->framebuffer; 

/* Translation matrix */
GzMatrix	scale = 
{ 
	3.25,	0.0,	0.0,	0.0, 
	0.0,	3.25,	0.0,	-3.25, 
	0.0,	0.0,	3.25,	3.5, 
	0.0,	0.0,	0.0,	1.0 
}; 
 
GzMatrix	rotateX = 
{ 
	1.0,	0.0,	0.0,	0.0, 
	0.0,	.7071,	.7071,	0.0, 
	0.0,	-.7071,	.7071,	0.0, 
	0.0,	0.0,	0.0,	1.0 
}; 
 
GzMatrix	rotateY = 
{ 
	.866,	0.0,	-0.5,	0.0, 
	0.0,	1.0,	0.0,	0.0, 
	0.5,	0.0,	.866,	0.0, 
	0.0,	0.0,	0.0,	1.0 
}; 

#if 1 	/* set up app-defined camera if desired, else use camera defaults */
	camera.position[X] = 0.0;
	camera.position[Y] = 5.0;
	camera.position[Z] = -20.0;

	camera.lookat[X] = 5.0;
	camera.lookat[Y] = 0.7;
	camera.lookat[Z] = 6.5;

	camera.worldup[X] = -0.2;
	camera.worldup[Y] = 1.0;
	camera.worldup[Z] = 0.0;

	camera.FOV = 35;           /* degrees *              /* degrees */

	status |= m_pRender->GzPutCamera(camera); 
#endif 

	/* Start Renderer */
	status |= m_pRender->GzBeginRender();

	/* Light */
	GzLight	light1 = { {-0.7071, 0.7071, 0}, {0.5, 0.5, 0.9} };
	GzLight	light2 = { {0, -0.7071, -0.7071}, {0.9, 0.2, 0.3} };
	GzLight	light3 = { {0.7071, 0.0, -0.7071}, {0.2, 0.7, 0.3} };
	GzLight	ambientlight = { {0, 0, 0}, {0.3, 0.3, 0.3} };

	/* Material property */
	GzColor specularCoefficient = { 0.3, 0.3, 0.3 };
	GzColor ambientCoefficient = { 0.1, 0.1, 0.1 };
	GzColor diffuseCoefficient = {0.7, 0.7, 0.7};

/* 
  renderer is ready for frame --- define lights and shader at start of frame 
*/

        /*
         * Tokens associated with light parameters
         */
        nameListLights[0] = GZ_DIRECTIONAL_LIGHT;
        valueListLights[0] = (GzPointer)&light1;
        nameListLights[1] = GZ_DIRECTIONAL_LIGHT;
        valueListLights[1] = (GzPointer)&light2;
        nameListLights[2] = GZ_DIRECTIONAL_LIGHT;
        valueListLights[2] = (GzPointer)&light3;
        status |= m_pRender->GzPutAttribute(3, nameListLights, valueListLights);

        nameListLights[0] = GZ_AMBIENT_LIGHT;
        valueListLights[0] = (GzPointer)&ambientlight;
        status |= m_pRender->GzPutAttribute(1, nameListLights, valueListLights);

        /*
         * Tokens associated with shading 
         */
        nameListShader[0]  = GZ_DIFFUSE_COEFFICIENT;
        valueListShader[0] = (GzPointer)diffuseCoefficient;

	/* 
	* Select either GZ_COLOR or GZ_NORMALS as interpolation mode  
	*/
        nameListShader[1]  = GZ_INTERPOLATE;
        interpStyle = GZ_COLOR;         /* Gouraud shading */
        //interpStyle = GZ_NORMALS;         /* Phong shading */
        valueListShader[1] = (GzPointer)&interpStyle;

        nameListShader[2]  = GZ_AMBIENT_COEFFICIENT;
        valueListShader[2] = (GzPointer)ambientCoefficient;
        nameListShader[3]  = GZ_SPECULAR_COEFFICIENT;
        valueListShader[3] = (GzPointer)specularCoefficient;
        nameListShader[4]  = GZ_DISTRIBUTION_COEFFICIENT;
        specpower = 32;
        valueListShader[4] = (GzPointer)&specpower;

        nameListShader[5]  = GZ_TEXTURE_MAP;
        valueListShader[5] = (GzPointer)(ptex_fun);

        status |= m_pRender->GzPutAttribute(6, nameListShader, valueListShader);

		AContainer *myContainer;

#if 0
		myContainer = new NaiveContainer();
#else
		myContainer = new KDTree();
#endif
		
		nameListContainer[0] = GZ_CONTAINER_TYPE;
		valueListContainer[0] = (GzPointer)myContainer;

		status |= m_pRender->GzPutAttribute(1, nameListContainer, valueListContainer);


	status |= m_pRender->GzPushMatrix(scale);  
	status |= m_pRender->GzPushMatrix(rotateY); 
	status |= m_pRender->GzPushMatrix(rotateX); 

	if (status) exit(GZ_FAILURE); 

	if (status) 
		return(GZ_FAILURE); 
	else 
		return(GZ_SUCCESS); 
}

int Application5::Render() 
{
	GzToken		nameListTriangle[3]; 	/* vertex attribute names */
	GzPointer	valueListTriangle[3]; 	/* vertex attribute pointers */
	GzCoord		vertexList[3];	/* vertex position coordinates */ 
	GzCoord		normalList[3];	/* vertex normals */ 
	GzTextureIndex  	uvList[3];		/* vertex texture map indices */
	GzPointer	valueListTriangle2[3]; 	/* vertex attribute pointers */
	GzCoord		vertexList2[3];	/* vertex position coordinates */
	GzCoord		normalList2[3];	/* vertex normals */
	GzTextureIndex  	uvList2[3];		/* vertex texture map indices */
	char		dummy[256]; 
	int			status;



	/* Initialize Display */
	status |= m_pRender->GzDefault();  /* init for new frame */
	
	/* 
	* Tokens associated with triangle vertex values 
	*/ 
	nameListTriangle[0] = GZ_POSITION; 
	nameListTriangle[1] = GZ_NORMAL; 
	nameListTriangle[2] = GZ_TEXTURE_INDEX;

	// I/O File open
	FILE* infile;
	if ((infile = fopen(INFILE, "r")) == NULL)
	{
		AfxMessageBox("The input file was not opened\n");
		return GZ_FAILURE;
	}
	FILE* infloor;
	if ((infloor = fopen(INFLOOR, "r")) == NULL)
	{
		AfxMessageBox("The input floor file was not opened\n");
		return GZ_FAILURE;
	}

	FILE *outfile;
	if ((outfile = fopen(OUTFILE, "wb")) == NULL)
	{
		AfxMessageBox("The output file was not opened\n");
		return GZ_FAILURE;
	}

	FILE* heatmapfile;
	if ((heatmapfile = fopen(HEATFILE, "wb")) == NULL)
	{
		AfxMessageBox("The heatmap file was not opened\n");
		return GZ_FAILURE;
	}

	FILE* quantfile;
	if ((quantfile = fopen(QUANTFILE, "wb")) == NULL) {
		AfxMessageBox("The quant file was not opened\n");
		return GZ_FAILURE;
	}

	FILE* scalefile;
	if ((scalefile = fopen(SCALEFILE, "wb")) == NULL)
	{
		AfxMessageBox("The scale file was not opened\n");
		return GZ_FAILURE;
	}

	m_pRender->texturing = false;

	/* 
	* Walk through the list of triangles, set color 
	* and render each triangle 
	*/ 
	while( fscanf(infile, "%s", dummy) == 1) { 	/* read in tri word */
	    fscanf(infile, "%f %f %f %f %f %f %f %f", 
		&(vertexList[0][0]), &(vertexList[0][1]),  
		&(vertexList[0][2]), 
		&(normalList[0][0]), &(normalList[0][1]), 	
		&(normalList[0][2]), 
		&(uvList[0][0]), &(uvList[0][1]) ); 
	    fscanf(infile, "%f %f %f %f %f %f %f %f", 
		&(vertexList[1][0]), &(vertexList[1][1]), 	
		&(vertexList[1][2]), 
		&(normalList[1][0]), &(normalList[1][1]), 	
		&(normalList[1][2]), 
		&(uvList[1][0]), &(uvList[1][1]) ); 
	    fscanf(infile, "%f %f %f %f %f %f %f %f", 
		&(vertexList[2][0]), &(vertexList[2][1]), 	
		&(vertexList[2][2]), 
		&(normalList[2][0]), &(normalList[2][1]), 	
		&(normalList[2][2]), 
		&(uvList[2][0]), &(uvList[2][1]) ); 

	    /* 
	     * Set the value pointers to the first vertex of the 	
	     * triangle, then feed it to the renderer 
	     * NOTE: this sequence matches the nameList token sequence
	     */ 
	     valueListTriangle[0] = (GzPointer)vertexList; 
		 valueListTriangle[1] = (GzPointer)normalList; 
		 valueListTriangle[2] = (GzPointer)uvList;
		 m_pRender->GzPutTriangle(3, nameListTriangle, valueListTriangle);

		 for (int i = 0; i < 3; i++)
		 {
			 for (int j = 0; j < 3; j++)
			 {
				 vertexList2[i][j] = vertexList[i][j];
			 }
		 }
		 vertexList2[0][2] = vertexList[0][2] - 4.2;
		 vertexList2[1][2] = vertexList[1][2] - 4.2;
		 vertexList2[2][2] = vertexList[2][2] - 4.2;

		 valueListTriangle2[0] = (GzPointer)vertexList2;
		 valueListTriangle2[1] = (GzPointer)normalList;
		 valueListTriangle2[2] = (GzPointer)uvList;
		 m_pRender->GzPutTriangle(3, nameListTriangle, valueListTriangle2);
	}

	m_pRender->texturing = true;
	while (fscanf(infloor, "%s", dummy) == 1) { 	/* read in tri word */
		fscanf(infloor, "%f %f %f %f %f %f %f %f",
			&(vertexList[0][0]), &(vertexList[0][1]),
			&(vertexList[0][2]),
			&(normalList[0][0]), &(normalList[0][1]),
			&(normalList[0][2]),
			&(uvList[0][0]), &(uvList[0][1]));
		fscanf(infloor, "%f %f %f %f %f %f %f %f",
			&(vertexList[1][0]), &(vertexList[1][1]),
			&(vertexList[1][2]),
			&(normalList[1][0]), &(normalList[1][1]),
			&(normalList[1][2]),
			&(uvList[1][0]), &(uvList[1][1]));
		fscanf(infloor, "%f %f %f %f %f %f %f %f",
			&(vertexList[2][0]), &(vertexList[2][1]),
			&(vertexList[2][2]),
			&(normalList[2][0]), &(normalList[2][1]),
			&(normalList[2][2]),
			&(uvList[2][0]), &(uvList[2][1]));

		/*
		 * Set the value pointers to the first vertex of the
		 * triangle, then feed it to the renderer
		 * NOTE: this sequence matches the nameList token sequence
		 */
		valueListTriangle[0] = (GzPointer)vertexList;
		valueListTriangle[1] = (GzPointer)normalList;
		valueListTriangle[2] = (GzPointer)uvList;
		m_pRender->GzPutTriangle(3, nameListTriangle, valueListTriangle);
	}

	m_pRender->Build();

	float *avgHeatBuffer = new float[this->m_nWidth * this->m_nHeight];
	memset(avgHeatBuffer, 0, this->m_nWidth * this->m_nHeight * sizeof(float));

#if 1
	const int NUM_TRIALS = 1; // Took about 20 seconds for Adam
#else
	const int NUM_TRIALS = 20; // Took about 6 minutes for Adam
#endif

	for (int n = 0; n < NUM_TRIALS; n++) {
		m_pRender->rayTrace();
		for (int i = 0; i < this->m_nWidth; i++) {
			for (int j = 0; j < this->m_nHeight; j++) {
				avgHeatBuffer[i + j * this->m_nWidth] += m_pRender->timebuffer[i + j * this->m_nWidth];
			}
		}
	}

	for (int i = 0; i < this->m_nWidth; i++) {
		for (int j = 0; j < this->m_nHeight; j++) {
			avgHeatBuffer[i + j * this->m_nWidth] /= NUM_TRIALS;
		}
	}

	// Not the best engineering, but the class's gz code base is not great engineering in the first place...
	char* header_buf = new char[2000];
	sprintf(header_buf, "P6 %d %d 255\r", this->m_nWidth, this->m_nHeight);
	fputs(header_buf, heatmapfile);
	fputs(header_buf, scalefile);
	delete[] header_buf;
	GzColor c;

	int buffersize = (3 * this->m_nWidth * this->m_nHeight);
	float totalpixeltime = 0;
	char* output = new char[buffersize + 1];
	char* scale = new char[buffersize + 1];
	for (int i = 0; i < buffersize; i += 3) {
		int pixel_num = i / 3;
		if (i + 2 > buffersize) {
			printf("Error: overflowing flush to ppm file");
			return GZ_FAILURE;
		}

		// 0 to 530000
		time2color(avgHeatBuffer[pixel_num], c);
		totalpixeltime += avgHeatBuffer[pixel_num] / 1000000;
		output[i] = intensityToChar2(c[RED]);
		output[i + 1] = intensityToChar2(c[GREEN]);
		output[i + 2] = intensityToChar2(c[BLUE]);

		time2color(530000 - 530000 * (pixel_num % this->m_nWidth) / this->m_nHeight, c);
		scale[i] = intensityToChar2(c[RED]);
		scale[i + 1] = intensityToChar2(c[GREEN]);
		scale[i + 2] = intensityToChar2(c[BLUE]);
	}
	output[buffersize] = '\0';
	scale[buffersize] = '\0';
	fwrite(output, sizeof(unsigned char), buffersize, heatmapfile);
	fwrite(scale, sizeof(unsigned char), buffersize, scalefile);
	fprintf(quantfile, "Average pixel render time was: %fms", totalpixeltime / this->m_nWidth / this->m_nHeight);

	m_pRender->GzFlushDisplay2File(outfile); 	/* write out or update display to file*/
	m_pRender->GzFlushDisplay2FrameBuffer();	// write out or update display to frame buffer

	delete[] avgHeatBuffer;

	/* 
	 * Close file
	 */ 

	if( fclose( infile ) )
      AfxMessageBox(_T( "The input file was not closed\n" ));

	if( fclose( outfile ) )
      AfxMessageBox(_T( "The output file was not closed\n" ));

	if (fclose(quantfile))
		AfxMessageBox(_T("The quant file was not closed\n"));

	if (fclose(heatmapfile))
		AfxMessageBox(_T("The quant file was not closed\n"));

	if (fclose(scalefile))
		AfxMessageBox(_T("The quant file was not closed\n"));
 
	if (status) 
		return(GZ_FAILURE); 
	else 
		return(GZ_SUCCESS); 
}

int Application5::Clean()
{
	/* 
	 * Clean up and exit 
	 */ 
	int	status = 0; 

	delete(m_pRender);
	status |= GzFreeTexture();
	
	if (status) 
		return(GZ_FAILURE); 
	else 
		return(GZ_SUCCESS);
}



