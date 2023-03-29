/*
 *	COMP559 Winter 2023 Final Project
 *	Submission by Edwin Pan 260843175
 *
 *	Simulator code reference is Matthias M�ller's online "Ten Minute Physics" resource, particularly their video on "How to write a FLIP water / fluid simulation" - youtube https://www.youtube.com/watch?v=XmzBREkK8kY and github https://github.com/matthias-research/pages/blob/master/tenMinutePhysics/18-flip.html .
 *	Instead of reverse engineering opengl code provided for assignments I figured I'd just take an opengl tutorial starting here https://www.youtube.com/watch?v=OR4fNpBjmq8&list=PLlrATfBNZ98foTJPJ_Ev03o2oq3-GGOS2&index=2
 *
 */


#include "includes.h"
#include "OpenGlHelpers.h"
#include "Scene.h"
#include "FlipFluid.h"



 /*
 * ========================================================================================================================
 *	Function Declarations
 * ========================================================================================================================
 */

#define ASSERT(x) if (!(x)) __debugbreak();
void setupSceneAndFlipFluid();
void setObstacle(float x, float y, bool reset);
void startDrag(float x, float y);
void drag(float x, float y);
void endDrag();

void toggleStart();
void simulate();
void update();

void draw();



 /*
 * ========================================================================================================================
 *	Main Constants and Variables
 * ========================================================================================================================
 */

int width = 1920;
int height = 1080;
float simHeight = 3.0f;
float cScale = height / simHeight;
float simWidth = width / cScale;

bool mouseDown = false;

GLFWwindow* window;
Scene scene;
FlipFluid fluid;

unsigned int pointShader = -1;
unsigned int meshShader = -1;
unsigned int pointVertexBuffer = -1;
unsigned int pointColorBuffer = -1;
unsigned int gridVertexBuffer = -1;
unsigned int gridColorBuffer = -1;
unsigned int diskVertBuffer = -1;
unsigned int diskIdBuffer = -1;




/*
* ========================================================================================================================
*	Main Function
* ========================================================================================================================
*/
int main() {


	/* Initialize GLFW */
	if (!glfwInit()) {return -1;}

	/* Create a windowed mode window and its OpenGL context */
	window = glfwCreateWindow(1024, 768, "COMP559 Winter 2023 Final Project 260843175", NULL, NULL);
	if (!window) {
		glfwTerminate();
		return -1;
	}
	// Mark current flgw context and initialize glew
	glfwMakeContextCurrent(window);
	if (glewInit() != GLEW_OK) {
		std::cout << "Error! Failed to initialize glew!" << std::endl;
		return -2;
	}
	std::cout << glGetString(GL_VERSION) << std::endl;

	// Rendering loop
	while (!glfwWindowShouldClose(window)) {
		update();
	}

	// Clean up
	//glDeleteProgram(shader);

	glfwTerminate();
	return 0;

}



/*
* ========================================================================================================================
*	Function Definitions
* ========================================================================================================================
*/

void setupSceneAndFlipFluid() {
	scene.obstacleRadius = 0.15;
	scene.overRelaxation = 1.9;

	scene.dt = 1.0 / 60.0;
	scene.numPressureIters = 50;
	scene.numParticleIters = 2;

	int res = 100;

	float tankHeight = 1.0 * simHeight;
	float tankWidth = 1.0 * simWidth;
	float h = tankHeight / res;
	float density = 1000.0;

	float relWaterHeight = 0.8;
	float relWaterWidth = 0.6;

	// dam break

	// compute number of particles

	float r = 0.3 * h;	// particle radius w.r.t. cell size
	float dx = 2.0 * r;
	float dy = sqrt(3.0) / 2.0 * dx;

	int numX = floor((relWaterWidth * tankWidth - 2.0 * h - 2.0 * r) / dx);
	int numY = floor((relWaterHeight * tankHeight - 2.0 * h - 2.0 * r) / dy);
	int maxParticles = numX * numY;

	// create fluid

	fluid = FlipFluid(density, tankWidth, tankHeight, h, r, maxParticles, &scene );

	// create particles

	fluid.numParticles = numX * numY;
	int p = 0;
	for (int i = 0; i < numX; i++) {
		for (int j = 0; j < numY; j++) {
			fluid.particlePos[p++] = h + r + dx * i + (j % 2 == 0 ? 0.0 : r);
			fluid.particlePos[p++] = h + r + dy * j;
		}
	}

	// setup grid cells for tank

	int n = fluid.fNumY;

	for (int i = 0; i < fluid.fNumX; i++) {
		for (int j = 0; j < fluid.fNumY; j++) {
			float s = 1.0;	// fluid
			if (i == 0 || i == fluid.fNumX - 1 || j == 0)
				s = 0.0;	// solid
			fluid.s[i * n + j] = s;
		}
	}

	setObstacle(3.0, 2.0, true);
}

void setObstacle(float x, float y, bool reset) {
	float vx = 0.0;
	float vy = 0.0;

	if (!reset) {
		vx = (x - scene.obstacleX) / scene.dt;
		vy = (y - scene.obstacleY) / scene.dt;
	}

	scene.obstacleX = x;
	scene.obstacleY = y;
	float r = scene.obstacleRadius;
	int n = fluid.fNumY;
	float cd = sqrt(2) * fluid.h;

	for (int i = 1; i < fluid.fNumX - 2; i++) {
		for (int j = 1; j < fluid.fNumY - 2; j++) {

			fluid.s[i * n + j] = 1.0;

			float dx = (i + 0.5) * fluid.h - x;
			float dy = (j + 0.5) * fluid.h - y;

			if (dx * dx + dy * dy < r * r) {
				fluid.s[i * n + j] = 0.0;
				fluid.u[i * n + j] = vx;
				fluid.u[(i + 1) * n + j] = vx;
				fluid.v[i * n + j] = vy;
				fluid.v[i * n + j + 1] = vy;
			}
		}
	}

	scene.showObstacle = true;
	scene.obstacleVelX = vx;
	scene.obstacleVelY = vy;
}

void startDrag(float x, float y) {

}

void drag(float x, float y) {

}

void endDrag() {

}

void toggleStart() {
	scene.paused = !scene.paused;
}

void simulate() {
	if (!scene.paused) {
		fluid.simulate(
			scene.dt, scene.gravity, scene.flipRatio, scene.numPressureIters, scene.numParticleIters,
			scene.overRelaxation, scene.compensateDrift, scene.separateParticles,
			scene.obstacleX, scene.obstacleY, scene.obstacleRadius);
		scene.frameNr++;
	}
}

void update() {
	simulate();
	draw();
}

void draw() {
	glClear(GL_COLOR_BUFFER_BIT); // Clear
	// DRAW CODE STARTS HERE

	//Prepare Shaders
	if (pointShader == -1) {
		ShaderProgramSource pointSources = ParseShader("res/shaders/point.shader");
		pointShader = CreateShader(pointSources.VertexSource, pointSources.FragmentSource);
	}
	if (meshShader == -1) {
		ShaderProgramSource meshSources = ParseShader("res/shaders/mesh.shader");
		meshShader = CreateShader(meshSources.VertexSource, meshSources.FragmentSource);
	}

	// Grid
	if (gridVertexBuffer == -1) {

		float* cellCenters = new float[2 * fluid.fNumCells];
		int p = 0;
		for (int i = 0; i < fluid.fNumX; i++) {
			for (int j = 0; j < fluid.fNumY; j++) {
				cellCenters[p++] = (i + 0.5) * fluid.h;
				cellCenters[p++] = (j + 0.5) * fluid.h;
			}
		}

		GLCall(glGenBuffers(1, &gridVertexBuffer));
		GLCall(glBindBuffer(GL_ARRAY_BUFFER, gridVertexBuffer));
		GLCall(glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 2 * fluid.fNumCells, cellCenters, GL_DYNAMIC_DRAW));
		GLCall(glBindBuffer(GL_ARRAY_BUFFER, 0)); //Assume '0' equivalent to null

	}
	if (gridColorBuffer == -1) {
		GLCall(glGenBuffers(1, &gridColorBuffer));
	}
	if (scene.showGrid) {
		float pointSize = 0.9 * fluid.h / simWidth * width;

		GLCall(glUseProgram(pointShader));
		GLCall(glUniform2f(glGetUniformLocation(pointShader, "domainSize"), simWidth, simHeight));
		GLCall(glUniform1f(glGetUniformLocation(pointShader, "pointSize"), pointSize));
		GLCall(glUniform1f(glGetUniformLocation(pointShader, "drawDisk"), 0.0f));

		GLCall(glBindBuffer(GL_ARRAY_BUFFER, gridVertexBuffer));
		GLCall(unsigned int posLoc = glGetAttribLocation(pointShader, "attrPosition"));
		GLCall(glEnableVertexAttribArray(posLoc));
		GLCall(glVertexAttribPointer(posLoc, 2, GL_FLOAT, false, 0, 0));

		GLCall(glBindBuffer(GL_ARRAY_BUFFER, gridColorBuffer));
		GLCall(glBufferData(GL_ARRAY_BUFFER, sizeof(float) * fluid.cellColor.size(), &fluid.cellColor.front(), GL_DYNAMIC_DRAW));

		GLCall(unsigned int colorLoc = glGetAttribLocation(pointShader, "attrColor"));
		GLCall(glEnableVertexAttribArray(colorLoc));
		GLCall(glVertexAttribPointer(colorLoc, 3, GL_FLOAT, false, 0, 0 ));

		GLCall(glDrawArrays(GL_POINTS, 0, fluid.fNumCells));

		GLCall(glDisableVertexAttribArray(posLoc));
		GLCall(glDisableVertexAttribArray(colorLoc));

		GLCall(glBindBuffer(GL_ARRAY_BUFFER, 0));
	}

	// Water
	if (scene.showParticles) {
		GLCall(glClear(GL_DEPTH_BUFFER_BIT));

		double pointSize = 2.0 * fluid.particleRadius / simWidth * width;

		GLCall(glUseProgram(pointShader));
		GLCall(glUniform2f(glGetUniformLocation(pointShader, "domainSize"), simWidth, simHeight));
		GLCall(glUniform1f(glGetUniformLocation(pointShader, "pointSize"), pointSize));
		GLCall(glUniform1f(glGetUniformLocation(pointShader, "drawDisk"), 1.0));

		if (pointVertexBuffer == -1) {
			GLCall(glGenBuffers(1, &pointVertexBuffer));
		}
		if (pointColorBuffer == -1) {
			GLCall(glGenBuffers(1, &pointColorBuffer));
		}

		GLCall(glBindBuffer(GL_ARRAY_BUFFER, pointVertexBuffer));
		GLCall(glBufferData(GL_ARRAY_BUFFER, sizeof(float) * fluid.particlePos.size(), &fluid.particlePos.front(), GL_DYNAMIC_DRAW));
		
		GLCall(unsigned int posLoc = glGetAttribLocation(pointShader, "attrPosition"));
		GLCall(glEnableVertexAttribArray(posLoc));
		GLCall(glVertexAttribPointer(posLoc, 2, GL_FLOAT, false, 0, 0));

		GLCall(glBindBuffer(GL_ARRAY_BUFFER, pointColorBuffer));
		GLCall(glBufferData(GL_ARRAY_BUFFER, sizeof(float) * fluid.particleColor.size(), &fluid.particleColor.front(), GL_DYNAMIC_DRAW));

		GLCall(unsigned int colorLoc = glGetAttribLocation(pointShader, "attrColor"));
		GLCall(glEnableVertexAttribArray(colorLoc));
		GLCall(glVertexAttribPointer(colorLoc, 3, GL_FLOAT, false, 0, 0));

		GLCall(glDrawArrays(GL_POINTS, 0, fluid.numParticles));

		GLCall(glDisableVertexAttribArray(posLoc));
		GLCall(glDisableVertexAttribArray(colorLoc));

		GLCall(glBindBuffer(GL_ARRAY_BUFFER, 0));
	}

	// Disk
	
	int numSegs = 50;
	if (diskVertBuffer == -1) {
		GLCall(glGenBuffers(1, &diskVertBuffer));
		double dphi = 2.0 * 3.1415926535898 / numSegs;
		float diskVerts[102];
		int p = 0;
		diskVerts[p++] = 0.0;
		diskVerts[p++] = 0.0;
		for (int i = 0; i < numSegs; i++) {
			diskVerts[p++] = cos(i * dphi);
			diskVerts[p++] = sin(i * dphi);
		}
		GLCall(glBindBuffer(GL_ARRAY_BUFFER, diskVertBuffer));
		GLCall(glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 102, diskVerts, GL_DYNAMIC_DRAW));
		GLCall(glBindBuffer(GL_ARRAY_BUFFER, 0));

		GLCall(glGenBuffers(1, &diskIdBuffer));
		unsigned short diskIds[150];
		p = 0;
		for (int i = 0; i < numSegs; i++) {
			diskIds[p++] = 0;
			diskIds[p++] = 1 + i;
			diskIds[p++] = 1 + (i + 1) % numSegs;
		}
		GLCall(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, diskIdBuffer));
		GLCall(glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned short) * 150, diskIds, GL_DYNAMIC_DRAW));
		GLCall(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0));
	}

	GLCall(glClear(GL_DEPTH_BUFFER_BIT));

	float diskColor[3] = { 1.0f, 0.0f, 0.0f };

	GLCall(glUseProgram(meshShader));
	GLCall(glUniform2f(glGetUniformLocation(meshShader, "domainSize"), simWidth, simHeight));
	GLCall(glUniform3f(glGetUniformLocation(meshShader, "color"), diskColor[0], diskColor[1], diskColor[2]));
	GLCall(glUniform2f(glGetUniformLocation(meshShader, "translation"), scene.obstacleX, scene.obstacleY));
	GLCall(glUniform1f(glGetUniformLocation(meshShader, "scale"), scene.obstacleRadius + fluid.particleRadius));

	GLCall(unsigned int posLoc = glGetAttribLocation(meshShader, "attrPosition"));
	GLCall(glEnableVertexAttribArray(posLoc));
	GLCall(glBindBuffer(GL_ARRAY_BUFFER, diskVertBuffer));
	GLCall(glVertexAttribPointer(posLoc, 2, GL_FLOAT, false, 0, 0));
	
	GLCall(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, diskIdBuffer));
	GLCall(glDrawElements(GL_TRIANGLES, 3 * numSegs, GL_UNSIGNED_SHORT, 0));

	GLCall(glDisableVertexAttribArray(posLoc));



	// DRAW CODE ENDS HERE
	glfwSwapBuffers(window); // Buffer Swap
	glfwPollEvents(); // Poll for Events
}
