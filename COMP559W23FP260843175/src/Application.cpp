/*
 *	COMP559 Winter 2023 Final Project
 *	Submission by Edwin Pan 260843175
 *
 *	Simulator code reference is Matthias M�ller's online "Ten Minute Physics" resource, particularly their video on "How to write a FLIP water / fluid simulation" - youtube https://www.youtube.com/watch?v=XmzBREkK8kY and github https://github.com/matthias-research/pages/blob/master/tenMinutePhysics/18-flip.html .
 *	Instead of reverse engineering opengl code provided for assignments I figured I'd just take an opengl tutorial starting here https://www.youtube.com/watch?v=OR4fNpBjmq8&list=PLlrATfBNZ98foTJPJ_Ev03o2oq3-GGOS2&index=2
 *	Used https://www.glfw.org/docs/3.3/input_guide.html as reference for how to get user inputs.
 *
 */


#include "imgui/imgui.h"
#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl3.h"
#include "includes.h"
#include "OpenGlHelpers.h"
#include "Scene.h"
#include "FlipFluid.cpp"



 /*
 * ========================================================================================================================
 *	Function Declarations
 * ========================================================================================================================
 */

// Initialization-Related
void finishSceneFluidSetup();
// Simulation Related
void setObstacle(float x, float y, bool reset);
void toggleStart();
void simulate();
// User-Input Related
void mousePositionCallback(GLFWwindow* window, double x, double y);
void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods);
void keyboardKeyCallback(GLFWwindow* window, int key, int scancode, int action, int mods);
// Graphics Related
void draw();
void drawFluids();
void drawUI();



/*
* ========================================================================================================================
*	Graphics Variables
* ========================================================================================================================
*/

// Shaders must be initialized after GLEW is initialized
unsigned int pointShader = -1;
unsigned int meshShader = -1;
unsigned int dotShader = -1; // My attempt at recreating Matthia's Muller's point shader which doesn't do shit here.
unsigned int pointVertexBuffer = -1;
unsigned int pointColorBuffer = -1;
unsigned int gridVertexBuffer = -1;
unsigned int gridColorBuffer = -1;
unsigned int diskVertBuffer = -1;
unsigned int diskIdBuffer = -1;
// TEMPORARY EXPERIMENTAL STUFF
unsigned int tempVertexBuffer = -1;
unsigned int tempColourBuffer = -1;



 /*
 * ========================================================================================================================
 *	Main Constants and Variables
 * ========================================================================================================================
 */

// Window and 2D Scene Properties
GLFWwindow* window;
Scene scene;
int width = 1920;
int height = 1080;
float simHeight = 3.0f;
float cScale = height / simHeight;
float simWidth = width / cScale;

// UI Variables
bool mouseDown = false;
bool stepForward = false;

// Statistics Variables
std::chrono::high_resolution_clock::time_point prevTick;
std::chrono::high_resolution_clock::time_point currTick;
bool prevTickWasPaused = true;

// Sim Physical Properties SET UP
int res = 64;
float tankHeight = 1.0 * simHeight;
float tankWidth = 1.0 * simWidth;
float h = tankHeight / res;
float density = 1000.0;
float relWaterHeight = 0.8;
float relWaterWidth = 0.6;
// Particle Properties SET UP
float r = 0.3 * h;	// particle radius w.r.t. cell size
float dx = 2.0 * r;
float dy = sqrt(3.0) / 2.0 * dx;
int numX = floor((relWaterWidth * tankWidth - 2.0 * h - 2.0 * r) / dx);
int numY = floor((relWaterHeight * tankHeight - 2.0 * h - 2.0 * r) / dy);
int maxParticles = numX * numY;
// Fluid-Encompassing Object Instantiation SET UP
FlipFluid fluid = FlipFluid(density, tankWidth, tankHeight, h, r, maxParticles, &scene);




/*
* ========================================================================================================================
*	Main Function
* ========================================================================================================================
*/
int main() {
	// Finish Setting Up Scene and Fluid objects
	finishSceneFluidSetup();

	// Reset Statistics Stuff
	prevTick = std::chrono::high_resolution_clock::now();
	currTick = std::chrono::high_resolution_clock::now();
	prevTickWasPaused = true;

	// Initialize Graphics Components
	if (!glfwInit()) {return -1;}
	// GLFW
	window = glfwCreateWindow(width, height, "COMP559 Winter 2023 Final Project 260843175", NULL, NULL);
	if (!window) {
		glfwTerminate();
		return -1;
	}
	// GLEW
	glfwMakeContextCurrent(window);
	if (glewInit() != GLEW_OK) {
		std::cout << "Error! Failed to initialize glew!" << std::endl;
		return -2;
	}
	std::cout << glGetString(GL_VERSION) << std::endl;

	// Assign User Input Behaviours
	glfwSetKeyCallback(window, keyboardKeyCallback);
	glfwSetCursorPosCallback(window, mousePositionCallback);
	glfwSetMouseButtonCallback(window, mouseButtonCallback);

	// Initialize ImGUI based on https://www.youtube.com/watch?v=VRwhNKoxUtk
	// Annoyingly, this MUST be after glfwSetXCallbacks are called, because
	// otherwise ImGUI will ignore all user inputs.
	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImGui::StyleColorsDark();
	ImGui_ImplGlfw_InitForOpenGL(window, true);
	ImGui_ImplOpenGL3_Init("#version 400");

	// Initialize Shaders
	pointShader = ParseAndCreateShader("res/shaders/point.shader");
	meshShader = ParseAndCreateShader("res/shaders/mesh.shader");

	// Rendering loop
	while (!glfwWindowShouldClose(window)) {
		simulate();
		draw();
	}

	// Clean up
	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();
	glDeleteProgram(pointShader);
	glDeleteProgram(meshShader);
	glfwTerminate();
	return 0;

}



/*
* ========================================================================================================================
*	Function Definitions
* ========================================================================================================================
*/

/// <summary>
/// Wraps up Scene and Fluid object setup.
/// Based on setupScene() in Matthia Muller's code, but that method
/// had to be split up in this C++ re-implementation C++ wants to
/// immediately initialize the FlipFluid object the moment it is called,
/// which caused compile-time errors in C++. Section of the original
/// setupScene() have thus been split up to before the main function
/// call among variable declarations and this method.
/// </summary>
void finishSceneFluidSetup() {
	// Scene Properties
	scene.obstacleRadius = 0.15;
	scene.overRelaxation = 1.9;
	scene.dt = 1.0 / 60.0;
	scene.numPressureIters = 50;
	scene.numParticleIters = 2;
	// Particle Initial Conditions
	fluid.numParticles = numX * numY;
	int p = 0;
	for (int i = 0; i < numX; i++) {
		for (int j = 0; j < numY; j++) {
			fluid.particlePos[p++] = h + r + dx * i + (j % 2 == 0 ? 0.0 : r);
			fluid.particlePos[p++] = h + r + dy * j;
		}
	}
	// Fluid Grid Setup
	int n = fluid.fNumY;
	for (int i = 0; i < fluid.fNumX; i++) {
		for (int j = 0; j < fluid.fNumY; j++) {
			float s = 1.0;	// fluid
			if (i == 0 || i == fluid.fNumX - 1 || j == 0)
				s = 0.0;	// solid: side of the tank
			fluid.s[i * n + j] = s;
		}
	}
	// Set Obstacle Initial Position
	setObstacle(3.0, 2.0, true);
}

/// <summary>
/// Updates the position of the obstacle to position x and y.
/// Updates the velocity of the obstacle unless reset is true.
/// If reset is true, then velocity of the obstacle is zeroed.
/// </summary>
/// <param name="x"></param>
/// <param name="y"></param>
/// <param name="reset"></param>
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

/// <summary>
/// Pauses or Resumes the scene.
/// </summary>
void toggleStart() {
	scene.paused = !scene.paused;
}

/// <summary>
/// Calls for the calculation of the next frame of the simulation.
/// </summary>
void simulate() {
	if (!scene.paused || stepForward) {
		fluid.simulate(
			scene.dt, scene.gravity, scene.flipRatio, scene.numPressureIters, scene.numParticleIters,
			scene.overRelaxation, scene.compensateDrift, scene.separateParticles,
			scene.obstacleX, scene.obstacleY, scene.obstacleRadius);
		scene.frameNr++;
		stepForward = false; // If keyboard asks for stepforward, step forward once.
	}
}

/// <summary>
/// Mouse Position Update callback. Updates the position of the obstacle
/// if mouse is currently pressed.
/// Based on the `function drag(x, y)` code in Matthia's code.
/// Adapted for C++ GLFW via https://www.glfw.org/docs/3.3/input_guide.html
/// </summary>
/// <param name="window"></param>
/// <param name="x"></param>
/// <param name="y"></param>
void mousePositionCallback(GLFWwindow* window, double x, double y) {
	if (mouseDown) {
		setObstacle(x / cScale, (height - y) / cScale, false);
	}
}

/// <summary>
/// Mouse Button Input callback. Starts or finishes the dragging of the obstacle
/// Based on `function startDrag(x, y)` and `function endDrag()` in Matthia's code.
/// Adapted for C++ GLFW via https://www.glfw.org/docs/3.3/input_guide.html
/// Also this stackoverflow page helps with getting xpos ypos since our
/// signature does not itself give this info: https://stackoverflow.com/questions/45130391/opengl-get-cursor-coordinate-on-mouse-click-in-c
/// Exclusion of ImGUI clicks based on code provided to us in assignment 3.
/// </summary>
/// <param name="window"></param>
/// <param name="x"></param>
/// <param name="y"></param>
void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods) {
	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS && !ImGui::GetIO().WantCaptureMouse) {
		// If left press
		double x, y;
		glfwGetCursorPos(window, &x, &y);
		y = height - y;
		setObstacle(x/cScale, y/cScale, true);
		mouseDown = true;

	}
	else if(button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {
		// If left release
		scene.obstacleVelX = 0.0f;
		scene.obstacleVelY = 0.0f;
		mouseDown = false;

	}
}

/// <summary>
/// Keyboard Key Input callback. Handles keyboard inputs.
/// </summary>
/// <param name="window"></param>
/// <param name="key"></param>
/// <param name="scancode"></param>
/// <param name="action"></param>
/// <param name="mods"></param>
void keyboardKeyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
	if (key == GLFW_KEY_SPACE && action == GLFW_PRESS) {
		toggleStart(); // Pause/Play
	}
	else if (key == GLFW_KEY_S && action == GLFW_PRESS) {
		stepForward = true; // Step
	}
	else if (key == GLFW_KEY_R && action == GLFW_PRESS) {
		finishSceneFluidSetup(); // Restart
		scene.frameNr = 0;
	}
}

/// <summary>
/// Handles the generation of a new frame.
/// </summary>
void draw() {
	glClear(GL_COLOR_BUFFER_BIT); // Clear
	ImGui_ImplOpenGL3_NewFrame();
	ImGui_ImplGlfw_NewFrame();
	ImGui::NewFrame();
	// DRAW CODE STARTS HERE


	drawFluids();
	drawUI();


	// DRAW CODE ENDS HERE
	ImGui::Render();
	ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
	glfwSwapBuffers(window); // Buffer Swap
	glfwPollEvents(); // Poll for Events
}

/// <summary>
/// OpenGL Draw Submethod. Synonymous with Matthias Muller's draw() method,
/// but rewritten in C++. Combines code from Matthia's repository as well as
/// youtube The Cherno's tutorial series on fundamentals of OpenGL in C++.
/// Does not clear nor swap buffers. Should be used within draw().
/// </summary>
void drawFluids() {

	//std::cout << "<Application.drawFluids()> Begun drawing." << std::endl;

	GLCall(glViewport(0, 0, width, height));

	/*
	* ====================================================================================================================
	*	Buffering
	*		If buffers not allocated, allocates them.
	* ====================================================================================================================
	*/

	// Grid Buffering
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
		GLCall(glBindBuffer(GL_ARRAY_BUFFER, 0)); 

		delete[] cellCenters;

		//std::cout << "<Application.drawFluids()> Grid Vertex Buffer initialized." << std::endl;
	}
	if (gridColorBuffer == -1) {
		GLCall(glGenBuffers(1, &gridColorBuffer));
		//std::cout << "<Application.drawFluids()> Grid Color Buffer initialized." << std::endl;
	}

	// Particles Buffering
	if (pointVertexBuffer == -1) {
		GLCall(glGenBuffers(1, &pointVertexBuffer));
		GLCall(glBindBuffer(GL_ARRAY_BUFFER, pointVertexBuffer));
		GLCall(glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 2 * fluid.numParticles, fluid.particlePos.data(), GL_DYNAMIC_DRAW));
		GLCall(glBindBuffer(GL_ARRAY_BUFFER, 0));
		//std::cout << "<Application.drawFluids()> Point Vertex Buffer initialized." << std::endl;
	}
	if (pointColorBuffer == -1) {
		GLCall(glGenBuffers(1, &pointColorBuffer));
		GLCall(glBindBuffer(GL_ARRAY_BUFFER, gridColorBuffer));
		GLCall(glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 3 * fluid.numParticles, fluid.particleColor.data(), GL_DYNAMIC_DRAW));
		GLCall(glBindBuffer(GL_ARRAY_BUFFER, 0));
		//std::cout << "<Application.drawFluids()> Point Colour Buffer initialized" << std::endl;
	}

	// Disk Buffering
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

		//std::cout << "<Application.drawFluids()> Disk Vertex Buffer initialized." << std::endl;
	}

	/*
	* ====================================================================================================================
	*	Binding and Drawing
	*		Actual drawing steps, if requested.
	* ====================================================================================================================
	*/

	// Grid Showing
	if (scene.showGrid) {
		float pointSize = 0.9 * fluid.h / simWidth * width;

		GLCall(glUseProgram(pointShader));
		GLCall(glUniform2f(glGetUniformLocation(pointShader, "domainSize"), simWidth, simHeight));
		GLCall(glUniform1f(glGetUniformLocation(pointShader, "pointSize"), pointSize));
		GLCall(glUniform1f(glGetUniformLocation(pointShader, "drawDisk"), 0.0f));

		GLCall(glBindBuffer(GL_ARRAY_BUFFER, gridVertexBuffer));
		GLCall(unsigned int posLoc = glGetAttribLocation(pointShader, "attrPosition"));
		GLCall(glEnableVertexAttribArray(posLoc));
		GLCall(glVertexAttribPointer(posLoc, 2, GL_FLOAT, GL_FALSE, 0, 0));

		GLCall(glBindBuffer(GL_ARRAY_BUFFER, gridColorBuffer));
		GLCall(glBufferData(GL_ARRAY_BUFFER, sizeof(float) * fluid.cellColor.size(), fluid.cellColor.data(), GL_DYNAMIC_DRAW));
		GLCall(unsigned int colorLoc = glGetAttribLocation(pointShader, "attrColor"));
		GLCall(glEnableVertexAttribArray(colorLoc));
		GLCall(glVertexAttribPointer(colorLoc, 3, GL_FLOAT, GL_FALSE, 0, 0));

		GLCall(glDrawArrays(GL_POINTS, 0, fluid.fNumCells));

		GLCall(glDisableVertexAttribArray(posLoc));
		GLCall(glDisableVertexAttribArray(colorLoc));

		GLCall(glBindBuffer(GL_ARRAY_BUFFER, 0));
		GLCall(glUseProgram(0));

		//std::cout << "<Application.drawFluids()> Grid shown." << std::endl;
	}

	// Water Showing
	if (scene.showParticles) {
		GLCall(glClear(GL_DEPTH_BUFFER_BIT));

		// Point Size
		float pointSize = 2.0 * fluid.particleRadius / simWidth * width; 
		// Vector Positions in Screen Space
		std::vector<float> vertexPositions = fluid.particlePos;
		for (int i = 0; i < fluid.numParticles; i++) {
			vertexPositions[i * 2 + 0] = vertexPositions[i * 2 + 0] / simWidth * 2 - 1;
			vertexPositions[i * 2 + 1] = vertexPositions[i * 2 + 1] / simHeight * 2 - 1;
		}

		GLCall(glUseProgram(0));
		GLCall(glEnable(GL_POINT_SMOOTH));
		GLCall(glEnableClientState(GL_VERTEX_ARRAY));
		GLCall(glEnableClientState(GL_COLOR_ARRAY));
		GLCall(glPointSize(pointSize));
		GLCall(glBindBuffer(GL_ARRAY_BUFFER, pointVertexBuffer));
		GLCall(glBufferData(GL_ARRAY_BUFFER, sizeof(float) * vertexPositions.size(), vertexPositions.data(), GL_DYNAMIC_DRAW));
		GLCall(glVertexPointer(2, GL_FLOAT, 0, 0));
		GLCall(glBindBuffer(GL_ARRAY_BUFFER, pointColorBuffer));
		GLCall(glBufferData(GL_ARRAY_BUFFER, sizeof(float) * fluid.particleColor.size(), fluid.particleColor.data(), GL_DYNAMIC_DRAW));
		GLCall(glColorPointer(3, GL_FLOAT, 0, 0));
		GLCall(glDrawArrays(GL_POINTS, 0, fluid.numParticles));
		GLCall(glBindBuffer(GL_ARRAY_BUFFER, 0));
		GLCall(glDisableClientState(GL_COLOR_ARRAY));
		GLCall(glDisableClientState(GL_VERTEX_ARRAY));
		GLCall(glDisable(GL_POINT_SMOOTH));
		GLCall(glUseProgram(0));

		//std::cout << "<Application.drawFluids()> Particles shown." << std::endl;
	}

	// Disk Drawing
	if (scene.showObstacle) {
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
		GLCall(glVertexAttribPointer(posLoc, 2, GL_FLOAT, GL_FALSE, 0, 0));
		GLCall(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, diskIdBuffer));
		GLCall(glDrawElements(GL_TRIANGLES, 3 * numSegs, GL_UNSIGNED_SHORT, 0));

		GLCall(glDisableVertexAttribArray(posLoc));
		GLCall(glUseProgram(0));

		//std::cout << "<Application.drawFluids()> Disk shown." << std::endl;
	}

	//std::cout << "<Application.drawFluids()> Finished drawing!" << std::endl;
}

/// <summary>
/// OpenGL Draw Submethod in charge of UI elements.
/// THIS FUNCTION IS BASED ON CODE PROVIDED TO US 
/// FOR ASSIGNEMNT 3 AS `fluid.cpp`.
/// </summary>
void drawUI() {

	prevTick = currTick;
	prevTickWasPaused = scene.paused;
	currTick = std::chrono::high_resolution_clock::now();
	double velocity = scene.dt / std::chrono::duration_cast<std::chrono::duration<double>>(currTick - prevTick).count();

	ImGui::Begin("Scene Settings");
	ImGui::Text("Shortkeys:\n\tSpace\t- pause/resume");
	ImGui::Checkbox("Paused", &scene.paused);
	ImGui::Text("Frame Number: %7i", scene.frameNr);
	if (!prevTickWasPaused) {
		ImGui::Text("Simulation Rate: %.2f times real speed", velocity);
	}
	else {
		ImGui::Text("Simulation Rate: PAUSED");
	}
	ImGui::SliderFloat("g", &scene.gravity, -9.81f, 9.81f);
	ImGui::SliderFloat("dt", &scene.dt, 0.001f, 0.1f);
	ImGui::SliderFloat("flipRatio", &scene.flipRatio, 0.0f, 1.0f);
	ImGui::SliderInt("numPressureIters", &scene.numPressureIters, 1, 200);
	ImGui::SliderInt("numParticleIters", &scene.numParticleIters, 1, 10);
	ImGui::SliderFloat("overRelaxation", &scene.overRelaxation, 1.0f, 2.0f);
	ImGui::Checkbox("Drift Compensation", &scene.compensateDrift);
	ImGui::Checkbox("Separate Particles", &scene.separateParticles);
	ImGui::Checkbox("Show Obstacle", &scene.showObstacle);
	ImGui::Checkbox("Show Particles", &scene.showParticles);
	ImGui::Checkbox("Show Grid", &scene.showGrid);
	ImGui::End();
}
