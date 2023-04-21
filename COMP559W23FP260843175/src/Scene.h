/*
* Scene code as a step in for the javascript scene data structure.
* Code here is adapted from Matthias Müller's tenMinutePhysics
* FLIP Simulator code available on his GitHub repository here=
* https=//github.com/matthias-research/pages/blob/master/tenMinutePhysics/18-flip.html
* 
* This is a pretty simpel file so it works as a header.
*
* I am adapting his code written in javascript into C++ code
* because I don't want to write in javascript.
*
* @author Edwin Pan (260843175) for COMP559 Winter 2023 Final Project
*/

#ifndef SCENE_H
#define SCENE_H

class Scene {
public:
	float gravity = -9.81f;
	float dt = 1.0f / 120.0f;
	float flipRatio = 0.9f;
	int numPressureIters = 100;
	int numParticleIters = 2;
	int frameNr = 0;
	float overRelaxation = 1.9f;
	bool compensateDrift = true;
	bool separateParticles = true;
	float obstacleX = 0.0f;
	float obstacleY = 0.0f;
	float obstacleRadius = 0.15f;
	bool paused = true;
	bool showObstacle = true;
	float obstacleVelX = 0.0f;
	float obstacleVelY = 0.0f;
	bool showParticles = true;
	bool showGrid = false;
	float cohesionMaxAccel;
	float cohesionFallOffRate; 
	float cohesionMaxDistance;

	Scene() {
		gravity = -9.81f;
		dt = 1.0f / 120.0f;
		flipRatio = 0.9f;
		numPressureIters = 100;
		numParticleIters = 2;
		frameNr = 0;
		overRelaxation = 1.9f;
		compensateDrift = true;
		separateParticles = true;
		obstacleX = 0.0f;
		obstacleY = 0.0f;
		obstacleRadius = 0.15f;
		paused = true;
		showObstacle = true;
		obstacleVelX = 0.0f;
		obstacleVelY = 0.0f;
		showParticles = true;
		showGrid = false;
		cohesionMaxAccel = 0.01f;
		cohesionFallOffRate = 10.0f;
		cohesionMaxDistance = 0.5f;
	}

};

#endif
