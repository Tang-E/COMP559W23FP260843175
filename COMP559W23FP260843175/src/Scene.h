/*
* Header file for scene properties datatype
*/

#ifndef SCENE_H
#define SCENE_H

class Scene {
	public:
		double gravity = -9.81;
		double dt = 1.0 / 120.0;
		double flipRatio = 0.9;
		int numPressureIters = 100;
		int numParticleIters = 2;
		int frameNr = 0;
		double overRelaxation = 1.9;
		bool compensateDrift = true;
		bool separateParticles = true;
		double obstacleX = 0.0;
		double obstacleY = 0.0;
		double obstacleRadius = 0.15;
		bool paused = true;
		bool showObstacle = true;
		double obstacleVelX = 0.0;
		double obstacleVelY = 0.0;
		bool showParticles = true;
		bool showGrid = false;
		//Simulator* fluid = nullptr;
		Scene();
};

#endif
