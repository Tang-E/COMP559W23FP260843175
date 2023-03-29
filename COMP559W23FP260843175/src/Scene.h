/*
* Header file for scene properties datatype
*/

#ifndef SCENE_H
#define SCENE_H

class Scene {
	public:
		float gravity = -9.81;
		float dt = 1.0 / 120.0;
		float flipRatio = 0.9;
		int numPressureIters = 100;
		int numParticleIters = 2;
		int frameNr = 0;
		float overRelaxation = 1.9;
		bool compensateDrift = true;
		bool separateParticles = true;
		float obstacleX = 0.0;
		float obstacleY = 0.0;
		float obstacleRadius = 0.15;
		bool paused = true;
		bool showObstacle = true;
		float obstacleVelX = 0.0;
		float obstacleVelY = 0.0;
		bool showParticles = true;
		bool showGrid = false;
		Scene();
};

#endif
