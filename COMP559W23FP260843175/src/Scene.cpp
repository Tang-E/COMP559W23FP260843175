/*
* Scene code as a step in for the javascript scene data structure.
* Code here is adapted from Matthias Müller's tenMinutePhysics
* FLIP Simulator code available on his GitHub repository here=
* https=//github.com/matthias-research/pages/blob/master/tenMinutePhysics/18-flip.html
*
* I am adapting his code written in javascript into C++ code
* because I don't want to write in javascript.
*
* @author Edwin Pan (260843175) for COMP559 Winter 2023 Final Project
*/

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

		Scene() {
			gravity = -9.81;
			dt = 1.0 / 120.0;
			flipRatio = 0.9;
			numPressureIters = 100;
			numParticleIters = 2;
			frameNr = 0;
			overRelaxation = 1.9;
			compensateDrift = true;
			separateParticles = true;
			obstacleX = 0.0;
			obstacleY = 0.0;
			obstacleRadius = 0.15;
			paused = true;
			showObstacle = true;
			obstacleVelX = 0.0;
			obstacleVelY = 0.0;
			showParticles = true;
			showGrid = false;
		}

};