/*
* Simulator Code header file
*/

#ifndef SIMULATOR_H
define SIMULATOR_H

class FlipFluid {

	public:

		FlipFluid(double density, int width, int height, double spacing, double particleRadius, int maxParticles);
		void integrateParticles(double dt, double gravity);
		void pushParticlesApart(int numIters);
		void handleParticleCollisions(double obstacleX, double obstacleY, double obstacleRadius);
		void updateParticleDensity();
		void transferVelocities(bool toGrid, float flipRatio);
		void solveIncompressibility(int numIters, double dt, double overRelaxation, bool compensateDrift = true);
		void updateParticleColors();
		void setSciColor(int cellNr, double val, double minVal, double maxVal);
		void updateCellColors();
		void simulate(double dt, double gravity, float flipRatio, int numPressureIters, int numParticleIters, double overRelaxation, bool compensateDrift, bool separateParticles, double obstacleX, double abstacleY, double obstacleRadius);

};

#endif
