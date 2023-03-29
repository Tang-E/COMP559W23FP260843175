/*
* Simulator Code header file
*/

#ifndef FLIPFLUID_H
#define FLIPFLUID_H

class FlipFluid {

	public:

		// Pointers

		Scene* scene;

		// Fluid Properties

		float density;
		int fNumX;
		int fNumY;
		float h;
		float fInvSpacing;
		int fNumCells;

		float* u;
		float* v;
		float* du;
		float* dv;
		float* prevU;
		float* prevV;
		float* p;
		float* s;
		int* cellType;
		float* cellColor;

		// Particle Properties

		int maxParticles;

		float* particlePos;
		float* particleColor;
		float* particleVel;
		float* particleDensity;
		float particleRestDensity;

		float particleRadius;
		float pInvSpacing;
		int pNumX;
		int pNumY;
		int pNumCells;

		int* numCellParticles;
		int* firstCellParticle;
		int* cellParticleIds;

		int numParticles;

		// Operations

		FlipFluid(float density, int width, int height, float spacing, float particleRadius, int maxParticles, Scene* scene);
		void integrateParticles(float dt, float gravity);
		void pushParticlesApart(int numIters);
		void handleParticleCollisions(float obstacleX, float obstacleY, float obstacleRadius);
		void updateParticleDensity();
		void transferVelocities(bool toGrid, float flipRatio);
		void solveIncompressibility(int numIters, float dt, float overRelaxation, bool compensateDrift = true);
		void updateParticleColors();
		void setSciColor(int cellNr, float val, float minVal, float maxVal);
		void updateCellColors();
		void simulate(float dt, float gravity, float flipRatio, int numPressureIters, int numParticleIters, float overRelaxation, bool compensateDrift, bool separateParticles, float obstacleX, float abstacleY, float obstacleRadius);

};

#endif
