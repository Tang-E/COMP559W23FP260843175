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

		std::vector<float> u;
		std::vector<float> v;
		std::vector<float> du;
		std::vector<float> dv;
		std::vector<float> prevU;
		std::vector<float> prevV;
		std::vector<float> p;
		std::vector<float> s;
		std::vector<int> cellType;
		std::vector<float> cellColor;

		// Particle Properties

		int maxParticles;

		std::vector<float> particlePos;
		std::vector<float> particleColor;
		std::vector<float> particleVel;
		std::vector<float> particleDensity;
		float particleRestDensity;

		float particleRadius;
		float pInvSpacing;
		int pNumX;
		int pNumY;
		int pNumCells;

		std::vector<int> numCellParticles;
		std::vector<int> firstCellParticle;
		std::vector<int> cellParticleIds;

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
