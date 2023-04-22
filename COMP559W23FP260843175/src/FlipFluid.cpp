/*
* Simulator Code in charge of the particle and grid simulations.
* Code here is adapted from Matthias Müller's tenMinutePhysics
* FLIP Simulator code available on his GitHub repository here:
* https://github.com/matthias-research/pages/blob/master/tenMinutePhysics/18-flip.html
*
* I am adapting his code written in javascript into C++ code
* because I don't want to write in javascript.
*
* @author Edwin Pan (260843175) for COMP559 Winter 2023 Final Project
*/


#include "includes.h"
#include "scene.h"
#include "HelperFuncs.h"


/// <summary>
/// FLIP Fluid Simulator object.
/// </summary>
class FlipFluid {

	private:
		glm::vec2 simDimensions;

	public:

		// Constants/Definitions

		const int U_FIELD = 0;
		const int V_FIELD = 1;
		const int FLUID_CELL = 0;
		const int AIR_CELL = 1;
		const int SOLID_CELL = 2;

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

		// Timing Statistics
		double previousFrameTime = -1;

		/// <summary>
		/// Constructor. Doesn't do everything though and will still require an external
		/// setupmethod which is located in Application.cpp.
		/// </summary>
		/// <param name="density"> Fluid Density </param>
		/// <param name="width"> Tank Width in Metres </param>
		/// <param name="height" Tank Width in Metres ></param>
		/// <param name="spacing"> Cell Width in Metres </param>
		/// <param name="particleRadius"> Particle Radius in Metres </param>
		/// <param name="maxParticles"></param>
		/// <param name="scene"></param>
		FlipFluid(float density, float width, float height, float spacing, float particleRadius, int maxParticles, Scene* scene) {

			// Domain Properties
			simDimensions = glm::vec2(width, height);

			// Pointer fields

			this->scene = scene;

			// Fluid Properties

			this->density = density;
			this->fNumX = floor(width / spacing) + 1;
			this->fNumY = floor(height / spacing) + 1;
			this->h = max(width / this->fNumX, height / this->fNumY);
			this->fInvSpacing = 1.0 / this->h;
			this->fNumCells = this->fNumX * this->fNumY;
			std::cout << "<FlipFluid> \n"
				<< "\tfNumX = " << this->fNumX << "\n"
				<< "\tfNumY = " << this->fNumY << "\n"
				<< "\th = " << this->h << "\n"
				<< "\tfInvSpacing = " << this->fInvSpacing << "\n"
				<< "\tfNumCells = " << this->fNumCells << std::endl;

			this->u.resize(fNumCells);
			this->v.resize(fNumCells);
			this->du.resize(fNumCells);
			this->dv.resize(fNumCells);
			this->prevU.resize(fNumCells);
			this->prevV.resize(fNumCells);
			this->p.resize(fNumCells);
			this->s.resize(fNumCells);
			this->cellType.resize(fNumCells);
			this->cellColor.resize(3 * fNumCells);

			// Particle Properties

			this->maxParticles = maxParticles;

			this->particlePos.resize(2 * this->maxParticles);
			this->particleColor.resize(3 * this->maxParticles);
			for (int i = 0; i < this->maxParticles; i++) {
				this->particleColor[3 * i + 2] = 1.0; // Make blue, not white!
			}
			this->particleVel.resize(2 * this->maxParticles);
			this->particleDensity.resize(this->fNumCells);
			this->particleRestDensity = 0.0;

			this->particleRadius = particleRadius;
			this->pInvSpacing = 1.0 / (2.2 * particleRadius);
			this->pNumX = floor((float)width * this->pInvSpacing) + 1;
			this->pNumY = floor((float)height * this->pInvSpacing) + 1;
			this->pNumCells = this->pNumX * this->pNumY;

			this->numCellParticles.resize(this->pNumCells);
			this->firstCellParticle.resize(this->pNumCells + 1);
			this->cellParticleIds.resize(this->maxParticles);

			this->numParticles = 0;

		}

		/// <summary>
		/// Simplecton integrates particles on step.
		/// Accelerates y velocity by gravity; then
		/// moves position by x and y velocity.
		/// </summary>
		/// <param name="dt"></param>
		/// <param name="gravity"></param>
		void integrateParticles(float dt, float gravity)
		{
			for (int i = 0; i < numParticles; i++) {
				particleVel[2 * i + 1] += dt * gravity;
				particlePos[2 * i] += particleVel[2 * i] * dt;
				particlePos[2 * i + 1] += particleVel[2 * i + 1] * dt;
			}
		}

		glm::vec2 calculateCohesiveAccelByPositionAndStrength(glm::vec2 p, glm::vec2 target, float strength, float maxAccel, float maxDistance, float fallOffCoeff) {
			glm::vec2 dp = target - p;
			float dist = glm::length(dp);
			float accel = maxAccel * std::exp(-fallOffCoeff * dist) * strength;
			return dp / dist * accel;
		}

		/// <summary>
		/// Calculates a acceleration by "cohesion" forces to a particle at
		/// position p toward every particle listed in particleIndices based
		/// on pull properties maxAccel, maxDistance, and fallOffCoeff.
		/// </summary>
		/// <param name="p"></param>
		/// <param name="particleIndices"></param>
		/// <param name="maxAccel"></param>
		/// <param name="maxDistance"></param>
		/// <param name="fallOffCoeff"></param>
		/// <returns></returns>
		glm::vec2 calculateCohesiveAccelForParticleIndices(glm::vec2 p, std::list<int> particleIndices, float maxAccel, float maxDistance, float fallOffCoeff) {
			// Start summing the acceleration a_i for particle at position p_i
			glm::vec2 a(0.0f, 0.0f);
			// For all nearby particles j
			for (int j : particleIndices) {
				glm::vec2 p_j(particlePos[2 * j + 0], particlePos[2 * j + 1]);
				glm::vec2 dp = p_j - p;
				float dist = glm::length(dp);
				if (dist > maxDistance) continue;
				float accel = maxAccel * std::exp(-fallOffCoeff * dist);
				a += dp / dist * accel;
			}
			// Return
			return a;
		}

		/// <summary>
		/// Simplecton integrate particles with a inter-particle cohesive force.
		/// We want a force that keeps nearby particles together and ignores 
		/// faraway particles. We thus cycle through each particle pair and
		/// give them a pull force toward each other based on negative exponential
		/// of their distance such that a maximum force can be set and a natural
		/// falloff occurs based on larger distance. Both those two terms are
		/// adjustable and we introduce a maxDistance that allows us to flat
		/// out ignore pairs of particles too far from each other.
		/// 
		/// This method is not in Matthias Muller's provided code.
		/// Custom implementation by Edwin Pan 260843175
		/// 
		/// </summary>
		/// <param name="dt"></param>
		/// <param name="maxAccel">Accecleration in m/s^2 if particles atop each other.</param>
		/// <param name="fallOffRate">Bigger Falloff Rate means greater sensitivity to far-ness</param>
		/// <param name="maxDistance">Max distance after which pair of particles' cohesion will be ignored</param>
		void addCohesionAccel(float dt, float cMaxAccel, float cMaxDistance, bool fastSimulate) {
			// Protect from invalid input
			if (cMaxAccel <= 0 || cMaxDistance <= 0) return;
			// Calculate the cohesive acceleration curve properties
			// accel(distance) = euler^( fallOffCoeff*distance )
			float epsilon = 1e-6; // Acceleration at cMaxDistance such that it is insignificant
			float fallOffCoeff = -std::log(epsilon) / cMaxDistance;

			// Subdivide the space. Divide space into squares put into sets of VERTICAL stripes.
			int numStripes = simDimensions.x / cMaxDistance + 1;
			int cellsPerStripe = simDimensions.y / cMaxDistance + 1;
			std::vector<std::vector<std::list<int>>> stripes; // The Grid
			stripes.resize(numStripes);
			for (int i = 0; i < numStripes; i++) {
				stripes[i].resize(cellsPerStripe);
				for (int j = 0; j < cellsPerStripe; j++) {
					stripes[i][j] = std::list<int>();
				}
			}

			// Put each particle into its respective grid block.
			for (int i = 0; i < numParticles; i++) {
				int stripeNo = clamp(particlePos[i * 2 + 0] / cMaxDistance, 0, numStripes);
				int subCellNo = clamp(particlePos[i * 2 + 1] / cMaxDistance, 0, cellsPerStripe);
				stripes[stripeNo][subCellNo].push_back(i);
			}
			
			// Fast Simulation Parametre
			std::vector<std::vector<glm::vec2>> meanPos;	// Mean of Positions per cell
			std::vector<std::vector<float>> numPos;			// Number of Positions per cell
			if (fastSimulate) {
				// Size the meanPos and numPos 2D arrays
				meanPos.resize(numStripes);
				numPos.resize(numStripes);
				for (int stripe = 0; stripe < numStripes; stripe++) {
					meanPos[stripe].resize(cellsPerStripe);
					numPos[stripe].resize(cellsPerStripe);
					for (int cell = 0; cell < cellsPerStripe; cell++) {
						meanPos[stripe][cell] = glm::vec2(0, 0);
						numPos[stripe][cell] = 1e-9;
					}
				}
				// Fill arrays with mean position within each cell and number of positions per cell
				for (int stripe = 0; stripe < numStripes; stripe++) {
					for (int cell = 0; cell < cellsPerStripe; cell++) {
						for (int i : stripes[stripe][cell]) {
							meanPos[stripe][cell].x += particlePos[2 * i + 0];
							meanPos[stripe][cell].y += particlePos[2 * i + 1];
							numPos[stripe][cell] += 1.0001f; // Code is FUCKING MORONIC. For some reason if I don't make use of decimals, the integer-ness makes vector freak the fuck out and crash. fucking stupid.
						}
						meanPos[stripe][cell] /= numPos[stripe][cell];
					}
				}

			}

			// For each vertical stripe
			for (int stripe = 0; stripe < numStripes; stripe++) {
				bool rightOK = stripe < numStripes - 1; // Right is +1
				bool leftOK = stripe > 0; // Left is -1

				// For each cell in stripe
				for (int cell = 0; cell < cellsPerStripe; cell++) {
					bool upOK = cell < cellsPerStripe - 1; // Up is +1
					bool downOK = cell > 0; // Down is -1
					// For each particle i in this cell

					for (int i : stripes[stripe][cell]) {
						// Position and Acceleration of this particle
						glm::vec2 p_i(particlePos[i * 2 + 0], particlePos[i * 2 + 1]);
						glm::vec2 a_i(0.0f, 0.0f);

						if (fastSimulate) { // Fast Simulate: Calculate Cohesion Force between particle i and regions
							// Calculate cohesive acceleration on this particle based on this cell and neighoubr cells
							if (rightOK && upOK)	a_i += calculateCohesiveAccelByPositionAndStrength(p_i, meanPos[stripe + 1][cell + 1], numPos[stripe + 1][cell + 1], cMaxAccel, cMaxDistance, fallOffCoeff);
							if (rightOK)			a_i += calculateCohesiveAccelByPositionAndStrength(p_i, meanPos[stripe + 1][cell + 0], numPos[stripe + 1][cell + 0], cMaxAccel, cMaxDistance, fallOffCoeff);
							if (rightOK && downOK)	a_i += calculateCohesiveAccelByPositionAndStrength(p_i, meanPos[stripe + 1][cell - 1], numPos[stripe + 1][cell - 1], cMaxAccel, cMaxDistance, fallOffCoeff);
							if (leftOK && upOK)		a_i += calculateCohesiveAccelByPositionAndStrength(p_i, meanPos[stripe - 1][cell + 1], numPos[stripe - 1][cell + 1], cMaxAccel, cMaxDistance, fallOffCoeff);
							if (leftOK)				a_i += calculateCohesiveAccelByPositionAndStrength(p_i, meanPos[stripe - 1][cell + 0], numPos[stripe - 1][cell + 0], cMaxAccel, cMaxDistance, fallOffCoeff);
							if (leftOK && downOK)	a_i += calculateCohesiveAccelByPositionAndStrength(p_i, meanPos[stripe - 1][cell - 1], numPos[stripe - 1][cell - 1], cMaxAccel, cMaxDistance, fallOffCoeff);
							if (upOK)				a_i += calculateCohesiveAccelByPositionAndStrength(p_i, meanPos[stripe + 0][cell + 1], numPos[stripe + 0][cell + 1], cMaxAccel, cMaxDistance, fallOffCoeff);
							if (downOK)				a_i += calculateCohesiveAccelByPositionAndStrength(p_i, meanPos[stripe + 0][cell - 1], numPos[stripe + 0][cell - 1], cMaxAccel, cMaxDistance, fallOffCoeff);
							a_i += calculateCohesiveAccelByPositionAndStrength(p_i, meanPos[stripe + 0][cell + 0], numPos[stripe + 0][cell + 0], cMaxAccel, cMaxDistance, fallOffCoeff);

						}
						else { // Slow Simulate: Calculate Cohesion Force between particle i and every single particle in regions
							// Calculate cohesive acceleration on this particle based on this cell and neighoubr cells
							if (rightOK && upOK)	a_i += calculateCohesiveAccelForParticleIndices(p_i, stripes[stripe + 1][cell + 1], cMaxAccel, cMaxDistance, fallOffCoeff);
							if (rightOK)			a_i += calculateCohesiveAccelForParticleIndices(p_i, stripes[stripe + 1][cell + 0], cMaxAccel, cMaxDistance, fallOffCoeff);
							if (rightOK && downOK)	a_i += calculateCohesiveAccelForParticleIndices(p_i, stripes[stripe + 1][cell - 1], cMaxAccel, cMaxDistance, fallOffCoeff);
							if (leftOK && upOK)		a_i += calculateCohesiveAccelForParticleIndices(p_i, stripes[stripe - 1][cell + 1], cMaxAccel, cMaxDistance, fallOffCoeff);
							if (leftOK)				a_i += calculateCohesiveAccelForParticleIndices(p_i, stripes[stripe - 1][cell + 0], cMaxAccel, cMaxDistance, fallOffCoeff);
							if (leftOK && downOK)	a_i += calculateCohesiveAccelForParticleIndices(p_i, stripes[stripe - 1][cell - 1], cMaxAccel, cMaxDistance, fallOffCoeff);
							if (upOK)				a_i += calculateCohesiveAccelForParticleIndices(p_i, stripes[stripe + 0][cell + 1], cMaxAccel, cMaxDistance, fallOffCoeff);
							if (downOK)				a_i += calculateCohesiveAccelForParticleIndices(p_i, stripes[stripe + 0][cell - 1], cMaxAccel, cMaxDistance, fallOffCoeff);
							std::list<int> particlesInThisCell = stripes[stripe][cell];
							particlesInThisCell.remove(i); // Ignore this own particle
							a_i += calculateCohesiveAccelForParticleIndices(p_i, particlesInThisCell, cMaxAccel, cMaxDistance, fallOffCoeff);

						}

						//// Integrate acceleration on particle i into velocity on particle i
						particleVel[i * 2 + 0] += a_i.x * dt;
						particleVel[i * 2 + 1] += a_i.y * dt;
					}
				}

			}
			
		}



		/// <summary>
		/// Applies repulsion between particles efficiently via Gauss-Seidel algorithm 
		/// that involves filtering out collision checks of particles that would be too 
		/// far from each other. Said filtering is something Matthias Müller has figured 
		/// out and I'm not going to bother trying to understand it.
		/// </summary>
		/// <param name="numIters"></param>
		void pushParticlesApart(int numIters)
		{
			float colorDiffusionCoeff = 0.001;

			// count particles per cell

			std::fill(numCellParticles.begin(), numCellParticles.end(), 0);

			for (int i = 0; i < numParticles; i++) {
				float x = particlePos[2 * i];
				float y = particlePos[2 * i + 1];

				int xi = clamp(floor(x * pInvSpacing), 0, pNumX - 1);
				int yi = clamp(floor(y * pInvSpacing), 0, pNumY - 1);
				int cellNr = xi * pNumY + yi;
				numCellParticles[cellNr]++;
			}

			// partial sums

			int first = 0;

			for (int i = 0; i < pNumCells; i++) {
				first += numCellParticles[i];
				firstCellParticle[i] = first;
			}
			firstCellParticle[pNumCells] = first;		// guard

			// fill particles into cells

			for (int i = 0; i < numParticles; i++) {
				float x = particlePos[2 * i];
				float y = particlePos[2 * i + 1];

				int xi = clamp(floor(x * pInvSpacing), 0, pNumX - 1);
				int yi = clamp(floor(y * pInvSpacing), 0, pNumY - 1);
				int cellNr = xi * pNumY + yi;
				firstCellParticle[cellNr]--;
				cellParticleIds[firstCellParticle[cellNr]] = i;
			}

			// push particles apart

			float minDist = 2.0 * particleRadius;
			float minDist2 = minDist * minDist;

			for (int iter = 0; iter < numIters; iter++) {

				for (int i = 0; i < numParticles; i++) {
					float px = particlePos[2 * i];
					float py = particlePos[2 * i + 1];

					float pxi = floor(px * pInvSpacing);
					float pyi = floor(py * pInvSpacing);
					int x0 = max(pxi - 1, 0);
					int y0 = max(pyi - 1, 0);
					int x1 = min(pxi + 1, pNumX - 1);
					int y1 = min(pyi + 1, pNumY - 1);

					for (int xi = x0; xi <= x1; xi++) {
						for (int yi = y0; yi <= y1; yi++) {
							int cellNr = xi * pNumY + yi;
							int first = firstCellParticle[cellNr];
							int last = firstCellParticle[cellNr + 1];
							for (int j = first; j < last; j++) {
								int id = cellParticleIds[j];
								if (id == i)
									continue;
								float qx = particlePos[2 * id];
								float qy = particlePos[2 * id + 1];

								float dx = qx - px;
								float dy = qy - py;
								float d2 = dx * dx + dy * dy;
								if (d2 > minDist2 || d2 == 0.0)
									continue;
								float d = sqrt(d2);
								float s = 0.5 * (minDist - d) / d;
								dx *= s;
								dy *= s;
								particlePos[2 * i] -= dx;
								particlePos[2 * i + 1] -= dy;
								particlePos[2 * id] += dx;
								particlePos[2 * id + 1] += dy;

								// diffuse colors

								for (int k = 0; k < 3; k++) {
									float color0 = particleColor[3 * i + k];
									float color1 = particleColor[3 * id + k];
									float color = (color0 + color1) * 0.5;
									particleColor[3 * i + k] = color0 + (color - color0) * colorDiffusionCoeff;
									particleColor[3 * id + k] = color1 + (color - color1) * colorDiffusionCoeff;
								}
							}
						}
					}
				}
			}
		}

		/// <summary>
		/// Handles particle collision with the obstacle and walls.
		/// </summary>
		/// <param name="obstacleX"></param>
		/// <param name="obstacleY"></param>
		/// <param name="obstacleRadius"></param>
		void handleParticleCollisions(float obstacleX, float obstacleY, float obstacleRadius)
		{
			float h = 1.0 / fInvSpacing;
			float r = particleRadius;
			float or1 = obstacleRadius;
			float or2 = or1 * or1;
			float minDist = obstacleRadius + r;
			float minDist2 = minDist * minDist;

			float minX = h + r;
			float maxX = (fNumX - 1) * h - r;
			float minY = h + r;
			float maxY = (fNumY - 1) * h - r;


			for (int i = 0; i < numParticles; i++) {
				float x = particlePos[2 * i];
				float y = particlePos[2 * i + 1];

				float dx = x - obstacleX;
				float dy = y - obstacleY;
				float d2 = dx * dx + dy * dy;

				// obstacle collision

				if (d2 < minDist2) {

					// var d = Math.sqrt(d2);
					// var s = (minDist - d) / d;
					// x += dx * s;
					// y += dy * s;

					particleVel[2 * i] = scene->obstacleVelX;
					particleVel[2 * i + 1] = scene->obstacleVelY;
				}

				// wall collisions

				if (x < minX) {
					x = minX;
					particleVel[2 * i] = 0.0;

				}
				if (x > maxX) {
					x = maxX;
					particleVel[2 * i] = 0.0;
				}
				if (y < minY) {
					y = minY;
					particleVel[2 * i + 1] = 0.0;
				}
				if (y > maxY) {
					y = maxY;
					particleVel[2 * i + 1] = 0.0;
				}
				particlePos[2 * i] = x;
				particlePos[2 * i + 1] = y;
			}
		}

		/// <summary>
		/// Updates the particle-density reading at each cell.
		/// Used as a countermeasure for drift - that is, the 
		/// loss of water volume as a consequence of particles
		/// that collide and are not accounted for. This accounts
		/// for them by representing high densities of particles
		/// which are then used in divergence calculations.
		/// </summary>
		void updateParticleDensity()
		{
			int n = fNumY;
			float h = this->h;
			float h1 = fInvSpacing;
			float h2 = 0.5 * h;

			std::fill(particleDensity.begin(), particleDensity.end(), 0.0f);

			for (int i = 0; i < numParticles; i++) {
				float x = particlePos[2 * i];
				float y = particlePos[2 * i + 1];

				x = clamp(x, h, (fNumX - 1) * h);
				y = clamp(y, h, (fNumY - 1) * h);

				int x0 = floor((x - h2) * h1);
				float tx = ((x - h2) - x0 * h) * h1;
				int x1 = min(x0 + 1, fNumX - 2);

				int y0 = floor((y - h2) * h1);
				float ty = ((y - h2) - y0 * h) * h1;
				int y1 = min(y0 + 1, fNumY - 2);

				float sx = 1.0 - tx;
				float sy = 1.0 - ty;

				if (x0 < fNumX && y0 < fNumY) particleDensity[x0 * n + y0] += sx * sy;
				if (x1 < fNumX && y0 < fNumY) particleDensity[x1 * n + y0] += tx * sy;
				if (x1 < fNumX && y1 < fNumY) particleDensity[x1 * n + y1] += tx * ty;
				if (x0 < fNumX && y1 < fNumY) particleDensity[x0 * n + y1] += sx * ty;
			}

			if (particleRestDensity == 0.0) {
				float sum = 0.0;
				int numFluidCells = 0;

				for (int i = 0; i < fNumCells; i++) {
					if (cellType[i] == FLUID_CELL) {
						sum += particleDensity[i];
						numFluidCells++;
					}
				}

				if (numFluidCells > 0)
					particleRestDensity = sum / numFluidCells;
			}

			// 			for (var xi = 1; xi < fNumX; xi++) {
			// 				for (var yi = 1; yi < fNumY; yi++) {
			// 					var cellNr = xi * n + yi;
			// 					if (cellType[cellNr] != FLUID_CELL)
			// 						continue;
			// 					var hx = h;
			// 					var hy = h;

			// 					if (cellType[(xi - 1) * n + yi] == SOLID_CELL || cellType[(xi + 1) * n + yi] == SOLID_CELL)
			// 						hx -= particleRadius;
			// 					if (cellType[xi * n + yi - 1] == SOLID_CELL || cellType[xi * n + yi + 1] == SOLID_CELL)
			// 						hy -= particleRadius;

			// 					var scale = h * h / (hx * hy)
			// 					d[cellNr] *= scale;
			// 				}
			// 			}
		}

		/// <summary>
		/// Function for transfering velocity data between the grid and the particles.
		/// If toGrid is true, then velocity data is transferred from the particles
		/// and into the grid. If toGrid is false, then velocity data is transferred
		/// from the grid and into the particles. flipRatio is useful for when toGrid
		/// is false as when velocity is being transferred from the grid and into the
		/// particles, we want to know whether to overried each particle's velocity
		/// by their grid values or whether to simply add to each particle's velocities.
		/// </summary>
		/// <param name="toGrid"> If True: To Grid; If False: To Particles</param>
		/// <param name="flipRatio"> 0%: Viscous; 100% Wild </param>
		void transferVelocities(bool toGrid, float flipRatio)
		{
			int n = fNumY;
			float h = this->h;
			float h1 = fInvSpacing;
			float h2 = 0.5 * h;

			if (toGrid) {

				prevU = u;
				prevV = v;

				std::fill(du.begin(), du.end(), 0.0f);
				std::fill(dv.begin(), dv.end(), 0.0f);
				std::fill(u.begin(), u.end(), 0.0f);
				std::fill(v.begin(), v.end(), 0.0f);

				for (int i = 0; i < fNumCells; i++)
					cellType[i] = s[i] == 0.0 ? SOLID_CELL : AIR_CELL;

				for (int i = 0; i < numParticles; i++) {
					float x = particlePos[2 * i];
					float y = particlePos[2 * i + 1];
					int xi = clamp(floor(x * h1), 0, fNumX - 1);
					int yi = clamp(floor(y * h1), 0, fNumY - 1);
					int cellNr = xi * n + yi;
					if (cellType[cellNr] == AIR_CELL)
						cellType[cellNr] = FLUID_CELL;
				}
			}

			for (int component = 0; component < 2; component++) {

				float dx = component == 0 ? 0.0 : h2;
				float dy = component == 0 ? h2 : 0.0;

				std::vector<float>* f = component == 0 ? &u : &v;
				std::vector<float>* prevF = component == 0 ? &prevU : &prevV;
				std::vector<float>* d = component == 0 ? &du : &dv;

				for (int i = 0; i < numParticles; i++) {
					float x = particlePos[2 * i];
					float y = particlePos[2 * i + 1];

					x = clamp(x, h, (fNumX - 1) * h);
					y = clamp(y, h, (fNumY - 1) * h);

					int x0 = min(floor((x - dx) * h1), fNumX - 2);
					float tx = ((x - dx) - x0 * h) * h1;
					int x1 = min(x0 + 1, fNumX - 2);

					int y0 = min(floor((y - dy) * h1), fNumY - 2);
					float ty = ((y - dy) - y0 * h) * h1;
					int y1 = min(y0 + 1, fNumY - 2);

					float sx = 1.0 - tx;
					float sy = 1.0 - ty;

					float d0 = sx * sy;
					float d1 = tx * sy;
					float d2 = tx * ty;
					float d3 = sx * ty;

					int nr0 = x0 * n + y0;
					int nr1 = x1 * n + y0;
					int nr2 = x1 * n + y1;
					int nr3 = x0 * n + y1;

					if (toGrid) {
						float  pv = particleVel[2 * i + component];
						(*f)[nr0] += pv * d0;  (*d)[nr0] += d0;
						(*f)[nr1] += pv * d1;  (*d)[nr1] += d1;
						(*f)[nr2] += pv * d2;  (*d)[nr2] += d2;
						(*f)[nr3] += pv * d3;  (*d)[nr3] += d3;
					}
					else {
						int offset = component == 0 ? n : 1;
						float valid0 = cellType[nr0] != AIR_CELL || cellType[nr0 - offset] != AIR_CELL ? 1.0 : 0.0;
						float valid1 = cellType[nr1] != AIR_CELL || cellType[nr1 - offset] != AIR_CELL ? 1.0 : 0.0;
						float valid2 = cellType[nr2] != AIR_CELL || cellType[nr2 - offset] != AIR_CELL ? 1.0 : 0.0;
						float valid3 = cellType[nr3] != AIR_CELL || cellType[nr3 - offset] != AIR_CELL ? 1.0 : 0.0;

						float v = particleVel[2 * i + component];
						float d = valid0 * d0 + valid1 * d1 + valid2 * d2 + valid3 * d3;

						if (d > 0.0) {

							float picV = (valid0 * d0 * (*f)[nr0] + valid1 * d1 * (*f)[nr1] + valid2 * d2 * (*f)[nr2] + valid3 * d3 * (*f)[nr3]) / d;
							float corr = (valid0 * d0 * ((*f)[nr0] - (*prevF)[nr0]) + valid1 * d1 * ((*f)[nr1] - (*prevF)[nr1])
								+ valid2 * d2 * ((*f)[nr2] - (*prevF)[nr2]) + valid3 * d3 * ((*f)[nr3] - (*prevF)[nr3])) / d;
							float flipV = v + corr;

							particleVel[2 * i + component] = (1.0 - flipRatio) * picV + flipRatio * flipV;
						}
					}
				}

				if (toGrid) {
					for (int i = 0; i < /*f.length*/ fNumCells; i++) {
						if ((*d)[i] > 0.0)
							(*f)[i] /= (*d)[i];
					}

					// restore solid cells

					for (int i = 0; i < fNumX; i++) {
						for (int j = 0; j < fNumY; j++) {
							bool solid = cellType[i * n + j] == SOLID_CELL;
							if (solid || (i > 0 && cellType[(i - 1) * n + j] == SOLID_CELL))
								u[i * n + j] = prevU[i * n + j];
							if (solid || (j > 0 && cellType[i * n + j - 1] == SOLID_CELL))
								v[i * n + j] = prevV[i * n + j];
						}
					}
				}
			}
		}

		/// <summary>
		/// Divergence Solver.
		/// Solves incompressibility through gauss-seidel by going over each 
		/// grid element and reducing divergence to 0 at every pass with the 
		/// ability to undershoot or overshoot with overRelaxation term.
		/// Can compensate for volume-drift (ie, high density particle areas).
		/// </summary>
		/// <param name="numIters"></param>
		/// <param name="dt"></param>
		/// <param name="overRelaxation"></param>
		/// <param name="compensateDrift"></param>
		void solveIncompressibility(int numIters, float dt, float overRelaxation, bool compensateDrift = true) {

			std::fill(p.begin(), p.end(), 0.0f);
			prevU = u;
			prevV = v;

			int n = fNumY;
			float cp = density * h / dt;

			//for (int i = 0; i < fNumCells; i++) { //wtf does this do for us??????
			//	float u = this->u[i];
			//	float v = this->v[i];
			//}

			for (int iter = 0; iter < numIters; iter++) {

				for (int i = 1; i < fNumX - 1; i++) {
					for (int j = 1; j < fNumY - 1; j++) {

						if (cellType[i * n + j] != FLUID_CELL)
							continue;

						int center = i * n + j;
						int left = (i - 1) * n + j;
						int right = (i + 1) * n + j;
						int bottom = i * n + j - 1;
						int top = i * n + j + 1;

						float s = this->s[center];
						float sx0 = this->s[left];
						float sx1 = this->s[right];
						float sy0 = this->s[bottom];
						float sy1 = this->s[top];
						s = sx0 + sx1 + sy0 + sy1;
						if (s == 0.0) // If none are blocks that could hold fluid
							continue;

						float div = u[right] - u[center] +
							v[top] - v[center];

						if (particleRestDensity > 0.0 && compensateDrift) {
							float k = 1.0;
							float compression = particleDensity[i * n + j] - particleRestDensity;
							if (compression > 0.0)
								div = div - k * compression;
						}

						float p = -div / s;
						p *= overRelaxation;
						this->p[center] += cp * p;

						this->u[center] -= sx0 * p;
						this->u[right] += sx1 * p;
						this->v[center] -= sy0 * p;
						this->v[top] += sy1 * p;
					}
				}
			}
		}

		/// <summary>
		/// Updates colours of particles based on density of local particle density
		/// </summary>
		void updateParticleColors()
		{
			// for (var i = 0; i < numParticles; i++) {
			// 	particleColor[3 * i] *= 0.99; 
			// 	particleColor[3 * i + 1] *= 0.99
			// 	particleColor[3 * i + 2] = 
			// 		clamp(particleColor[3 * i + 2] + 0.001, 0.0, 1.0)
			// }

			// return;

			float h1 = fInvSpacing;

			for (int i = 0; i < numParticles; i++) {

				float s = 0.01;

				particleColor[3 * i] = clamp(particleColor[3 * i] - s, 0.0, 1.0);
				particleColor[3 * i + 1] = clamp(particleColor[3 * i + 1] - s, 0.0, 1.0);
				particleColor[3 * i + 2] = clamp(particleColor[3 * i + 2] + s, 0.0, 1.0);

				float x = particlePos[2 * i];
				float y = particlePos[2 * i + 1];
				int xi = clamp(floor(x * h1), 1, fNumX - 1);
				int yi = clamp(floor(y * h1), 1, fNumY - 1);
				int cellNr = xi * fNumY + yi;

				float d0 = particleRestDensity;

				if (d0 > 0.0) {
					float relDensity = particleDensity[cellNr] / d0;
					if (relDensity < 0.7) {
						float s = 0.8;
						particleColor[3 * i] = s;
						particleColor[3 * i + 1] = s;
						particleColor[3 * i + 2] = 1.0;
					}
				}
			}
		}

		/// <summary>
		/// Helper Function: Updates the colour of the specific cell based on val
		/// </summary>
		/// <param name="cellNr"></param>
		/// <param name="val"></param>
		/// <param name="minVal"></param>
		/// <param name="maxVal"></param>
		void setSciColor(int cellNr, float val, float minVal, float maxVal)
		{
			val = min(max(val, minVal), maxVal - 0.0001);
			float d = maxVal - minVal;
			val = d == 0.0 ? 0.5 : (val - minVal) / d;
			float m = 0.25;
			int num = floor(val / m);
			float s = (val - num * m) / m;
			float r, g, b;

			switch (num) {
			case 0: r = 0.0; g = s; b = 1.0; break;
			case 1: r = 0.0; g = 1.0; b = 1.0 - s; break;
			case 2: r = s; g = 1.0; b = 0.0; break;
			case 3: r = 1.0; g = 1.0 - s; b = 0.0; break;
			}

			cellColor[3 * cellNr] = r;
			cellColor[3 * cellNr + 1] = g;
			cellColor[3 * cellNr + 2] = b;
		}

		/// <summary>
		/// Updates the colour of the cells
		/// </summary>
		void updateCellColors()
		{
			std::fill(cellColor.begin(), cellColor.end(), 0.0f);

			for (int i = 0; i < fNumCells; i++) {

				if (cellType[i] == SOLID_CELL) {
					cellColor[3 * i] = 0.5;
					cellColor[3 * i + 1] = 0.5;
					cellColor[3 * i + 2] = 0.5;
				}
				else if (cellType[i] == FLUID_CELL) {
					float d = particleDensity[i];
					if (particleRestDensity > 0.0)
						d /= particleRestDensity;
					setSciColor(i, d, 0.0, 2.0);
				}
			}
		}

		// Simulates the next simulation timestep
		void simulate(float dt, float gravity, float flipRatio, int numPressureIters, int numParticleIters, float overRelaxation, bool compensateDrift, bool separateParticles, float obstacleX, float abstacleY, float obstacleRadius, 
			bool cohesionOn, float cMaxAccel, float cMaxDistance, bool fastCohesionSimulate)
		{

			std::chrono::high_resolution_clock::time_point startTime = std::chrono::high_resolution_clock::now();


			int numSubSteps = 1;
			float sdt = dt / numSubSteps;

			for (int step = 0; step < numSubSteps; step++) {

				if (cohesionOn) addCohesionAccel(sdt, cMaxAccel, cMaxDistance, fastCohesionSimulate);
				integrateParticles(sdt, gravity);
				if (separateParticles) pushParticlesApart(numParticleIters);
				handleParticleCollisions(obstacleX, abstacleY, obstacleRadius);
				transferVelocities(true, 999.9f);
				updateParticleDensity();
				solveIncompressibility(numPressureIters, sdt, overRelaxation, compensateDrift);
				transferVelocities(false, flipRatio);
			}

			updateParticleColors();
			updateCellColors();


			std::chrono::high_resolution_clock::time_point endTime = std::chrono::high_resolution_clock::now();
			previousFrameTime = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime).count();

		}





};