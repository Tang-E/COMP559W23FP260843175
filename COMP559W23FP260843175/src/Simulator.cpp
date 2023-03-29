/*
* Simulator Code in charge of the particle and grid simulations.
* Code here is adapted from Matthias M�ller's tenMinutePhysics
* FLIP Simulator code available on his GitHub repository here:
* https://github.com/matthias-research/pages/blob/master/tenMinutePhysics/18-flip.html
* 
* I am adapting his code written in javascript into C++ code
* because I don't want to write in javascript.
* 
* @author Edwin Pan (260843175) for COMP559 Winter 2023 Final Project
* 
* Extra notes:
*	I use memset() to fill arrays of double with zeros. It is not guaranteed that
*	the float and double formatting produces 0 when all bits are false - and according
*	to here https://stackoverflow.com/questions/4629853/is-it-legal-to-use-memset-0-on-an-array-of-doubles
*	it is only true of "__STDC_IEC_559__" is defined. So we start with that.
*/


#include <iostream>
#include <cmath>
#include <cstring>
#include <algorithm>

#include "scene.cpp"
#include "HelperFuncs.cpp"



class FlipFluid {

	// Constants
	const int U_FIELD = 0;
	const int V_FIELD = 1;
	const int FLUID_CELL = 0;
	const int AIR_CELL = 1;
	const int SOLID_CELL = 2;

	// Fluid Properties

	double density;
	int fNumX;
	int fNumY;
	double h;
	double fInvSpacing;
	int fNumCells;

	double* u;
	double* v;
	double* du;
	double* dv;
	double* prevU;
	double* prevV;
	double* p;
	double* s;
	int* cellType;
	float* cellColor;

	// Particle Properties

	int maxParticles;

	double* particlePos;
	float* particleColor;
	double* particleVel;
	double* particleDensity;
	double particleRestDensity;

	double particleRadius;
	double pInvSpacing;
	int pNumX;
	int pNumY;
	int pNumCells;

	int* numCellParticles;
	int* firstCellParticle;
	int* cellParticleIds;

	int numParticles;

	FlipFluid(double density, int width, int height, double spacing, double particleRadius, int maxParticles) {

		// Fluid Properties

		this->density = density;
		this->fNumX = floor(width / spacing) + 1;
		this->fNumY = floor(height / spacing) + 1;
		this->h = std::fmax(width/this->fNumX, height/this->fNumY);
		this->fInvSpacing = 1.0 / this->h;
		this->fNumCells = this->fNumX * this->fNumY;

		this->u = new double[fNumCells];
		this->v = new double[fNumCells];
		this->du = new double[fNumCells];
		this->dv = new double[fNumCells];
		this->prevU = new double[fNumCells];
		this->prevV = new double[fNumCells];
		this->p = new double[fNumCells];
		this->s = new double[fNumCells];
		this->cellType = new int[fNumCells];
		this->cellColor = new float[3 * fNumCells];

		// Particle Properties

		this->maxParticles = maxParticles;

		this->particlePos = new double[this->maxParticles];
		this->particleColor = new float[3 * this->maxParticles];
		for (int i = 0; i < this->maxParticles; i++) {
			this->particleColor[3 * i + 2] = 1.0; // Make blue, not white!
		}
		this->particleVel = new double[2 * this->maxParticles];
		this->particleDensity = new double[this->fNumCells];
		this->particleRestDensity = 0.0;

		this->particleRadius = particleRadius;
		this->pInvSpacing = 1.0 / (2.2 * particleRadius);
		this->pNumX = floor(width * this->pInvSpacing) + 1;
		this->pNumY = floor(height * this->pInvSpacing) + 1;
		this->pNumCells = this->pNumX * this->pNumY;

		this->numCellParticles = new int[this->pNumCells];
		this->firstCellParticle = new int[this->pNumCells + 1];
		this->cellParticleIds = new int[this->maxParticles];

		this->numParticles = 0;

		// Check for dependency regarding memset(), used in functions
		#ifndef __STDC_IEC_559__ 
		std::cout << "ALERT! float and double types do not match IEC60599 single format and memset() operations will fail!" << std::endl;
		#endif

	}

	~FlipFluid() {
		delete[] u;
		delete[] v;
		delete[] du;
		delete[] dv;
		delete[] prevU;
		delete[] prevV;
		delete[] p;
		delete[] s;
		delete[] cellType;
		delete[] cellColor;
		delete[] particlePos;
		delete[] particleColor;
		delete[] particleVel;
		delete[] particleDensity;
		delete[] numCellParticles;
		delete[] firstCellParticle;
		delete[] cellParticleIds;
	}

	/*
	* Simplecton integrate particles with only gravity as force.
	*/
	void integrateParticles(double dt, double gravity)
	{
		for (int i = 0; i < numParticles; i++) {
			particleVel[2 * i + 1] += dt * gravity;
			particlePos[2 * i] += particleVel[2 * i] * dt;
			particlePos[2 * i + 1] += particleVel[2 * i + 1] * dt;
		}
	}



	/*
	* Applies repulsion between particles efficiently via Gauss-Seidel
	* Complicated algorithm that involves filtering out collision checks
	* of particles that would be too far from each other. Said filtering is
	* something Matthias M�ller has figured out and I'm not going to bother
	* trying to reconstruct it.
	*/
	void pushParticlesApart(int numIters)
	{
		double colorDiffusionCoeff = 0.001;

		// count particles per cell

		std::memset(numCellParticles, 0, pNumCells * sizeof(double));

		for (int i = 0; i < numParticles; i++) {
			double x = particlePos[2 * i];
			double y = particlePos[2 * i + 1];

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
			double x = particlePos[2 * i];
			double y = particlePos[2 * i + 1];

			int xi = clamp(floor(x * pInvSpacing), 0, pNumX - 1);
			int yi = clamp(floor(y * pInvSpacing), 0, pNumY - 1);
			int cellNr = xi * pNumY + yi;
			firstCellParticle[cellNr]--;
			cellParticleIds[firstCellParticle[cellNr]] = i;
		}

		// push particles apart

		double minDist = 2.0 * particleRadius;
		double minDist2 = minDist * minDist;

		for (int iter = 0; iter < numIters; iter++) {

			for (int i = 0; i < numParticles; i++) {
				double px = particlePos[2 * i];
				double py = particlePos[2 * i + 1];

				double pxi = floor(px * pInvSpacing);
				double pyi = floor(py * pInvSpacing);
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
							double qx = particlePos[2 * id];
							double qy = particlePos[2 * id + 1];

							double dx = qx - px;
							double dy = qy - py;
							double d2 = dx * dx + dy * dy;
							if (d2 > minDist2 || d2 == 0.0)
								continue;
							double d = sqrt(d2);
							double s = 0.5 * (minDist - d) / d;
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

	/*
	* Handles particle collisions with the obstacle
	*/
	void handleParticleCollisions(double obstacleX, double obstacleY, double obstacleRadius)
	{
		double h = 1.0 / fInvSpacing;
		double r = particleRadius;
		double or1 = obstacleRadius;
		double or2 = or1 * or1 ;
		double minDist = obstacleRadius + r;
		double minDist2 = minDist * minDist;

		double minX = h + r;
		double maxX = (fNumX - 1) * h - r;
		double minY = h + r;
		double maxY = (fNumY - 1) * h - r;


		for (int i = 0; i < numParticles; i++) {
			double x = particlePos[2 * i];
			double y = particlePos[2 * i + 1];

			double dx = x - obstacleX;
			double dy = y - obstacleY;
			double d2 = dx * dx + dy * dy;

			// obstacle collision

			if (d2 < minDist2) {

				// var d = Math.sqrt(d2);
				// var s = (minDist - d) / d;
				// x += dx * s;
				// y += dy * s;

				particleVel[2 * i] = scene.obstacleVelX;
				particleVel[2 * i + 1] = scene.obstacleVelY;
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

	/*
	* TODO: Investigate commented out section of this code.
	* Recalculates particle rest density and then commentedly
	* would recalculate particle densities?
	*/
	void updateParticleDensity()
	{
		int n = fNumY;
		double h = h;
		double h1 = fInvSpacing;
		double h2 = 0.5 * h;

		//double d = particleDensity;
		//d.fill(0.0);
		double* d = new double[fNumCells];
		memset(d, 0, sizeof(double) * fNumCells);

		for (int i = 0; i < numParticles; i++) {
			double x = particlePos[2 * i];
			double y = particlePos[2 * i + 1];

			x = clamp(x, h, (fNumX - 1) * h);
			y = clamp(y, h, (fNumY - 1) * h);

			int x0 = floor((x - h2) * h1);
			double tx = ((x - h2) - x0 * h) * h1;
			int x1 = min(x0 + 1, fNumX - 2);

			int y0 = floor((y - h2) * h1);
			double ty = ((y - h2) - y0 * h) * h1;
			int y1 = min(y0 + 1, fNumY - 2);

			double sx = 1.0 - tx;
			double sy = 1.0 - ty;

			if (x0 < fNumX && y0 < fNumY) d[x0 * n + y0] += sx * sy;
			if (x1 < fNumX && y0 < fNumY) d[x1 * n + y0] += tx * sy;
			if (x1 < fNumX && y1 < fNumY) d[x1 * n + y1] += tx * ty;
			if (x0 < fNumX && y1 < fNumY) d[x0 * n + y1] += sx * ty;
		}

		if (particleRestDensity == 0.0) {
			double sum = 0.0;
			int numFluidCells = 0;

			for (int i = 0; i < fNumCells; i++) {
				if (cellType[i] == FLUID_CELL) {
					sum += d[i];
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

	/*
	* Velocity Transfer Stuff between Particles and Grid
	*/
	void transferVelocities(bool toGrid, float flipRatio)
	{
		int n = fNumY;
		double h = h;
		double h1 = fInvSpacing;
		double h2 = 0.5 * h;

		if (toGrid) {

			copy(u, prevU, fNumCells);
			copy(v, prevV, fNumCells);

			memset(du, 0, sizeof(double) * fNumCells);
			memset(dv, 0, sizeof(double) * fNumCells);
			memset(u, 0, sizeof(double) * fNumCells);
			memset(v, 0, sizeof(double) * fNumCells);

			for (int i = 0; i < fNumCells; i++)
				cellType[i] = s[i] == 0.0 ? SOLID_CELL : AIR_CELL;

			for (int i = 0; i < numParticles; i++) {
				double x = particlePos[2 * i];
				double y = particlePos[2 * i + 1];
				int xi = clamp(floor(x * h1), 0, fNumX - 1);
				int yi = clamp(floor(y * h1), 0, fNumY - 1);
				int cellNr = xi * n + yi;
				if (cellType[cellNr] == AIR_CELL)
					cellType[cellNr] = FLUID_CELL;
			}
		}

		for (int component = 0; component < 2; component++) {

			double dx = component == 0 ? 0.0 : h2;
			double dy = component == 0 ? h2 : 0.0;

			double* f = component == 0 ? u : v;
			double* prevF = component == 0 ? prevU : prevV;
			double* d = component == 0 ? du : dv;

			for (int i = 0; i < numParticles; i++) {
				double x = particlePos[2 * i];
				double y = particlePos[2 * i + 1];

				x = clamp(x, h, (fNumX - 1) * h);
				y = clamp(y, h, (fNumY - 1) * h);

				int x0 = min(floor((x - dx) * h1), fNumX - 2);
				double tx = ((x - dx) - x0 * h) * h1;
				int x1 = min(x0 + 1, fNumX - 2);

				int y0 = min(floor((y - dy) * h1), fNumY - 2);
				double ty = ((y - dy) - y0 * h) * h1;
				int y1 = min(y0 + 1, fNumY - 2);

				double sx = 1.0 - tx;
				double sy = 1.0 - ty;

				double d0 = sx * sy;
				double d1 = tx * sy;
				double d2 = tx * ty;
				double d3 = sx * ty;

				int nr0 = x0 * n + y0;
				int nr1 = x1 * n + y0;
				int nr2 = x1 * n + y1;
				int nr3 = x0 * n + y1;

				if (toGrid) {
					double  pv = particleVel[2 * i + component];
					f[nr0] += pv * d0;  d[nr0] += d0;
					f[nr1] += pv * d1;  d[nr1] += d1;
					f[nr2] += pv * d2;  d[nr2] += d2;
					f[nr3] += pv * d3;  d[nr3] += d3;
				}
				else {
					int offset = component == 0 ? n : 1;
					double valid0 = cellType[nr0] != AIR_CELL || cellType[nr0 - offset] != AIR_CELL ? 1.0 : 0.0;
					double valid1 = cellType[nr1] != AIR_CELL || cellType[nr1 - offset] != AIR_CELL ? 1.0 : 0.0;
					double valid2 = cellType[nr2] != AIR_CELL || cellType[nr2 - offset] != AIR_CELL ? 1.0 : 0.0;
					double valid3 = cellType[nr3] != AIR_CELL || cellType[nr3 - offset] != AIR_CELL ? 1.0 : 0.0;

					double v = particleVel[2 * i + component];
					double d = valid0 * d0 + valid1 * d1 + valid2 * d2 + valid3 * d3;

					if (d > 0.0) {

						double picV = (valid0 * d0 * f[nr0] + valid1 * d1 * f[nr1] + valid2 * d2 * f[nr2] + valid3 * d3 * f[nr3]) / d;
						double corr = (valid0 * d0 * (f[nr0] - prevF[nr0]) + valid1 * d1 * (f[nr1] - prevF[nr1])
							+ valid2 * d2 * (f[nr2] - prevF[nr2]) + valid3 * d3 * (f[nr3] - prevF[nr3])) / d;
						double flipV = v + corr;

						particleVel[2 * i + component] = (1.0 - flipRatio) * picV + flipRatio * flipV;
					}
				}
			}

			if (toGrid) {
				for (int i = 0; i < /*f.length*/ fNumCells; i++) {
					if (d[i] > 0.0)
						f[i] /= d[i];
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

	/*
	* Solves incompressibility through gauss-seidel by going over each
	* grid element and reducing divergence to 0 at every pass with the
	* ability to undershoot or overshoot with overRelaxation term.
	*/
	void solveIncompressibility(int numIters, double dt, double overRelaxation, bool compensateDrift = true) {

		memset(p, 0, sizeof(double) * fNumCells);
		copy(u, prevU, fNumCells);
		copy(v, prevV, fNumCells);

		int n = fNumY;
		double cp = density * h / dt;

		for (int i = 0; i < fNumCells; i++) { //wtf does this do for us??????
			double u = this->u[i];
			double v = this->v[i];
		}

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

					double s = this->s[center];
					double sx0 = this->s[left];
					double sx1 = this->s[right];
					double sy0 = this->s[bottom];
					double sy1 = this->s[top];
					s = sx0 + sx1 + sy0 + sy1;
					if (s == 0.0)
						continue;

					double div = u[right] - u[center] +
						v[top] - v[center];

					if (particleRestDensity > 0.0 && compensateDrift) {
						double k = 1.0;
						double compression = particleDensity[i * n + j] - particleRestDensity;
						if (compression > 0.0)
							div = div - k * compression;
					}

					double p = -div / s;
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

	/*
	* Updates colours of particles based on density of local particle density
	*/
	void updateParticleColors()
	{
		// for (var i = 0; i < numParticles; i++) {
		// 	particleColor[3 * i] *= 0.99; 
		// 	particleColor[3 * i + 1] *= 0.99
		// 	particleColor[3 * i + 2] = 
		// 		clamp(particleColor[3 * i + 2] + 0.001, 0.0, 1.0)
		// }

		// return;

		double h1 = fInvSpacing;

		for (int i = 0; i < numParticles; i++) {

			double s = 0.01;

			particleColor[3 * i] = clamp(particleColor[3 * i] - s, 0.0, 1.0);
			particleColor[3 * i + 1] = clamp(particleColor[3 * i + 1] - s, 0.0, 1.0);
			particleColor[3 * i + 2] = clamp(particleColor[3 * i + 2] + s, 0.0, 1.0);

			double x = particlePos[2 * i];
			double y = particlePos[2 * i + 1];
			int xi = clamp(floor(x * h1), 1, fNumX - 1);
			int yi = clamp(floor(y * h1), 1, fNumY - 1);
			int cellNr = xi * fNumY + yi;

			double d0 = particleRestDensity;

			if (d0 > 0.0) {
				double relDensity = particleDensity[cellNr] / d0;
				if (relDensity < 0.7) {
					double s = 0.8;
					particleColor[3 * i] = s;
					particleColor[3 * i + 1] = s;
					particleColor[3 * i + 2] = 1.0;
				}
			}
		}
	}

	/*
	* Helper Function: Updates the colour of the specific cell based on val
	*/
	void setSciColor(int cellNr, double val, double minVal, double maxVal)
	{
		val = min(max(val, minVal), maxVal - 0.0001);
		double d = maxVal - minVal;
		val = d == 0.0 ? 0.5 : (val - minVal) / d;
		double m = 0.25;
		int num = floor(val / m);
		double s = (val - num * m) / m;
		double r, g, b;

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

	/*
	* Updates the colour of the cells
	*/
	void updateCellColors()
	{
		memset(cellColor, 0, sizeof(float) * 3 * fNumCells);

		for (int i = 0; i < fNumCells; i++) {

			if (cellType[i] == SOLID_CELL) {
				cellColor[3 * i] = 0.5;
				cellColor[3 * i + 1] = 0.5;
				cellColor[3 * i + 2] = 0.5;
			}
			else if (cellType[i] == FLUID_CELL) {
				double d = particleDensity[i];
				if (particleRestDensity > 0.0)
					d /= particleRestDensity;
				setSciColor(i, d, 0.0, 2.0);
			}
		}
	}

	/*
	* Moves the FLIP water simulation one dt-seconds forward.
	*/
	void simulate(double dt, double gravity, float flipRatio, int numPressureIters, int numParticleIters, double overRelaxation, bool compensateDrift, bool separateParticles, double obstacleX, double abstacleY, double obstacleRadius)
	{
		int numSubSteps = 1;
		double sdt = dt / numSubSteps;

		for (int step = 0; step < numSubSteps; step++) {
			integrateParticles(sdt, gravity);
			if (separateParticles)
				pushParticlesApart(numParticleIters);
			handleParticleCollisions(obstacleX, abstacleY, obstacleRadius);
			transferVelocities(true, flipRatio);
			updateParticleDensity();
			solveIncompressibility(numPressureIters, sdt, overRelaxation, compensateDrift);
			transferVelocities(false, flipRatio);
		}

		updateParticleColors();
		updateCellColors();

	}





};