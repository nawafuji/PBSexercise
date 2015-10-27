#include "gauss_seidel.h"
#include "fluid2d.h"

#define F_IDX(x,y) ((x) + (y) * _xRes)

// Problem 1
void ExSolvePoisson(int _xRes, int _yRes, int _iterations, double _accuracy, double* _field, double* _b)
{
	// Note that the boundaries are handles by the framework, so you iterations should be similar to:
	for (int y = 1; y < _yRes - 1; y++)
	{
		for (int x = 1; x < _xRes - 1; x++)
		{
			// Use the F_IDX macro to access the input fields _field and _b
			// _field[F_IDX(x, y)] ... _b[F_IDX(x, y)]
		}
	}
	
	// For your debugging, and ours, please add these prints after every iteration
	//if(it == _iterations - 1) 
	//	printf("Pressure solver: it=%d , res=%f \n", it, residual);
	//if(residual < _accuracy)
	//	printf("Pressure solver: it=%d , res=%f, converged \n", it, residual);
}

// Problem 2
void ExCorrectVelocities(int _xRes, int _yRes, double _dt, const double *_pressure, double *_xVelocity, double *_yVelocity)
{
	// Note: velocity u_{i+1/2} is practically stored at i+1
}

// Problem 3
void ExAdvectWithSemiLagrange(int _xRes, int _yRes, double _dt, double *_xVelocity, double *_yVelocity, double *_density, double *_densityTemp, double *_xVelocityTemp, double *_yVelocityTemp)
{
	// Note: velocity u_{i+1/2} is practically stored at i+1
}
