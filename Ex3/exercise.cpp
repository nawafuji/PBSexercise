#include "gauss_seidel.h"
#include "fluid2d.h"

#define F_IDX(x,y) ((x) + (y) * _xRes)

// Problem 1
void ExSolvePoisson(int _xRes, int _yRes, int _iterations, double _accuracy, double* _field, double* _b)
{
	// Note that the boundaries are handles by the framework, so you iterations should be similar to:
  // Use the F_IDX macro to access the input fields _field and _b
  // _field[F_IDX(x, y)] ... _b[F_IDX(x, y)]

  for (int it = 0; it < _iterations; it++) {
    for (int y = 1; y < _yRes - 1; y++)
    {
      for (int x = 1; x < _xRes - 1; x++)
      {
        _field[F_IDX(x,y)] = (_b[F_IDX(x,y)] + _field[F_IDX(x+1,y)] + _field[F_IDX(x-1,y)] + _field[F_IDX(x,y+1)] + _field[F_IDX(x,y-1)])/4.0;
      }
    }

    float residual = 0;
    for (int y = 1; y < _yRes - 1; y++)
    {
      for (int x = 1; x < _xRes - 1; x++)
      {
        residual += abs(_b[F_IDX(x,y)] - (4.0*_field[F_IDX(x,y)] - _field[F_IDX(x+1,y)] - _field[F_IDX(x-1,y)] - _field[F_IDX(x,y+1)] - _field[F_IDX(x,y-1)]));
      }
    }
    residual = sqrt(residual)/(_xRes-2)/(_yRes-2);

	// For your debugging, and ours, please add these prints after every iteration
    if(it == _iterations - 1) 
    {
      printf("Pressure solver: it=%d , res=%f \n", it, residual);
    }

    if(residual < _accuracy)
    {
      printf("Pressure solver: it=%d , res=%f, converged \n", it, residual);
      break;
    }
  }
}

// Problem 2
void ExCorrectVelocities(int _xRes, int _yRes, double _dt, const double *_pressure, double *_xVelocity, double *_yVelocity)
{
	// Note: velocity u_{i+1/2} is practically stored at i+1
  for (int y = 1; y < _yRes-1; y++)
  {
    for (int x = 1; x < _xRes-1; x++)
    {
      _xVelocity[F_IDX(x,y)] -= _dt*(_pressure[F_IDX(x,y)] - _pressure[F_IDX(x-1,y)])*_xRes;
      _yVelocity[F_IDX(x,y)] -= _dt*(_pressure[F_IDX(x,y)] - _pressure[F_IDX(x,y-1)])*_yRes;
    }
  }
}

// Problem 3
void ExAdvectWithSemiLagrange(int _xRes, int _yRes, double _dt, double *_xVelocity, double *_yVelocity, double *_density, double *_densityTemp, double *_xVelocityTemp, double *_yVelocityTemp)
{
	// Note: velocity u_{i+1/2} is practically stored at i+1
  for (int y = 1; y < _yRes-1; y++)
  {
    for (int x = 1; x < _xRes-1; x++)
    {
      //for xVelocity
      {
        double bx = x - _xVelocity[F_IDX(x,y)]*_dt*_xRes;
        double by = y - (_yVelocity[F_IDX(x,y)]+_yVelocity[F_IDX(x,y+1)]+_yVelocity[F_IDX(x-1,y)]+_yVelocity[F_IDX(x-1,y+1)])/4*_dt*_yRes;
        if(bx<1)bx=1;
        if(bx>_xRes-2.5)bx=_xRes-2;
        if(by<1)by=1;
        if(by>_yRes-2)by=_yRes-2;
        // if(bx<1.5)bx=1.5;
        // if(bx>_xRes-2.5)bx=_xRes-2.5;
        // if(by<1.5)by=1.5;
        // if(by>_yRes-2.5)by=_yRes-2.5;
        int ax = (int)bx;
        int ay = (int)by;
        bx -= ax;
        by -= ay;
        _xVelocityTemp[F_IDX(x,y)] = (1-by)*((1-bx)*_xVelocity[F_IDX(ax,ay)] + bx*_xVelocity[F_IDX(ax+1,ay)]) + by*((1-bx)*_xVelocity[F_IDX(ax,ay+1)] + bx*_xVelocity[F_IDX(ax+1,ay+1)]);

      }
      //for yVelocity
      {
        double bx = x - (_xVelocity[F_IDX(x,y)]+_xVelocity[F_IDX(x+1,y)]+_xVelocity[F_IDX(x,y-1)]+_xVelocity[F_IDX(x+1,y-1)])/4*_dt*_xRes;
        double by = y - _yVelocity[F_IDX(x,y)]*_dt*_yRes;
        if(bx<1)bx=1;
        if(bx>_xRes-2.5)bx=_xRes-2;
        if(by<1)by=1;
        if(by>_yRes-2)by=_yRes-2;
        // if(bx<1.5)bx=1.5;
        // if(bx>_xRes-2.5)bx=_xRes-2.5;
        // if(by<1.5)by=1.5;
        // if(by>_yRes-2.5)by=_yRes-2.5;
        int ax = (int)bx;
        int ay = (int)by;
        bx -= ax;
        by -= ay;
        _yVelocityTemp[F_IDX(x,y)] = (1-by)*((1-bx)*_yVelocity[F_IDX(ax,ay)] + bx*_yVelocity[F_IDX(ax+1,ay)]) + by*((1-bx)*_yVelocity[F_IDX(ax,ay+1)] + bx*_yVelocity[F_IDX(ax+1,ay+1)]);
      }
      //for density
      {
        double bx = x - (_xVelocity[F_IDX(x,y)]+_xVelocity[F_IDX(x+1,y)])/2*_dt*_xRes;
        double by = y - (_yVelocity[F_IDX(x,y)]+_yVelocity[F_IDX(x,y+1)])/2*_dt*_yRes;
        if(bx<1)bx=1;
        if(bx>_xRes-2)bx=_xRes-2;
        if(by<1)by=1;
        if(by>_yRes-2)by=_yRes-2;
        int ax = (int)bx;
        int ay = (int)by;
        bx -= ax;
        by -= ay;
      _densityTemp[F_IDX(x,y)] = (1-by)*((1-bx)*_density[F_IDX(ax,ay)] + bx*_density[F_IDX(ax+1,ay)]) + by*((1-bx)*_density[F_IDX(ax,ay+1)] + bx*_density[F_IDX(ax+1,ay+1)]);
      }
    }
  }
  for (int y = 1; y < _yRes-1; y++)
  {
    for (int x = 1; x < _xRes-1; x++)
    {
      _xVelocity[F_IDX(x,y)] = _xVelocityTemp[F_IDX(x,y)];
      _yVelocity[F_IDX(x,y)] = _yVelocityTemp[F_IDX(x,y)];
      _density[F_IDX(x,y)] = _densityTemp[F_IDX(x,y)];
    }
  }
}
