#include "gauss_seidel.h"
#include "fluid2d.h"

#define F_IDX(x,y) ((x) + (y) * _xRes)

// Problem 1
void ExSolvePoisson(int _xRes, int _yRes, int _iterations, double _accuracy, double* _field, double* _b)
{
    double residual,err;
    
    for(int it=0; it<_iterations;++it){
    
	// Note that the boundaries are handles by the framework, so you iterations should be similar to:
	for (int y = 1; y < _yRes - 1; y++)
	{
		for (int x = 1; x < _xRes - 1; x++)
		{
			// Use the F_IDX macro to access the input fields _field and _b
			// _field[F_IDX(x, y)] ... _b[F_IDX(x, y)]
            _field[F_IDX(x, y)]= 1./4.* ( _b[F_IDX(x, y)]+_field[F_IDX(x-1, y)]+ _field[F_IDX(x, y-1)]+ _field[F_IDX(x+1, y)]+ _field[ F_IDX(x, y+1)] );
            
		}
	}
        
  /*
        
        if(it==_iterations-1)
            printf("Pressure solver: iter=%d , res=%f \n",it, residual);
        if(residual<_accuracy) {
            printf("Pressure solver: iter=%d , converged \n",it);
            break; // optional
  
        }*/
    }
	
}

// Problem 2
void ExCorrectVelocities(int _xRes, int _yRes, double _dt, const double *_pressure, double *_xVelocity, double *_yVelocity)
{
	// Note: velocity u_{i+1/2} is practically stored at i+1
    for (int y = 1; y < _yRes - 1; y++)
    {
        for (int x = 1; x < _xRes - 1; x++)
        {
            _xVelocity[F_IDX(x, y)] = _xVelocity[F_IDX(x, y)] - _dt*_xRes*(_pressure[F_IDX(x, y)]-_pressure[F_IDX(x-1, y)]);
            _yVelocity[F_IDX(x, y)] = _yVelocity[F_IDX(x, y)] - _dt*_yRes*(_pressure[F_IDX(x, y)]-_pressure[F_IDX(x, y-1)]);
        }
    }
    
}

// Problem 3
void ExAdvectWithSemiLagrange(int _xRes, int _yRes, double _dt, double *_xVelocity, double *_yVelocity, double *_density, double *_densityTemp, double *_xVelocityTemp, double *_yVelocityTemp)
{
	// Note: velocity u_{i+1/2} is practically stored at i+1
    double x1,x2,y1,y2;
    int i1,i2,j1,j2;
    
    double f1,f2,f3,f4;
    
    double  xCentral, yCentral;
    
    for (int y = 1; y< _yRes-1; ++y) {
        for (int x = 1 ; x< _xRes-1; ++x) {
            
            
            //backtracing velocities
            xCentral = (double)1./_xRes*(x+0.5)-_dt*0.5*(_xVelocity[F_IDX(x, y)]+_xVelocity[F_IDX(x+1, y)]);
            
            yCentral = (double)1./_yRes*(y+0.5)-_dt*0.5*(_yVelocity[F_IDX(x, y)]+_yVelocity[F_IDX(x, y+1)]);
            
            
            i1=xCentral*_xRes+0.5;
            i2=xCentral*_xRes-0.5;
            j1=yCentral*_yRes+0.5;
            j2=yCentral*_yRes-0.5;
            
            x1=(double)(i1+0.5)/_xRes;
            y1=(double)(j1+0.5)/_yRes;
            x2=(double)(i2+0.5)/_xRes;;
            y2=(double)(j2+0.5)/_yRes;;
            
            
            f1=_density[F_IDX(i1, j1)]*(x2-xCentral)*(y2-yCentral);
            f2=_density[F_IDX(i2, j1)]*(xCentral-x1)*(y2-yCentral);
            f3=_density[F_IDX(i1, j2)]*(x2-xCentral)*(yCentral-y1);
            f4=_density[F_IDX(i2, j2)]*(xCentral-x1)*(yCentral-y1);
            
            _densityTemp[F_IDX(x, y)]= (f1+f2+f3+f4)*(_yRes*_xRes);
            
            
            //************************************************
            
            //find x_p with backtracing and interpolation of the velocities at the center
            xCentral=(double) (x+0.5)/_xRes-_dt*(_xVelocity[F_IDX(x, y)]+_xVelocity[F_IDX(x+1, y)]+_xVelocity[F_IDX(x, y-1)]+_xVelocity[F_IDX(x+1, y-1)])/4.;
            yCentral=(double) (y)/_yRes-_dt*_yVelocity[F_IDX(x, y)];
            
            
            i1=xCentral*_xRes+0.5;
            i2=xCentral*_xRes-0.5;
            j1=yCentral*_yRes+1.;
            j2=yCentral*_yRes;
            
            x1=(double)(i1+0.5)/_xRes;
            y1=(double)(j1)/_yRes;
            x2=(double)(i2+0.5)/_xRes;;
            y2=(double)(j2)/_yRes;;
            
            f1=_yVelocity[F_IDX(i1, j1)]*(x2-xCentral)*(y2-yCentral);
            f2=_yVelocity[F_IDX(i2, j1)]*(xCentral-x1)*(y2-yCentral);
            f3=_yVelocity[F_IDX(i1, j2)]*(x2-xCentral)*(yCentral-y1);
            f4=_yVelocity[F_IDX(i2, j2)]*(xCentral-x1)*(yCentral-y1);
            
            _yVelocityTemp[F_IDX(x, y)]= (f1+f2+f3+f4)*(_xRes*_yRes);//((x2-x1)*(y2-y1));
            
            //************************************************

            yCentral=(double) (y+0.5)/_yRes-_dt*(_yVelocity[F_IDX(x, y)]+_yVelocity[F_IDX(x-1, y)]+_yVelocity[F_IDX(x, y+1)]+_yVelocity[F_IDX(x-1, y+1)])/4.;
            xCentral= (double) (x)/_xRes-_dt*_xVelocity[F_IDX(x, y)];
            
            i1=xCentral*_xRes+1;
            i2=xCentral*_xRes;
            j1=yCentral*_yRes+0.5;
            j2=yCentral*_yRes-0.5;
 
            x1=(double)(i1)/_xRes;
            y1=(double)(j1+0.5)/_yRes;
            x2=(double)(i2)/_xRes;;
            y2=(double)(j2+0.5)/_yRes;;
            
            f1=_xVelocity[F_IDX(i1, j1)]*(x2-xCentral)*(y2-yCentral);
            f2=_xVelocity[F_IDX(i2, j1)]*(xCentral-x1)*(y2-yCentral);
            f3=_xVelocity[F_IDX(i1, j2)]*(x2-xCentral)*(yCentral-y1);
            f4=_xVelocity[F_IDX(i2, j2)]*(xCentral-x1)*(yCentral-y1);
            
            _xVelocityTemp[F_IDX(x, y)]= (f1+f2+f3+f4)*(_yRes*_xRes);
            
            //************************************************

  
            
        }
    }
    
    
  
   for (int y = 1; y< _yRes-1; ++y) {
        for (int x = 1 ; x< _xRes-1; ++x) {
            _xVelocity[F_IDX(x, y)]=_xVelocityTemp[F_IDX(x, y)];
            _yVelocity[F_IDX(x, y)]=_yVelocityTemp[F_IDX(x, y)];
            _density[F_IDX(x, y)]=_densityTemp[F_IDX(x, y)];
        }
    }
  

}
