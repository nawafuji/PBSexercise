#include "Vec2.h"

// gravitational acceleration (9.81)
static const double g = 9.81;


// Exercise 1
// hanging mass point
void AdvanceTimeStep1(double k, double m, double d, double L, double dt, int method, double p1, double v1, double& p2, double& v2)
{
  double F=0.,k1=0.,k2=0.,l1=0.;
    //Euler
    if (method==1)
    {
        
        F=-k*(abs(p2-p1)-L)*(p2-p1)/abs(p2-p1);        //not sure if this is correct
        p2 = p2+dt*v2;
        v2= v2+ dt*(F-m*g-v2*d)/m;
       
        
        
    }else if(method==3) //midpoint
    {
        k1=v2;
        F=-k*(abs(p2-p1)-L)*(p2-p1)/abs(p2-p1);     //not sure if this is correct
        l1=(F-m*g-v2*d)/m;
        k2=v2+dt*l1;
        
        p2=p2+dt/2*(k1+k2);
        
        v2=v2+dt/2*(l1+(-k*(p2+v2*dt)-d*k2)/m);
    }
}


// Exercise 3
// falling triangle
void AdvanceTimeStep3(double k, double m, double d, double L, double dt,
                      Vec2& p1, Vec2& v1, Vec2& p2, Vec2& v2, Vec2& p3, Vec2& v3)
{
	p1 += Vec2(1,1);
}
