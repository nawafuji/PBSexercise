#include "Vec2.h"
#include <cmath>
using namespace std;
// gravitational acceleration (9.81)
static const double g = 9.81;


// Exercise 1
// hanging mass point
void AdvanceTimeStep1(double k, double m, double d, double L, double dt, int method, double p1, double v1, double& p2, double& v2)
{
  double F=0.,k1=0.,k2=0.,l1=0.;
    F=-k*(abs(p2-p1)-L)*(p2-p1)/abs(p2-p1);
    
    //Euler
    if (method==1)
    {
        p2 = p2+dt*v2;
        v2= v2+ dt*(F-m*g-v2*d)/m;
        
    }else if(method==2)     //symplectic euler
    {
        
        v2=v2+ dt*(F-m*g -v2*d)/m;
        p2= p2 +v2*dt;
    }
    else if(method==3) //midpoint
    {
        k1=v2;
        
        l1=(F-m*g-v2*d)/m;
        k2=v2+dt*l1;
        
        p2=p2+dt/2*(k1+k2);
        
        v2=v2+dt/2*(l1+(-k*(p2+v2*dt)-d*k2)/m);
    }else if (method==4)        //semi-implicit Euler
    {
        v2= ((m-dt*(-d))*v2 +(F-m*g-v2*d)*dt)/(m-dt*(-d)-dt*dt*(-k));
        p2=p2+dt*v2;
    }else if(method==5)		//analytic
    {
        double a=-d/2*m,b=sqrt(4*k*m-d*d)/(2*m);
        double c1=m*g/k,c2=-a/b*c1;
        
        double p2_0=p2;
        p2=c1*exp(a*dt)*cos(b*dt)+c2*exp(a*dt)*sin(b*dt)-L-m*g/k;
        v2=c1*exp(a*dt)*(a*cos(b*dt)-sin(b*dt)*b)+c2*exp(a*dt)*(a*sin(b*dt)+cos(b*dt)*b);
    }
}


// Exercise 3
// falling triangle
void AdvanceTimeStep3(double k, double m, double d, double L, double dt,
                      Vec2& p1, Vec2& v1, Vec2& p2, Vec2& v2, Vec2& p3, Vec2& v3)
{
	p1 += Vec2(1,1);
}
