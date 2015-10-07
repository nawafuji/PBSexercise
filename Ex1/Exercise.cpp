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

  switch(method){
    case 1:
        p2 = p2+dt*v2;
        v2= v2+ dt*(F-m*g-v2*d)/m;
        break;
        
    case 2:
        v2=v2+ dt*(F-m*g -v2*d)/m;
        p2= p2 +v2*dt;
        break;

    case 3:
        k1=v2;
        
        l1=(F-m*g-v2*d)/m;
        k2=v2+dt*l1;
        
        p2=p2+dt/2*(k1+k2);
        
        v2=v2+dt/2*(l1+(-k*(p2+v2*dt)-d*k2)/m);
        break;

    case 4:
        v2= ((m-dt*(-d))*v2 +(F-m*g-v2*d)*dt)/(m-dt*(-d)-dt*dt*(-k));
        p2=p2+dt*v2;
        break;

    case 5:
        double a=-d/2*m,b=sqrt(4*k*m-d*d)/(2*m);
        double c1=m*g/k,c2=-a/b*c1;
        
        double p2_0=p2;
        p2=c1*exp(a*dt)*cos(b*dt)+c2*exp(a*dt)*sin(b*dt)-L-m*g/k;
        v2=c1*exp(a*dt)*(a*cos(b*dt)-sin(b*dt)*b)+c2*exp(a*dt)*(a*sin(b*dt)+cos(b*dt)*b);
        break;

    }

}


// Exercise 3
// falling triangle
void AdvanceTimeStep3(double k, double m, double d, double L, double dt,
                      Vec2& p1, Vec2& v1, Vec2& p2, Vec2& v2, Vec2& p3, Vec2& v3)
{
	//p1 += Vec2(1,1);
  int kr = 100;
  Vec2 F1 = Vec2::ZERO;
  Vec2 F2 = Vec2::ZERO;
  Vec2 F3 = Vec2::ZERO;

  F1=(k*((p2-p1).length()-L)/(p2-p1).length())*(p2-p1)+(k*((p3-p1).length()-L)/(p3-p1).length())*(p3-p1);
  F2=(k*((p1-p2).length()-L)/(p1-p2).length())*(p1-p2)+(k*((p3-p2).length()-L)/(p3-p2).length())*(p3-p2);
  F3=(k*((p2-p3).length()-L)/(p2-p3).length())*(p2-p3)+(k*((p1-p3).length()-L)/(p1-p3).length())*(p1-p3);
  Vec2 G = Vec2(0, -g);

  Vec2 Fr1=kr*Vec2(0,-1-p1.y);
  Vec2 Fr2=kr*Vec2(0,-1-p2.y);
  Vec2 Fr3=kr*Vec2(0,-1-p3.y);

  if(p1.y<-1){
    v1= v1+ dt/m*(F1+m*G-d*v1+Fr1);
  }
  else{
    v1= v1+ dt/m*(F1+m*G-d*v1);
  }
  if(p2.y<-1){
    v2= v2+ dt/m*(F2+m*G-d*v2+Fr2);
  }
  else{
    v2= v2+ dt/m*(F2+m*G-d*v2);
  }
  if(p3.y<-1){
    v3= v3+ dt/m*(F3+m*G-d*v3+Fr3);
  }
  else{
    v3= v3+ dt/m*(F3+m*G-d*v3);
  }
  p1 = p1+dt*v1;
  p2 = p2+dt*v2;
  p3 = p3+dt*v3;

}
