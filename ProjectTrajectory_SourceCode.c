#include <math.h> 
#include <stdio.h>
#define Msun 1.989e30 
#define Rsun 695700000.0 
#define G 6.67408e-11 
#define Mearth 5.972e24 
#define Rearth 6378000.0

double ax(double x, double y) //acceleration of earth in x direction
{
	double r=x*x+y*y; 
	double r1=r*sqrt(r);
	return (G*Msun*(-1.0)*x)/r1;
}

double ay(double x, double y) //acceleration of earth in y direction
{
	double r=x*x+y*y; 
	double r1=r*sqrt(r);
	return G*Msun*(-1.0)*y/r1;
}

double axm(double x, double y) //acceleration of moon in x direction
{
	double r=x*x+y*y; 
	double r1=r*sqrt(r);
	return G*Mearth*(-1.0)*x/r1;
}
double aym(double x, double y) //acceleration of moon in y direction
{
	double r=x*x+y*y; double r1=r*sqrt(r);
	return G*Mearth*(-1.0)*y/r1;
}
//_______________________________________________________________________							

double f(double x)
{
	return (-1.0)*G*Msun*Mearth/x;
}

double simpson13(double func(double x),double a,double b,double n)
{
	double x,h,I,S=0; 
	int i;
	h=fabs(b-a)/n;
	for(i=1;i<n;++i)
	{
		x=a+i*h; 
		if(i%2==0) 
			S=S+2*func(x);
		else
			S=S+4*func(x);
	}
	I=(h/3)*(func(a)+func(b)+S); 
	return I;
}

double trapezoidal(double func(double x),double a,double b,int n)
{
	int i;
	double x,h,I,S=0; 
	h=fabs(b-a)/n; 
	for(i=1;i<n;++i)
	{
		x=a+i*h; 
		S=S+func(x);
	}
	I=(h/2)*(func(a)+func(b)+2*S); 
	return I;
}
//_______________________________________________________________________							

double fn(double r,double E) //function of central potential
{
	double l=Mearth*Rearth*30300.0; //angular momentum
	//printf("%f	%f	%f\n",l,E,r);
	return (1/(r*r*sqrt(2*Mearth*E/l+2*Mearth*G*Msun/(l*r)-1/(r*r))));
}

double gauss(double f(double r,double E),double r,double E,double a, double b) //gaussian quadrature
{
	double x1, x2;
	x1=((b-a)/2.0)*(1/1.73)+((b+a)/2);
	x2=((b-a)/2.0)*(-1/1.73)+((b+a)/2);
	return (b-a)/a*(f(x1,E)+f(x2,E));
}
//_______________________________________________________________________							

void central() //evaluating the integral to get R vs theta graph
{ 
	FILE*fp=NULL;
	fp=fopen("central.txt","w"); 
	double E,r,r0,rm;
	printf("Enter the value of E:"); 
	scanf("%lf",&E);
	printf("Enter the lower and upper limit r0 & rm:\n"); 
	scanf("%lf%lf",&r0,&rm);
	for (r=r0;r<=rm;r=r+1.0) // varying R from r0 to rm
	{
		printf("%lf\t%lf\n",r,gauss(fn,r,E,r0,r));
		fprintf(fp,"%lf\t%lf\n",r,gauss(fn,r,E,r0,r));
	}
}
//_______________________________________________________________________							

double energy (double x,double y,double vx, double vy)
{
	double r=sqrt(x*x+y*y);
	return 0.5*(vx*vx+vy*vy)-G*Msun/r;
}
//_______________________________________________________________________							

double rk()
{ 
	FILE*fp=NULL;
	fp=fopen("trajectory.txt","w");
	FILE*fpm=NULL;
	fpm=fopen("t_moon.txt","w");
	FILE*fpe=NULL;
	fpe=fopen("energy.txt","w"); double xm=-384400000.0;
	double ym=0.0; double vxm=0.0; double vym=-1022.0;
	double k1xm,k2xm,k3xm,k4xm,k1ym,k2ym,k3ym,k4ym,k1vxm,k2vxm,k3vxm,k4vxm,k1vym,k2vym,k3vym,k4vym;
	double x=-147095000000.0;
	double y=0.0; double vx=0.0; double vy=-30300.0; double h =86400.0; int m=0;
	double k1x,k2x,k3x,k4x,k1y,k2y,k3y,k4y,k1vx,k2vx,k3vx,k4vx,k1vy,k2vy,k3vy,k4vy; 
	for(int n=0;n<=365;n++)
	{
	fprintf(fp,"%f	%f\n",x,y);
	{
	k1x=vx; k1y=vy; 
	k1vx=ax(x,y);
	k1vy=ay(x,y);
	k2x=vx+k1vx*h/2.0; 
	k2y=vy+k1vy*h/2.0;
	k2vx=ax(x+k1x*h/2.0,y+k1y*h/2.0); 
	k2vy=ay(x+k1x*h/2.0,y+k1y*h/2.0);
	k3x=vx+k2vx*h/2.0; k3y=vy+k2vy*h/2.0; 
	k3vx=ax(x+k2x*h/2.0,y+k2y*h/2.0); 
	k3vy=ay(x+k2x*h/2.0,y+k2y*h/2.0);
	k4x=vx+k3vx*h; k4y=vy+k3vy*h; 
	k4vx=ax(x+k3x*h,y+k3y*h); 
	k4vy=ay(x+k3x*h,y+k3y*h);
	x=x+(h/6.0)*(k1x+2.0*k2x+2.0*k3x+k4x); 
	y=y+(h/6.0)*(k1y+2.0*k2y+2.0*k3y+k4y);
	vx=vx+(h/6.0)*(k1vx+2.0*k2vx+2.0*k3vx+k4vx); 
	vy=vy+(h/6.0)*(k1vy+2.0*k2vy+2.0*k3vy+k4vy);
	//___________________________________________________________________							
	
	fprintf(fpm,"%f	%f\n",xm+x,ym+y); 
	k1xm=vxm;
	k1ym=vym; 
	k1vxm=axm(xm,ym); 
	k1vym=aym(xm,ym);
	k2xm=vxm+k1vxm*h/2.0; 
	k2ym=vym+k1vym*h/2.0; 
	k2vxm=axm(xm+k1xm*h/2.0,ym+k1ym*h/2.0); 
	k2vym=aym(xm+k1xm*h/2.0,ym+k1ym*h/2.0);
	k3xm=vxm+k2vxm*h/2.0; 
	k3ym=vym+k2vym*h/2.0; 
	k3vxm=axm(xm+k2xm*h/2.0,ym+k2ym*h/2.0); 
	k3vym=aym(xm+k2xm*h/2.0,ym+k2ym*h/2.0);
	k4xm=vxm+k3vxm*h; k4ym=vym+k3vym*h; 
	k4vxm=axm(xm+k3xm*h,ym+k3ym*h); 
	k4vym=aym(xm+k3xm*h,ym+k3ym*h);
	xm=xm+(h/6.0)*(k1xm+2.0*k2xm+2.0*k3xm+k4xm); 
	ym=ym+(h/6.0)*(k1ym+2.0*k2ym+2.0*k3ym+k4ym); 
	vxm=vxm+(h/6.0)*(k1vxm+2.0*k2vxm+2.0*k3vxm+k4vxm); 
	vym=vym+(h/6.0)*(k1vym+2.0*k2vym+2.0*k3vym+k4vym);
	}
	//n++;
	double e=energy(x,y,vx,vy);
	//printf("%d	%f\n",n,e);
	fprintf(fpe,"%d	%f\n",n,e);
	//fprintf(fp2,"%f	%d\n",y,n);
	//fprintf(fp3,"%f	%d\n",vx,n);
	//fprintf(fp4,”%f	%d\n",vy,n);
	}//while(n<365);
}
//_______________________________________________________________________						

double rk_mars()
{ 
	FILE*fp=NULL;
	fp=fopen("t_mars.txt","w"); double x=-224490000000.0;
	double y=0.0; double vx=0.0; double vy=-24070.0;
	//double a[4];
	//int n=0;
	double h =88620.0;
	double k1x,k2x,k3x,k4x,k1y,k2y,k3y,k4y,k1vx,k2vx,k3vx,k4vx,k1vy,k2vy,k3vy,k4vy; 
	for(int n=0;n<=687;n++)
	{
		fprintf(fp,"%f	%f\n",x,y); k1x=vx;
		k1y=vy; k1vx=ax(x,y);
		k1vy=ay(x,y);
		k2x=vx+k1vx*h/2.0;
	 	k2y=vy+k1vy*h/2.0; 
		k2vx=ax(x+k1x*h/2.0,y+k1y*h/2.0);
		k2vy=ay(x+k1x*h/2.0,y+k1y*h/2.0);
		k3x=vx+k2vx*h/2.0;
		k3y=vy+k2vy*h/2.0; 
		k3vx=ax(x+k2x*h/2.0,y+k2y*h/2.0); 
		k3vy=ay(x+k2x*h/2.0,y+k2y*h/2.0);
		k4x=vx+k3vx*h;
		k4y=vy+k3vy*h; 
		k4vx=ax(x+k3x*h,y+k3y*h); 
		k4vy=ay(x+k3x*h,y+k3y*h);
		x=x+(h/6.0)*(k1x+2.0*k2x+2.0*k3x+k4x); 
		y=y+(h/6.0)*(k1y+2.0*k2y+2.0*k3y+k4y); 
		vx=vx+(h/6.0)*(k1vx+2.0*k2vx+2.0*k3vx+k4vx); 
		vy=vy+(h/6.0)*(k1vy+2.0*k2vy+2.0*k3vy+k4vy);
	}
}
//_______________________________________________________________________							

void sun()
{ 
	FILE*fp=NULL;
	fp=fopen("sun.txt","w"); 
	for(int i=1;i<=365;i++)
	{
		double x=Rsun*cos(i); double y=Rsun*sin(i);
		//printf("%f	%f\n",x,y);
		fprintf(fp,"%f	%f\n",x,y);
	}
}
//========================================================================

int main()
{
	rk();
	sun(); 
	rk_mars(); 
	central(); 
	int n=2;
	double a=sqrt(146000000000),b=sqrt(152000000000),I,A;
	//printf("%f %f",a,b);
	do
	{
		I=A;
		n=n+2; 
		A=simpson13(f,a,b,n);
	}while(fabs(A-I)>=0.00001);
	printf("\nThe energy integral using Simpson's Rule is : %lf\n",A);
	n=2;
	do
	{
		I=A; n++;
		A=trapezoidal(f,a,b,n);
	}while(fabs(A-I)>=0.00001);
	printf("The energy integral using Trapezoidal Rule is: %lf",A);
}

/* OUTPUT
The energy integral using Simpson's Rule is :
-15963924776896031000000000000000000000000000.000000
The energy integral using Trapezoidal Rule is:
-15963924777020315000000000000000000000000000.000000
 
Process exited after 0.2448 seconds with return value 0 Press any key to continue . . .
*/
