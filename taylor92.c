// functions used to derive phase shifts, according to Taylor 1992  
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
#include "ptime_time.h"
#include "T2toolkit.h"
#include "tempo2pred.h"
//#include "nrutil.h"
#define ITMAX 100000  // Maximum allowed number of iterations.
#define EPS 1.0e-16 // Machine double floating-point precision.
//#define EPS 3.0e-8 // Machine floating-point precision.
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
double pi=3.1415926;

double A7 (double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn)
// calculate function A7 in Taylor 92
//double A7 (int n, double *amp_s, double *amp_p, double *phi_s, double *phi_p, double phase)
{
	double A7=0;
	int i,j;

	for (i = 0; i < nchn; i++)
	{
		for (j = 0; j < num; j++)
	    {
			A7+=(j+1)*a_s[i][j]*a_p[i][j]*sin(p_s[i][j]-p_p[i][j]+(j+1)*phase);
		//printf ("%lf %lf\n", a_s[i], p_s[i]);
		}
	}
	
	return A7;
}

/*
double A7_multi (double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn, double *rms, double *b)
// calculate function A7 in Taylor 92, for multi-frequency channel
//double A7 (int n, double *amp_s, double *amp_p, double *phi_s, double *phi_p, double phase)
{
	double A7=0;
	int i,j;

	for (i = 0; i < nchn; i++)
	{
		for (j = 0; j < num; j++)
	    {
			A7+=(b[i]*(j+1)*a_s[i][j]*a_p[i][j]*sin(p_s[i][j]-p_p[i][j]+(j+1)*phase))/(rms[i]*rms[i]);
		//printf ("%lf %lf\n", a_s[i], p_s[i]);
		}
	}
	
	return A7;
}
*/

double A7_multi (double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn, double *rms)
// calculate function A7 in Taylor 92, for multi-frequency channel
//double A7 (int n, double *amp_s, double *amp_p, double *phi_s, double *phi_p, double phase)
{
	double A7=0;
	int i,j;
	double c1, c2, s ;

	for (i = 0; i < nchn; i++)
	{
		c1 = 0.0;
		c2 = 0.0;
		s = 0.0;
		for (j = 0; j < num; j++)
		{
			c1 += (j+1)*a_s[i][j]*a_p[i][j]*sin(p_s[i][j]-p_p[i][j]+(j+1)*phase);
			c2 += a_s[i][j]*a_p[i][j]*cos(p_s[i][j]-p_p[i][j]+(j+1)*phase);
			s += a_s[i][j]*a_s[i][j];
		//printf ("%lf %lf\n", a_s[i], p_s[i]);
		}
		A7 += (c1*c2)/(s*rms[i]*rms[i]);
	}
	
	return A7;
}

double A9 (double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn)
// calculate function A9 in Taylor 92
{
	double A9=0.0, sum=0.0;
	int i,j;

	for (i = 0; i < nchn; i++)
	{
	    for (j = 0; j < num; j++)
	    {
		    A9+=a_s[i][j]*a_p[i][j]*cos(p_s[i][j]-p_p[i][j]+(j+1)*phase);
		    sum+=a_s[i][j]*a_s[i][j];
		    //printf ("%lf %lf\n", a_s[i], p_s[i]);
		}
	}
	
	A9=A9/sum;

	return A9;
}

/*
double A9_multi (double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn, double *rms)
{
// calculate function A9 in Taylor 92, for multi-frequency channel
	double A9=0.0, sum=0.0;
	int i,j;

	for (i = 0; i < nchn; i++)
	{
	    for (j = 0; j < num; j++)
	    {
		    A9+=(a_s[i][j]*a_p[i][j]*cos(p_s[i][j]-p_p[i][j]+(j+1)*phase))/(rms[i]*rms[i]);
		    sum+=(a_s[i][j]*a_s[i][j])/(rms[i]*rms[i]);
		    //printf ("%lf %lf\n", a_s[i], p_s[i]);
		}
	}
	
	A9=A9/sum;

	return A9;
}
*/

int dft_profiles (int N, double *in, fftw_complex *out)
// dft of profiles
{
	//  dft of profiles 
	///////////////////////////////////////////////////////////////////////
	
	//printf ("%lf\n", in[0]);
	//double *in;
	//fftw_complex *out;
	fftw_plan p;
	
	//in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	//out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_MEASURE);

	fftw_execute(p);

	fftw_destroy_plan(p);
	//fftw_free(in); 
	//fftw_free(out);
  
	return 0;
}

//int error (double phase, double b, double a,  double *errphase, double *errb)
int error (double phase, double b, double *errphase, double *errb, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn)
// calculate the errors of phase, a and b according to Talyor 1992  
{
	double rms,gk,s1,s2;
	int i,j,n;

	gk=0.0;
	s1=0.0;
	s2=0.0;
	n=0;

	for (i = 0; i < nchn; i++)
	{
	    for (j = 0; j < num; j++)
	    {
		    //gk+=a_p[i]*a_p[i]+b*b*a_s[i]*a_s[i]-2.0*b*a_s[i]*a_p[i]*cos(p_p[i]-p_s[i]-(i+1)*phase)+a*a*1024.0*1024.0-2.0*a*1024.0*a_p[i]*cos(p_p[i])+2.0*a*b*1024.0*a_s[i]*cos(p_s[i]+(i+1)*phase);
		    gk+=a_p[i][j]*a_p[i][j]+b*b*a_s[i][j]*a_s[i][j]-2.0*b*a_s[i][j]*a_p[i][j]*cos(p_p[i][j]-p_s[i][j]-(j+1)*phase);
		//printf ("%lf %lf\n", a_s[i], p_s[i]);
		s1+=(j+1)*(j+1)*a_p[i][j]*a_s[i][j]*cos(p_p[i][j]-p_s[i][j]-(j+1)*phase);
		s2+=a_s[i][j]*a_s[i][j];
		n++;
		}
	}
	
	rms=sqrt(gk/n);

	(*errphase)=rms/sqrt(2.0*fabs(b*s1));
	(*errb)=rms/sqrt(2.0*fabs(s2));

	return 0;
}

/*
int error_multi (double phase, double *errphase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn, double *rms, double *bx)
// calculate the errors of phase, a and b according to Talyor 1992  
{
	double s1;
	int i,j,n;

	n=0;
	s1=0.0;

	for (i = 0; i < nchn; i++)
	{
	    for (j = 0; j < num; j++)
	    {
		    s1+=(bx[i]*(j+1)*(j+1)*a_p[i][j]*a_s[i][j]*cos(p_p[i][j]-p_s[i][j]-(j+1)*phase))/(rms[i]*rms[i]);
		    //s2+=(a_s[i][j]*a_s[i][j])/(rms[i]*rms[i]);
		n++;
		}
	}
	
	(*errphase)=1.0/sqrt(2.0*fabs(s1));
	//(*errb)=1.0/sqrt(2.0*s2);

	return 0;
}
*/

int error_multi (double phase, double *errphase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn, double *rms)
// calculate the errors of phase, a and b according to Talyor 1992  
{
	double s1;
	int i,j;

	s1=0.0;

	double c1, c2, c3, s;

	for (i = 0; i < nchn; i++)
	{
		c1 = 0.0;
		c2 = 0.0;
		c3 = 0.0;
		s = 0.0;
		for (j = 0; j < num; j++)
		{
			c1 += (j+1)*a_s[i][j]*a_p[i][j]*sin(p_s[i][j]-p_p[i][j]+(j+1)*phase);
			c2 += a_s[i][j]*a_p[i][j]*cos(p_s[i][j]-p_p[i][j]+(j+1)*phase);
			c3 += (j+1)*(j+1)*a_s[i][j]*a_p[i][j]*cos(p_s[i][j]-p_p[i][j]+(j+1)*phase);
			s += a_s[i][j]*a_s[i][j];
			//s2+=(a_s[i][j]*a_s[i][j])/(rms[i]*rms[i]);
		}
		s1 += (c2*c3-c1*c1)/(s*rms[i]*rms[i]);
	}
	
	(*errphase)=1.0/sqrt(2.0*fabs(s1));
	//(*errb)=1.0/sqrt(2.0*s2);

	return 0;
}

double zbrent(double (*func)(double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn), double x1, double x2, double tol, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn)
//	Using Brent’s method, find the root of a function func known to lie between x1 and x2. The root, returned as zbrent, will be refined until its accuracy is tol.
{
	int iter;
	double a=x1,b=x2,c=x2,d,e,min1,min2;
	double fa=(*func)(a, a_s, a_p, p_s, p_p, num, nchn),fb=(*func)(b, a_s, a_p, p_s, p_p, num, nchn),fc,p,q,r,s,tol1,xm;

	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
		puts ("Root must be bracketed in zbrent\n");

	fc=fb;
	for (iter=1;iter<=ITMAX;iter++) 
	{
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) 
		{
			c=a;   // Rename a, b, c and adjust bounding interval d.
			fc=fa;
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) 
		{
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}

		tol1=2.0*EPS*fabs(b)+0.5*tol;   // Convergence check.
		xm=0.5*(c-b);

		if (fabs(xm) <= tol1 || fb == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) 
		{
			s=fb/fa;  // Attempt inverse quadratic interpolation.

			if (a == c) 
			{
				p=2.0*xm*s;
				q=1.0-s;
			} 
			else 
			{
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;  // Check whether in bounds.

			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);

			if (2.0*p < (min1 < min2 ? min1 : min2)) 
			{
				e=d;  // Accept interpolation.
				d=p/q;
			} 
			else 
			{
				d=xm; // Interpolation failed, use bisection.
				e=d;
			}
		} 
		else  // Bounds decreasing too slowly, use bisection.
		{
			d=xm;
			e=d;
		}
		a=b;  //  Move last best guess to a.
		fa=fb;
		if (fabs(d) > tol1)     //  Evaluate new trial root.
			b += d;
		else
			b += SIGN(tol1,xm);

		fb=(*func)(b, a_s, a_p, p_s, p_p, num, nchn);
	}

	puts ("Maximum number of iterations exceeded in zbrent\n");

	return 0.0;
}

double zbrent_multi(double (*func)(double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn, double *rms), double x1, double x2, double tol, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn, double *rms)
//double zbrent_multi(double (*func)(double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn, double *rms, double *bx), double x1, double x2, double tol, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn, double *rms, double *bx)
{
	int iter;
	double a=x1,b=x2,c=x2,d,e,min1,min2;
	double fa=(*func)(a, a_s, a_p, p_s, p_p, num, nchn, rms),fb=(*func)(b, a_s, a_p, p_s, p_p, num, nchn, rms),fc,p,q,r,s,tol1,xm;
	//double fa=(*func)(a, a_s, a_p, p_s, p_p, num, nchn, rms, bx),fb=(*func)(b, a_s, a_p, p_s, p_p, num, nchn, rms, bx),fc,p,q,r,s,tol1,xm;

	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
		puts ("Root must be bracketed in zbrent\n");

	fc=fb;
	for (iter=1;iter<=ITMAX;iter++) 
	{
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) 
		{
			c=a;   // Rename a, b, c and adjust bounding interval d.
			fc=fa;
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) 
		{
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}

		tol1=2.0*EPS*fabs(b)+0.5*tol;   // Convergence check.
		xm=0.5*(c-b);

		if (fabs(xm) <= tol1 || fb == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) 
		{
			s=fb/fa;  // Attempt inverse quadratic interpolation.

			if (a == c) 
			{
				p=2.0*xm*s;
				q=1.0-s;
			} 
			else 
			{
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;  // Check whether in bounds.

			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);

			if (2.0*p < (min1 < min2 ? min1 : min2)) 
			{
				e=d;  // Accept interpolation.
				d=p/q;
			} 
			else 
			{
				d=xm; // Interpolation failed, use bisection.
				e=d;
			}
		} 
		else  // Bounds decreasing too slowly, use bisection.
		{
			d=xm;
			e=d;
		}
		a=b;  //  Move last best guess to a.
		fa=fb;
		if (fabs(d) > tol1)     //  Evaluate new trial root.
			b += d;
		else
			b += SIGN(tol1,xm);

		fb=(*func)(b, a_s, a_p, p_s, p_p, num, nchn, rms);
		//fb=(*func)(b, a_s, a_p, p_s, p_p, num, nchn, rms, bx);
	}

	puts ("Maximum number of iterations exceeded in zbrent\n");

	return 0.0;
}

/*
int find_peak (int n, double *s, int *position)
{
	int i;
	double temp[n];
	double peak;

	for (i=0;i<n;i++)
	{
		temp[i]=s[i];
	}

	double a,b,c;
	for (i=0;i<n-1;i++)
	{
		a=temp[i];
		b=temp[i+1];
		c=(a>=b ? a : b);

		temp[i+1]=c;
	}
	peak=temp[n-1];

	for (i=0;i<n;i++)
	{
		if (fabs(peak-s[i])<1.0e-3)
		{
			(*position)=i;
		}
	}

	return 0;
}

double find_peak_value (int n, double *s)
{
	int i;
	double temp[n];

	for (i=0;i<n;i++)
	{
		temp[i]=s[i];
	}

	double a,b,c;
	for (i=0;i<n-1;i++)
	{
		a=temp[i];
		b=temp[i+1];
		c=(a>=b ? a : b);

		temp[i+1]=c;
	}

	return temp[n-1];
}
*/

int get_toa (double *s, double *p, double *phasex, double *errphasex, double psrfreq, int nphase, double *rms, double *bx)
// calculate the phase shift between template and simulated (or real) data 
// error of phase can be calculated
// initial guess of phase shift is added
// try to do two-dim template matching
{
    //int nphase=1024;
    int nchn=1;

	// read a std
	
	//puts(argv[1]);
	//puts(argv[2]);
	//double t[nphase*nchn],s[nphase*nchn];
	//int n;

	//readfile(name,&n,t,s);
	//printf ("%d\n", n);
	//puts(argv[1]);

	//////////////////////////////////////////////////////////////////////////
	// simulate data

	//double p[nphase*nchn];
	//double SNR=atof(argv[2]);
	//simulate(n,SNR,s,p);//*/
	
	/////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////
	// dft profile and template
	
	//nchn = n/nphase;
	//printf ("%d\n", nchn);
	int k;  // k=nphase/2

	//double amp_s[nchn][nphase/2],amp_p[nchn][nphase/2];  // elements for calculating A7
	//double phi_s[nchn][nphase/2],phi_p[nchn][nphase/2];
	double amp_s[nchn][NP],amp_p[nchn][NP];  // elements for calculating A7
	double phi_s[nchn][NP],phi_p[nchn][NP];  // the second dim should be NP, which is large enough for different observations

	preA7(&k, amp_s, amp_p, phi_s, phi_p, s, p, nphase, nchn);
	//printf ("%d\n", nchn);
	
	// initial guess of the phase
	int peak_s, peak_p;	

	find_peak(nphase,s,&peak_s);
	find_peak(nphase,p,&peak_p);

	int d;
	double step;
	double ini_phase,up_phase,low_phase;

	d = InitialGuess (s, p, nphase);
	//d=peak_p-peak_s;
	step=2.0*3.1415926/(10.0*nphase);
	//step=2.0*3.1415926/10240.0;

	if (d>=nphase/2)
	{
		ini_phase=2.0*3.1415926*(nphase-1-d)/nphase;
		//ini_phase=2.0*3.1415926*(1023-d)/1024.0;
		up_phase=ini_phase+step;
		low_phase=ini_phase-step;
		while (A7(up_phase, amp_s, amp_p, phi_s, phi_p, k, nchn)*A7(low_phase, amp_s, amp_p, phi_s, phi_p, k, nchn)>0.0)
		{
		    up_phase+=step;
		    low_phase-=step;
		}
	}
	else
	{
		ini_phase=-2.0*3.1415926*d/nphase;
		//ini_phase=-2.0*3.1415926*d/1024.0;
		up_phase=ini_phase+step;
		low_phase=ini_phase-step;
		while (A7(up_phase, amp_s, amp_p, phi_s, phi_p, k, nchn)*A7(low_phase, amp_s, amp_p, phi_s, phi_p, k, nchn)>0.0)
		{
		    up_phase+=step;
		    low_phase-=step;
		}
	}

    // calculate phase shift, a and b
    double phase,b;
    phase=zbrent(A7, low_phase, up_phase, 1.0e-16, amp_s, amp_p, phi_s, phi_p, k, nchn);
    //phase=zbrent(A7, -1.0, 1.0, 1.0e-16);
    //phase=zbrent(A7, -0.005, 0.005, 1.0e-16);
    b=A9(phase, amp_s, amp_p, phi_s, phi_p, k, nchn);
    //a=A4(b);
	(*bx) = b;

		
	//printf ("%.10lf %.10lf\n", phase, A7(phase));
	//printf ("%.10lf \n", ((phase/3.1415926)*5.75/2.0)*1.0e+3);  // microseconds
	//printf ("%.10lf \n", b);
	//printf ("%.10lf \n", a);
	//printf ("///////////////////////// \n");
		
	
	// calculate the errors of phase and b
    double errphase, errb;	

	error(phase,b,&errphase,&errb, amp_s, amp_p, phi_s, phi_p, k,nchn);
	printf ("%.10lf %.10lf\n", ((phase/3.1415926)/(psrfreq*2.0))*1.0e+6, ((errphase/3.1415926)/(psrfreq*2.0))*1.0e+6);  // microseconds
	//printf ("%.10lf %.10lf\n", ((phase/3.1415926)*4.569651/2.0)*1.0e+3, ((errphase/3.1415926)*4.569651/2.0)*1.0e+3);  // microseconds
	//printf ("errphase %.10lf \n", ((errphase/3.1415926)*5.75/2.0)*1.0e+6);
	//printf ("errb %.10lf \n", errb);
	
	(*phasex) = phase;
	(*errphasex) = errphase;

	// calculate the rms
	cal_rms(phase,b,rms, amp_s, amp_p, phi_s, phi_p, k,nchn);

	return 0;
}

int get_toa_multi (double *s, double *p, double *rms, int nchn, double *phasex, double *errphasex, double psrfreq, int nphase)
//int get_toa_multi (double *s, double *p, double *rms, double *bx, int nchn, double *phasex, double *errphasex, double psrfreq, int nphase)
{
    //int nphase=1024;

	// read a std
	
	//puts(argv[1]);
	//puts(argv[2]);
	//double t[nphase*nchn],s[nphase*nchn];
	//int n;

	//readfile(name,&n,t,s);
	//printf ("%d\n", n);
	//puts(argv[1]);

	//////////////////////////////////////////////////////////////////////////
	// simulate data

	//double p[nphase*nchn];
	//double SNR=atof(argv[2]);
	//simulate(n,SNR,s,p);//*/
	
	/////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////
	// dft profile and template
	
	//nchn = n/nphase;
	//printf ("%d\n", nchn);
	int k;  // k=nphase/2

	//double amp_s[nchn][nphase/2],amp_p[nchn][nphase/2];  // elements for calculating A7
	//double phi_s[nchn][nphase/2],phi_p[nchn][nphase/2];
	double amp_s[nchn][NP],amp_p[nchn][NP];  // elements for calculating A7
	double phi_s[nchn][NP],phi_p[nchn][NP];  // the second dim should be NP, which is large enough for different observations

	preA7(&k, amp_s, amp_p, phi_s, phi_p, s, p, nphase, nchn);
	
	// initial guess of the phase
  int peak_s, peak_p;	

	find_peak(nphase,s,&peak_s);
	find_peak(nphase,p,&peak_p);

	int d;
	double step;
	double ini_phase,up_phase,low_phase;

	d = InitialGuess (s, p, nphase);
	//d=peak_p-peak_s;
	step=2.0*3.1415926/(10.0*nphase);
	//step=2.0*3.1415926/10240.0;

	if (d>=nphase/2)
	{
		ini_phase=2.0*3.1415926*(nphase-1-d)/nphase;
		up_phase=ini_phase+step;
		low_phase=ini_phase-step;
		while (A7_multi(up_phase, amp_s, amp_p, phi_s, phi_p, k, nchn, rms)*A7_multi(low_phase, amp_s, amp_p, phi_s, phi_p, k, nchn, rms)>0.0)
		//while (A7_multi(up_phase, amp_s, amp_p, phi_s, phi_p, k, nchn, rms, bx)*A7_multi(low_phase, amp_s, amp_p, phi_s, phi_p, k, nchn, rms,bx)>0.0)
		{
		    up_phase+=step;
		    low_phase-=step;
		}
	}
	else
	{
		ini_phase=-2.0*3.1415926*d/nphase;
		up_phase=ini_phase+step;
		low_phase=ini_phase-step;
		while (A7_multi(up_phase, amp_s, amp_p, phi_s, phi_p, k, nchn, rms)*A7_multi(low_phase, amp_s, amp_p, phi_s, phi_p, k, nchn, rms)>0.0)
		//while (A7_multi(up_phase, amp_s, amp_p, phi_s, phi_p, k, nchn, rms, bx)*A7_multi(low_phase, amp_s, amp_p, phi_s, phi_p, k, nchn, rms, bx)>0.0)
		{
		    up_phase+=step;
		    low_phase-=step;
		}
	}

    // calculate phase shift, a and b
    double phase;
    phase=zbrent_multi(A7_multi, low_phase, up_phase, 1.0e-16, amp_s, amp_p, phi_s, phi_p, k, nchn, rms);
    //phase=zbrent_multi(A7_multi, low_phase, up_phase, 1.0e-16, amp_s, amp_p, phi_s, phi_p, k, nchn, rms, bx);
    //phase=zbrent(A7, -1.0, 1.0, 1.0e-16);
    //phase=zbrent(A7, -0.005, 0.005, 1.0e-16);
    //b=A9_multi(phase, amp_s, amp_p, phi_s, phi_p, k, nchn, rms);
    //a=A4(b);

		
	//printf ("%.10lf %.10lf\n", phase, A7(phase));
	//printf ("%.10lf \n", ((phase/3.1415926)*5.75/2.0)*1.0e+3);  // microseconds
	//printf ("%.10lf \n", b);
	//printf ("%.10lf \n", a);
	//printf ("///////////////////////// \n");
		
	
	// calculate the errors of phase and b
  double errphase;	

	error_multi(phase, &errphase, amp_s, amp_p, phi_s, phi_p, k, nchn, rms);
	//error_multi(phase, &errphase, amp_s, amp_p, phi_s, phi_p, k, nchn, rms, bx);
	printf ("multi-template\n");
	printf ("%.10lf %.10lf\n", ((phase/3.1415926)/(psrfreq*2.0))*1.0e+6, ((errphase/3.1415926)/(psrfreq*2.0))*1.0e+6);  // microseconds
	//printf ("%.10lf %.10lf\n", ((phase/3.1415926)*4.569651/2.0)*1.0e+3, ((errphase/3.1415926)*4.569651/2.0)*1.0e+3);  // microseconds
	//printf ("errphase %.10lf \n", ((errphase/3.1415926)*5.75/2.0)*1.0e+6);
	//printf ("errb %.10lf \n", errb);
	
	(*phasex) = phase;
	(*errphasex) = errphase;

	return 0;
}

int preA7 (int *k, double amp_s[][NP], double amp_p[][NP], double phi_s[][NP], double phi_p[][NP], double *s, double *p, int nphase, int nchn)
// preparation for calculating A7 of Talyor 1992  
{
	// nphase is the dimention of one profile, nchn is number of profiles
	// k is the dimention of amp of one profile 
	int i,j;
	
	/////////////////////////////////////////////////////////////////////////////////
	double test[nphase];  // initialize the system, don't know why....

	for (i=0;i<nphase;i++)
	{
		test[i]=s[i];
	}
	fftw_complex *out_t;
	out_t = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nphase);
	dft_profiles(nphase,test,out_t);
	//////////////////////////////////////////////////////////////////////////////

    fftw_complex *out_s;
	fftw_complex *out_p;
	
	out_s = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nphase);
	out_p = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nphase);
	
	double s_temp[nphase];  // store one template and profile
	double p_temp[nphase];  

	int n;
	double r_s[nphase/2],im_s[nphase/2];
	double r_p[nphase/2],im_p[nphase/2];
	for (i = 0; i < nchn; i++)
	{
	    for (j=0;j<nphase;j++)
	    {
		    s_temp[j]=s[i*nphase + j];
		    p_temp[j]=p[i*nphase + j];
	    }

	    dft_profiles(nphase,s_temp,out_s);
	    //printf ("%lf %lf\n", out_s[1][0], out_s[1][1]);

	    dft_profiles(nphase,p_temp,out_p);

	    //double amp_s[N/2],phi_s[N/2];
	    //double amp_p[N/2],phi_p[N/2];

		n = 0;
	    for (j = 0; j <= nphase/2-1; j++)
	    {
		    r_s[j]=out_s[j+1][0];
		    im_s[j]=out_s[j+1][1];
		    r_p[j]=out_p[j+1][0];
		    im_p[j]=out_p[j+1][1];
		    //printf ("%lf %lf\n", r_p[i], im_p[i]);
		    //printf ("%lf %lf\n", out_s[i][0], out_s[i][1]);
		    n++;
	    }
	    //printf ("%d\n", n);
	    //printf ("%d %d\n", nphase, nchn);

	    for (j = 0; j < n; j++)
	    {
		    amp_s[i][j]=sqrt(r_s[j]*r_s[j]+im_s[j]*im_s[j]);
		    amp_p[i][j]=sqrt(r_p[j]*r_p[j]+im_p[j]*im_p[j]);
		    phi_s[i][j]=atan2(im_s[j],r_s[j]);
		    phi_p[i][j]=atan2(im_p[j],r_p[j]);
		    //printf ("%lf %lf %lf\n", r_s[i], im_s[i], amp_s[i]);
		    //printf ("%lf %lf %lf\n", r_p[i], im_p[i], amp_p[i]);
		    //printf ("%lf\n", amp_s[i]);
		    //printf ("%lf\n", amp_p[i]);
	    }
	}
	(*k)=n;

	fftw_free(out_s); 
	fftw_free(out_p); 
	fftw_free(out_t); 

	return 0;
}

int cal_rms (double phase, double b, double *rms, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn)
// calculate the rms of each subchannel  
{
	double gk;
	int i,j,n;

	gk=0.0;
	n=0;

	for (i = 0; i < nchn; i++)
	{
	    for (j = 0; j < num; j++)
	    {
		    //gk+=a_p[i]*a_p[i]+b*b*a_s[i]*a_s[i]-2.0*b*a_s[i]*a_p[i]*cos(p_p[i]-p_s[i]-(i+1)*phase)+a*a*1024.0*1024.0-2.0*a*1024.0*a_p[i]*cos(p_p[i])+2.0*a*b*1024.0*a_s[i]*cos(p_s[i]+(i+1)*phase);
		    gk+=a_p[i][j]*a_p[i][j]+b*b*a_s[i][j]*a_s[i][j]-2.0*b*a_s[i][j]*a_p[i][j]*cos(p_p[i][j]-p_s[i][j]-(j+1)*phase);
		//printf ("%lf %lf\n", a_s[i], p_s[i]);
		n++;
		}
	}
	
	(*rms)=sqrt(gk/n);

	return 0;
}

/*
int simulate (int n, double SNR, double *s, double *p)
// simulate pulse profiles, adding white noise; return simulated profiles
{
	// simulate a profile with white noise
	///////////////////////////////////////////////////////////////////////
	// initialize gsl 
	
	int i;
	const gsl_rng_type * T;
	gsl_rng * r;

	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
  
	////////////////////////////////////////////////////////////////////////
	//  determine the amplitude of white noise according to SNR
	
	double scale;   // the scale multiply to white noise to get certain SNR
	double amp_noise, noise[n];
	double peak_s;

	for (i=0;i<n;i++)
	{
		noise[i]=gsl_ran_gaussian(r,1.0);
		amp_noise+=noise[i]*noise[i];
	}
	
	amp_noise=sqrt(amp_noise/n);

	peak_s=find_peak_value(n,s);   // find the peak flux of the std
    //printf ("peak of std: %g\n", peak_s);

	scale=peak_s/(SNR*amp_noise);
    //printf ("%g\n", scale);
	
	//////////////////////////////////////////////////////////////////////////
	//  add noise to std ==> p

	for (i=0;i<n;i++)
	{
		p[i]=(s[i]+scale*noise[i]);
	}

	
	double peak_p;
	peak_p=find_peak(n,p);  // find the peak flux of the profile
    //printf ("peak of profile: %g\n", peak_p);

	//  normalize the std and profile
	
	for (i=0;i<n;i++)
	{
		p[i]=p[i]/peak_p;
		s[i]=s[i]/peak_s;
		//printf ("%g %g\n", s[i], p[i]);
	}

	//double rms=0.0;
	//int m=0;

	//for (i=300;i<700;i++)
	//{
    //	rms+=p[i]*p[i];
    //	m++;
	//}

	//rms=sqrt(rms/m);
	//printf("rms is: %f\n", rms);

	gsl_rng_free (r);
  
	return 0;
}
*/

int form_toa_multi (char *name_data, char *name_predict, int subint, int nchn, long int imjd, long int smjd, double offs, double phase, double e_phase, long double *t, long double *e_dt, double *freqout)
{
	int h;
	h = subint;

	long double dt;  
	double offset;   // offset of each subint
	long double mjd0;  // the mjd of each subint
	T2Predictor pred;
	int ret;
	double period, weight, frequency;
	double freq[nchn], wts[nchn];

	// transform phase shift to TOAs
	// get the freq of the subint
	read_freq(name_data, h, freq, nchn);
	read_wts(name_data, h, wts, nchn);
	frequency = 0.0;
	weight = 0.0;

	int z;
	for (z = 0; z < nchn; z++)
	{
		frequency += freq[z]*wts[z];
		weight += wts[z];
	}
	frequency = frequency/weight;
	(*freqout) = frequency;
    printf ("Frequency is %lf\n", frequency);

	// get the period
    print_t2pred(name_predict);   // output t2pred.dat

	T2Predictor_Init(&pred);  // prepare the predictor
	if (ret=T2Predictor_Read(&pred,(char *)"t2pred.dat"))
    {
		printf("Error: unable to read predictor\n");
		exit(1);
	}

	// get the offset of each subint
	offset = read_offs(name_data, h);

	// get the period at mjd0
	mjd0 = (long double)(imjd) + ((long double)(smjd) + (long double)(offs) + (long double)(offset))/86400.0L;
	printf ("imjd is %ld \n", imjd);
	printf ("mjd0 is %.15Lf \n", mjd0);

	period = 1.0/T2Predictor_GetFrequency(&pred,mjd0,frequency);
    printf ("Period is %.15lf\n", period);
	
	// transform phase shift to time shift
    //dt = (phase/PI)*period/2.0;
    //e_dt = (e_phase/PI)*period/2.0;
    dt = ((long double)(phase)/PI)*((long double)(period))/2.0L;
    (*e_dt) = ((long double)(e_phase)/PI)*((long double)(period))/2.0L;
    printf ("dt is %.10Lf +/- %.10Lf\n", dt, e_dt);

	// calculate the TOA
    (*t) = (long double)(imjd) + ((long double)(smjd) + (long double)(offs) - (long double)(dt) + (long double)(offset))/86400.0L;
    //t = imjd;
		
    printf ("offset is %lf\n", offset);
	//fprintf (fp, "%s  %lf  %.15Lf  %Lf  7\n", fname, frequency, t, e_dt*1e+6);

	return 0;
}

int form_toa (char *name_data, char *name_predict, int subint, int chn, int nchn, long int imjd, long int smjd, double offs, double phase, double e_phase, long double *t, long double *e_dt, double *freqout)
// chn is the channel to form toa
// nchn is the total number of subchn
{
	int h;
	h = subint;

	long double dt;  
	double offset;   // offset of each subint
	long double mjd0;  // the mjd of each subint
	T2Predictor pred;
	int ret;
	double period, frequency;
	double freq[nchn];

	// transform phase shift to TOAs
	// get the freq of the subint
	read_freq(name_data, h, freq, nchn);

	frequency = freq[chn];
	(*freqout) = frequency;
    printf ("Frequency is %lf\n", frequency);

	// get the period
    print_t2pred(name_predict);   // output t2pred.dat

	T2Predictor_Init(&pred);  // prepare the predictor
	if (ret=T2Predictor_Read(&pred,(char *)"t2pred.dat"))
    {
		printf("Error: unable to read predictor\n");
		exit(1);
	}

	// get the offset of each subint
	offset = read_offs(name_data, h);

	// get the period at mjd0
	mjd0 = (long double)(imjd) + ((long double)(smjd) + (long double)(offs) + (long double)(offset))/86400.0L;
	printf ("imjd is %ld \n", imjd);
	printf ("mjd0 is %.15Lf \n", mjd0);

	period = 1.0/T2Predictor_GetFrequency(&pred,mjd0,frequency);
    printf ("Period is %.15lf\n", period);
	
	// transform phase shift to time shift
    //dt = (phase/PI)*period/2.0;
    //e_dt = (e_phase/PI)*period/2.0;
    dt = ((long double)(phase)/PI)*((long double)(period))/2.0L;
    (*e_dt) = ((long double)(e_phase)/PI)*((long double)(period))/2.0L;
    printf ("dt is %.10Lf +/- %.10Lf\n", dt, e_dt);

	// calculate the TOA
    (*t) = (long double)(imjd) + ((long double)(smjd) + (long double)(offs) - (long double)(dt) + (long double)(offset))/86400.0L;
    //t = imjd;
		
    printf ("offset is %lf\n", offset);
	//fprintf (fp, "%s  %lf  %.15Lf  %Lf  7\n", fname, frequency, t, e_dt*1e+6);

	return 0;
}

int find_peak (int n, double *s, int *position)
{
	int i;
	double temp[n];
	double peak;

	for (i = 0; i < n; i++)
	{
		temp[i] = s[i];
	}

	double a, b, c;
	for (i = 0; i < n-1; i++)
	{
		//a = fabs(temp[i]);
		//b = fabs(temp[i+1]);
		a = temp[i];
		b = temp[i+1];
		c = (a >= b ? a : b);

		temp[i+1] = c;
	}
	peak = temp[n-1];

	for (i = 0; i < n; i++)
	{
		if (fabs(peak-fabs(s[i])) < 1.0e-3)
		{
			(*position) = i;
		}
	}

	return 0;
}

double find_peak_value (int n, double *s)
{
	int i;
	double temp[n];

	for (i = 0; i < n; i++)
	{
		temp[i] = s[i];
	}

	double a, b, c;
	for (i = 0; i < n-1; i++)
	{
		a = temp[i];
		b = temp[i+1];
		//a = fabs(temp[i]);
		//b = fabs(temp[i+1]);
		c = (fabs(a) >= fabs(b) ? a : b);

		temp[i+1] = c;
	}

	return temp[n-1];
}

int corr (double *s, double *p, int nphase)
{
	/*
	FILE *fp1, *fp2;

	if ((fp1 = fopen(argv[1], "r")) == NULL)
	{
		fprintf (stdout, "Can't open file\n");
		exit(1);
	}

	if ((fp2 = fopen(argv[2], "r")) == NULL)
	{
		fprintf (stdout, "Can't open file\n");
		exit(1);
	}

	float x1[1024], x2[1024];
	int i = 0;
	while (fscanf (fp1, "%f", &x1[i]) == 1)
	{
		i++;
	}

	i = 0;
	while (fscanf (fp2, "%f", &x2[i]) == 1)
	{
		i++;
	}
	*/

	double r[nphase];
	int i, j;
	for (j = 0; j < nphase; j++)
	{
		r[j] = 0.0;
		for (i = 0; i < nphase; i++)
		{
			if ((i+j) > (nphase-1))
			{
				r[j] += p[i]*s[i+j-(nphase-1)];
			}
			else
			{
				r[j] += p[i]*s[i+j];
			}
			//printf ("%f %f\n", x1[i], x2[i]);
		}
	}

	int shift;
	find_peak (nphase, r,  &shift);
	/*
	for (j = 0; j < 1024; j++)
	{
		printf ("%f\n", r[j]);
	}
	*/

	return -shift;
}

int def_off_pulse (int nphase, double *in, double frac_off)
// define the off pulse region based on I, return the starting index of off pulse region
// using frac_off to calculate the off pulse region
{
	int n = nphase;
	int num_off = (int)(n*frac_off);
	int i,j;
	double small;
	double temp;
	int index = 0;

	for (i = 0; i < n; i++)
	{
		if (i == 0)
		{
			small = 0.0;
			for(j = 0; j < num_off; j++)
			{
				small += (in[j])*(in[j]);
			}
			small = sqrt(small/num_off);
		}
			
		temp = 0.0;
		for(j = 0; j < num_off; j++)
		{
			if ((i+j) > n-1)
			{
				temp += (in[(i+j)-(n-1)])*(in[(i+j)-(n-1)]);
			}
			else 
			{
				temp += (in[i+j])*(in[i+j]);
			}
		}
		temp = sqrt(temp/num_off);

		small = (temp <= small ? temp : small);
		index = (temp <= small ? i : index);
		//printf ("%d %lf %lf\n", index, small, ave);
	}

	return index;
}

int off_pulse (int nphase, int index, double *in, double *out, double frac_off)
// get the off_pulse region
{
	int n = nphase;
	int num_off = (int)(n*frac_off);
	int i;

	for (i = 0; i < num_off; i++)
	{
		if ((index+i) > n-1)
		{
			out[i] = in[(index+i)-(n-1)];
		}
		else 
		{
			out[i] = in[index+i];
		}
	}

	return 0;
}

int remove_baseline (double *in, int index, double frac_off, int n, double *out)
{
	// define the off_pulse range, frac_off is the fraction of the phase
	// index is the starting point of the off_pulse range
	int num_off = (int)(n*frac_off);
	double off_0[num_off];

	off_pulse (n, index, in, off_0, frac_off);

	int i;
	double baseline = 0.0;
    for (i = 0; i < num_off; i++)
    {
        baseline += off_0[i];
        //average_s += s_off[i];
    }
	baseline = baseline/num_off;

    //printf ("the baseline of std is: %lf \n", baseline);
    //printf ("average is: %lf %lf\n", average, average_s);

	for (i = 0; i < n; i++)
	{
		out[i] = (in[i]-baseline);
		//s_norm[i] = (s[i]-baseline)/(s_peak-baseline);
	}
	
	return 0;
}

int pre_diff (double *s, int nphase, int index, double frac_off, double *s_out)
{
	int n = nphase;
	
	// remove the baseline
	remove_baseline (s, index, frac_off, n, s_out);

    return 0;
}


int InitialGuess (double *s, double *p, int nphase) 
{
	int index;
	double frac_off = 0.05;  // set to be 0.1

	// remove the baseline of template
	index = def_off_pulse (nphase, s, frac_off);

	double s_out[nphase];
	pre_diff (s, nphase, index, frac_off, s_out);

	// remove the baseline of profile
	index = def_off_pulse (nphase, p, frac_off);

	double p_out[nphase];
	pre_diff (p, nphase, index, frac_off, p_out);

	// Guess the phase shift
	int d;
	d = corr (s_out, p_out, nphase);

	return d;
}

