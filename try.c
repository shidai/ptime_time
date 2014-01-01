#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "T2toolkit.h"

void convertStokesAABBCRCI(double *stokes_act,double *obsact);

int main()
{
  long seed = TKsetSeed();
  double AA,BB;
  double stokes[4];
  double Ex2,Sx2,Sy2,Ey2,Nx2,ExEys,EyExs,Ny2;
  double ax,c,ay,d;
  double ax_r,ax_i,ax2,c2,ay_r,ay_i,ay2,d2;
  double c_r,c_i,d_r,d_i;
  double noise_axex;
  double noise_axsx;
  double noise_cey;
  double noise_csy;
  double noise_nx;
  double noise_ct1;
  double noise_ct2;
  double noise_ct3;
  double noise_ct4;
  double noise_ct5;
  double noise_ct6;
  double noise_ct7;
  double noise_ct8;
  double noise_ct9;
  double noise_ct10;
  double ndof;
  int nbin;
  int i;
  double c_bw = 256.0e6/1024.0;
  double t_sub = 60;
  double rho2;

  nbin = 256;
  Nx2 = 0.1;  // Noise term for LNA 1
  Ny2 = 0.1;  // Noise term for LNA 2
  Sx2 = 0.1;  // Sky noise
  Sy2 = 0.1;  // Sky noise

  ax_r = 0.01;
  ax_i = 0.04;
  c_r  = 0.02;
  c_i  = 0.005;

  ay_r = 0.02;
  ay_i = 0.05;
  d_r  = 0.02;
  d_i = 0.005;

  ax2 = ax_r*ax_r+ax_i*ax_i;
  c2 = c_r*c_r+c_i*c_i;
  ax =sqrt(ax2);
  c = sqrt(c2);

  ay2 = ay_r*ay_r+ay_i*ay_i;
  d2 = d_r*d_r+d_i*d_i;
  ay =sqrt(ay2);
  d = sqrt(d2);

  for (i=0;i<nbin;i++)
    {
      if (i > 50 && i < 200) // Cal ON
	{
	  Ex2 = 1;
	  Ey2 = 1;
	  stokes[0] = 1;
	  stokes[1] = 0;
	  stokes[2] = 1;
	  stokes[3] = 0;
	}
      else // Cal OFF
	{ 
	  Ex2 = 0;
	  Ey2 = 0;
	  stokes[0] = 1;
	  stokes[1] = 0;
	  stokes[2] = 0;
	  stokes[3] = 0;
	}

      // Calculate mean value of AA
      AA = ax2*(Ex2+Sx2) + c2*(Ey2 + Sy2) + Nx2 + stokes[2]*(ax_r*c_r + ax_i*c_i)
	+ stokes[3]*(ax_i*c_r-ax_r*c_i);

      // Calculate auto noise terms
      ndof = c_bw*t_sub/(double)nbin;
      noise_axex = ax2*Ex2/sqrt(ndof);
      noise_axsx = ax2*Sx2/sqrt(ndof);
      noise_cey  = c2*Ey2/sqrt(ndof);
      noise_csy  = c2*Sy2/sqrt(ndof);
      noise_nx   = Nx2/sqrt(ndof);

      AA += TKgaussDev(&seed)*noise_axex;
      AA += TKgaussDev(&seed)*noise_axsx;
      AA += TKgaussDev(&seed)*noise_cey;
      AA += TKgaussDev(&seed)*noise_csy;
      AA += TKgaussDev(&seed)*noise_nx;

      // Calculate cross noise terms
      
      rho2 = pow(stokes[2]/2.0,2)+pow(stokes[3]/2.0,2);
      if (Ex2 > 0 && Ey2 > 0)
	noise_ct1 = 2*ax*c*sqrt(Ex2)*sqrt(Ey2)/sqrt(ndof)*sqrt(0.5*(1+rho2/Ex2/Ey2));
      else noise_ct1 = 0;

      noise_ct2 = 2*ax*c*sqrt(Sx2)*sqrt(Sy2)/sqrt(ndof);
      noise_ct3 = 2*ax*c*sqrt(Ex2)*sqrt(Sy2)/sqrt(ndof);
      noise_ct4 = 2*ax*c*sqrt(Sx2)*sqrt(Ey2)/sqrt(ndof);
      noise_ct5 = 2*c*c*sqrt(Ey2)*sqrt(Sy2)/sqrt(ndof);
      noise_ct6 = 2*ax*ax*sqrt(Ex2)*sqrt(Sx2)/sqrt(ndof);
      noise_ct7 = 2*ax*sqrt(Ex2)*sqrt(Nx2)/sqrt(ndof);
      noise_ct8 = 2*ax*sqrt(Sx2)*sqrt(Nx2)/sqrt(ndof);
      noise_ct9 = 2*c*sqrt(Ey2)*sqrt(Nx2)/sqrt(ndof);
      noise_ct10 = 2*c*sqrt(Sy2)*sqrt(Nx2)/sqrt(ndof);

      AA += TKgaussDev(&seed)*noise_ct1;
      AA += TKgaussDev(&seed)*noise_ct2;
      AA += TKgaussDev(&seed)*noise_ct3;
      AA += TKgaussDev(&seed)*noise_ct4;
      AA += TKgaussDev(&seed)*noise_ct5;
      AA += TKgaussDev(&seed)*noise_ct6;
      AA += TKgaussDev(&seed)*noise_ct7;
      AA += TKgaussDev(&seed)*noise_ct8;


      // NOW CALCULATE BB
      // Calculate mean value of BB
      BB = ay2*(Ey2+Sy2) + d2*(Ex2 + Sx2) + Ny2 + stokes[2]*(ay_r*d_r + ay_i*d_i)
	+ stokes[3]*(ay_i*d_r-ay_r*d_i);

      // Calculate auto noise terms
      noise_axex = ay2*Ey2/sqrt(ndof);  // This naming scheme is misleading - just copied from AA
      noise_axsx = ay2*Sy2/sqrt(ndof);
      noise_cey  = d2*Ex2/sqrt(ndof);
      noise_csy  = d2*Sx2/sqrt(ndof);
      noise_nx   = Ny2/sqrt(ndof);

      BB += TKgaussDev(&seed)*noise_axex;
      BB += TKgaussDev(&seed)*noise_axsx;
      BB += TKgaussDev(&seed)*noise_cey;
      BB += TKgaussDev(&seed)*noise_csy;
      BB += TKgaussDev(&seed)*noise_nx;

      // Calculate cross noise terms
      
      rho2 = pow(stokes[2]/2.0,2)+pow(stokes[3]/2.0,2);
      if (Ey2 > 0 && Ex2 > 0)
	noise_ct1 = 2*ay*d*sqrt(Ey2)*sqrt(Ex2)/sqrt(ndof)*sqrt(0.5*(1+rho2/Ex2/Ey2));
      else noise_ct1 = 0;

      noise_ct2 = 2*ay*d*sqrt(Sy2)*sqrt(Sx2)/sqrt(ndof);
      noise_ct3 = 2*ay*d*sqrt(Ey2)*sqrt(Sx2)/sqrt(ndof);
      noise_ct4 = 2*ay*d*sqrt(Sy2)*sqrt(Ex2)/sqrt(ndof);
      noise_ct5 = 2*d*d*sqrt(Ex2)*sqrt(Sx2)/sqrt(ndof);
      noise_ct6 = 2*ay*ay*sqrt(Ey2)*sqrt(Sy2)/sqrt(ndof);
      noise_ct7 = 2*ay*sqrt(Ey2)*sqrt(Ny2)/sqrt(ndof);
      noise_ct8 = 2*ay*sqrt(Sy2)*sqrt(Ny2)/sqrt(ndof);
      noise_ct9 = 2*d*sqrt(Ex2)*sqrt(Ny2)/sqrt(ndof);
      noise_ct10 = 2*d*sqrt(Sx2)*sqrt(Ny2)/sqrt(ndof);

      BB += TKgaussDev(&seed)*noise_ct1;
      BB += TKgaussDev(&seed)*noise_ct2;
      BB += TKgaussDev(&seed)*noise_ct3;
      BB += TKgaussDev(&seed)*noise_ct4;
      BB += TKgaussDev(&seed)*noise_ct5;
      BB += TKgaussDev(&seed)*noise_ct6;
      BB += TKgaussDev(&seed)*noise_ct7;
      BB += TKgaussDev(&seed)*noise_ct8;

      
      printf("%d %g %g\n",i,AA,BB);
      //  double noise_axex;
      //  double noise_axsx;
      //  double noise_cey;
      //  double noise_csy;
      //  double noise_nx;


    }
}



void convertStokesAABBCRCI(double *stokes_act,double *obsact)
{
  double I,Q,U,V,A,B,C,D;

  I = stokes_act[0];
  Q = stokes_act[1];
  U = stokes_act[2];
  V = stokes_act[3];

  A = 0.5*(I+Q);
  B = I-A;
  C = U/2.0;
  D = V/2.0;

  obsact[0] = A;
  obsact[1] = B;
  obsact[2] = C;
  obsact[3] = D;

}
