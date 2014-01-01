#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
#include <fftw3.h>
#include "fitsio.h"
//#include "ptime.h"

#define NP 2048
#define PI 3.14159265359

long int stt_imjd ( char *name );
long int stt_smjd ( char *name );
double stt_offs ( char *name );

int get_nchan ( char *name );
int get_npol ( char *name );
int get_nphase ( char *name );
int get_subint ( char *name );


int read_std ( char *name, int subint, double *profile, int nphase, int mode, int nchn, double *temp);
int read_prof ( char *name, int subint, double *profile, int phase );


int print_t2pred ( char *name );
double read_offs ( char *name, int subint);
int read_freq ( char *name, int subint, double *freq, int nchan );
int read_wts ( char *name, int subint, double *wts, int nchan );

double read_psrfreq ( char *name );

int readfile ( char *filename, int *ntxt, double *x, double *y );

int dft_profiles (int N, double *in, fftw_complex *out);

//int simulate (int n, double SNR, double *s, double *p);

double find_peak_value (int n, double *s);

int preA7 (int *k, double amp_s[][NP], double amp_p[][NP], double phi_s[][NP], double phi_p[][NP], double *s, double *p, int nphase, int nchn);

//double A7 (double phase);
double A7 (double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn);
double A7_multi (double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn, double *rms);

//double A9 (double phase);
double A9 (double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn);
double A9_multi (double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn, double *rms);

//double zbrent(double (*func)(double), double x1, double x2, double tol);
double zbrent(double (*func)(double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn), double x1, double x2, double tol, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn);

double zbrent_multi(double (*func)(double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn, double *rms), double x1, double x2, double tol, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn, double *rms);

int find_peak (int n, double *s, int *position);

//int error (double phase, double b, double *errphase, double *errb);
int error (double phase, double b, double *errphase, double *errb, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn);

int error_multi (double phase, double b, double *errphase, double *errb, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn, double *rms);

// calculate the rms of each profile
int cal_rms (double phase, double b, double *rms, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn);

double get_toa (double *s, double *p, double psrfreq, int nphase);

int get_toa_multi (double *s, double *p, double *rms, int nchn, double *phasex, double *errphasex, double psrfreq, int nphase);
