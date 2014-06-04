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


int check_std ( char *name, int subint, int mode, int nchn, int nphase);
int read_std ( char *name, int subint, double *profile, int nphase, int mode, int nchn);
int read_prof ( char *name, int subint, double *profile, int phase );


int print_t2pred ( char *name );
double read_offs ( char *name, int subint);
int read_freq ( char *name, int subint, double *freq, int nchan );
int read_wts ( char *name, int subint, double *wts, int nchan );

double read_psrfreq ( char *name );

int readfile ( char *filename, int *ntxt, double *x, double *y );

int dft_profiles (int N, double *in, fftw_complex *out);

//int simulate (int n, double SNR, double *s, double *p);


int preA7 (int *k, double amp_s[][NP], double amp_p[][NP], double phi_s[][NP], double phi_p[][NP], double *s, double *p, int nphase, int nchn);

//double A7 (double phase);
double A7 (double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn);
//double A7_multi (double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn, double *rms, double *b);
double A7_multi (double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn, double *rms);

//double A9 (double phase);
double A9 (double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn);

int A9_multi (double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn, double *b);

//double zbrent(double (*func)(double), double x1, double x2, double tol);
double zbrent(double (*func)(double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn), double x1, double x2, double tol, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn);

double zbrent_multi(double (*func)(double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn, double *rms), double x1, double x2, double tol, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn, double *rms);
//double zbrent_multi(double (*func)(double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn, double *rms, double *bx), double x1, double x2, double tol, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn, double *rms, double *bx);


//int error (double phase, double b, double *errphase, double *errb);
int error (double phase, double b, double *errphase, double *errb, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn);

int error_multi (double phase, double *errphase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn, double *rms);
//int error_multi (double phase, double *errphase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn, double *rms, double *bx);

// calculate the rms of each profile
int cal_rms (double phase, double b, double *rms, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn);

int get_toa (double *s, double *p, double *phasex, double *errphasex, double psrfreq, int nphase, double *rms, double *bx);

int get_toa_multi (double *s, double *p, double *rms, int nchn, double *phasex, double *errphasex, double psrfreq, int nphase);
//int get_toa_multi (double *s, double *p, double *rms, double *bx, int nchn, double *phasex, double *errphasex, double psrfreq, int nphase);

// transform phase shifts to MJD TOAs
int form_toa_multi (char *name_data, char *name_predict, int subint, int nchn, long int imjd, long int smjd, double offs, double phase, double e_phase, long double *t, long double *e_dt, double *frequency);

int form_toa (char *name_data, char *name_predict, int subint, int chn, int nchn, long int imjd, long int smjd, double offs, double phase, double e_phase, long double *t, long double *e_dt, double *frequency);

// initial guess
int find_peak (int n, double *s, int *position);

double find_peak_value (int n, double *s);

int corr (double *s, double *p, int nphase);

int def_off_pulse (int nphase, double *in, double frac_off);

int off_pulse (int nphase, int index, double *in, double *out, double frac_off);

int remove_baseline (double *in, int index, double frac_off, int n, double *out);

int pre_diff (double *s, int nphase, int index, double frac_off, double *s_out);

int InitialGuess (double *s, double *p, int nphase);

int preA7_QUV (double *p, int nphase, int nchn, double *real_p, double *ima_p);

int rotate (int N, double *real_p, double *real_p_rotate, double *ima_p, double *ima_p_rotate, double rot);

int align (int N, double phase, double b, double a, double *real_p, double *real_p_align, double *ima_p, double *ima_p_align, double rotate);

int inverse_dft (double *real_p, double *ima_p, int ncount, double *p_new);
