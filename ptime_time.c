// calculate the phase shift between template and simulated (or real) data 
// error of phase can be calculated
// initial guess of phase shift is added
// try to do two-dim template matching
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "ptime_time.h"
#include "T2toolkit.h"
#include "tempo2pred.h"
#include "simulatePseudoBB.h"

int main (int argc, char *argv[])
{
	int h,i,j,z;

	if (argc != 7)
	{
		printf ("Usage: ptime_time -f fname -std tname (-pt tname) -o oname\n"
	            "Derive the TOAs "
	            "fname: data file; tname: templates; oname: output .tim; -std: standard template format; -pt: ptime template;\n");
	    exit (0);
	}


	//////////////////////////////////////////////////////
	char fname[128];   // name of data file
	char tname[128];   // name of template
	char oname[128];   // name of output .tim
	int mode; // to distinguish different type of templates

	for (i=0;i<argc;i++)
    {
		if (strcmp(argv[i],"-f") == 0)
		{
			strcpy(fname,argv[++i]);
		}
		else if (strcmp(argv[i],"-std")==0)
		{
			strcpy(tname,argv[++i]);
			mode = 0; // standard template format
			printf ("standard template format\n");
			//sscanf(argv[++i],"%d",&nbin);
		}
		else if (strcmp(argv[i],"-pt")==0)
		{
			strcpy(tname,argv[++i]);
			mode = 1; // ptime template
			printf ("ptime template format\n");
		}
		else if (strcmp(argv[i],"-o")==0)
		{
			strcpy(oname,argv[++i]);
		}
    }

	// name of different extension
	char name_data[50]; 
	char name_predict[50]; 
	char name_psrparam[50]; 

	strcpy(name_data,fname);
	strcpy(name_predict,fname);
	strcpy(name_psrparam,fname);

	char data[] = "[SUBINT]";
	char predict[] = "[T2PREDICT]";
	char psrparam[] = "[PSRPARAM]";

	strcat(name_data, data);
	strcat(name_predict, predict);
	strcat(name_psrparam, psrparam);

	//puts(name_data);
	//puts(name_predict);
	////////////////////////////////////////////////////
	
	double psrfreq;
	psrfreq = read_psrfreq(name_psrparam);
	printf ("%.15lf\n", psrfreq);
	
	////////////////////////////////////////////////
	long int imjd, smjd;
	double offs;
    int nphase;
	int nchn;
	int nsub;
	int npol;
	
	imjd = stt_imjd(fname);
	smjd = stt_smjd(fname);
	offs = stt_offs(fname);

    nchn = get_nchan(name_data);	
    npol = get_npol(name_data);	
    nsub = get_subint(name_data);	
    nphase = get_nphase(name_data);	

	//printf ("%d\n", nchn);
	////////////////////////////////////////////////

	// read a std
	char std[50];
	strcpy(std,tname);
	if ( mode == 0)
	{
		strcat(std, data);
	}
	//puts(argv[1]);
	//puts(argv[2]);
	double s_multi[nphase*nchn*npol];
	double s_temp[nphase];
	//double tt[nphase];
	//int n;

	//readfile(argv[1],&n,tt,s);
	//read_prof(std,1,s_multi,nphase);
	
	// check the channel number of template
	check_std(std,1,mode,nchn);

	read_std(std,1,s_multi,nphase,mode,nchn);

	/*
	for (i = 0; i < nphase*nchn*npol; i++)
	{
		printf ("%d %lf\n", i, s_multi[i]);
	}
	exit (0);

	int i;
	for (i = 0; i < nphase; i++)
	{
	    printf ("%lf\n", s[i]);
	}
	//puts(argv[1]);
	*/
	/////////////////////////////////////////////////////////////////////////////////
	FILE *fp;
	if ((fp = fopen(oname, "w+")) == NULL)
	{
        fprintf (stdout, "Can't open file\n");
		exit(1);
	}
	
	fprintf (fp, "FORMAT 1\n");

	////////////////////////////////////////////////////////////////////////////////

	double p_multi[nchn*npol*nphase];
	double p_temp[nphase];
    //double SNR; 

	double rms[nchn];  // rms for each profile
	double phase, e_phase;
	long double dt, e_dt;  
	long double t;     // TOA
	double offset;   // offset of each subint
	long double mjd0;  // the mjd of each subint
	T2Predictor pred;
	int ret;
	double period, frequency, weight;
	double freq[nchn], wts[nchn];
	for (h = 1; h <= nsub; h++)
	{
	    //////////////////////////////////////////////////////////////////////////
	    // simulate data

		//SNR = 500.0 + 200.0*i;
	    //simulate(n,SNR,s,p_temp);

	    read_prof(name_data,h,p_multi,nphase);
	    //readfile(argv[2],&n,tt,p_multi);

		for (i = 0; i < nchn; i++)
		{
			for (j = 0; j < nphase; j++)
			{
				//printf ("%lf %lf\n", p_multi[j], s[j]);
				//s_multi[i*nphase + j] = s[j];
				p_temp[j] = p_multi[i*nphase + j];
				s_temp[j] = s_multi[i*nphase + j];
			}

			// calculate toa, rms for each profile
			rms[i] = get_toa(s_temp, p_temp, psrfreq, nphase);
		}

		// do template matching, get the phase shift
		get_toa_multi(s_multi, p_multi, rms, nchn, &phase, &e_phase, psrfreq, nphase);

		////////////////////////////////////////////////////////////////////////////////////////

		// transform phase shift to TOAs
		// get the freq of the subint
		read_freq(name_data, h, freq, nchn);
		read_wts(name_data, h, wts, nchn);
		frequency = 0.0;
		weight = 0.0;
		for (z = 0; z < nchn; z++)
		{
			frequency += freq[z]*wts[z];
			weight += wts[z];
		}
		frequency = frequency/weight;
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
        e_dt = ((long double)(e_phase)/PI)*((long double)(period))/2.0L;
	    printf ("dt is %.10Lf +/- %.10Lf\n", dt, e_dt);

		// calculate the TOA
        t = (long double)(imjd) + ((long double)(smjd) + (long double)(offs) - (long double)(dt) + (long double)(offset))/86400.0L;
        //t = imjd;
		
	    printf ("offset is %lf\n", offset);
		fprintf (fp, "%s  %lf  %.15Lf  %Lf  7\n", fname, frequency, t, e_dt*1e+6);
	}

    if (fclose (fp) != 0)
		fprintf (stderr, "Error closing\n");


	return 0;
}
