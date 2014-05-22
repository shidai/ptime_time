// read PSRFITS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fitsio.h"
#include "ptime_time.h"
#include "ptime.h"

long int stt_imjd ( char *name )
{  
//double *read_arrival_time( char *input, long *nrows )
    fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
    int status;

    status = 0;

    if ( fits_open_file(&fptr, name, READONLY, &status) )          // open the file
    {
        printf( "error while openning file\n" );
    }

	long int imjd;
    if ( fits_read_key(fptr, TLONG, (char *)"STT_IMJD", &imjd, NULL, &status) )           // get the row number
    {
        printf( "error while getting the npol number\n" );
		//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
	}
    printf ("STT_IMJD: %ld\n", imjd);
	///////////////////////////////////////////////////////////////////////////

    if ( fits_close_file(fptr, &status) )
    {
        printf( " error while closing the file " );
    }

    return imjd;
}

long int stt_smjd ( char *name )
{  
//double *read_arrival_time( char *input, long *nrows )
    fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
    int status;

    status = 0;

    if ( fits_open_file(&fptr, name, READONLY, &status) )          // open the file
    {
        printf( "error while openning file\n" );
    }

	long int smjd;
    if ( fits_read_key(fptr, TLONG, (char *)"STT_SMJD", &smjd, NULL, &status) )           // get the row number
    {
        printf( "error while getting the npol number\n" );
		//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
	}
    printf ("STT_SMJD: %ld\n", smjd);
	///////////////////////////////////////////////////////////////////////////

    if ( fits_close_file(fptr, &status) )
    {
        printf( " error while closing the file " );
    }

    return smjd;
}

double stt_offs ( char *name )
{  
//double *read_arrival_time( char *input, long *nrows )
    fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
    int status;

    status = 0;

    if ( fits_open_file(&fptr, name, READONLY, &status) )          // open the file
    {
        printf( "error while openning file\n" );
    }

    double stt_offs;
    if ( fits_read_key(fptr, TDOUBLE, (char *)"STT_OFFS", &stt_offs, NULL, &status) )           // get the row number
    {
        printf( "error while getting the npol number\n" );
		//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
	}
    printf ("STT_OFFS: %lf\n", stt_offs);
	///////////////////////////////////////////////////////////////////////////

    if ( fits_close_file(fptr, &status) )
    {
        printf( " error while closing the file " );
    }

    return stt_offs;
}

int get_nchan ( char *name )
{  
//double *read_arrival_time( char *input, long *nrows )
    fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
    int status;

    status = 0;

    if ( fits_open_file(&fptr, name, READONLY, &status) )          // open the file
    {
        printf( "error while openning file\n" );
    }

	int nchan;
    if ( fits_read_key(fptr, TINT, (char *)"NCHAN", &nchan, NULL, &status) )           // get the row number
    {
        printf( "error while getting the npol number\n" );
		//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
	}
    printf ("number of nchan: %d\n", nchan);
	///////////////////////////////////////////////////////////////////////////

    if ( fits_close_file(fptr, &status) )
    {
        printf( " error while closing the file " );
    }

    return nchan;
}

int get_npol ( char *name )
{  
//double *read_arrival_time( char *input, long *nrows )
    fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
    int status;

    status = 0;

    if ( fits_open_file(&fptr, name, READONLY, &status) )          // open the file
    {
        printf( "error while openning file\n" );
    }

	//////////////////////////////////////////////////////////////////////////
	int npol;
    if ( fits_read_key(fptr, TINT, (char *)"NPOL", &npol, NULL, &status) )           // get the row number
    {
        printf( "error while getting the npol number\n" );
		//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
	}
    printf ("number of npol: %d\n", npol);

	///////////////////////////////////////////////////////////////////////////

    if ( fits_close_file(fptr, &status) )
    {
        printf( " error while closing the file " );
    }

    return npol;
}

int get_nphase ( char *name )
{  
//double *read_arrival_time( char *input, long *nrows )
    fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
    int status;

    status = 0;

    if ( fits_open_file(&fptr, name, READONLY, &status) )          // open the file
    {
        printf( "error while openning file\n" );
    }

	//////////////////////////////////////////////////////////////////////////
	int nbin;
    if ( fits_read_key(fptr, TINT, (char *)"NBIN", &nbin, NULL, &status) )           // get the row number
    {
        printf( "error while getting the nbin number\n" );
		//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
	}
    printf ("number of nbin: %d\n", nbin);

	///////////////////////////////////////////////////////////////////////////

    if ( fits_close_file(fptr, &status) )
    {
        printf( " error while closing the file " );
    }

    return nbin;
}

int get_subint ( char *name )
{  
//double *read_arrival_time( char *input, long *nrows )
    fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
    int status;
    long int nrows;

    status = 0;

    if ( fits_open_file(&fptr, name, READONLY, &status) )          // open the file
    {
        printf( "error while openning file\n" );
    }

    if ( fits_get_num_rows(fptr, &nrows, &status) )           // get the row number
    {
        printf( "error while getting the row number\n" );
    }
    printf ("number of subint: %ld\n", nrows);
    
	///////////////////////////////////////////////////////////////////////////

    if ( fits_close_file(fptr, &status) )
    {
        printf( " error while closing the file " );
    }

    return nrows;
}

int print_t2pred ( char *name )
{  
    fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
    int status;
    int colnum;
    long int nrows;

    status = 0;

    if ( fits_open_file(&fptr, name, READONLY, &status) )          // open the file
    {
        printf( "error while openning file\n" );
    }

    if ( fits_get_num_rows(fptr, &nrows, &status) )           // get the row number
    {
        printf( "error while getting the row number\n" );
    }
    //printf ("%ld\n", nrows);
    
    //if ( fits_get_colnum(fptr, CASEINSEN, "TSUBINT", &colnum, &status) )           // get the row number
    if ( fits_get_colnum(fptr, CASEINSEN, "PREDICT", &colnum, &status) )           // get the row number
    {
        printf( "error while getting the colnum number\n" );
		//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
	}
    //printf ("%d\n", colnum);

    int frow;
    int	felem = 1;
    int nelem = 1;
    int	anynull = 0;
    char nval[]="NULL";

	char **line;
	line = (char **)malloc(sizeof(char *));
	line[0] = (char *)malloc(sizeof(char)*1024);

	// print predictor to t2pred.dat
	FILE *fp;
	char filename[]="t2pred.dat";
		    
	if ((fp = fopen(filename, "w+")) == NULL)
	{
		fprintf (stdout, "Can't open file\n");
		exit(1);
	}

	int i;
    for (i = 1; i <= nrows; i++)                             // print the results
	{
		frow = i;

		fits_read_col(fptr, TSTRING, colnum, frow, felem, nelem, nval, line, &anynull, &status);           // read the column
	//fits_read_col_str(fptr,colnum,frow,felem,nelem,nval,line,&anynull,&status);

		//puts(line[0]);
        //printf("%s\n", line[0]);
		fprintf (fp, "%s\n", line[0]);
	}

    if (fclose (fp) != 0)
	    fprintf (stderr, "Error closing\n");


    if ( fits_close_file(fptr, &status) )
    {
        printf( " error while closing the file \n" );
    }

    return 0;
}

int read_freq ( char *name, int subint, double *freq, int nchan)
//int main (int argc, char *argv[] )
{  
    fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
    int status;
    int colnum;

    status = 0;

    if ( fits_open_file(&fptr, name, READONLY, &status) )          // open the file
    //if ( fits_open_file(&fptr, argv[1], READONLY, &status) )          // open the file
    {
        printf( "error while openning file\n" );
    }

    if ( fits_get_colnum(fptr, CASEINSEN, "DAT_FREQ", &colnum, &status) )           // get the colnum number
    {
        printf( "error while getting the colnum number\n" );
		//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
	}
    //printf ("%d\n", colnum);

    int frow;
    int	felem = 1;
    int nelem = nchan;
    int null = 0;
    int	anynull = 0;

	//int subint = 1;
	//int nchan = 8;
	//double freq[nchan];
	frow = subint;

    fits_read_col(fptr, TDOUBLE, colnum, frow, felem, nelem, &null, freq, &anynull, &status);           // read the column

	//int i;
    //for (i = 0; i < nchan; i++)                             // print the results
	//{
		//puts(line[0]);
    //    printf("%lf\n", freq[i]);
		//fprintf (fp, "%s\n", line[0]);
	//}

    if ( fits_close_file(fptr, &status) )
    {
        printf( " error while closing the file \n" );
    }

    return 0;
}

int read_wts ( char *name, int subint, double *wts, int nchan)
//int main (int argc, char *argv[] )
{  
    fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
    int status;
    int colnum;

    status = 0;

    if ( fits_open_file(&fptr, name, READONLY, &status) )          // open the file
    //if ( fits_open_file(&fptr, argv[1], READONLY, &status) )          // open the file
    {
        printf( "error while openning file\n" );
    }

    if ( fits_get_colnum(fptr, CASEINSEN, "DAT_WTS", &colnum, &status) )           // get the colnum number
    {
        printf( "error while getting the colnum number\n" );
		//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
	}
    //printf ("%d\n", colnum);

    int frow;
    int	felem = 1;
    int nelem = nchan;
    int null = 0;
    int	anynull = 0;

	//int subint = 1;
	//int nchan = 8;
	//double wts[nchan];
	frow = subint;

    fits_read_col(fptr, TDOUBLE, colnum, frow, felem, nelem, &null, wts, &anynull, &status);           // read the column

	//int i;
    //for (i = 0; i < nchan; i++)                             // print the results
	//{
		//puts(line[0]);
    //    printf("%lf\n", wts[i]);
		//fprintf (fp, "%s\n", line[0]);
	//}

    if ( fits_close_file(fptr, &status) )
    {
        printf( " error while closing the file \n" );
    }

    return 0;
}

double read_offs ( char *name, int subint)
//int main (int argc, char *argv[] )
{  
    fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
    int status;
    int colnum;
    long int nrows;

    status = 0;

    if ( fits_open_file(&fptr, name, READONLY, &status) )          // open the file
    {
        printf( "error while openning file\n" );
    }

    if ( fits_get_num_rows(fptr, &nrows, &status) )           // get the row number
    {
        printf( "error while getting the row number\n" );
    }
    //printf ("%ld\n", nrows);
    
    if ( fits_get_colnum(fptr, CASEINSEN, "OFFS_SUB", &colnum, &status) )           // get the colnum number
    {
        printf( "error while getting the colnum number\n" );
		//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
	}
    //printf ("%d\n", colnum);

    int frow;
    int	felem = 1;
    int nelem = 1;
    int null = 0;
    int	anynull = 0;

    double offs;

	//int i;
    //for (i = 1; i <= nrows; i++)                             // print the results
	//{
		frow = subint;

        fits_read_col(fptr, TDOUBLE, colnum, frow, felem, nelem, &null, &offs, &anynull, &status);           // read the column

		//puts(line[0]);
        //printf("subint offset is %Lf\n", offs);
		//fprintf (fp, "%s\n", line[0]);
	//}

    if ( fits_close_file(fptr, &status) )
    {
        printf( " error while closing the file \n" );
    }

    return offs;
}

//int main ( int argc, char *argv[] )
int read_prof ( char *name, int subint, double *profile, int nphase)
{  
//double *read_arrival_time( char *input, long *nrows )
    //int subint = 1;
	//double profile[8*1024];
    fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
    int status;
    int colnum;
    long int nrows;

    status = 0;

    //if ( fits_open_file(&fptr, argv[1], READONLY, &status) )          // open the file
    if ( fits_open_file(&fptr, name, READONLY, &status) )          // open the file
    {
        printf( "error while openning file\n" );
    }

    if ( fits_get_num_rows(fptr, &nrows, &status) )           // get the row number
    {
        printf( "error while getting the row number\n" );
    }
    //printf ("%ld\n", nrows);
    
    //if ( fits_get_colnum(fptr, CASEINSEN, "TSUBINT", &colnum, &status) )           // get the row number
    if ( fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status) )           // get the row number
    {
        printf( "error while getting the colnum number\n" );
		//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
	}
    //printf ("%d\n", colnum);

	//////////////////////////////////////////////////////////////////////////
	int npol;
    if ( fits_read_key(fptr, TINT, (char *)"NPOL", &npol, NULL, &status) )           // get the row number
    {
        printf( "error while getting the npol number\n" );
		//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
	}
    //printf ("%d\n", npol);

	int nchan;
    if ( fits_read_key(fptr, TINT, (char *)"NCHAN", &nchan, NULL, &status) )           // get the row number
    {
        printf( "error while getting the npol number\n" );
		//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
	}
    //printf ("%d\n", nchan);
	///////////////////////////////////////////////////////////////////////////

	int nbin;
    int frow;
    int felem;
    int nelem;
    int null;
    int anynull;
    //double *profile;     // the array to store the profile   

	nbin = nphase;
    //profile = ( double *)malloc( (nchan*npol*nbin) * sizeof( double ) );               // allocate space for column value
    frow = subint;
    felem = 1;
    nelem = nbin*nchan*npol;
    //nelem = 1024;
    null = 0;
    anynull = 0;

    fits_read_col(fptr, TDOUBLE, colnum, frow, felem, nelem, &null, profile, &anynull, &status);           // read the column

	//int i;
    //for (i = 0; i < (nchan*npol*nbin); i++)                             // print the results
    //    printf("%d %lf \n", i, profile[i]);

    if ( fits_close_file(fptr, &status) )
    {
        printf( " error while closing the file " );
    }

    return 0;
}


double read_psrfreq ( char *name )
//int main(int argc, char *argv[])
{
    fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
    int status;

    status = 0;

    //if ( fits_open_file(&fptr, argv[1], READONLY, &status) )          // open the file
    if ( fits_open_file(&fptr, name, READONLY, &status) )          // open the file
    {
        printf( "error while openning file\n" );
    }

	//////////////////////////////////////////////////////////////////////////
	
	char freq[100];
	char F0[100];

	int colnum = 1;
    int frow = 7;
    int	felem = 1;
    int nelem = 1;
    int	anynull = 0;
    char nval[]="NULL";

	char **line;
	line = (char **)malloc(sizeof(char *));
	line[0] = (char *)malloc(sizeof(char)*1024);

	fits_read_col(fptr, TSTRING, colnum, frow, felem, nelem, nval, line, &anynull, &status);           // read the column

	//puts(line[0]);

	int nchar = strlen(line[0]);
	//printf("strlen %d\n", nchar);

	int i;
	for (i = 0; i < nchar; i++)
	{
		F0[i] = line[0][i];
	}
	//printf("F0 %s\n", F0);

	double a[10];
	int l = 0;
	int j = 0;
	for (i = 0; i < nchar; i++)
	{
		if( (F0[i] >= '0' && F0[i] <= '9') || F0[i] =='.' ) 
		{ 
			freq[l] = F0[i];
			l++;
		}
		else if(l > 0)
		{
			freq[l] = '\0';
			a[j]=atof(freq);
			j++;
			l=0;
		}
	}

	/*
	for(i = 0; i < j; i++)
	{
        printf("%.15lf\n", a[i]);
	}
	*/
 
	double psrfreq;
	psrfreq = a[1];

    if ( fits_close_file(fptr, &status) )
    {
        printf( " error while closing the file\n " );
    }

	return psrfreq;
}

int check_std ( char *name, int subint, int mode, int nchn, int nphase)
// check if the template has right number of channel, and return the number of channel
{  
	// currently, npol should be 1
	// nchn is the number of channel of the data profile
	if (mode == 0)
	{
		fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
		int status;
		int colnum;
		status = 0;

		//if ( fits_open_file(&fptr, argv[1], READONLY, &status) )          // open the file
		if ( fits_open_file(&fptr, name, READONLY, &status) )          // open the file
		{
			printf( "error while openning file\n" );
		}

		if ( fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status) )           // get the colnum number
		{
			printf( "error while getting the colnum number\n" );
			//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
		}
		//printf ("%d\n", colnum);

		//////////////////////////////////////////////////////////////////////////
		int npol;
		if ( fits_read_key(fptr, TINT, (char *)"NPOL", &npol, NULL, &status) )           // get the pol number
		{
			printf( "error while getting the npol number\n" );
			//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
		}
		//printf ("%d\n", npol);

		int nchan;
		if ( fits_read_key(fptr, TINT, (char *)"NCHAN", &nchan, NULL, &status) )           // get the chan number
		{
			printf( "error while getting the npol number\n" );
			//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
		}
		printf ("STD %d\n", nchan);
		///////////////////////////////////////////////////////////////////////////

        //////////////////////////////////////////////////////////////////////////
		int nbin;
		if ( fits_read_key(fptr, TINT, (char *)"NBIN", &nbin, NULL, &status) )           // get the row number
		{
			printf( "error while getting the nbin number\n" );
		   //fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
		}
		printf ("number of nbin of std: %d\n", nbin);
		
	    ///////////////////////////////////////////////////////////////////////////
		//
		if ( fits_close_file(fptr, &status) )
		{
			printf( " error while closing the file\n" );
		}

		// check if nbin = nphase
		if ( nbin != nphase )
		{
			printf ("The phase bin number of template != the phase bin number of data\n");
	        exit (0);
	    }
				
		// output the template
		// check if nchan = nchn
		if ( nchan == nchn )
		{
			printf ("The channel number of template = the channel number of data\n");
		}
		else 
		{
			printf ("The channel number of template != the channel number of data\n");

			if ( nchan != 1 )
			{
				printf ("Can't not do template matching! The channel number of template should be 1 or equal to that of data.\n");
				exit (0);
			}
			else 
			{
				printf ("The template has only one channel. The same template will be used for all channels.\n");
			}
		}


	}
	else if (mode == 1)
	{
		tmplStruct tmpl;
		initialiseTemplate(&tmpl);
		printf("Reading template\n");
		readTemplate_ptime(name,&tmpl);
	    printf("Complete reading template\n");

		if ( tmpl.nchan == nchn )
		{
			printf ("The channel number of template = the channel number of data\n");
		}
		else 
		{
			printf ("The channel number of template != the channel number of data\n");

			if ( tmpl.nchan != 1 )
			{
				printf ("Can't not do template matching! The channel number of template should be 1 or equal to that of data.\n");
				exit (0);
			}
			else 
			{
				printf ("The template has only one channel. The same template will be used for all channels.\n");
			}
		}


	}

	return 0;
}

int read_std ( char *name, int subint, double *profile, int nphase, int mode, int nchn)
{  
	int i,j;
	// currently, npol should be 1
	// nchn is the number of channel of the data profile
	if (mode == 0)
	{
		fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
		int status;
		int colnum;
		long int nrows;

		status = 0;

		//if ( fits_open_file(&fptr, argv[1], READONLY, &status) )          // open the file
		if ( fits_open_file(&fptr, name, READONLY, &status) )          // open the file
		{
			printf( "error while openning file\n" );
		}

		if ( fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status) )           // get the colnum number
		{
			printf( "error while getting the colnum number\n" );
			//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
		}
		//printf ("%d\n", colnum);

		//////////////////////////////////////////////////////////////////////////
		int npol;
		if ( fits_read_key(fptr, TINT, (char *)"NPOL", &npol, NULL, &status) )           // get the pol number
		{
			printf( "error while getting the npol number\n" );
			//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
		}
		//printf ("%d\n", npol);

		int nchan;
		if ( fits_read_key(fptr, TINT, (char *)"NCHAN", &nchan, NULL, &status) )           // get the chan number
		{
			printf( "error while getting the npol number\n" );
			//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
		}
		//printf ("%d\n", nchan);
		///////////////////////////////////////////////////////////////////////////

		// read the data
		int nbin;
		int frow;
		int felem;
		int nelem;
		int null;
		int anynull;

		//double *temp;     // the array to store the data   
		//temp = ( double *)malloc( (nchan*npol*nbin) * sizeof( double ) );               // allocate space for column value

		nbin = nphase;
		frow = subint;
		felem = 1;
		nelem = nbin*nchan*npol;
		//nelem = 1024;
		null = 0;
		anynull = 0;

		printf ("%d\n", nbin);
		fits_read_col(fptr, TDOUBLE, colnum, frow, felem, nelem, &null, profile, &anynull, &status);           // read the column
		//fits_read_col(fptr, TDOUBLE, colnum, frow, felem, nelem, &null, temp, &anynull, &status);           // read the column

		if ( fits_close_file(fptr, &status) )
		{
			printf( " error while closing the file\n" );
		}

		if ( nchan == 1 )
		{
			for (i = 1; i < nchn; i++)
			{
				for (j = 0; j < nphase; j++)
				{
					profile[i*nphase+j] = profile[j];
				}
			}
		}
	}
	else if (mode == 1)
	{
		tmplStruct tmpl;
		initialiseTemplate(&tmpl);
		printf("Reading template\n");
		readTemplate_ptime(name,&tmpl);
	    printf("Complete reading template\n");

		int i,j;
		double phi;
		for (i = 0; i < tmpl.nchan; i++)
		{
			for (j = 0; j < nphase; j++)
			{
				phi = j/(double)nphase;
				profile[i*nphase+j] = (double)evaluateTemplateChannel(&tmpl,phi,i,0,0);
			}
		}

		if ( tmpl.nchan == 1 )
		{
			for (i = 1; i < nchn; i++)
			{
				for (j = 0; j < nphase; j++)
				{
					profile[i*nphase+j] = profile[j];
				}
			}
		}
	}

	return 0;
}

