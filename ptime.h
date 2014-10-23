// ptime header file
// G. Hobbs, S. Dai
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

typedef struct vMises {
  double concentration; // Concentration of each component
  double height;        // Height of each component
  double centroid;      // Centroid of each component
  double concentration_err; // error of Concentration 
  double height_err;        // error of Height
  double centroid_err;      // error of Centroid
} vMises;

typedef struct component {
  int Comp;            
  int nVm;             // Number of components for each channel for each Stokes
  vMises *vonMises;
  int vmMemoryAllocated;
  int nVmAllocated;
} component;

/*
typedef struct component {
  double concentration; // Concentration of each component
  double height;        // Height of each component
  double centroid;      // Centroid of each component
  double concentration_err; // error of Concentration 
  double height_err;        // error of Height
  double centroid_err;      // error of Centroid
} component;
*/

typedef struct polStruct {
  int stokes;            // 1 = I, 2 = Q, 3 = U, 4 = V
  int nComp;             // Number of components for each channel for each Stokes
  int allVm;
  component *comp;
  int compMemoryAllocated;
  int nCompAllocated;
} polStruct;

typedef struct channelStruct {
  int nstokes;             // 1 = I, 4 = I,Q,U,V
  double freqLow;
  double freqHigh;
  polStruct *pol; 
  int polMemoryAllocated;
  int nPolAllocated;
} channelStruct;

typedef struct tmplStruct {
  // Common to all templates
  char dte[1024];      // Date template made
  char user[1024];   // Person who made the template
  float templateVersion; // Version of template header
  char source[128]; // Source name
  char profileFile[1024]; // Profile file name
  char units[1024];   // Unit definition
  double dedispersed; // = 0 by default

  int nchan; // Number of channels
  channelStruct *channel; // Channels
  int channelMemoryAllocated;
  int nChannelAllocated;
} tmplStruct;


// Initialises the template, but does not allocate memory for the profiles
// This routine should be called at the start of the program
void initialiseTemplate(tmplStruct *tmpl)
{
  strcpy(tmpl->dte,"UNSET");
  strcpy(tmpl->user,"UNSET");
  tmpl->templateVersion = 0;
  strcpy(tmpl->source,"UNSET");
  strcpy(tmpl->profileFile,"UNSET");
  strcpy(tmpl->units,"UNSET");
  tmpl->dedispersed = 0;
  tmpl->nchan = 0;
  tmpl->channelMemoryAllocated=0;
  tmpl->nChannelAllocated=0;
}

// Reads a template from disk
void readTemplate_ptime(char *file,tmplStruct *tmpl)
{
  FILE *fin;
  int i;
  char line[4096];
  char firstword[4096];
  char dummy[4096];
  int nchan=-1;
  int nstokes=1;
  int chan,stokes,comp,ivm,icomp;
  char stokesStr[1024];
  double f1,f2;

  // Read primary header
  if (!(fin = fopen(file,"r"))){
    printf("Unable to open file: >%s<\n",file);
    exit(1);
  }
  while (!feof(fin))
    {
      if (fgets(line,4096,fin) != NULL){
	if (line[0] == '#') // Comment line
	  {
	    // Do nothing
	  }
	else {
	  sscanf(line,"%s",firstword);
	  if (strcasecmp(firstword,"TEMPLATE_VERSION:")==0)
		{
	    sscanf(line,"%s %f",dummy,&(tmpl->templateVersion));
			//printf("%f\n", dummy);
		}
	  else if (strcasecmp(firstword,"SOURCE:")==0)
	    sscanf(line,"%s %s",dummy,(tmpl->source));
	  else if (strcasecmp(firstword,"PROFILE_FILE:")==0)
	    sscanf(line,"%s %s",dummy,(tmpl->profileFile));
	  else if (strcasecmp(firstword,"DATE:")==0)
	    sscanf(line,"%s %s",dummy,(tmpl->dte));
	  else if (strcasecmp(firstword,"UNITS:")==0)
	    sscanf(line,"%s %s",dummy,(tmpl->units));
	  else if (strcasecmp(firstword,"ID:")==0)
	    sscanf(line,"%s %s",dummy,(tmpl->user));
	  else if (strcasecmp(firstword,"DM_CORRECTION:")==0)
	    sscanf(line,"%s %lf",dummy,&(tmpl->dedispersed));
	  else if (strcasecmp(firstword,"NCHAN:")==0)
	    sscanf(line,"%s %d",dummy,&nchan);
	  else if (strcasecmp(firstword,"STOKES:")==0)
	    {
	      sscanf(line,"%s %s",dummy,stokesStr);
	      if (strcmp(stokesStr,"I")==0)
		nstokes=1;
	      else if (strcmp(stokesStr,"Q")==0 || strcmp(stokesStr,"U")==0 || strcmp(stokesStr,"V")==0)
		nstokes=4;
	    }
	}
      }
    }
  fclose(fin);
  // Do some checks
  if (nchan < 0){
    printf("Have not defined any channels. Unable to continue\n");
    exit(1);
  }
  // Allocate memory for these channels
  tmpl->nchan = nchan;
  if (tmpl->channelMemoryAllocated == 0){
    if (!(tmpl->channel = (channelStruct *)malloc(sizeof(channelStruct)*nchan))){
      printf("ERROR in allocated memory for channels\n");
      exit(1);
    }
    tmpl->channelMemoryAllocated = 1;
    tmpl->nChannelAllocated = nchan;
    for (i=0;i<nchan;i++)
      tmpl->channel[i].polMemoryAllocated = 0;
  }

  chan = -1;
  stokes = -1;
  comp=-1;
  // Now read the data
  if (!(fin = fopen(file,"r"))){
    printf("Unable to open file: >%s<\n",file);
    exit(1);
  }
  while (!feof(fin))
    {
      if (fgets(line,4096,fin)!=NULL)
			{
				if (line[0] == '#') // Comment line
				{
					// Do nothing
				}
				else 
				{
					sscanf(line,"%s",firstword);
					if (strcasecmp(firstword,"STOKES:")==0)
					{
						sscanf(line,"%s %s",dummy,stokesStr);
						if (strcmp(stokesStr,"I")==0)
							stokes=0;
						else if (strcmp(stokesStr,"Q")==0)
							stokes = 1;
						else if (strcmp(stokesStr,"U")==0)
							stokes = 2;
						else if (strcmp(stokesStr,"V")==0)
							stokes = 3;
					} 
					else if (strcasecmp(firstword,"FREQUENCY_RANGE:")==0)
					{
						if (stokes==0)
							chan++;
						sscanf(line,"%s %lf %lf",dummy,&f1,&f2);
						if (f1 < f2)
						{
							tmpl->channel[chan].freqLow = f1;
							tmpl->channel[chan].freqHigh = f2;
							//printf("%lf %lf\n",tmpl->channel[chan].freqLow,tmpl->channel[chan].freqHigh);
						} 
						else 
						{
							tmpl->channel[chan].freqLow = f2;
							tmpl->channel[chan].freqHigh = f1;
						}
						tmpl->channel[chan].nstokes = nstokes;
						// Allocate memory
						if (tmpl->channel[chan].polMemoryAllocated==0)
						{
							if (!(tmpl->channel[chan].pol = (polStruct *)malloc(sizeof(polStruct)*nstokes)))
							{
								printf("Error in allocated memory for Stokes\n");
								exit(1);
							}
		
							tmpl->channel[chan].polMemoryAllocated = 1;
							for (i=0;i<nstokes;i++)
								tmpl->channel[chan].pol[i].compMemoryAllocated = 0;
						}
	      //	    }
						tmpl->channel[chan].nPolAllocated = nstokes;
				}
				else if (strcasecmp(firstword,"NCOMP:")==0)
				{
					int ncomp;
					sscanf(line,"%s %d",dummy,&ncomp);
					tmpl->channel[chan].pol[stokes].nComp = ncomp;
					tmpl->channel[chan].pol[stokes].stokes = stokes;
					if (tmpl->channel[chan].pol[stokes].compMemoryAllocated==0)
					{
						if (!(tmpl->channel[chan].pol[stokes].comp = (component *)malloc(sizeof(component)*ncomp)))
						{
							printf("Error in allocated memory for components\n");
							exit(1);
						}
		
						tmpl->channel[chan].pol[stokes].compMemoryAllocated = 1;
					}
					tmpl->channel[chan].pol[stokes].nCompAllocated = ncomp;
					//printf ("%d\n",ncomp);
				}
				else if (strcasecmp(firstword,"NVonMises:")==0)
				{
					int nAllVm;
					sscanf(line,"%s %d",dummy,&nAllVm);
					tmpl->channel[chan].pol[stokes].allVm = nAllVm;
					//printf ("%d\n",nAllVm);
	      
					comp=0;
					icomp=0;
					ivm=0;
				}
				else // Look for the number of Von Mises for each component
				{
					char substr[4096];
					strcpy(substr,firstword);
					substr[4]='\0';
	      
					if (strcasecmp(substr,"VonM")==0)
					{
						//printf ("%d\n",comp);
						tmpl->channel[chan].pol[stokes].comp[comp].vmMemoryAllocated=0;
						sscanf(line,"%s %d",dummy, &(tmpl->channel[chan].pol[stokes].comp[comp].nVm));
						//printf ("COMP%d has %d Von Mises functions\n",comp+1, tmpl->channel[chan].pol[stokes].comp[comp].nVm);
						if (tmpl->channel[chan].pol[stokes].comp[comp].vmMemoryAllocated==0)
						{
							if (!(tmpl->channel[chan].pol[stokes].comp[comp].vonMises = (vMises *)malloc(sizeof(vMises)*tmpl->channel[chan].pol[stokes].comp[comp].nVm)))
							{
								printf("Error in allocated memory for components\n");
								exit(1);
							}
							tmpl->channel[chan].pol[stokes].comp[comp].vmMemoryAllocated = 1;
							//printf ("%d\n",tmpl->channel[chan].pol[stokes].comp[comp].vmMemoryAllocated);
						}
						tmpl->channel[chan].pol[stokes].comp[comp].nVmAllocated = tmpl->channel[chan].pol[stokes].comp[comp].nVm;
						comp++;
					}
					else if (strcasecmp(substr,"COMP")==0)
					{
						//printf ("%d\n",tmpl->channel[chan].pol[stokes].comp[icomp].nVm);
						if (ivm != tmpl->channel[chan].pol[stokes].comp[icomp].nVm-1)
						{
							sscanf(line,"%s %lf %lf %lf %lf %lf %lf",dummy,
									&(tmpl->channel[chan].pol[stokes].comp[icomp].vonMises[ivm].height),
									&(tmpl->channel[chan].pol[stokes].comp[icomp].vonMises[ivm].height_err),
									&(tmpl->channel[chan].pol[stokes].comp[icomp].vonMises[ivm].concentration),
									&(tmpl->channel[chan].pol[stokes].comp[icomp].vonMises[ivm].concentration_err),
									&(tmpl->channel[chan].pol[stokes].comp[icomp].vonMises[ivm].centroid),
									&(tmpl->channel[chan].pol[stokes].comp[icomp].vonMises[ivm].centroid_err));
							//printf("%lf %lf %lf\n",tmpl->channel[chan].pol[stokes].comp[icomp].vonMises[ivm].height,tmpl->channel[chan].pol[stokes].comp[icomp].vonMises[ivm].concentration,tmpl->channel[chan].pol[stokes].comp[icomp].vonMises[ivm].centroid);
							//printf ("%d %d\n",ivm, icomp);
							ivm++;
						}
						else 
						{
							sscanf(line,"%s %lf %lf %lf %lf %lf %lf",dummy,
									&(tmpl->channel[chan].pol[stokes].comp[icomp].vonMises[ivm].height),
									&(tmpl->channel[chan].pol[stokes].comp[icomp].vonMises[ivm].height_err),
									&(tmpl->channel[chan].pol[stokes].comp[icomp].vonMises[ivm].concentration),
									&(tmpl->channel[chan].pol[stokes].comp[icomp].vonMises[ivm].concentration_err),
									&(tmpl->channel[chan].pol[stokes].comp[icomp].vonMises[ivm].centroid),
									&(tmpl->channel[chan].pol[stokes].comp[icomp].vonMises[ivm].centroid_err));
							//printf ("%d %d\n",ivm,icomp);
							icomp++;
							ivm = 0;
						}
		  
						//&(tmpl->channel[chan].pol[stokes].comp[comp].height),
						//&(tmpl->channel[chan].pol[stokes].comp[comp].height_err),
						//&(tmpl->channel[chan].pol[stokes].comp[comp].concentration),
						//&(tmpl->channel[chan].pol[stokes].comp[comp].concentration_err),
						//&(tmpl->channel[chan].pol[stokes].comp[comp].centroid),
						//&(tmpl->channel[chan].pol[stokes].comp[comp].centroid_err));
					}
				}
			}
		}
  }
  fclose(fin);
}

// Evaluate a single template component
double evaluateTemplateComponent(tmplStruct *tmpl,double phi,int chan,int stokes,int comp,double phiRot)
{
  double result=0;
  //result = tmpl->channel[chan].pol[stokes].comp[comp].height *
  int k;
  for (k=0;k<tmpl->channel[chan].pol[stokes].comp[comp].nVm;k++)
  	result += fabs(tmpl->channel[chan].pol[stokes].comp[comp].vonMises[k].height) *
    	exp(tmpl->channel[chan].pol[stokes].comp[comp].vonMises[k].concentration*
	(cos((phi - tmpl->channel[chan].pol[stokes].comp[comp].vonMises[k].centroid - phiRot)*2*M_PI)-1));
  return result;
}

// Evaluate a given frequency channel and polarisation
double evaluateTemplateChannel(tmplStruct *tmpl,double phi,int chan,int stokes,double phiRot)
{
  double result=0;
  int k;
  for (k=0;k<tmpl->channel[chan].pol[stokes].nComp;k++)
    result += evaluateTemplateComponent(tmpl,phi,chan,stokes,k,phiRot);
  return result;
}

// Allocate specific amount of memory
void allocateMemoryTemplateDefault(tmplStruct *tmpl,int nchan,int npol,int ncomp,int nvm)
{
  int i,j,k;
  tmpl->channel = (channelStruct *)malloc(sizeof(channelStruct)*nchan);
  tmpl->channelMemoryAllocated = 1;
  tmpl->nChannelAllocated = nchan;

  for (i=0;i<nchan;i++)
    {
      tmpl->channel[i].pol = (polStruct *)malloc(sizeof(polStruct)*npol);
      tmpl->channel[i].polMemoryAllocated = 1;
      tmpl->channel[i].nPolAllocated = npol;

			for (j=0;j<npol;j++)
			{
				tmpl->channel[i].pol[j].comp = (component *)malloc(sizeof(component)*ncomp);
				tmpl->channel[i].pol[j].compMemoryAllocated = 1;
				tmpl->channel[i].pol[j].nCompAllocated = ncomp;

				for (k=0;k<ncomp;k++)
				{
					tmpl->channel[i].pol[j].comp[k].vonMises = (vMises *)malloc(sizeof(vMises)*nvm);
					tmpl->channel[i].pol[j].comp[k].vmMemoryAllocated = 1;
					tmpl->channel[i].pol[j].comp[k].nVmAllocated = nvm;
				}
			}
    }
}


void saveTemplate(char *fname,tmplStruct *tmpl)
{
  FILE *fout;
  char version[128] = "1.0";
  int i,j,k,h;
  double centre;

  if (!(fout = fopen(fname,"w")))
    {
      printf("Unable to open file %s\n",fname);
      return;
    }
  fprintf(fout,"TEMPLATE_VERSION: %s\n",version);
  fprintf(fout,"PROFILE_FILE: %s\n",tmpl->profileFile);
  fprintf(fout,"SOURCE: %s\n",tmpl->source);
  fprintf(fout,"DATE: %s\n",tmpl->dte);
  fprintf(fout,"ID: %s\n",tmpl->user);
  fprintf(fout,"UNITS: %s\n",tmpl->units);
  fprintf(fout,"NCHAN: %d\n",tmpl->nchan);

  for (i=0;i<tmpl->nchan;i++)
    {
      fprintf(fout,"#\n");
      for (j=0;j<tmpl->channel[i].nstokes;j++)
	{
	  fprintf(fout,"#\n");
	  if (tmpl->channel[i].pol[j].stokes == 0)
	    fprintf(fout,"STOKES: I\n");
	  else if (tmpl->channel[i].pol[j].stokes == 1)
	    fprintf(fout,"STOKES: Q\n");
	  else if (tmpl->channel[i].pol[j].stokes == 2)
	    fprintf(fout,"STOKES: U\n");
	  else if (tmpl->channel[i].pol[j].stokes == 3)
	    fprintf(fout,"STOKES: V\n");
	  fprintf(fout,"FREQUENCY_RANGE: %f %f\n",tmpl->channel[i].freqLow,tmpl->channel[i].freqHigh);
	  fprintf(fout,"NCOMP: %d\n",tmpl->channel[i].pol[j].nComp);
		fprintf(fout,"NVonMises: %d\n", tmpl->channel[i].pol[j].allVm);
  
		for (k=0;k<tmpl->channel[i].pol[j].nComp;k++)
		{
			fprintf(fout,"VonM: %d\n", tmpl->channel[i].pol[j].comp[k].nVm);
		}

	  for (k=0;k<tmpl->channel[i].pol[j].nComp;k++)
	  {
	  	for (h=0;h<tmpl->channel[i].pol[j].comp[k].nVm;h++)
		{
	    		centre = tmpl->channel[i].pol[j].comp[k].vonMises[h].centroid;
	    		if (centre > 1) centre-=1;
	    		if (centre < 0) centre+=1;
	    		//fprintf(fout,"COMP%d%d: %g %g %g %g %g %g\n",k+1,h+1,tmpl->channel[i].pol[j].comp[k].vonMises[h].height,
							fprintf(fout,"COMP%d: %.4lf %.6lf %.4lf %.6lf %.8lf %.10lf\n",k+1,tmpl->channel[i].pol[j].comp[k].vonMises[h].height,
		    tmpl->channel[i].pol[j].comp[k].vonMises[h].height_err,
		    tmpl->channel[i].pol[j].comp[k].vonMises[h].concentration,
		    tmpl->channel[i].pol[j].comp[k].vonMises[h].concentration_err,
		    centre,
		    tmpl->channel[i].pol[j].comp[k].vonMises[h].centroid_err);
	    //fprintf(fout,"COMP%d: %g %g %g\n",k+1,tmpl->channel[i].pol[j].comp[k].height,
	//	    tmpl->channel[i].pol[j].comp[k].concentration,
	//	    centre);
		}
	  }
	}
    }
  fclose(fout);
}
