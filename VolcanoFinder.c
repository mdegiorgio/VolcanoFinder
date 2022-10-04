#include <stdarg.h>
#include <math.h>
#include "VolcanoFinder.h"
#include "my_rand.h"
#include "freq.h"
#include "sort.h"
//#include "matrix.h"

double *log_fact = NULL;
double **log_binom = NULL;

//stores the probabilities over a grid of values of alpha*d
//indexed by prob[n][x][gridsize_dvalue][2][gridsize]
//gridsize_dvalue is the number of values such that theta, 2theta, ..., theta*Do
//gridsize is hard-coded, not same as argument to program
//prob[n][x][d][0][g] holds probabilities, 
//prob[n][x][d][1][g] holds vector that spline function needs
double *****prob=NULL;
double *ads=NULL; //holds vector of alpha*d values used to get prob

int isblock = 0, block = 1, nblock=1, block_start = 0, block_stop = 0; // To break into blocks for trivial parallelization
int introgression_model = 1;
int invar=0;
int datasize=0, sweep_width, nmax, nmin, xmax, closest_data;
double lowbound[MAXPAR], upbound[MAXPAR];
double minFreqPos=-1.0, maxFreqPos=-1.0, est_theta, est_thetaL;
struct datatype *data=NULL;

int gridsize_dvalue = 1; // Sets size of Divergence-value grid
double *dvalue_grid; // Divergence-value grid

double gen_dist = 0.0; // DIVERGENCE WITH INTROGRESSED POPULATION
int curr_d = 0;
double D0 = 0.0; // DIVERGENCE WITH OUTGROUP POPULATION

int alpha_grid_flag = 1;
double *vals_int;
double *likes_int;


void InitFact(void) 
{	
	int n = 0;
	
	log_fact = (double*)malloc((nmax+1)*sizeof(double));
	log_fact[0] = 0.0;
	for(n = 1; n <= nmax; n++) {
		log_fact[n] = log_fact[n-1] + log((double)n);
	}
}

void InitBinom(void)
{
	int n = 0;
	int k = 0;
 
	printf("Initializing binomial coefficients\n");

	log_binom = (double**)malloc((nmax+1)*sizeof(double*));
	for(n = 0; n <= nmax; n++) {
		log_binom[n] = (double*)malloc((nmax+1)*sizeof(double));
		
		for(k = 0; k <= n; k++) {
			log_binom[n][k] = log_fact[n] - log_fact[n-k] - log_fact[k];
		}
	}
}

void getlims(void) {
  int i;
  invar=0;
  if (datasize==0) {
    nmin=nmax=0;
    return;
  }
  nmin = nmax = data[0].n;
  minFreqPos = data[0].loc;
  maxFreqPos = minFreqPos;
  for (i=0; i<datasize; i++) {
    if(data[i].loc < minFreqPos) minFreqPos = data[i].loc;
    if(maxFreqPos < data[i].loc) maxFreqPos = data[i].loc;
    if (data[i].n>nmax) nmax=data[i].n;
    if (data[i].n<nmin) nmin=data[i].n;
    if (invar==0 && data[i].x==data[i].n) invar=1;
    if (invar!=2 && data[i].x==0 && data[i].n > 0) invar=2;
    if (invar!=2 && data[i].x==data[i].n && data[i].n > 0 && data[i].folded) invar=2;
  }
  if (invar) xmax=nmax+1;
  else xmax=nmax;
}

// variables for storing user-defined grid
int flagUserDefinedGrid = 0;
int datasize_grid=0;
double gridposmin, gridposmax;
double *data_grid=NULL;

// obtains the limits for the user-defined grid
void getlims_grid()
{
	int i;

	if(datasize_grid==0) {
		gridposmin=gridposmax=0.0;

		return;
	}

	gridposmin = data_grid[0];
	gridposmax = data_grid[0];

	for(i=0; i<datasize_grid; i++) {
		if(data_grid[i] > gridposmax) {
			gridposmax=data_grid[i];
		}

		if(data_grid[i]<gridposmin) {
			gridposmin=data_grid[i];
		}
	}
}


FILE *my_fopen(char *fn, char *mode) {
  FILE *rv = fopen(fn, mode);
  if (rv==NULL) {
    fprintf(stderr, "error opening %s mode \"%s\"\n", fn, mode);
    exit(-1);
  }
  return rv;
}


void get_nlim() {
  int i;
  if (datasize==0) {
    nmax=nmin=0;
    return;
  }
  nmax=nmin=data[0].n;
  for (i=1; i<datasize; i++) {
    if (data[i].n > nmax) nmax=data[i].n;
    if (data[i].n < nmin) nmin=data[i].n;
  }
}


void get_lims() {
  fflush(stdout);
  get_nlim();
  if (invar) xmax = nmax+1;
  else xmax = nmax;
}

#define NCOLUMN 4
#define LOC 0
#define X 1
#define N 2
#define FOLDED 3

//all column headers except folded are required.
//if folded is not present, all SNPs will be assumed unfolded
char colname[NCOLUMN][10]={"position", "x", "n", "folded"};
char outcolname[4][10]={"location", "LR", "alpha", "D"};

//get rid of some data so we can test missing data
void losedata(double per) {
  struct datatype *newdata = malloc(datasize*sizeof(struct datatype));
  int newdatasize=0;
  int i,j , x, n;
  for (i=0; i<datasize; i++) {
    x = data[i].x;
    n = data[i].n;
    for (j=0; j<data[i].n; j++) {
      if (uniform()  < per) {
	if (j < data[i].x) x--;
	n--;
      }
    }
    if (x!=0 && x!=n) {
      newdata[newdatasize].x = x;
      newdata[newdatasize].n = n;
      newdata[newdatasize].folded=data[i].folded;
      newdata[newdatasize].loc = data[i].loc;
      newdatasize++;
    }
  }
  if (datasize > 0) free(data);
  datasize = newdatasize;
  data = newdata;
  data = realloc(newdata, datasize*sizeof(struct datatype));
  printf("done losedata datasize=%i\n", datasize);
}

void readms_error(char *infn) {
  printf("error reading data from %s\n", infn);
  exit(-1);
}


int readsnps_ms(char *infn) {
  static FILE *infile=NULL;
  static char *curr_infile=NULL;
  static int numchr=0, rep=0;
  char str[1000], c;
  int i, j;

  if (infn!=curr_infile) {
    if (infile!=NULL) 
      fclose(infile);
    printf("opening %s\n", infn);
    rep=0;
    infile=fopen(infn, "r");
    curr_infile=infn;
    if (EOF==fscanf(infile, "%s %i", str, &numchr)) readms_error(infn);
    if (strstr(str, "ms")==NULL || numchr <=0) readms_error(infn);
    printf("numchr=%i\n", numchr);
  }
  while (1) {
    if (EOF==fscanf(infile, "%s", str)) {
      fclose(infile);
      infile=NULL;
      curr_infile=NULL;
      return 0;
    }
    if (strcmp(str, "segsites:")==0) break;
  }
  if(EOF==fscanf(infile, "%i",&datasize)) readms_error(infn);
  printf("reading %s rep %i\n", infn, ++rep);
  printf("datasize=%i\n", datasize);
  if (datasize<=0) {
    printf("readsnps_ms: segsites must be greater than 0!\n");
    readms_error(infn);
  }
  data = malloc(datasize*sizeof(struct datatype));
  if(EOF==fscanf(infile, "%s", str)) readms_error(infn);
  if(strcmp(str, "positions:")!=0) readms_error(infn);
  for (i=0; i<datasize; i++) {
    if (EOF==fscanf(infile, "%lf", &data[i].loc))
      readms_error(infn);
    data[i].folded=0;
    data[i].x=0;
    data[i].n=numchr;
  }
  for (i=0; i<numchr; i++) {
    while ('\n'!=(c=fgetc(infile))) if (isspace(c)==0) readms_error(infn);
    for (j=0; j<datasize; j++) {
      c=fgetc(infile);
      if (c!='0' && c!='1') readms_error(infn);
      data[j].x += (c=='1');
    }
  }
  for (i=0; i<datasize; i++) {
    if (data[i].x==0 || data[i].x>=data[i].n) {
      printf("invalid data in msfile %s data[%i].x=%i n=%i\n", infn,
	     i, data[i].x, data[i].n);
      readms_error(infn);
    }
  }
  getlims();
  return 1;
}

//this is the format the old version of yuseob's program outputs
//it doesn't output n anywhere so you have to know what it is
void readsnps_yuseob(char *infn, int rep, int n) {
  FILE *infile = my_fopen(infn, "r");
  char str[1000], c;
  int count=0, nsites, i;
  if (datasize > 0) free(data);
  datasize=0;
  while (1) {
    while (EOF!=fscanf(infile, "%s", str))
      if (strcmp(str, "sites")==0) break;
    if (strcmp(str, "sites")!=0) break;
    while (':'!=(c=fgetc(infile))) assert(isspace(c));
    fscanf(infile, "%i", &nsites);
    if (rep < 0 || count==rep) {
      if (datasize==0) data = malloc(nsites*sizeof(struct datatype));
      else data = realloc(data, (datasize+nsites)*sizeof(struct datatype));
      for (i=0; i<nsites; i++) {
	assert(EOF!=fscanf(infile, "%lf %i", &data[datasize+i].loc,
			   &data[datasize+i].x));
	data[datasize+i].folded=0;
	data[datasize+i].n=n;
      }
      datasize += nsites;
      if (count==rep) break;
    }
    count++;
  }
  fclose(infile);
  losedata(0.1);
  getlims();
  printf("done getlims_yuseob datasize=%i nmax=%i nmin=%i xmax=%i invar=%i\n",
	 datasize, nmax, nmin, xmax, invar);
}


void printdata(char *outfn) {
  int i;
  FILE *outfile = fopen(outfn, "w");
  fprintf(outfile, "position\tx\tn\tfolded\n");
  for (i=0; i<datasize; i++)
    fprintf(outfile, "%f\t%i\t%i\t%i\n", data[i].loc,
	    data[i].x, data[i].n, data[i].folded);
  fclose(outfile);
}

void readsnps_error(char *infn) {
  printf("error reading %s\n", infn);
  exit(-1);
}

void readsnps(char *infn) {
  FILE *infile=my_fopen(infn, "r");
  int colpos[NCOLUMN], pos=0, col=0, i, j;
  char c, str[1000];
  if (datasize > 0) free(data);
  datasize=0;
  c=fgetc(infile);
  for (i=0; i<NCOLUMN; i++) colpos[i]=-1;
  while (c!='\n' && c!=EOF) {
    str[0]=c;
    pos=1;
    while ('\t'!=(c=fgetc(infile)) && c!='\n' && c!=EOF)
      str[pos++]=c;
    str[pos]='\0';
    for (i=0; i<NCOLUMN; i++)
      if (strcmp(str, colname[i])==0) {
	colpos[i]=col;
	break;
      }
    col++;
    if (c=='\n' || c==EOF) break;
    c=fgetc(infile);
  }
  if (colpos[LOC]<0 || colpos[X]<0 || colpos[N]<0) {
    fprintf(stderr, "readsnps: infile should have columns named position, x, and n (and optionally folded\n");
    exit(-1);
  }
  while (EOF!=(c=fgetc(infile)))
    if (c=='\n') datasize++;
  fclose(infile);
  infile=my_fopen(infn, "r");
  while ('\n'!=(c=fgetc(infile)) && c!=EOF);
  data = malloc(datasize*sizeof(struct datatype));
  for (i=0; i<datasize; i++) {
    for (j=0; j<col; j++) {
      if (colpos[LOC]==j) 
	{if(EOF==fscanf(infile, "%lf", &data[i].loc)) readsnps_error(infn);}
      else if (colpos[X]==j) 
	{if(EOF==fscanf(infile, "%i", &data[i].x)) readsnps_error(infn);}
      else if (colpos[N]==j) 
	{if(EOF==fscanf(infile, "%i", &data[i].n)) readsnps_error(infn);}
      else if (colpos[FOLDED]==j) 
	{if(EOF==fscanf(infile, "%i", &data[i].folded)) readsnps_error(infn);}
      else {
	while ('\t'!=(c=fgetc(infile)) && c!='\n' && c!=EOF);
	if (c=='\n' || c==EOF) if(j!=col-1) readsnps_error(infn);
      }
    }
  }
  fclose(infile);
  if (colpos[FOLDED] < 0)
    for (i=0; i<datasize; i++) 
      data[i].folded=0;

  
  for(i = 1; i < datasize; i++) {
    if(data[i].loc < data[i-1].loc) {
      fprintf(stderr, "readsnps: infile should have loci sorted in non-descending ordder\n");
      exit(-1);
    }
  }

  getlims();
  printf("done readsnps datasize=%i nmax=%i nmin=%i xmax=%i invar=%i\n", datasize, nmax, nmin, xmax, invar);
}

void readgrid(char *infn) { // Okay, we have read in the corret data file
  FILE *infile=my_fopen(infn, "r");
  int i;
  char c;
  if (datasize_grid > 0) free(data_grid);
  datasize_grid=0;
  c=fgetc(infile);
  while (EOF!=(c=fgetc(infile)))
    if (c=='\n') datasize_grid++;
  fclose(infile);
  infile=my_fopen(infn, "r");
  data_grid = malloc(datasize_grid*sizeof(double));
  for (i=0; i<datasize_grid; i++) {
    fscanf(infile, "%lf", &data_grid[i]);
  }
  fclose(infile);
  getlims_grid();

  printf("done readgrid datasize_grid=%i gridposmin=%le gridposmax=%le\n", datasize_grid, gridposmin, gridposmax);
}

void writeblocks(FILE *outfile, char *outfn, int currBlock, int numBlock) {
    double block_loc, block_LR, block_alpha, block_D;
    char infn[1000];
    sprintf(infn, "%s_%d_%d", outfn, currBlock, numBlock);
    FILE *infile=fopen(infn, "r");  
    int num_sites=0;  
    int i;
    char c;

	c=fgetc(infile);
    
	while(EOF!=(c=fgetc(infile))) {
		if(c=='\n') {
			num_sites++;
		}
	}

	fclose(infile);

    infile=fopen(infn, "r");

	for(i=0; i<num_sites; i++) {
        fscanf(infile, "%lf\t%lf\t%lf\t%lf\n", &block_loc, &block_LR, &block_alpha, &block_D);
        fprintf(outfile, "%f\t%f\t%e\t%f\n", block_loc, block_LR, block_alpha, block_D);
    }

	fclose(infile);
   //fclose(outfile);
}

void writeprobs(char *outfn, int gridsize) {
    FILE *outfile=NULL;  
    int n, x, d, j, k, minx, maxx;
    char curr_outfn[1000];
    
    sprintf(curr_outfn, "%s_dvalues", outfn);
    outfile = fopen(curr_outfn, "w");
    fprintf(outfile, "%d\n", gridsize_dvalue); // Write number of grid points
    for(d = 0; d < gridsize_dvalue; d++) {
        fprintf(outfile, "%lf\n", dvalue_grid[d]);
    }
    fclose(outfile);

    sprintf(curr_outfn, "%s_lookuptable", outfn);
    outfile = fopen(curr_outfn, "w");
    fprintf(outfile, "%d\t%d\t%d\t%d\n", invar, nmin, nmax, gridsize); // Write necessary information to reconstruct probability grid
    for(n = nmin; n <= nmax; n++) {
        if(invar == 0) {
            minx=1;
            maxx=n;
        }
        else if(invar == 1) {
            minx = 1;
            maxx = n+1;
        }
        else {
            minx = 0;
            maxx = n+1;
        }

        for(x = minx; x < maxx; x++) {                
            for(d = 0; d < gridsize_dvalue; d++) {
                for(j = 0; j < 2; j++) {
                    for(k = 0; k < gridsize; k++) {
                        fprintf(outfile, "%d\t%d\t%d\t%d\t%d\t%e\n", n, x, d, j, k, prob[n][x][d][j][k]);
                    }
                }
            }
        }
    }
	fclose(outfile);
}

void readprobs(char *infn) {
    FILE *infile=NULL;  
    int gridsize, n, x, d, j, k, minx, maxx, tempn, tempx, tempd, tempj, tempk, i, numlog, numrest;
    double tempprob, maxval, minval, interval;
    char curr_infn[1000];
    
    sprintf(curr_infn, "%s_dvalues", infn);
    infile = fopen(curr_infn, "r");
    fscanf(infile, "%d\n", &gridsize_dvalue); // Read number of grid points
    dvalue_grid = (double*)malloc(gridsize_dvalue*sizeof(double)); // Divergence-value grid
    for(d = 0; d < gridsize_dvalue; d++) {
        fscanf(infile, "%lf\n", &dvalue_grid[d]); // Read grid points
    }
    fclose(infile);

    sprintf(curr_infn, "%s_lookuptable", infn);
    infile = fopen(curr_infn, "r");
    fscanf(infile, "%d\t%d\t%d\t%d\n", &invar, &nmin, &nmax, &gridsize); // Read necessary information to reconstruct probability grid

    prob = malloc((nmax+1)*sizeof(double****));  // Added component for divergence D
    for(n = nmin; n <= nmax; n++) {
        if(invar == 0) {
            minx=1;
            maxx=n;
        }
        else if(invar == 1) {
            minx = 1;
            maxx = n+1;
        }
        else {
            minx = 0;
            maxx = n+1;
        }
        prob[n] = malloc(maxx*sizeof(double***));

        for(x = minx; x < maxx; x++) {      
            prob[n][x] = malloc(gridsize_dvalue*sizeof(double**));

            for(d = 0; d < gridsize_dvalue; d++) {
                prob[n][x][d] = malloc(2*sizeof(double*));

                for(j = 0; j < 2; j++) {
                    prob[n][x][d][j] = malloc(gridsize*sizeof(double));

                    for(k = 0; k < gridsize; k++) {
                        fscanf(infile, "%d\t%d\t%d\t%d\t%d\t%lf\n", &tempn, &tempx, &tempd, &tempj, &tempk, &tempprob);
                        
                        if(tempn != n || tempx != x || tempd != d || tempj != j || tempk != k) {
                            printf("ERROR: File name generated properly\n");
                            exit(-1);                         
                        }
        
                        prob[tempn][tempx][tempd][tempj][tempk] = tempprob;
                    }
                }
            }
        }
    }
	fclose(infile);

    // Sets ads grid structure
    ads = malloc(gridsize*sizeof(double));
    numlog = gridsize*0.5;
    numrest = gridsize-numlog;
    //this is a bit ad-hoc but i'm choosing values of ad that cover
    //the entire range of minval to maxval, but concentrate on the higher
    //values which are usually more relevant
    //so first 50% of values are chosen on a log scale, other rest on linear
    minval=1.0e-8;
    maxval=BIG_ENOUGH;
    interval = log(maxval / minval) / (double)(numlog - 1);
    for(i = 0; i < numlog; i++)
        ads[i] = exp(log(minval)+i*interval);
    minval=0.1;
    maxval=3.0;
    interval = (maxval - minval) / (double)(numrest - 1);
    for(i = 0; i < numrest; i++)
        ads[i + numlog] = minval + interval * i;
    dsort(gridsize, ads - 1);
    for(i = 1; i < gridsize; i++)
        if(ads[i] == ads[i - 1]) {
            printf("ads[%i]=ads[%i]=%e\n\n\n", i, i-1, ads[i]);
            exit(-1);
        }
}


void readGridDValues(char *infn) {
    FILE *infile=NULL;  
    int d;
    
    infile = fopen(infn, "r");
    fscanf(infile, "%d\n", &gridsize_dvalue); // Read number of grid points
    dvalue_grid = (double*)malloc(gridsize_dvalue*sizeof(double)); // Divergence-value grid
    for(d = 0; d < gridsize_dvalue; d++) {
        fscanf(infile, "%lf\n", &dvalue_grid[d]); // Read grid points
    }
    fclose(infile);
}


//probability of choosing j from sample size s given frequency 
//spectrum p
double p_j_s(int j, int s, double *p) { // Calculation of pj,H
  int i;
  double rv=0.0;
  if (j < 0 || j > s) return 0.0;
  for (i=j; i<=nmax; i++) // ADVISED: Formula in paper goes to n-1, this formula goes to n
    if (s-j >=0 && i<xmax) {
     // rv += p[i]*xchoosey(i,j)*xchoosey(nmax-i,s-j)/xchoosey(nmax,s);
      if(0 <= j && j <= i && 0 <= s-j && s-j <= nmax-i && 0 <= s && s <= nmax) {
          rv += p[i]*exp(log_binom[i][j] + log_binom[nmax-i][s-j] - log_binom[nmax][s]); // USING THIS INSTEAD TO SPEED UP CALCULATION
      }
    }
  return rv;
}
 
double p_kescapes(int k, double ad) { // Calculation of Pe(k)
  double pe = exp(-ad); // Strangely, this is just the opposite of what is given in the paper
  return exp( log_binom[nmax][k] ) *pow((1.0-pe),k)*pow(pe,(nmax-k)); // But this reverses the opposite and makes it correct
  //return xchoosey(nmax, k)*pow((1.0-pe),k)*pow(pe,(nmax-k)); // But this reverses the opposite and makes it correct
}

double get_pstar_sub_introgression(int i, double* p_un, double ad) { // Calculation of normalized S_i^*(n|d,s)
  int k;
  double muTMRCA = 0.0;
//  static double lastad=-1.0;
  static double *rv=NULL;
 //   D0=p_j_s(1, 1, p_un); // S_1(1)
    
    //printf("D0=%lf\n",D0);    
//  if (lastad!=ad) { // WHAT IS THIS?
    if (rv==NULL) rv = malloc((nmax+1)*sizeof(double));
    if(invar==0 && introgression_model == 1) { // MODEL 1
          rv[i] = (gen_dist/2.0)*( p_kescapes(i, ad) + p_kescapes(nmax-i, ad) ); // Last term of S_i^*(n|d,s) calculation
          for (k=i+1; k<=nmax; k++) { // First set of terms for S_i^*(n|d,s) calculation
	        rv[i] += p_kescapes(k, ad)*p_j_s(i, k, p_un);
          }
    }
    else if(invar==1 && introgression_model == 1) { // MODEL 1a (including invariant sites)
          if(D0 < gen_dist/2.0) {
            printf("Unfortunately D0 < D/2, and so model 2a cannot be used\n");
            exit(-1);
          }

          rv[i] = 0.0;
          if (i < nmax) {
   	        rv[i] = (gen_dist/2.0)*( p_kescapes(i, ad) + p_kescapes(nmax-i, ad) ); // Last term of S_i^*(n|d,s) calculation
            for (k=i+1; k<=nmax; k++) { // First set of terms for S_i^*(n|d,s) calculation
	            rv[i] += p_kescapes(k, ad)*p_j_s(i, k, p_un);
            }          
          }
          else if (i == nmax) {
            for (k=1; k<nmax; k++) { // First set of terms for S_n^*(n|d,s) calculation
	            rv[i] += p_kescapes(k, ad);
            }
            rv[i] *= (D0 - gen_dist/2.0); // Multiple by (D_0 - D/2)
            rv[i] += D0 * (p_kescapes(0, ad) + p_kescapes(nmax, ad));   // Last term of S_n^*(n|d,s) calculation
          }
    }
    else if(invar==0 && introgression_model == 2) { // MODEL 2 (accounting for coalescence time in recipient population)
          // Compute muTMRCA(n-1)
          for(k = 1; k <= nmax-2; k++) {
            muTMRCA += k*p_j_s(k, nmax-1, p_un);
            muTMRCA *= (1.0/(nmax-1.0));
          }          
        
          if(gen_dist/2.0 <= muTMRCA) {
            printf("Unfortunately D/2 <= muTMRCA(n-1), and so model 2 cannot be used\n");
            exit(-1);
          }

          // Compute muTMRCA(i)
          muTMRCA = 0.0;
          for(k = 1; k <= i-1; k++) {
            muTMRCA += k*p_j_s(k, i, p_un);
            muTMRCA *= (1.0/i);
          }

          rv[i] = (gen_dist/2.0 - muTMRCA)* p_kescapes(i, ad); // Last term of S_i^*(n|d,s) calculation
          rv[i] += (gen_dist/2.0) * p_kescapes(nmax-i, ad); // Middle term of S_i^*(n|d,s) calculation 
          for (k=i+1; k<=nmax; k++) { // First set of terms for S_i^*(n|d,s) calculation
	        rv[i] += p_kescapes(k, ad)*p_j_s(i, k, p_un);
          }
    }
    else if(invar==1 && introgression_model == 2) { // MODEL 2a (including invariant sites and accounting for coalesacence time in recipient population)
          // Compute muTMRCA(n)
          for(k = 1; k <= nmax-1; k++) {
            muTMRCA += k*p_j_s(k, nmax, p_un);
            muTMRCA *= (1.0/nmax);
          }          
        
          if(D0 < gen_dist/2.0) {
            printf("Unfortunately D0 < D/2, and so model 2a cannot be used\n");
            exit(-1);
          }

          if(gen_dist/2.0 <= muTMRCA) {
            printf("Unfortunately D/2 <= muTMRCA(n), and so model 2a cannot be used\n");
            exit(-1);
          }
          
          rv[i] = 0.0;
          if (i < nmax) {
            if(gen_dist/2.0 - D0 + p_j_s(i,i,p_un) >= 0) {
                rv[i] = ( gen_dist/2.0 - D0 + p_j_s(i, i, p_un) )* p_kescapes(i, ad); // Last term of S_i^*(n|d,s) calculation
            }
            rv[i] += (gen_dist/2.0) * p_kescapes(nmax-i, ad); // Middle term of S_i^*(n|d,s) calculation 
            for (k=i+1; k<=nmax; k++) { // First set of terms for S_i^*(n|d,s) calculation
	          rv[i] += p_kescapes(k, ad)*p_j_s(i, k, p_un);
            }          
          }
          else if (i == nmax) {
            for (k=1; k<nmax; k++) { // First set of terms for S_n^*(n|d,s) calculation
	            rv[i] += p_kescapes(k, ad);
            }
            rv[i] *= (D0 - gen_dist/2.0); // Multiple by (D_0 - D/2)
            rv[i] += D0 * p_kescapes(0, ad) + p_j_s(nmax, nmax, p_un) * p_kescapes(nmax, ad);   // Last term of S_n^*(n|d,s) calculation
          }
          
    }
    else {
        printf("ERROR: Only include polymorphisms or substitutions\n");
        exit(-1);
    }
//    lastad = ad;
 // }
//  printf("rv[%d]=%lf  norm_rv[%d]=%lf\n",i,rv[i],i,rv[i]/sum);
 //   printf("rv[%d]=%lf",i,rv[i]);
    return rv[i];
}

double get_pstar_introgression(int i, double* p, double ad) {
  double pr;//, sum;

  ad = ad + 1E-7; // Calculation becomes unstable when too close to the test site

//  sum=1.0; // If sum < 1, then conditioning on only variable sites
//  if (invar==0 || invar==1)
//    sum -= get_pstar_sub_introgression(0, p, ad); // Subtract on invariable sites (except fixed differences)
//  if (invar==0) 
//    sum -= get_pstar_sub_introgression(nmax, p, ad); // Subtract out fixed differences
//  assert(sum <= 1.0 && sum > 0.0);
//  pr = get_pstar_sub_introgression(i, p, ad)/sum; // Calculate the pB*, with B=i, conditioning on the types of sites observed (variable or invariable)
  pr = get_pstar_sub_introgression(i, p, ad); // Calculate the pB*, with B=i, conditioning on the types of sites observed (variable or invariable)
  return pr;
}

double splint(double* xvec, double *yvec, double *yvec2, int n, 
	      double x)
{
  int lowpos,hipos,pos;
  double diff,b,a;
  double y;
  
  if (x < xvec[0]) return xvec[0];
  lowpos=0;
  hipos=n-1;
  while (hipos-lowpos > 1) {
    pos=(hipos+lowpos)/2;
    if (xvec[pos] > x) hipos=pos;
    else lowpos=pos;
  }
  diff=xvec[hipos]-xvec[lowpos];
  if (diff == 0.0) {fprintf(stderr, "splint error\n"); exit(-1);}
  a=(xvec[hipos]-x)/diff;
  b=(x-xvec[lowpos])/diff;
  y=a*yvec[lowpos]+b*yvec[hipos]+((a*a*a-a)*yvec2[lowpos]+(b*b*b-b)*yvec2[hipos])*(diff*diff)/6.0;
  if (y<=0.0 || y>=1.0) {
    y=(yvec[lowpos]*a+b*yvec[hipos]);
  }
  if (y <0.0 || y>1.0) {
    printf("y=%e lowpos=%i hipos=%i a=%e b=%e yvec=%e, %e\n", 
	   y, lowpos, hipos, a, b, yvec[lowpos], yvec[hipos]);
    y = (yvec[lowpos]+yvec[hipos])/2.0;
  }
  return y;
}


void spline(double *x, double *y ,int n, double yp1, double ypn,
	    double *y2)
{
  int i,k;
  double p,qn,sig,un,*u;
  
  /*	for (i=1; i<=n; i++) {
	
  printf("%i\t%e\t%e\t%e\n", i, x[i], y[i], y2[i]);
  }
  exit(-1);*/
  
  u=malloc(n*sizeof(double));
  if (yp1 > 0.99e30)
    y2[0]=u[0]=0.0;
  else {
    y2[0] = -0.5;
    u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
    if (x[1]-x[0]==0.0) {
      printf("spline: x[1]=x[0]=%e\n", x[1]); exit(-1);}
    
  }
  for (i=1;i<=n-2;i++) {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    if (x[i+1]==x[i-1]) {
      printf("spline2: x[%i]=x[%i]=%e\n", i+1, i-1, x[i+1]);
      exit(-1);}
    p=sig*y2[i-1]+2.0;
    if (p==0.0) {
      printf("spline: p=%e\n", p); exit(-1);}
    y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    if (x[i+1]==x[i]) {
      printf("spline3: x[%i]=x[%i]=%e\n", i+1, i, x[i+1]); exit(-1);}
    if (x[i]==x[i-1]) {
      printf("spline4: x[%i]=x[%i]=%e\n", i, i-1, x[i]); exit(-1);}
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }
  if (ypn > 0.99e30)
    qn=un=0.0;
  else {
    printf("spline: Else!\n"); exit(-1);
    qn=0.5;
    un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
  }
  y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
  if (qn*y2[n-2]+1.0==0.0) {
    printf("spline: qn=%e y2[%i]=%e\n", qn, n-2, y2[n-2]); exit(-1);}
  for (k=n-2;k>=0;k--)
    y2[k]=y2[k]*y2[k+1]+u[k];
  free(u);
}

void calcprobs_introgression(double *p_un, int gridsize) {
  double pr, maxval, minval, interval, sum;
  int i, n, x, d, j, k, minx, maxx, maxxi, xi, numlog, numrest;
  double *spectra_values;
  printf("calcprobs nmax=%i nmin=%i xmax=%i invar=%i\n",
	 nmax, nmin, xmax, invar);
  fflush(stdout);

  ads = malloc(gridsize*sizeof(double));
  numlog = gridsize*0.5;
  numrest = gridsize-numlog;
  //this is a bit ad-hoc but i'm choosing values of ad that cover
  //the entire range of minval to maxval, but concentrate on the higher
  //values which are usually more relevant
  //so first 50% of values are chosen on a log scale, other rest on linear
  minval=1.0e-8;
  maxval=BIG_ENOUGH;
  interval=log(maxval/minval)/(double)(numlog-1);
  for (i=0; i<numlog; i++)
    ads[i]=exp(log(minval)+i*interval);
  minval=0.1;
  maxval=3.0;
  interval=(maxval-minval)/(double)(numrest-1);
  for (i=0; i<numrest; i++)
    ads[i+numlog]=minval+interval*i;
  dsort(gridsize, ads-1);
  for (i=1; i<gridsize; i++)
    if (ads[i]==ads[i-1]) {
      printf("ads[%i]=ads[%i]=%e\n\n\n", i, i-1, ads[i]);
      exit(-1);
    }
  prob = malloc((nmax+1)*sizeof(double****));  // Added component for divergence D
  for (n=nmin; n<=nmax; n++) {
    if (invar==0) {minx=1; maxx=n;}
    else if (invar==1) {minx=1; maxx=n+1;}
    else {minx=0; maxx=n+1;}
    prob[n] = malloc(maxx*sizeof(double***));
    for (x=minx; x<maxx; x++) {        
      prob[n][x] = malloc(gridsize_dvalue*sizeof(double**));
      for(d = 0; d < gridsize_dvalue; d++) {
        prob[n][x][d] = malloc(2*sizeof(double*));
        for (j=0; j<2; j++) {
	      prob[n][x][d][j] = malloc(gridsize*sizeof(double));
	      for (k=0; k<gridsize; k++)
	        prob[n][x][d][j][k] = 0.0;
        }
      }
    }
  }

  spectra_values = malloc((nmax+1)*sizeof(double)); // STORE unnormalized spectra after sweep
  for(d = 0; d < gridsize_dvalue; d++) { // Added component for divergence D
    gen_dist = dvalue_grid[d]; // Set the genetic distance D
    for (k=0; k<gridsize; k++) { 
      sum=0.0; // ADDED
      for(x = 1; x < xmax; x++) { // ADDED
        spectra_values[x] = get_pstar_introgression(x, p_un, ads[k]); // Get unnoramlized spectra after sweep
        sum += spectra_values[x]; // ADDED
      }    // ADDED
      for (n=nmin; n<=nmax; n++) {
        if (invar==0) {minx=1; maxx=n;}
        else if (invar==1) {minx=1; maxx=n+1;}
        else {minx=0; maxx=n+1;}
        for (x=minx; x<maxx; x++) {
	      maxxi = x + nmax - n;
	      for (xi=x; xi<=maxxi; xi++) {
            pr = exp(log_binom[xi][x] + log_binom[nmax-xi][n-x] - log_binom[nmax][n]);
	        //pr = exp(ln_xchoosey(xi, x) + ln_xchoosey(nmax-xi, n-x) - ln_xchoosey(nmax, n));
            prob[n][x][d][0][k] += pr*spectra_values[xi]/sum; // DO NORMALIZATION ONLY ONCE (ADDED)  
	        //prob[n][x][d][0][k] += pr*get_pstar_introgression(xi, p_un, ads[k])/sum; // DO NORMALIZATION ONLY ONCE (ADDED)
          }
        }
      }
      
    }
  }
  free(spectra_values);

  for (n=nmin; n<=nmax; n++) {
    if (invar==0) {minx=1; maxx=n;}
    else if (invar==1) {minx=1; maxx=n+1;}
    else {minx=0; maxx=n+1;}
    for (x=minx; x<maxx; x++) {
      for(d = 0; d < gridsize_dvalue; d++) { // Added this to loop across divergence values D
        spline(ads, prob[n][x][d][0], gridsize, 1.0e30, 
	       1.0e30, prob[n][x][d][1]);
      }
    }
  }
  printf("done calcprob\n");
  fflush(stdout);
}

//actually returns LR now
// Uses genetic distance with introgressed population
// Uses unnormalized spectrum
double ln_likelihood_introgression(double *p_un, double alpha, double sweep) {
  int s, n, x, gridsize=500, first;//, start=closest_data, stop=closest_data;
  double like=0.0, dist, ad, pr;

  if (prob==NULL)
    calcprobs_introgression(p_un, gridsize); // Changed to calcprobs_introgression and uses unnormalized spectrum

//  for(s = closest_data; s >= 0; s--) { // Go backward from closest data point to find start of calculation
//    dist = fabs(data[s].loc-sweep);
//    if (dist < 0.01) dist = 0.01;  //don't want to be right on top of sweep
//    ad = alpha*dist;
//    if (ad >= BIG_ENOUGH) {
//      break;    
//    }
//    else {  
//      start = s;
//    }
//  }
  
//  for(s = closest_data; s < datasize; s++) { // Go forward from closest data point to find end of calculation
//    dist = fabs(data[s].loc-sweep);
//    if (dist < 0.01) dist = 0.01;  //don't want to be right on top of sweep
//    ad = alpha*dist;
//    if (ad >= BIG_ENOUGH) {
//      break;   
//    }
//    else { 
//      stop = s;
//    }
//  }
//  printf("(closest_data=%d, start=%d, stop=%d, datasize=%d),  (%lf, %lf, %lf),  alpha=%e,  dist=%e,  ad=%e,  BIG_ENOUGH=%e\n",closest_data,start,stop,datasize,data[closest_data].loc,data[start].loc,data[stop].loc,alpha,dist,ad,BIG_ENOUGH);
//  start=0;
//  stop=datasize-1;
  //for (s=0; s<datasize; s++) {
  sweep_width=0;
//  for (s=start; s<=stop; s++) { // Changed to only consider relevant data points
  for(s = closest_data; s >= 0; s--) { // Go backward from closest data point to find start of calculation
    dist = fabs(data[s].loc-sweep);
    if (dist < 0.01) dist = 0.01;  //don't want to be right on top of sweep
    ad = alpha*dist;
    if (ad >= BIG_ENOUGH) break;
    sweep_width++;
    n=data[s].n;
    x = data[s].x;
    pr =0.0;
   
  calc_backward:
    if (ad < ads[0]) 
      pr += prob[n][x][curr_d][0][0];
    else 
      pr += splint(ads, prob[n][x][curr_d][0],
		   prob[n][x][curr_d][1], 
		   gridsize, ad);
 
    if (pr >= 0.0 && pr < 1.0);
    else {
      //printf("1: val=%e\n", pr);  // FOR DEBUGGING
      //printf("%i %i %e\n", nmax-n, x, ad);  // FOR DEBUGGING
      //printf("%e %e\t%e %e\n", ads[gridsize-2], 
	  //   prob[n][x][curr_d][0][gridsize-2],
	  //   ads[gridsize-1],
	  //   prob[n][x][curr_d][0][gridsize-1]); // FOR DEBUGGING
      //      debug_splint=1;
      splint(ads, prob[n][x][curr_d][0],
	     prob[n][x][curr_d][1],
	     gridsize, ad);
    }
    
    if (pr < 0.0 || pr >= 1.00000001) {
        if(pr >= 1.00000001) { // Mod for errors in numerical precision
          pr = 1.0;
        } 
        else {
          printf("val=%e\n", pr); exit(-1);
        }
    }
    assert(pr >=0.0 && pr<1.00000001);
    
    if (first && data[s].folded) {
      x = n-x;
      first=0;
      if (x!=n-x)
	goto calc_backward;
    }
    like += log(pr) - data[s].baselike;
  }

  for(s = closest_data; s < datasize; s++) { // Go forward from closest data point to find end of calculation
    dist = fabs(data[s].loc-sweep);
    if (dist < 0.01) dist = 0.01;  //don't want to be right on top of sweep
    ad = alpha*dist;
    if (ad >= BIG_ENOUGH) break;
    sweep_width++;
    n=data[s].n;
    x = data[s].x;
    pr =0.0;
    first=1;
  calc_forward:
    if (ad < ads[0]) 
      pr += prob[n][x][curr_d][0][0];
    else 
      pr += splint(ads, prob[n][x][curr_d][0],
		   prob[n][x][curr_d][1], 
		   gridsize, ad);
    
    if (pr >= 0.0 && pr < 1.0);
    else {
      //printf("1: val=%e\n", pr);  // FOR DEBUGGING
      //printf("%i %i %e\n", nmax-n, x, ad);  // FOR DEBUGGING
      //printf("%e %e\t%e %e\n", ads[gridsize-2], 
	  //   prob[n][x][curr_d][0][gridsize-2],
	  //   ads[gridsize-1],
	  //   prob[n][x][curr_d][0][gridsize-1]); // FOR DEBUGGING
      //      debug_splint=1;
      splint(ads, prob[n][x][curr_d][0],
	     prob[n][x][curr_d][1],
	     gridsize, ad);
      
    }
    if (pr < 0.0 || pr >= 1.00000001) {
        if(pr >= 1.00000001) { // Mod for errors in numerical precision
          pr = 1.0;
        } 
        else {
          printf("val=%e\n", pr); exit(-1);
        }
    }
    assert(pr >=0.0 && pr<1.00000001);
    
    if (first && data[s].folded) {
      x = n-x;
      first=0;
      if (x!=n-x)
	goto calc_forward;
    }
    like += log(pr) - data[s].baselike;
  }

  return like;
}

//returns maximum likelihood for sweep at given location (sweep) 
//and also sets MLE for alpha
// assumes model with adaptive introgression
double findmax_alpha_introgression(double *alpha, double *pi, double sweep) {
  int startgridsize=100, gridsize, max=-1, i, sw[100], count=0;
  double tol=1.0e-6, minalpha, maxalpha, mind=-1, maxd=-1, 
    interval, like, totd=0;
  gridsize=startgridsize;
  if(alpha_grid_flag) { // If first time running findmax_alpha_introgression
    vals_int = malloc(gridsize*sizeof(double));
    likes_int = malloc(gridsize*sizeof(double));
    alpha_grid_flag = 0;
  }
  for (i=0; i<datasize; i++) {
    if (i==0 || fabs(data[i].loc-sweep) < mind)
      mind = fabs(data[i].loc-sweep);
    if (i==0 || fabs(data[i].loc-sweep) > maxd)
      maxd = fabs(data[i].loc-sweep);
    totd += fabs(data[i].loc-sweep);
  }
  if (mind<0.01) mind=0.01;
  maxalpha=(BIG_ENOUGH+1)/mind;
//  minalpha = BIG_ENOUGH/(totd/datasize);
  minalpha = BIG_ENOUGH/(totd/datasize)/10.0; // MODIFIED AND DIVIDED BY AN ORDER OF MAGNITUDE 
  //  printf("minalpha=%e maxalpha=%e\n", minalpha, maxalpha);
//int firstflag = 1;
//int redoflag = 0;
  // Ensures that the search is fine enough for a strong sweep
  while ((maxalpha-minalpha)/((maxalpha+minalpha)/2.0) > tol) {
    interval = log(maxalpha/minalpha)/(gridsize-1);
    for (i=0; i<gridsize; i++) {
      vals_int[i] = exp(log(minalpha)+i*interval);
      if (i!=0 && sw[i-1]==0) likes_int[i] = likes_int[i-1];
      else likes_int[i] = ln_likelihood_introgression(pi, vals_int[i], sweep);
      //printf("%i %f %e like=%f sweep_width=%i\n", i, sweep, vals_int[i], likes_int[i], sweep_width);
      sw[i] = sweep_width;
      if (i==0 || likes_int[i] > likes_int[max])
	    max = i;
    }
//    if(firstflag==1){
//        firstflag=0;
//        if(max < startgridsize*0.9){
           // printf("\tHELLO max is %d\n", max);
//            redoflag=1;
//            break;
//        }
//    }
    if (max==0) 
      minalpha = exp(log(minalpha)-gridsize*interval);
    else minalpha = vals_int[max-1];
    if (max==gridsize-1) maxalpha += (vals_int[max-1]-vals_int[max-2]);
    else maxalpha = vals_int[max+1];
    gridsize = 5;
    count++;
    /*    if (count > 10000) {
      printf("%i %e %e\n", count, minalpha, maxalpha);
      fflush(stdout);
      }*/
  }

//  double tempmaxalpha=(BIG_ENOUGH+1)/mind;
//  double tempminalpha = BIG_ENOUGH/(totd/datasize);
//  if(vals_int[max]<tempminalpha)
//    printf("NEGATIVE\tD=%lf\tlike=%lf\n", gen_dist, likes_int[max]);
//  else
//    printf("POSITIVE:    D=%lf\tlike=%lf\n", gen_dist, likes_int[max]);
    
if(0){
 //   double tempmaxalpha=(BIG_ENOUGH+1)/mind;
    double tempminalpha = BIG_ENOUGH/(totd/datasize);
 //   double tempinterval = log(maxalpha/minalpha)/(gridsize-1);
//  printf("\t\tMAX = %d,   COUNT = %d\n",max, count);
//    if(redoflag) {

    if(vals_int[max] <= tempminalpha) {
        printf("redoing: %e <= %e\n", vals_int[max], tempminalpha);
    //if(max <= -1) { // Alpha grid was not large enough scale (adaptively expand)
        gridsize=startgridsize;
        max = -1;
        count = 0;
        maxalpha=(BIG_ENOUGH+1)/mind;
        minalpha = BIG_ENOUGH/(totd/datasize)/10.0; // Adaptively expand by an order of magnitude smaller
  
        // Ensures that the search is fine enough for a strong sweep
        while ((maxalpha-minalpha)/((maxalpha+minalpha)/2.0) > tol) {
            interval = log(maxalpha/minalpha)/(gridsize-1);
            for (i=0; i<gridsize; i++) {
                vals_int[i] = exp(log(minalpha)+i*interval);
                if (i!=0 && sw[i-1]==0) likes_int[i] = likes_int[i-1];
                else likes_int[i] = ln_likelihood_introgression(pi, vals_int[i], sweep);
                //printf("%i %f %e like=%f sweep_width=%i\n", i, sweep, vals_int[i], likes_int[i], sweep_width);
                sw[i] = sweep_width;
                if (i==0 || likes_int[i] > likes_int[max])
                    max = i;
            }
            if (max==0) 
                minalpha = exp(log(minalpha)-gridsize*interval);
            else minalpha = vals_int[max-1];
            if (max==gridsize-1) maxalpha += (vals_int[max-1]-vals_int[max-2]);
            else maxalpha = vals_int[max+1];
            gridsize = 5;
            count++;
            /*    if (count > 10000) {
            printf("%i %e %e\n", count, minalpha, maxalpha);
            fflush(stdout);
            }*/
        } 
    }    
}
  like = likes_int[max];
  sweep_width = sw[max];
  *alpha = vals_int[max];
//  free(vals_int);
//  free(likes_int);
  return like;
}

// Normalizes the frequency spectrum
// Takes p_un and normalizes to p
void NormalizeSpectrum(double *p, double *p_un)
{
	int i;
    double sum=0.0;

    for(i = 0; i < xmax; i++) {
        sum += p_un[i];
    }

    assert(sum > 0.0);    

    // Normalize
    for(i = 0; i < xmax; i++) {
        p[i] = p_un[i]/sum;
        assert(p[i] <= 1.0);
    }  
}

// Find introgression sweeps
// Uses unnormalizeds spectrum (p_un)
double find_introgression_sweeps(char *outfn, double *p_un, int gridsize, int noisy, int msrep) {
  double *p = malloc(nmax*sizeof(double)); // make new allele frequency spectrum
  double alpha, alpha_d, lr, lr_d, maxlr=0.0, minlike=0.0, smax=0, smin=1e100, sweep;
  int i, rep, maxd=0;
  FILE* outfile;
  double maxalpha=-1.0, maxsweep=-1.0;

  NormalizeSpectrum(p,p_un); // Computed the normalized spectrum for use when computing base likelihoods

  for (i=0; i<datasize; i++) {
    minlike += (data[i].baselike = 
		likelihood_freq_onesnp(data[i].x, data[i].n,
				       data[i].folded, 1, p)); // Calculate base likelihood
  }

  for (i=0; i<datasize; i++) {
    if (data[i].loc < smin || i==0) smin = data[i].loc;
    if (data[i].loc > smax || i==0) smax = data[i].loc;
  }
  printf("findsweeps smin=%e smax=%e gridsize=%i minlike=%f\n", 
	 smin, smax, gridsize, minlike);

  outfile=my_fopen(outfn, "w");
  fclose(outfile);

  if (noisy==1 && isblock==0) {
    outfile=my_fopen(outfn, "a");
    fprintf(outfile, "location\tLR\talpha\tD\n");
    fclose(outfile);
  }
  if (noisy==0 && msrep==1 && isblock==0) {
    outfile = my_fopen(outfn, "a");
    fprintf(outfile, "rep\tLR\tlocation\talpha\n");
    fclose(outfile);
  }

  closest_data = 0; // Closest data point to the test site is the first data point
//  closest_data = block_start; // Closest data point to the test site is the first data point of the block
//  for (rep=0; rep<gridsize; rep++) {
  for (rep=block_start; rep<=block_stop; rep++) { // Cycle through block of test sites  
    if (gridsize==1) { //gridsize=1 is not recommended
      if(flagUserDefinedGrid) { sweep = data_grid[0]; }        
      else { sweep=(smin+smax)/2.0; }
    }
    else {
      if(flagUserDefinedGrid) { sweep = data_grid[rep]; }        
      else { sweep = smin + (smax-smin)*rep/(gridsize-1); }
    }
    
    while( (closest_data+1 < datasize ) && ( fabs(data[closest_data+1].loc - sweep) <= fabs(data[closest_data].loc - sweep) ) ) { // Find closest data point to current sweep
      closest_data++;
    } 

    curr_d = 0;
    gen_dist = dvalue_grid[0];
    lr = findmax_alpha_introgression(&alpha, p_un, sweep); // Change to findmax_alpha_introgression    
    lr_d = lr;
    maxd = 0;
    alpha_d = alpha;    
    for(curr_d = 1; curr_d < gridsize_dvalue; curr_d++) { // Cycle through all D values
        gen_dist = dvalue_grid[curr_d];

        lr = findmax_alpha_introgression(&alpha, p_un, sweep); // Change to findmax_alpha_introgression    
 
        if (lr > lr_d) { // Changed to rep=block_start from rep=0
          lr_d = lr;
          maxd = curr_d;
          alpha_d = alpha;
        }    
    }
    if (rep==block_start || lr_d > maxlr) { // Changed to rep=block_start from rep=0
      maxlr = lr_d;
      maxalpha = alpha_d;
      maxsweep = sweep;
    }
    
    
    if (noisy) {
      outfile = my_fopen(outfn, "a");
      fprintf(outfile, "%f\t%f\t%e\t%f\n", sweep, lr_d, alpha_d, dvalue_grid[maxd]);
      fclose(outfile);
      printf("pos %f\tLR=%f\talpha=%e\tD=%f\n", sweep, lr_d, alpha_d, dvalue_grid[maxd]);
    }
  }
  if (noisy==0) {
    outfile=my_fopen(outfn, "a");
    fprintf(outfile, "%i\t%f\t%f\t%e\n", msrep, maxlr, maxsweep, maxalpha);
    fclose(outfile);
  }
  printf("maxsweep LR=%f loc=%f alpha=%e D=%f", maxlr, maxsweep, maxalpha, dvalue_grid[maxd]);
  if (msrep!=-1) printf("\tmsrep=%i\n", msrep);
  else fputc('\n', stdout);
 
  return maxlr;
}

void usage() {
  printf("usage:\n");
  printf("\tto pre-compute probability lookup table:\n");
  printf("\t\t./VolcanoFinder -p SpectFile D P MODEL nmin nmax xmin xmax LookupPrefix\n");
  printf("\tto find introgression sweeps:\n");
  printf("\t\t./VolcanoFinder -i G FreqFile SpectFile D P MODEL OutFile\n");
  printf("\t\t./VolcanoFinder -ig g FreqFile SpectFile D P MODEL OutFile\n");
  printf("\t\t./VolcanoFinder -iu GridFile FreqFile SpectFile D P MODEL OutFile\n");
  printf("\tto find introgression sweeps using pre-computed lookup table:\n");
  printf("\t\t./VolcanoFinder -pi G FreqFile SpectFile LookupPrefix OutFile\n");
  printf("\t\t./VolcanoFinder -pig g FreqFile SpectFile LookupPrefix OutFile\n");
  printf("\t\t./VolcanoFinder -piu GridFile FreqFile SpectFile LookupPrefix OutFile\n");
  printf("\tto find introgression sweeps in a given genomic block:\n");
  printf("\t\t./VolcanoFinder -bi G FreqFile SpectFile D P MODEL OutFile BLOCK NBLOCK\n");
  printf("\t\t./VolcanoFinder -big g FreqFile SpectFile D P MODEL OutFile BLOCK NBLOCK\n");
  printf("\tto find introgression sweeps in a given genomic block with a pre-computed lookup table:\n");
  printf("\t\t./VolcanoFinder -pbi G FreqFile SpectFile LookupPrefix OutFile BLOCK NBLOCK\n");
  printf("\t\t./VolcanoFinder -pbig g FreqFile SpectFile LookupPrefix OutFile BLOCK NBLOCK\n");
  printf("\tto merge scans across genomic blocks:\n");
  printf("\t\t./VolcanoFinder -m OutFile NBLOCK\n");
  exit(-1);
}

int main(int argc, char *argv[]) {
  double *p_un;
  float inignore;
  double space_between_grid_points=-1.0, block_interval_start, block_interval_stop;
  int gridsize=-1, isnoisy=1, i, my_xmin=0, inlen, isPol=1;
  char snpfn[1000], outfn[1000], freqfn[1000], gridfn[1000], lookupfn[1000], dvalue_gridfn[1000];
  FILE *outfile;
  if (argc < 3) usage();
//  if (strlen(argv[1])!=2 || argv[1][0]!='-')
 //   usage();
  if (strcmp(argv[1], "-p") == 0) {
    if (argc!=11) {
        usage();
    }
    else {
        printf("You have chosen to pre-compute the probability lookup table for VolcanoFinder\n");
    }

    sprintf(freqfn, "%s", argv[2]);

    isPol = atoi(argv[4]);    	
    introgression_model = atoi(argv[5]);
    if(introgression_model != 1 && introgression_model != 2) {
        printf("The MODEL should take values 1 or 2\n");
        usage();
    }

    nmin = atoi(argv[6]);
    nmax = atoi(argv[7]);
    my_xmin = atoi(argv[8]);
    xmax = atoi(argv[9]);    
    if(my_xmin > 0 && xmax < nmax) {
        xmax = nmax;
        invar = 0;
    }    
    else if(my_xmin > 0 && xmax == nmax) {
        xmax = nmax + 1;
        invar = 1;
    }    
    else if(my_xmin == 0 && xmax == nmax) {
        xmax = nmax + 1;
        invar = 2;
    }
   else {
        printf("ERROR: Invalid set constraints on xmin, xmax, nmin, and nmax\n");
        exit(-1);
    }

    sprintf(outfn, "%s", argv[10]);
    
    p_un = loadfreq_unnormalized(freqfn); // Loading S_j(n) for introgression sweeps
   
    InitFact();
	InitBinom();

    D0=p_j_s(1, 1, p_un); // S_1(1)
    est_theta=p_j_s(1, 2, p_un); // S_1(2)
    est_thetaL = 0.0;
    for(i = 1; i < nmax; i++) {
        est_thetaL += i * p_j_s(i, nmax, p_un); // i*S_i(n)
    } 
    est_thetaL =  est_thetaL / (nmax - 1.0);


    sprintf(dvalue_gridfn, "%s", argv[3]);
    if(sscanf(dvalue_gridfn, "%f %n", &inignore, &inlen) == 1 && !dvalue_gridfn[inlen]) { // Determine if number. If not, then read D from files
        gen_dist = atof(argv[3]);
        if(gen_dist >= 0.0) { // If positive D, then use the D inputted
            gridsize_dvalue = 1;
            dvalue_grid = (double*)malloc(sizeof(double)); // Divergence-value grid
            dvalue_grid[0] = gen_dist;

            // Check to see if user-defined D value is compatible with minimum
            if(isPol == 0) { // If fixed differences are not polarized                 
                if(introgression_model == 1 && (gen_dist < est_theta || gen_dist > D0)) {
                    printf("ERROR: D cannot be smaller than %e or greater than %e for Model 1 when not polarized\n", est_theta, D0);
                    exit(-1);
                }            
                else if(introgression_model == 2 && (gen_dist < 2.0*(nmax-1.0)*est_thetaL/nmax || gen_dist > D0)) {
                    printf("ERROR: D cannot be smaller than %e or greater than %e for Model 2 when not polarized\n", 2.0*(nmax-1.0)*est_thetaL/nmax, D0);
                    exit(-1);
                }
            }
            else { // If fixed differences are polarized 
                if(introgression_model == 1 && (gen_dist < est_theta || gen_dist > 2.0*D0)) {
                    printf("ERROR: D cannot be smaller than %e or greater than %e for Model 1 when polarized\n", est_theta, 2.0*D0);
                    exit(-1);
                }            
                else if(introgression_model == 2 && (gen_dist < 2.0*(nmax-1.0)*est_thetaL/nmax || gen_dist > 2.0*D0)) {
                    printf("ERROR: D cannot be smaller than %e or greater than %e for Model 2 when polarized\n", 2.0*(nmax-1.0)*est_thetaL/nmax, 2.0*D0);
                    exit(-1);
                }
            }
        }
        else { // If negative D, then cycle through set of positive D values        
            if(introgression_model == 1) {                
                if(isPol == 0) { // If fixed differences are not polarized            
                    gridsize_dvalue = (int)(D0/est_theta); // Set grid of divergence times D with introgressing population (from theta to D0)
                }
                else { // If fixed differences are polarized
                    gridsize_dvalue = (int)(2.0*D0/est_theta); // Set grid of divergence times D with introgressing population (from theta to 2D0)
                }                    
                    
                dvalue_grid = (double*)malloc(gridsize_dvalue*sizeof(double)); // Divergence-value grid
                for(i = 0; i < gridsize_dvalue; i++) {
                    dvalue_grid[i] = (i+1.0)*est_theta;
                }
            }
            else if(introgression_model == 2) {  
                if(isPol == 0) { // If fixed differences are not polarized                
                    gridsize_dvalue = (int)((D0 - 2.0*(nmax-1.0)*est_thetaL/nmax)/est_thetaL); // Set grid of divergence times D with introgressing population (2(n-1)theta_L/n < D < D0)
                }
                else { // If fixed differences are polarized
                    gridsize_dvalue = (int)(2.0*(D0 - (nmax-1.0)*est_thetaL/nmax)/est_thetaL); // Set grid of divergence times D with introgressing population (2(n-1)theta_L/n < D < 2D0)
                }                 
                
                dvalue_grid = (double*)malloc(gridsize_dvalue*sizeof(double)); // Divergence-value grid
                for(i = 0; i < gridsize_dvalue; i++) {
                //    dvalue_grid[i] = (2.0 + i)*(nmax-1.0)*est_thetaL/nmax; // Start at 2*(n-1)*theta_L/n
                    dvalue_grid[i] = (2.0 + i)*est_thetaL; // Start at 2*theta_L
                }
            }
            else {
                printf("ERROR: Did not set model to 1 or 2\n");
                usage();
            }
        }
    }
    else { // Read grid values from file
        readGridDValues(dvalue_gridfn); 

        for(i = 0; i < gridsize_dvalue; i++) { // Check to see if user-defined grid is compatible with minimum
            if(isPol == 0) { // If fixed differences are not polarized                 
                if(introgression_model == 1 && (dvalue_grid[i] < est_theta || dvalue_grid[i] > D0)) {
                    printf("ERROR: D cannot be smaller than %e or greater than %e for Model 1 when not polarized\n", est_theta, D0);
                    exit(-1);
                }            
                else if(introgression_model == 2 && (dvalue_grid[i] < 2.0*(nmax-1.0)*est_thetaL/nmax || dvalue_grid[i] > D0)) {
                    printf("ERROR: D cannot be smaller than %e or greater than %e for Model 2 when not polarized\n", 2.0*(nmax-1.0)*est_thetaL/nmax, D0);
                    exit(-1);
                }
            }
            else { // If fixed differences are polarized 
                if(introgression_model == 1 && (dvalue_grid[i] < est_theta || dvalue_grid[i] > 2.0*D0)) {
                    printf("ERROR: D cannot be smaller than %e or greater than %e for Model 1 when polarized\n", est_theta, 2.0*D0);
                    exit(-1);
                }            
                else if(introgression_model == 2 && (dvalue_grid[i] < 2.0*(nmax-1.0)*est_thetaL/nmax || dvalue_grid[i] > 2.0*D0)) {
                    printf("ERROR: D cannot be smaller than %e or greater than %e for Model 2 when polarized\n", 2.0*(nmax-1.0)*est_thetaL/nmax, 2.0*D0);
                    exit(-1);
                }
            }
        }
    }

    calcprobs_introgression(p_un, 500);
    writeprobs(outfn, 500);
    free(p_un);
      
  }
  else if(strcmp(argv[1], "-i") == 0 || strcmp(argv[1], "-ig") == 0  || strcmp(argv[1], "-iu") == 0) {
    if(argc!=9) {
      usage();
    }
    else {
      printf("You have chosen to get introgression sweeps using pre-computed frequency spectra\n");
    }

    if(strcmp(argv[1], "-i") == 0) { // Given gridsize G
      gridsize = atoi(argv[2]);
      if (gridsize<=0) {
        printf("gridsize should be > 0\n"); usage();
      }
    }
    else if(strcmp(argv[1], "-ig") == 0) { // Given space between grid point in g nucleotides
      space_between_grid_points = atoi(argv[2]);
      if (space_between_grid_points<=0) {
        printf("The space between grid points (g) should be > 0\n"); usage();
      } 
    }
    else if(strcmp(argv[1], "-iu") == 0) { // Given file of gridpoints GridFile
      sprintf(gridfn, "%s", argv[2]);
      readgrid(gridfn);
      flagUserDefinedGrid = 1;
      gridsize=datasize_grid;
    }
    

    sprintf(snpfn, "%s", argv[3]);
    sprintf(freqfn, "%s", argv[4]);

    isPol = atoi(argv[6]);    	
    introgression_model = atoi(argv[7]);
    if(introgression_model != 1 && introgression_model != 2) {
        printf("The MODEL should take values 1 or 2\n");
        usage();
    }

    sprintf(outfn, "%s", argv[8]);
    readsnps(snpfn);
    
    if(strcmp(argv[1], "-ig") == 0) {
	    gridsize = (int)((maxFreqPos-minFreqPos + space_between_grid_points/2.0 )/space_between_grid_points);
    	printf("(%lf,%lf,%d)\n", minFreqPos, maxFreqPos, gridsize);
    }   

    block_start = 0; // Assume 1 block, and so pick first test site
    block_stop = gridsize - 1; // Assume 1 block, and so pick last test site

    p_un = loadfreq_unnormalized(freqfn); // Loading S_j(n) for introgression sweeps
   
    InitFact();
	InitBinom();

    D0=p_j_s(1, 1, p_un); // S_1(1)
    est_theta=p_j_s(1, 2, p_un); // S_1(2)
    est_thetaL = 0.0;
    for(i = 1; i < nmax; i++) {
        est_thetaL += i * p_j_s(i, nmax, p_un); // i*S_i(n)
    } 
    est_thetaL =  est_thetaL / (nmax - 1.0);

    sprintf(dvalue_gridfn, "%s", argv[5]);
    if(sscanf(dvalue_gridfn, "%f %n", &inignore, &inlen) == 1 && !dvalue_gridfn[inlen]) { // Determine if number. If not, then read D from files
        gen_dist = atof(argv[5]);
        if(gen_dist >= 0.0) { // If positive D, then use the D inputted
            gridsize_dvalue = 1;
            dvalue_grid = (double*)malloc(sizeof(double)); // Divergence-value grid
            dvalue_grid[0] = gen_dist;

            // Check to see if user-defined D value is compatible with minimum
            if(isPol == 0) { // If fixed differences are not polarized                 
                if(introgression_model == 1 && (gen_dist < est_theta || gen_dist > D0)) {
                    printf("ERROR: D cannot be smaller than %e or greater than %e for Model 1 when not polarized\n", est_theta, D0);
                    exit(-1);
                }            
                else if(introgression_model == 2 && (gen_dist < 2.0*(nmax-1.0)*est_thetaL/nmax || gen_dist > D0)) {
                    printf("ERROR: D cannot be smaller than %e or greater than %e for Model 2 when not polarized\n", 2.0*(nmax-1.0)*est_thetaL/nmax, D0);
                    exit(-1);
                }
            }
            else { // If fixed differences are polarized 
                if(introgression_model == 1 && (gen_dist < est_theta || gen_dist > 2.0*D0)) {
                    printf("ERROR: D cannot be smaller than %e or greater than %e for Model 1 when polarized\n", est_theta, 2.0*D0);
                    exit(-1);
                }            
                else if(introgression_model == 2 && (gen_dist < 2.0*(nmax-1.0)*est_thetaL/nmax || gen_dist > 2.0*D0)) {
                    printf("ERROR: D cannot be smaller than %e or greater than %e for Model 2 when polarized\n", 2.0*(nmax-1.0)*est_thetaL/nmax, 2.0*D0);
                    exit(-1);
                }
            }
        }
        else { // If negative D, then cycle through set of positive D values        
            if(introgression_model == 1) {                
                if(isPol == 0) { // If fixed differences are not polarized            
                    gridsize_dvalue = (int)(D0/est_theta); // Set grid of divergence times D with introgressing population (from theta to D0)
                }
                else { // If fixed differences are polarized
                    gridsize_dvalue = (int)(2.0*D0/est_theta); // Set grid of divergence times D with introgressing population (from theta to 2D0)
                }                    
                    
                dvalue_grid = (double*)malloc(gridsize_dvalue*sizeof(double)); // Divergence-value grid
                for(i = 0; i < gridsize_dvalue; i++) {
                    dvalue_grid[i] = (i+1.0)*est_theta;
                }
            }
            else if(introgression_model == 2) {  
                if(isPol == 0) { // If fixed differences are not polarized                
                    gridsize_dvalue = (int)((D0 - 2.0*(nmax-1.0)*est_thetaL/nmax)/est_thetaL); // Set grid of divergence times D with introgressing population (2(n-1)theta_L/n < D < D0)
                }
                else { // If fixed differences are polarized
                    gridsize_dvalue = (int)(2.0*(D0 - (nmax-1.0)*est_thetaL/nmax)/est_thetaL); // Set grid of divergence times D with introgressing population (2(n-1)theta_L/n < D < 2D0)
                }                 
                
                dvalue_grid = (double*)malloc(gridsize_dvalue*sizeof(double)); // Divergence-value grid
                for(i = 0; i < gridsize_dvalue; i++) {
                //    dvalue_grid[i] = (2.0 + i)*(nmax-1.0)*est_thetaL/nmax; // Start at 2*(n-1)*theta_L/n
                    dvalue_grid[i] = (2.0 + i)*est_thetaL; // Start at 2*theta_L
                }
            }
            else {
                printf("ERROR: Did not set model to 1 or 2\n");
                usage();
            }
        }
    }
    else { // Read grid values from file
        readGridDValues(dvalue_gridfn); 

        for(i = 0; i < gridsize_dvalue; i++) { // Check to see if user-defined grid is compatible with minimum
            if(isPol == 0) { // If fixed differences are not polarized                 
                if(introgression_model == 1 && (dvalue_grid[i] < est_theta || dvalue_grid[i] > D0)) {
                    printf("ERROR: D cannot be smaller than %e or greater than %e for Model 1 when not polarized\n", est_theta, D0);
                    exit(-1);
                }            
                else if(introgression_model == 2 && (dvalue_grid[i] < 2.0*(nmax-1.0)*est_thetaL/nmax || dvalue_grid[i] > D0)) {
                    printf("ERROR: D cannot be smaller than %e or greater than %e for Model 2 when not polarized\n", 2.0*(nmax-1.0)*est_thetaL/nmax, D0);
                    exit(-1);
                }
            }
            else { // If fixed differences are polarized 
                if(introgression_model == 1 && (dvalue_grid[i] < est_theta || dvalue_grid[i] > 2.0*D0)) {
                    printf("ERROR: D cannot be smaller than %e or greater than %e for Model 1 when polarized\n", est_theta, 2.0*D0);
                    exit(-1);
                }            
                else if(introgression_model == 2 && (dvalue_grid[i] < 2.0*(nmax-1.0)*est_thetaL/nmax || dvalue_grid[i] > 2.0*D0)) {
                    printf("ERROR: D cannot be smaller than %e or greater than %e for Model 2 when polarized\n", 2.0*(nmax-1.0)*est_thetaL/nmax, 2.0*D0);
                    exit(-1);
                }
            }
        }
    }

    /*double h =0.0;
    for(i = 0; i < gridsize;i++) {
        h=0.0; 
        for(j=1; j<=29; j++) {
            h += 2.0*(j/30.0)*(1.0-j/30.0);      
            h += 2.0*(j/30.0)*(1.0-j/30.0)*prob[30][j][curr_d][0][i];        
        }
        printf("%d %lf %lf\n",i,ads[i],h);
    } */
    find_introgression_sweeps(outfn, p_un, gridsize, isnoisy, -1);
    free(p_un);
    free(vals_int);
    free(likes_int);
   /* int i=0;
    int j=0;
    for(i=0; i < 500; i++) {
        printf("%d \%lf", i, ads[i]);
        for(j=1; j<xmax; j++) {
            printf(" %lf",prob[50][j][curr_d][0][i]);         
        }
        printf("\n");
    }
    exit(-1);*/
  }
  else if(strcmp(argv[1], "-pi") == 0 || strcmp(argv[1], "-pig") == 0  || strcmp(argv[1], "-piu") == 0) {
    if(argc!=7) {
      usage();
    }
    else {
      printf("You have chosen to get introgression sweeps using pre-computed frequency spectra and pre-computed lookup table\n");
    }

    if(strcmp(argv[1], "-pi") == 0) { // Given gridsize G
      gridsize = atoi(argv[2]);
      if (gridsize<=0) {
        printf("gridsize should be > 0\n"); usage();
      }
    }
    else if(strcmp(argv[1], "-pig") == 0) { // Given space between grid point in g nucleotides
      space_between_grid_points = atoi(argv[2]);
      if (space_between_grid_points<=0) {
        printf("The space between grid points (g) should be > 0\n"); usage();
      } 
    }
    else if(strcmp(argv[1], "-piu") == 0) { // Given file of gridpoints GridFile
      sprintf(gridfn, "%s", argv[2]);
      readgrid(gridfn);
      flagUserDefinedGrid = 1;
      gridsize=datasize_grid;
    }
    
    sprintf(snpfn, "%s", argv[3]);
    sprintf(freqfn, "%s", argv[4]);
    sprintf(lookupfn, "%s", argv[5]);
    sprintf(outfn, "%s", argv[6]);
    readsnps(snpfn);
    
    if(strcmp(argv[1], "-pig") == 0) {
	    gridsize = (int)((maxFreqPos-minFreqPos + space_between_grid_points/2.0 )/space_between_grid_points);
    	printf("(%lf,%lf,%d)\n", minFreqPos, maxFreqPos, gridsize);
    }   

    block_start = 0; // Assume 1 block, and so pick first test site
    block_stop = gridsize - 1; // Assume 1 block, and so pick last test site

    p_un = loadfreq_unnormalized(freqfn); // Loading S_j(n) for introgression sweeps
   
    InitFact();
	InitBinom();

    D0=p_j_s(1, 1, p_un); // S_1(1)
    readprobs(lookupfn);    

    /*double h =0.0;
    for(i = 0; i < gridsize;i++) {
        h=0.0; 
        for(j=1; j<=29; j++) {
            h += 2.0*(j/30.0)*(1.0-j/30.0);      
            h += 2.0*(j/30.0)*(1.0-j/30.0)*prob[30][j][curr_d][0][i];        
        }
        printf("%d %lf %lf\n",i,ads[i],h);
    } */
    find_introgression_sweeps(outfn, p_un, gridsize, isnoisy, -1);
    free(p_un);
    free(vals_int);
    free(likes_int);
   /* int i=0;
    int j=0;
    for(i=0; i < 500; i++) {
        printf("%d \%lf", i, ads[i]);
        for(j=1; j<xmax; j++) {
            printf(" %lf",prob[50][j][curr_d][0][i]);         
        }
        printf("\n");
    }
    exit(-1);*/
  }
  //else if(strcmp(argv[1], "-bi") == 0 || strcmp(argv[1], "-big") == 0  || strcmp(argv[1], "-biu") == 0) {
  else if(strcmp(argv[1], "-bi") == 0 || strcmp(argv[1], "-big") == 0) {
    if(argc!=11) {
      usage();
    }
    else {
      block = atoi(argv[9]);
      nblock = atoi(argv[10]);
      if(block < 1 || nblock < block) {
        printf("ERROR: Block number (%d) is out of range\n", block);
        exit(-1); 
      }

      isblock = 1;
      printf("You have chosen to get introgression sweeps using pre-computed frequency spectra in block %d of %d\n", block, nblock);
    }

    if(strcmp(argv[1], "-bi") == 0) { // Given gridsize G
      gridsize = atoi(argv[2]);
      if (gridsize<=0) {
        printf("gridsize should be > 0\n"); usage();
      }
    }
    else if(strcmp(argv[1], "-big") == 0) { // Given space between grid point in g nucleotides
      space_between_grid_points = atoi(argv[2]);
      if (space_between_grid_points<=0) {
        printf("The space between grid points (g) should be > 0\n"); usage();
      } 
    }
//    else if(strcmp(argv[1], "-biu") == 0) { // Given file of gridpoints GridFile
//      sprintf(gridfn, "%s", argv[2]);
//      readgrid(gridfn);
//      flagUserDefinedGrid = 1;
//      gridsize=datasize_grid;
//    }
    

    sprintf(snpfn, "%s", argv[3]);
    sprintf(freqfn, "%s", argv[4]);

    isPol = atoi(argv[6]);    	
    introgression_model = atoi(argv[7]);
    if(introgression_model != 1 && introgression_model != 2) {
        printf("The MODEL should take values 1 or 2\n");
        usage();
    }

    sprintf(outfn, "%s_%d_%d", argv[8], block, nblock);
    readsnps(snpfn);
    
    if(strcmp(argv[1], "-bi") == 0) { // Added this to fix block issues
        space_between_grid_points = (maxFreqPos - minFreqPos) / gridsize;
    }

    if(strcmp(argv[1], "-big") == 0) {
        gridsize = 1 + (int)((maxFreqPos-minFreqPos)/space_between_grid_points);
	    //gridsize = (int)((maxFreqPos-minFreqPos + space_between_grid_points/2.0 )/space_between_grid_points);
    }   

    // Modified this so that the blocks are of unequal numbers of SNPs, but equal physical lengths
    block_interval_start = minFreqPos + (block - 1.0) * (maxFreqPos + 1.0 - minFreqPos) / nblock;
    block_interval_stop = minFreqPos + block * (maxFreqPos + 1.0 - minFreqPos) / nblock;
    for(i = 0; i < gridsize; i++) {
        if(block_interval_start <= minFreqPos + i*space_between_grid_points && minFreqPos + i*space_between_grid_points < block_interval_stop) {
            block_start = i;
            break;
        }      
    }
    for(i = gridsize - 1; i >= block_start; i--) {
        if(block_interval_start <= minFreqPos + i*space_between_grid_points && minFreqPos + i*space_between_grid_points < block_interval_stop) {
            block_stop = i;
            break;
        }      
    }
    //block_start = (block-1)*((int)(gridsize / nblock)); // First test site for block
    //block_stop = block_start + ((int)(gridsize / nblock)) - 1; // Last test site for block
    //if(block_stop >= gridsize) { // Check to ensure final block does not exceed range
      //block_stop = gridsize - 1;
    //}

    p_un = loadfreq_unnormalized(freqfn); // Loading S_j(n) for introgression sweeps

    InitFact();
	InitBinom();
  
    D0=p_j_s(1, 1, p_un); // S_1(1)
    est_theta=p_j_s(1, 2, p_un); // S_1(2)
    est_thetaL = 0.0;
    for(i = 1; i < nmax; i++) {
        est_thetaL += i * p_j_s(i, nmax, p_un); // i*S_i(n)
    } 
    est_thetaL =  est_thetaL / (nmax - 1.0);
    
    sprintf(dvalue_gridfn, "%s", argv[5]);
    if(sscanf(dvalue_gridfn, "%f %n", &inignore, &inlen) == 1 && !dvalue_gridfn[inlen]) { // Determine if number. If not, then read D from files
        gen_dist = atof(argv[5]);
        if(gen_dist >= 0.0) { // If positive D, then use the D inputted
            gridsize_dvalue = 1;
            dvalue_grid = (double*)malloc(sizeof(double)); // Divergence-value grid
            dvalue_grid[0] = gen_dist;

            // Check to see if user-defined D value is compatible with minimum
            if(isPol == 0) { // If fixed differences are not polarized                 
                if(introgression_model == 1 && (gen_dist < est_theta || gen_dist > D0)) {
                    printf("ERROR: D cannot be smaller than %e or greater than %e for Model 1 when not polarized\n", est_theta, D0);
                    exit(-1);
                }            
                else if(introgression_model == 2 && (gen_dist < 2.0*(nmax-1.0)*est_thetaL/nmax || gen_dist > D0)) {
                    printf("ERROR: D cannot be smaller than %e or greater than %e for Model 2 when not polarized\n", 2.0*(nmax-1.0)*est_thetaL/nmax, D0);
                    exit(-1);
                }
            }
            else { // If fixed differences are polarized 
                if(introgression_model == 1 && (gen_dist < est_theta || gen_dist > 2.0*D0)) {
                    printf("ERROR: D cannot be smaller than %e or greater than %e for Model 1 when polarized\n", est_theta, 2.0*D0);
                    exit(-1);
                }            
                else if(introgression_model == 2 && (gen_dist < 2.0*(nmax-1.0)*est_thetaL/nmax || gen_dist > 2.0*D0)) {
                    printf("ERROR: D cannot be smaller than %e or greater than %e for Model 2 when polarized\n", 2.0*(nmax-1.0)*est_thetaL/nmax, 2.0*D0);
                    exit(-1);
                }
            }
        }
        else { // If negative D, then cycle through set of positive D values        
            if(introgression_model == 1) {                
                if(isPol == 0) { // If fixed differences are not polarized            
                    gridsize_dvalue = (int)(D0/est_theta); // Set grid of divergence times D with introgressing population (from theta to D0)
                }
                else { // If fixed differences are polarized
                    gridsize_dvalue = (int)(2.0*D0/est_theta); // Set grid of divergence times D with introgressing population (from theta to 2D0)
                }                    
                    
                dvalue_grid = (double*)malloc(gridsize_dvalue*sizeof(double)); // Divergence-value grid
                for(i = 0; i < gridsize_dvalue; i++) {
                    dvalue_grid[i] = (i+1.0)*est_theta;
                }
            }
            else if(introgression_model == 2) {  
                if(isPol == 0) { // If fixed differences are not polarized                
                    gridsize_dvalue = (int)((D0 - 2.0*(nmax-1.0)*est_thetaL/nmax)/est_thetaL); // Set grid of divergence times D with introgressing population (2(n-1)theta_L/n < D < D0)
                }
                else { // If fixed differences are polarized
                    gridsize_dvalue = (int)(2.0*(D0 - (nmax-1.0)*est_thetaL/nmax)/est_thetaL); // Set grid of divergence times D with introgressing population (2(n-1)theta_L/n < D < 2D0)
                }                 
                
                dvalue_grid = (double*)malloc(gridsize_dvalue*sizeof(double)); // Divergence-value grid
                for(i = 0; i < gridsize_dvalue; i++) {
                //    dvalue_grid[i] = (2.0 + i)*(nmax-1.0)*est_thetaL/nmax; // Start at 2*(n-1)*theta_L/n
                    dvalue_grid[i] = (2.0 + i)*est_thetaL; // Start at 2*theta_L
                }
            }
            else {
                printf("ERROR: Did not set model to 1 or 2\n");
                usage();
            }
        }
    }
    else { // Read grid values from file
        readGridDValues(dvalue_gridfn); 

        for(i = 0; i < gridsize_dvalue; i++) { // Check to see if user-defined grid is compatible with minimum
            if(isPol == 0) { // If fixed differences are not polarized                 
                if(introgression_model == 1 && (dvalue_grid[i] < est_theta || dvalue_grid[i] > D0)) {
                    printf("ERROR: D cannot be smaller than %e or greater than %e for Model 1 when not polarized\n", est_theta, D0);
                    exit(-1);
                }            
                else if(introgression_model == 2 && (dvalue_grid[i] < 2.0*(nmax-1.0)*est_thetaL/nmax || dvalue_grid[i] > D0)) {
                    printf("ERROR: D cannot be smaller than %e or greater than %e for Model 2 when not polarized\n", 2.0*(nmax-1.0)*est_thetaL/nmax, D0);
                    exit(-1);
                }
            }
            else { // If fixed differences are polarized 
                if(introgression_model == 1 && (dvalue_grid[i] < est_theta || dvalue_grid[i] > 2.0*D0)) {
                    printf("ERROR: D cannot be smaller than %e or greater than %e for Model 1 when polarized\n", est_theta, 2.0*D0);
                    exit(-1);
                }            
                else if(introgression_model == 2 && (dvalue_grid[i] < 2.0*(nmax-1.0)*est_thetaL/nmax || dvalue_grid[i] > 2.0*D0)) {
                    printf("ERROR: D cannot be smaller than %e or greater than %e for Model 2 when polarized\n", 2.0*(nmax-1.0)*est_thetaL/nmax, 2.0*D0);
                    exit(-1);
                }
            }
        }
    }
    
    /*double h =0.0;
    for(i = 0; i < gridsize;i++) {
        h=0.0; 
        for(j=1; j<=29; j++) {
            h += 2.0*(j/30.0)*(1.0-j/30.0);      
            h += 2.0*(j/30.0)*(1.0-j/30.0)*prob[30][j][curr_d][0][i];        
        }
        printf("%d %lf %lf\n",i,ads[i],h);
    } */
    find_introgression_sweeps(outfn, p_un, gridsize, isnoisy, -1);
    free(p_un);
    free(vals_int);
    free(likes_int);
   /* int i=0;
    int j=0;
    for(i=0; i < 500; i++) {
        printf("%d \%lf", i, ads[i]);
        for(j=1; j<xmax; j++) {
            printf(" %lf",prob[50][j][curr_d][0][i]);         
        }
        printf("\n");
    }
    exit(-1);*/
  }
//  else if(strcmp(argv[1], "-pbi") == 0 || strcmp(argv[1], "-pbig") == 0  || strcmp(argv[1], "-pbiu") == 0) {
  else if(strcmp(argv[1], "-pbi") == 0 || strcmp(argv[1], "-pbig") == 0) {
    if(argc!=9) {
      usage();
    }
    else {
      block = atoi(argv[7]);
      nblock = atoi(argv[8]);
      if(block < 1 || nblock < block) {
        printf("ERROR: Block number (%d) is out of range\n", block);
        exit(-1); 
      }

      isblock = 1;
      printf("You have chosen to get introgression sweeps using pre-computed frequency spectra and pre-computed lookup table in block %d of %d\n", block, nblock);
    }

    if(strcmp(argv[1], "-pbi") == 0) { // Given gridsize G
      gridsize = atoi(argv[2]);
      if (gridsize<=0) {
        printf("gridsize should be > 0\n"); usage();
      }
    }
    else if(strcmp(argv[1], "-pbig") == 0) { // Given space between grid point in g nucleotides
      space_between_grid_points = atoi(argv[2]);
      if (space_between_grid_points<=0) {
        printf("The space between grid points (g) should be > 0\n"); usage();
      } 
    }
//    else if(strcmp(argv[1], "-pbiu") == 0) { // Given file of gridpoints GridFile
//      sprintf(gridfn, "%s", argv[2]);
//      readgrid(gridfn);
//      flagUserDefinedGrid = 1;
//      gridsize=datasize_grid;
//    }
    

    sprintf(snpfn, "%s", argv[3]);
    sprintf(freqfn, "%s", argv[4]);
    sprintf(lookupfn, "%s", argv[5]);
    sprintf(outfn, "%s_%d_%d", argv[6], block, nblock);
    readsnps(snpfn);
    
    if(strcmp(argv[1], "-pbi") == 0) { // Added this to fix block issues
        space_between_grid_points = (maxFreqPos - minFreqPos) / gridsize;
    }

    if(strcmp(argv[1], "-pbig") == 0) {
        gridsize = 1 + (int)((maxFreqPos-minFreqPos)/space_between_grid_points);
	    //gridsize = (int)((maxFreqPos-minFreqPos + space_between_grid_points/2.0 )/space_between_grid_points);
    }   

    // Modified this so that the blocks are of unequal numbers of SNPs, but equal physical lengths
    block_interval_start = minFreqPos + (block - 1.0) * (maxFreqPos + 1.0 - minFreqPos) / nblock;
    block_interval_stop = minFreqPos + block * (maxFreqPos + 1.0 - minFreqPos) / nblock;
    for(i = 0; i < gridsize; i++) {
        if(block_interval_start <= minFreqPos + i*space_between_grid_points && minFreqPos + i*space_between_grid_points < block_interval_stop) {
            block_start = i;
            break;
        }      
    }
    for(i = gridsize - 1; i >= block_start; i--) {
        if(block_interval_start <= minFreqPos + i*space_between_grid_points && minFreqPos + i*space_between_grid_points < block_interval_stop) {
            block_stop = i;
            break;
        }      
    }
    //block_start = (block-1)*((int)(gridsize / nblock)); // First test site for block
    //block_stop = block_start + ((int)(gridsize / nblock)) - 1; // Last test site for block
    //if(block_stop >= gridsize) { // Check to ensure final block does not exceed range
      //block_stop = gridsize - 1;
    //}
    
    p_un = loadfreq_unnormalized(freqfn); // Loading S_j(n) for introgression sweeps

    InitFact();
	InitBinom();
  
    D0=p_j_s(1, 1, p_un); // S_1(1)
    readprobs(lookupfn);
    
    /*double h =0.0;
    for(i = 0; i < gridsize;i++) {
        h=0.0; 
        for(j=1; j<=29; j++) {
            h += 2.0*(j/30.0)*(1.0-j/30.0);      
            h += 2.0*(j/30.0)*(1.0-j/30.0)*prob[30][j][curr_d][0][i];        
        }
        printf("%d %lf %lf\n",i,ads[i],h);
    } */
    find_introgression_sweeps(outfn, p_un, gridsize, isnoisy, -1);
    free(p_un);
    free(vals_int);
    free(likes_int);
   /* int i=0;
    int j=0;
    for(i=0; i < 500; i++) {
        printf("%d \%lf", i, ads[i]);
        for(j=1; j<xmax; j++) {
            printf(" %lf",prob[50][j][curr_d][0][i]);         
        }
        printf("\n");
    }
    exit(-1);*/
  }
  else if(strcmp(argv[1], "-m") == 0) {
    if(argc!=4) {
      usage();
    }
    else {
      nblock = atoi(argv[3]);
      printf("You have chosen to merge introgression sweeps across %d blocks\n", nblock);
    }

    sprintf(outfn, "%s", argv[2]);
    
    outfile=my_fopen(outfn, "w");
    fprintf(outfile, "location\tLR\talpha\tD\n");
    for(block = 1; block <= nblock; block++) {
         writeblocks(outfile, outfn, block, nblock);
    }
    fclose(outfile);
  }
  else usage();
  return 0;
}
