void getfreq(char *outfn);
void getfreq_unnormalized(char *outfn);
double *loadfreq(char *infn);
double *loadfreq_unnormalized(char *infn);
double likelihood_freq_onesnp(int x, int n, int folded, int count, double *p);
