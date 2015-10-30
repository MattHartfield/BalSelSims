/* BalSelSims.c 

Simulation of balancing selection simulation with sampling (finite population effects)
Selection; reproduction; mutation based on deterministic recursions, then sampling
Repeats for 'Length' generations

Simulation uses routines found with the GNU Scientific Library (GSL)
(http://www.gnu.org/software/gsl/)
Since GSL is distributed under the GNU General Public License 
(http://www.gnu.org/copyleft/gpl.html), you must download it 
separately from this file.

This program can be compiled in e.g. GCC using a command like:
gcc BalSelSims -lm -lgsl -lgslcblas -I/usr/local/include -L/usr/local/lib BalSelSims.c

Then run by executing:
./BalSelSims N s rec sex self gc reps
Where:
- N is the population size
- s is the fitness disadvantage of homozygotes
- rec is recombination rate
- sex is rate of sex (a value between 0 = obligate asex, and 1 = obligate sex)
- self is selfing rate
- gc is gene conversion
- reps is how many times to introduce linked neutral allele

Note that haplotypes are defined as:
x1 = ab
x2 = Ab
x3 = aB
x4 = AB

Genotypes defined as:
g11 = g1 = ab/ab
g12 = g2 = Ab/ab
g13 = g3 = aB/ab
g14 = g4 = AB/ab
g22 = g5 = Ab/Ab
g23 = g6 = Ab/aB
g24 = g7 = Ab/AB
g33 = g8 = aB/aB
g34 = g9 = aB/AB
g44 = g10 = AB/AB

*/

/* Preprocessor statements */
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stddef.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/* Function prototypes */
void geninit(double *geninit);
void selection(double *geninit);
void reproduction(double *geninit);
void gconv(double *geninit);
void neutinit(double *geninit,const gsl_rng *r);
double ncheck(double *geninit);
double pcheck(double *geninit);

/* Global variable declaration */
unsigned int N = 0;		/* Pop size */
double s = 0;			/* Fitness disadvantage of homozyogtes */
double rec = 0;			/* Recombination rate */
double sex = 0;			/* Rate of sexual reproduction */
double self = 0;		/* Rate of self-fertilisation */
double gc = 0;			/* Rate of gene conversion */

/* Main program */
int main(int argc, char *argv[]){
	unsigned int g, i; 				/* Counters. Reps counter, geno counter */
	unsigned int reps;				/* Length of simulation (no. of introductions of neutral site) */
	double Bcheck = 0;				/* Frequency of B after each reproduction */
	double Acheck = 0;				/* Frequency of polymorphism */
	double Hsum = 0;				/* Summed heterozygosity over transit time of neutral allele */
	
	/* GSL random number definitions */
	const gsl_rng_type * T; 
	gsl_rng * r;
	
	/* This reads in data from command line. */
	if(argc != 8){
		fprintf(stderr,"Invalid number of input values.\n");
		exit(1);
	}
	N = strtod(argv[1],NULL);
	s = strtod(argv[2],NULL);
	rec = strtod(argv[3],NULL);
	sex = strtod(argv[4],NULL);
	self = strtod(argv[5],NULL);
	gc = strtod(argv[6],NULL);
	reps = strtod(argv[7],NULL);
	
	/* Arrays definition and memory assignment */
	double *genotype = calloc(10,sizeof(double));				/* Genotype frequencies */
	unsigned int *gensamp = calloc(10,sizeof(unsigned int));	/* New population samples */
	  
	/* create a generator chosen by the 
    environment variable GSL_RNG_TYPE */
     
	gsl_rng_env_setup();
	if (!getenv("GSL_RNG_SEED")) gsl_rng_default_seed = time(0);
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	
	/* Initialising genotypes */
	geninit(genotype);
    
    /* Run simulation for 2000 generations to create a burn in */
    for(g = 0; g < 2000; g++){

    	/* Selection routine */
    	selection(genotype);
    	
    	/* Reproduction routine */
    	reproduction(genotype);
    	
       	/* Gene conversion routine */
       	gconv(genotype);
       	
       	/* Sampling based on new frequencies */
       	gsl_ran_multinomial(r,10,N,genotype,gensamp);
       	for(i = 0; i < 10; i++){
			*(genotype + i) = (*(gensamp + i))/(1.0*N);
       	}
       	
       	/* Printing out results (for testing) */
		/*
	   	for(i = 0; i < 10; i++){
			printf("%.10lf ", *(genotype + i));
		}
		printf("\n");
		*/
    }
    
	/* Reintroducing neutral genotype, resetting hap sum */	
	neutinit(genotype,r);
    Bcheck = ncheck(genotype);
    Hsum = Bcheck*(1-Bcheck);
    /* printf("%.10lf %.10lf\n",Bcheck,Hsum); */
    
    /* Introduce and track neutral mutations 'reps' times */
    g = 0;
    while(g < reps){

    	/* Selection routine */
    	selection(genotype);
    	
    	/* Reproduction routine */
    	reproduction(genotype);
    	
    	/* Gene conversion routine */
       	gconv(genotype);
       	
    	/* Sampling based on new frequencies */
       	gsl_ran_multinomial(r,10,N,genotype,gensamp);
       	for(i = 0; i < 10; i++){
       		*(genotype + i) = (*(gensamp + i))/(1.0*N);
       	}
       	
       	/* Checking state of haplotypes: if B fixed reset so can start fresh next time */
		Bcheck = ncheck(genotype);
		Hsum += Bcheck*(1-Bcheck);
		/* printf("%.10lf %.10lf\n",Bcheck,Hsum); */
		
		/* If polymorphism fixed then abandon simulation */
       	Acheck = pcheck(genotype);
       	if(Acheck == 0){
       		g = reps;
       	}
       	
       	if(Bcheck == 0 || Bcheck == 1){
       		printf("%.10lf\n",Hsum);
       		g++;
       		/* printf("Rep Number %d\n",g); */
       		
       		if(Bcheck == 1){
       			/* Reset genotypes so B becomes ancestral allele */
       			*(genotype + 0) = *(genotype + 7);
       			*(genotype + 1) = *(genotype + 8);
       			*(genotype + 4) = *(genotype + 9);
       			*(genotype + 7) = 0;
       			*(genotype + 8) = 0;
       			*(genotype + 9) = 0;      			
       		}
     		 
     		/* Reintroducing neutral genotype, resetting hap sum */     		
	    	neutinit(genotype,r);
			Bcheck = ncheck(genotype);
       		Hsum = Bcheck*(1-Bcheck);
       	}
    
	}	/* End of simulation */
	
	
	/* Freeing memory and wrapping up */
 	gsl_rng_free(r);
 	free(gensamp);
	free(genotype);
	/* printf("The End!\n"); */
	return 0;
}

/* Initialising genotypes */
void geninit(double *geninit){
		
	/* Basic idea: before neutral mutation introduced, g0 = 1/4; g1 = 1/2; g4 = 1/4. Deterministic frequencies. */

	*(geninit + 0) = 0.25;
	*(geninit + 1) = 0.50;
	*(geninit + 2) = 0;
	*(geninit + 3) = 0;
	*(geninit + 4) = 0.25;
	*(geninit + 5) = 0;
	*(geninit + 6) = 0;
	*(geninit + 7) = 0;
	*(geninit + 8) = 0;
	*(geninit + 9) = 0;

}	/* End of gen initiation routine */

/* Initialising NEUTRAL allele */
void neutinit(double *geninit, const gsl_rng *r){

	double *probin = calloc(4,sizeof(double));					/* Probability inputs, determine location of neutral allele */
	unsigned int *probout = calloc(4,sizeof(unsigned int));		/* Output from multinomial sampling */
	
	/* Basic idea: g11 at freq p^2; g12 at freq 2pq; g22 at freq q^2.
	So weighted sampling to determine which genotype the neutral allele arises on
	*/
	
	/* Prob definitions */
	*(probin + 0) = *(geninit + 0);
	*(probin + 1) = (*(geninit + 1))/2.0;
	*(probin + 2) = (*(geninit + 1))/2.0;
	*(probin + 3) = *(geninit + 4);
	
	gsl_ran_multinomial(r,4,1,probin,probout);
	
	/* Redefining genotypes depending on outcome */
	if(*(probout + 0) == 1){
		*(geninit + 0) = (*(geninit + 0)) - 1/(1.0*N);
		*(geninit + 2) = (*(geninit + 2)) + 1/(1.0*N);
	}
	else if(*(probout + 1) == 1){
		*(geninit + 1) = (*(geninit + 1)) - 1/(1.0*N);
		*(geninit + 3) = (*(geninit + 3)) + 1/(1.0*N);	
	}
	else if(*(probout + 2) == 1){
		*(geninit + 1) = (*(geninit + 1)) - 1/(1.0*N);
		*(geninit + 5) = (*(geninit + 5)) + 1/(1.0*N);	
	}
	else if(*(probout + 3) == 1){
		*(geninit + 4) = (*(geninit + 4)) - 1/(1.0*N);
		*(geninit + 6) = (*(geninit + 6)) + 1/(1.0*N);	
	}
	
 	free(probout);
	free(probin);
	
}	/* End of gen initiation routine */

/* Selection routine */
void selection(double *geninit){
	double Waa, WAa, WAA;		/* Fitness of locus A (bal sel locus) */
	double Wmean;				/* Mean fitness */
	
	Waa = 1-s;
	WAa = 1;
	WAA = 1-s;
	
	/* Mean fitness calculation */
	Wmean = ((*(geninit + 0))*Waa) + ((*(geninit + 1))*WAa) + ((*(geninit + 2))*Waa) + ((*(geninit + 3))*WAa) + ((*(geninit + 4))*WAA) + ((*(geninit + 5))*WAa) + ((*(geninit + 6))*WAA) + ((*(geninit + 7))*Waa) + ((*(geninit + 8))*WAa) + ((*(geninit + 9))*WAA);
	
	/* Changing frequencies by selection */
	*(geninit + 0) = ((*(geninit + 0))*Waa)/Wmean;
	*(geninit + 1) = ((*(geninit + 1))*WAa)/Wmean;
	*(geninit + 2) = ((*(geninit + 2))*Waa)/Wmean;
	*(geninit + 3) = ((*(geninit + 3))*WAa)/Wmean;
	*(geninit + 4) = ((*(geninit + 4))*WAA)/Wmean;
	*(geninit + 5) = ((*(geninit + 5))*WAa)/Wmean;
	*(geninit + 6) = ((*(geninit + 6))*WAA)/Wmean;
	*(geninit + 7) = ((*(geninit + 7))*Waa)/Wmean;
	*(geninit + 8) = ((*(geninit + 8))*WAa)/Wmean;
	*(geninit + 9) = ((*(geninit + 9))*WAA)/Wmean;
	
}	/* End of selection routine */

/* Reproduction routine */
void reproduction(double *geninit){
	/* Fed-in genotype frequencies (for ease of programming) */
	double g11s, g12s, g13s, g14s, g22s, g23s, g24s, g33s, g34s, g44s;
	/* Genotype frequencies after sex (outcross and selfing) */	
	double g11SX, g12SX, g13SX, g14SX, g22SX, g23SX, g24SX, g33SX, g34SX, g44SX;
	/* Genotype frequencies after ASEX */	
	double g11AS, g12AS, g13AS, g14AS, g22AS, g23AS, g24AS, g33AS, g34AS, g44AS;	
	/* Haplotypes */
	double x1, x2, x3, x4;
	
	/* Initial definition of genotypes */
	g11s = *(geninit + 0);
	g12s = *(geninit + 1);
	g13s = *(geninit + 2);
	g14s = *(geninit + 3);
	g22s = *(geninit + 4);
	g23s = *(geninit + 5);
	g24s = *(geninit + 6);
	g33s = *(geninit + 7);
	g34s = *(geninit + 8);
	g44s = *(geninit + 9);
	
	/* Baseline change in haplotype frequencies */
	x1 = g11s + (g12s + g13s + g14s)/2.0 - ((g14s - g23s)*rec)/2.0;
	x2 = g22s + (g12s + g23s + g24s)/2.0 + ((g14s - g23s)*rec)/2.0;
	x3 = g33s + (g13s + g23s + g34s)/2.0 + ((g14s - g23s)*rec)/2.0;
	x4 = g44s + (g14s + g24s + g34s)/2.0 - ((g14s - g23s)*rec)/2.0;
	
	/* Change in SEXUAL frequencies (both outcrossing and selfing) */
	g11SX = (g11s + (g12s + g13s + g14s*pow((1 - rec),2) + g23s*pow(rec,2))/4.0)*self*sex + (1 - self)*pow(x1,2)*sex;
	g22SX = (g22s + (g12s + g24s + g23s*pow((1 - rec),2) + g14s*pow(rec,2))/4.0)*self*sex + (1 - self)*pow(x2,2)*sex;
	g33SX = (g33s + (g13s + g34s + g23s*pow((1 - rec),2) + g14s*pow(rec,2))/4.0)*self*sex + (1 - self)*pow(x3,2)*sex;
	g44SX = (g44s + (g24s + g34s + g14s*pow((1 - rec),2) + g23s*pow(rec,2))/4.0)*self*sex + (1 - self)*pow(x4,2)*sex;
	g12SX = ((g12s + (g14s + g23s)*(1 - rec)*rec)*self*sex)/2.0 + 2.0*(1 - self)*x1*x2*sex;
	g13SX = ((g13s + (g14s + g23s)*(1 - rec)*rec)*self*sex)/2.0 + 2.0*(1 - self)*x1*x3*sex;
	g14SX = ((g14s*pow((1 - rec),2) + g23s*pow(rec,2))*self*sex)/2.0 + 2.0*(1 - self)*x1*x4*sex;
	g23SX = ((g23s*pow((1 - rec),2) + g14s*pow(rec,2))*self*sex)/2.0 + 2.0*(1 - self)*x2*x3*sex;
	g24SX = ((g24s + (g14s + g23s)*(1 - rec)*rec)*self*sex)/2.0 + 2.0*(1 - self)*x2*x4*sex;
	g34SX = ((g34s + (g14s + g23s)*(1 - rec)*rec)*self*sex)/2.0 + 2.0*(1 - self)*x3*x4*sex;
	
	/* Change in ASEXUAL frequencies */
	g11AS = g11s*(1 - sex);
	g12AS = g12s*(1 - sex);
	g13AS = g13s*(1 - sex);
	g14AS = g14s*(1 - sex);
	g22AS = g22s*(1 - sex);
	g23AS = g23s*(1 - sex);
	g24AS = g24s*(1 - sex);	
	g33AS = g33s*(1 - sex);	
	g34AS = g34s*(1 - sex);	
	g44AS = g44s*(1 - sex);
	
	/* Combining to give overall frequency change following reproduction */
	*(geninit + 0) = g11AS + g11SX;
	*(geninit + 1) = g12AS + g12SX;
	*(geninit + 2) = g13AS + g13SX;
	*(geninit + 3) = g14AS + g14SX;
	*(geninit + 4) = g22AS + g22SX;
	*(geninit + 5) = g23AS + g23SX;
	*(geninit + 6) = g24AS + g24SX;
	*(geninit + 7) = g33AS + g33SX;
	*(geninit + 8) = g34AS + g34SX;
	*(geninit + 9) = g44AS + g44SX;
	
}	/* End of reproduction routine */

/* Gene conversion routine */
void gconv(double *geninit){
	
	/* Fed-in genotype frequencies (for ease of programming) */
	double g11r, g12r, g13r, g14r, g22r, g23r, g24r, g33r, g34r, g44r;
	/* Frequencies after gene conversion */
	double g11gc, g12gc, g13gc, g14gc, g22gc, g23gc, g24gc, g33gc, g34gc, g44gc;
	
	/* Initial definition of genotypes */
	g11r = *(geninit + 0);
	g12r = *(geninit + 1);
	g13r = *(geninit + 2);
	g14r = *(geninit + 3);
	g22r = *(geninit + 4);
	g23r = *(geninit + 5);
	g24r = *(geninit + 6);
	g33r = *(geninit + 7);
	g34r = *(geninit + 8);
	g44r = *(geninit + 9);
	
	/* Gene conversion equations */
	g11gc = g11r + (gc*g12r)/4.0 + (gc*g13r)/4.0;
	g12gc = g12r*(1 - gc/2.0) + ((g14r + g23r)*gc)/4.0;
	g13gc = g13r*(1 - gc/2.0) + ((g14r + g23r)*gc)/4.0;
	g14gc = (1-gc)*g14r;
	g22gc = g22r + (g12r*gc)/4.0 + (g24r*gc)/4.0;
	g23gc = (1-gc)*g23r;
	g24gc = g24r*(1 - gc/2.0) + ((g14r + g23r)*gc)/4.0;
	g33gc = g33r + (g13r*gc)/4.0 + (g34r*gc)/4.0;
	g34gc = g34r*(1 - gc/2.0) + ((g14r + g23r)*gc)/4.0;
	g44gc = g44r + (g24r*gc)/4.0 + (g34r*gc)/4.0;
	
	/* Output */
	*(geninit + 0) = g11gc;
	*(geninit + 1) = g12gc;
	*(geninit + 2) = g13gc;
	*(geninit + 3) = g14gc;
	*(geninit + 4) = g22gc;
	*(geninit + 5) = g23gc;
	*(geninit + 6) = g24gc;
	*(geninit + 7) = g33gc;
	*(geninit + 8) = g34gc;
	*(geninit + 9) = g44gc;
}		/* End of gene conversion routine */

/* Has neutral allele fixed or not?  Measuring freq of B */
double ncheck(double *geninit){
	/* Fed-in genotype frequencies (for ease of programming) */
	double g11s, g12s, g13s, g14s, g22s, g23s, g24s, g33s, g34s, g44s;
	/* Haplotypes ONLY CONTAINING B */
	double x3, x4;
	double Btot = 0;        /* Total frequency of B */
	
	/* Initial definition of genotypes */
	g11s = *(geninit + 0);
	g12s = *(geninit + 1);
	g13s = *(geninit + 2);
	g14s = *(geninit + 3);
	g22s = *(geninit + 4);
	g23s = *(geninit + 5);
	g24s = *(geninit + 6);
	g33s = *(geninit + 7);
	g34s = *(geninit + 8);
	g44s = *(geninit + 9);
	
	/* Calculation of haplotypes containing B */
	
	x3 = g33s + (g13s + g23s + g34s)/2.0;
	x4 = g44s + (g14s + g24s + g34s)/2.0;
	
	/* Checking */
	Btot = x3 + x4;
	return Btot;
	
}	/* End of B check routine */

/* Checking if balancing polymorphism lost or not */
double pcheck(double *geninit){
	/* Fed-in genotype frequencies (for ease of programming) */
	double g11s, g12s, g13s, g14s, g22s, g23s, g24s, g33s, g34s, g44s;
	double x2, x4;
	double Atot = 0;        /* Total frequency of A */
	
	/* Initial definition of genotypes */
	g11s = *(geninit + 0);
	g12s = *(geninit + 1);
	g13s = *(geninit + 2);
	g14s = *(geninit + 3);
	g22s = *(geninit + 4);
	g23s = *(geninit + 5);
	g24s = *(geninit + 6);
	g33s = *(geninit + 7);
	g34s = *(geninit + 8);
	g44s = *(geninit + 9);
	
	/* Calculation of haplotypes containing B */
	
	x2 = g22s + (g12s + g23s + g24s)/2.0;
	x4 = g44s + (g14s + g24s + g34s)/2.0;
	
	/* Checking */
	Atot = x2 + x4;
	return Atot;
	
}	/* End of A check routine */

/* End of program */
