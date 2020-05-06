#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

int main (void)
{
	int N=128;
	int i;
	double data[2*N];
	float x_min=-10;
	float x_max=10;
	float Dx=(x_max-x_min)/(N-1);
	float k[N];
	float f_k[N];

	//Define the sinc function.
	for (i = 0; i < N; i++)
	{
		IMAG(data,i) =0;
		if ((x_min+Dx*i)!=0)
		  	REAL(data,i) = sin(x_min+Dx*i)/(x_min+Dx*i);
		else
			REAL(data,i) =1;
	}

	for (i = 0; i < N; i++)
	{
	   k[i]=-M_PI/Dx+i/(N*Dx)*2*M_PI;
	}

	gsl_fft_complex_radix2_forward (data, 1, N);

	FILE *Fil;
	Fil = fopen("Q_3.txt","w");

	/*The next two for loops are to get the transformed 
	values corresponding to increasing order of k values.*/
	for (i = 0; i <= N/2; i++)				
	{
		f_k[i+N/2-1]=Dx*fabs(REAL(data,i))/sqrt(2*M_PI);
	}
	for (i = N/2+1; i < N; i++)
	{
	    f_k[i-N/2-1]=Dx*fabs(REAL(data,i))/sqrt(2*M_PI);
	}
	
	for (i = 0; i < N; i++)
	{
		fprintf (Fil, "%e  %e\n",k[i], f_k[i]);
	}

	return 0;
}
//To compile: gcc -Wall Q3.c -lm -lgsl -lgslcblas -o Q3
