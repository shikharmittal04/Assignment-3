//Question 4 (First part): Fourier transform of Gaussian using FFT algorithm.
//The plotting part is in a python file Q4.py
#include <stdio.h>    
#include <complex.h>
#include <math.h>

double _Complex DFT(double F[],int n, float Q)	//FFT algorithm
{
	if(n!=1)
	{
		double ev[n/2];
		double od[n/2];
		for(int j=0;j<n/2;j++)
		{
			ev[j]=F[2*j];
		}
		for(int j=0;j<n/2;j++)
		{
			od[j]=F[2*j+1];
		}
		return DFT(ev,n/2,Q)+cexp(-2*M_PI*I*Q/n)* DFT(od,n/2,Q);
	}
	else
	{
		return F[0];
	}
}

int main() 
{
	int N=64;
	double f[N];
	float x_min=-5;
	float x_max=5;
	float Dx=(x_max-x_min)/(N-1);
	for(int i=0;i<N;i++)
	{
		f[i]=exp(-pow(x_min+i*Dx,2));	//Function whose Fourier transform is required.
	}
	
	double _Complex Fourier[N];
	double k[N];
	float q;
	FILE *data;
	data = fopen("Q_4.txt","w");
	if (!data)
	perror("fopen");
	
	for(int i=0;i<N;i++)	//Print data to a text file.
	{
		q=-N/2+i;			//Following the convention used by python for k values.
		k[i]=2*M_PI*q/(N*Dx);
		Fourier[i]=Dx*sqrt(1/(2*M_PI))*cexp(-I*k[i]*x_min)*DFT(f,N,q);
		fprintf(data,"%f   %f\n",k[i],creal(Fourier[i]));
	}
	fclose(data);
}
//To compile: gcc -Wall Q4.c -lm -lgsl -lgslcblas -o Q4
//To run	: ./Q4
