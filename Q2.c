//Question 2: Fourier transform of f(x)=sinc(x), using FFT.
#include <stdio.h>    
#include <complex.h>
#include <math.h>

double _Complex DFT(double F[],int n, float Q)  //FFT algorithm.
//Performs FFT recursively without the sqrt(N) normalisation!
{
	if(n!=1)
	{
		//Divide the array in 2 parts.
		double ev[n/2]; //To store even indexed elements
		double od[n/2]; //To store odd  indexed elements
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

int main() //Main function.
{
	int N=128;
	double f[N];
	float x_min=-10;
	float x_max=10;
	float Dx=(x_max-x_min)/(N-1);
	for(int i=0;i<N;i++)//Define the sinc function
	{
		if(x_min+i*Dx != 0)
		{
			f[i]=sin(x_min+i*Dx)/(x_min+i*Dx);
		}
		else
		{
			f[i]=1;
		}	
	}
	
	double _Complex Fourier[N];
	double k[N];
	float q;
	FILE *data;
	data = fopen("Q_2.txt","w");
	for(int i=0;i<N;i++)
	{
		q=-N/2+i;
		k[i]=2*M_PI*q/(N*Dx);	//Following the convention used by python for k values.
		Fourier[i]=Dx*sqrt(1/(2*M_PI))*cexp(-I*k[i]*x_min)*DFT(f,N,q);
		fprintf(data,"%f\t%f\n",k[i],creal(Fourier[i])); /*Send data to a file to be
		used later in python code.*/
	}
	fclose(data);	
}
//To compile: gcc -Wall Q2.c -lm -lgsl -lgslcblas -o Q2
//To run	: ./Q2
