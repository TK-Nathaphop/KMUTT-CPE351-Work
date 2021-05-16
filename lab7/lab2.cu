#include <stdio.h> 
#include <stdlib.h> 

__global__ void calPi(double * pi_D)
	{
	int t_rank = (blockIdx.x*blockDim.x) + threadIdx.x ;
	int x = 2 + (2 * t_rank);
	double y = ((double)4/x) * (double)1/(x+1) * (double)1/(x+2);

	if(t_rank % 2 != 0)
		y = -y;

	pi_D[t_rank] = y;
	}

int main()
	{
	printf("pi calculate...\n");

	int thread_size = 500, block_size = 5;
	double *pi_D;
	double *pi_H;

	pi_H = (double*) malloc(sizeof(double) *thread_size*block_size);
	cudaMalloc( (void **)&pi_D, sizeof(double)*thread_size*block_size);

	calPi<<<block_size,thread_size>>>(pi_D);

	cudaMemcpy(pi_H, pi_D, thread_size*block_size * sizeof(double), cudaMemcpyDeviceToHost);
	cudaFree(pi_D);

	double pi = 3;
	for(int i = 0 ; i < thread_size*block_size ; i++)
		pi = (double) pi + pi_H[i];

	/* Change double to string, prevent it round the decimal */
	char result[12];
	sprintf(result, "%.11lf", pi);
	result[strlen(result)-1] = '\0';
	printf("calculated pi = %s \n",result);
	}