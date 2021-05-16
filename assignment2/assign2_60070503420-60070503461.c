#include <omp.h> 
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#pragma GCC optimize("O2")
// A utility function to swap two elements 
static inline void swap(double* a, double* b) 
{ 
	//printf("swap\n");
    double temp = *a; 
    *a = *b; 
    *b = temp; 
} 

static inline void quickSort(double *arr, int elements)
{

	#define  MAX_LEVELS  1000

	double  piv; 
  	int beg[MAX_LEVELS], end[MAX_LEVELS], i=0, L, R, swap ;

  	beg[0]=0; end[0]=elements;
  	while (i>=0)
  	{
    	L=beg[i]; R=end[i]-1;
    	if (L<R)
    	{
      		piv=arr[L];
      		while (L<R)
      		{
        		while (arr[R]>=piv && L<R)
        			R--;
        		if (L<R)
        			arr[L++]=arr[R];
        		while (arr[L]<=piv && L<R)
        			L++;
        		if (L<R)
        			arr[R--]=arr[L];
        	}
	      	arr[L]=piv;
	      	beg[i+1]=L+1;
	      	end[i+1]=end[i];
	      	end[i++]=L;
	      	if (end[i]-beg[i]>end[i-1]-beg[i-1])
	      	{
	        	swap=beg[i];
	        	beg[i]=beg[i-1];
	        	beg[i-1]=swap;
	        	swap=end[i];
	        	end[i]=end[i-1];
	        	end[i-1]=swap;
	    	}
		}
    	else
    	{
 	     	i--;
  		}
	}
}

static inline int onceQsort( double *data ,int start ,int end)
{
	//printf("Sorting %d - %d\n",start,end);
	if(start < end)
	{
		//swap(&data[(end - start)/2 ],&data[end]);
		int sep_point = start - 1; 
		int pivot_point = (end+start)/2;
		// printf("Start:%7d End:%7d PivotPoint%d\n",start,end,pivot_point);
		swap(&data[pivot_point],&data[end]);
		register double pivot_value = data[end];
		for(register int i = start ; i < end  ; i++)
		{
			if(data[i] < pivot_value)
			{
				sep_point++;
				swap(&data[sep_point],&data[i]);	
			}
				
		}
		sep_point++;
		swap(&data[sep_point],&data[end]);
		return sep_point;
	}
	return -1;
}

void printAll(double data[],int start ,int end)
{
	for(int i = start; i < end; i++)
		printf("%.01lf ",data[i]);
	printf("\n");
}


void initCoreQsort(double *data, int lastPos, int posSep[], int rankSize)
{
	int rank = 0;
	if(rankSize == 1)
	{
		posSep[0] = 0;
		posSep[1] = lastPos;
	}
	else if(rankSize == 2)
	{
		posSep[0] = 0;
		posSep[1] = onceQsort(data, 0, lastPos);	
		posSep[2] = lastPos;
	}
	else if(rankSize == 4)
	{
		posSep[2] = onceQsort(data, 0, lastPos);
		#pragma omp parallel shared(data) private(rank)
		{
			rank = omp_get_thread_num();
			/* Second & Third separation */
			if(rank == 0)
				posSep[1] = onceQsort(data, 0, posSep[2]-1);
			else if(rank==1)
				posSep[3] = onceQsort(data, posSep[2]+1, lastPos);
			#pragma omp barrier
		}
		/* Set first and last pos */
		posSep[4] = lastPos;
		posSep[0] = 0;
	}
	else if(rankSize == 8)
	{
		/* First separation */
		posSep[4] = onceQsort(data, 0, lastPos);

		#pragma omp parallel shared(data) private(rank)
		{
			rank = omp_get_thread_num();
			/* Second & Third separation */
			if(rank == 0)
				posSep[2] = onceQsort(data, 0, posSep[4]-1);
			else if(rank==1)
				posSep[6] = onceQsort(data, posSep[4]+1, lastPos);
			#pragma omp barrier

			if(rank == 0)
				posSep[1] = onceQsort(data, 0, posSep[2]-1);
			else if(rank==1)
				posSep[3] = onceQsort(data, posSep[2]+1, posSep[4]-1);
			else if(rank==2)
				posSep[5] = onceQsort(data, posSep[4]+1, posSep[6]-1);
			else if(rank==3)
				posSep[7] = onceQsort(data, posSep[6]+1, lastPos);
			#pragma omp barrier
		}
		/* Set first and last pos */
		posSep[8] = lastPos;
		posSep[0] = 0;
	}
	else
	{
		/* First separation */
		posSep[8] = onceQsort(data, 0, lastPos);
		#pragma omp parallel shared(data) private(rank)
		{
			rank = omp_get_thread_num();
			/* Second & Third separation */
			if(rank == 0)
				posSep[4] = onceQsort(data, 0, posSep[8]-1);
			else if(rank==1)
				posSep[12] = onceQsort(data, posSep[8]+1, lastPos);
			#pragma omp barrier

			if(rank == 0)
				posSep[2] = onceQsort(data, 0, posSep[4]-1);
			else if(rank==1)
				posSep[6] = onceQsort(data, posSep[4]+1, posSep[8]-1);
			else if(rank==2)
				posSep[10] = onceQsort(data, posSep[8]+1, posSep[12]-1);
			else if(rank==3)
				posSep[14] = onceQsort(data, posSep[12]+1, lastPos);
			#pragma omp barrier

			if(rank == 0)
			{
				posSep[1] = onceQsort(data, 0, posSep[2]-1);
				posSep[3] = onceQsort(data, posSep[2]+1, posSep[4]-1);
			}
			else if(rank==1)
			{
				posSep[5] = onceQsort(data, posSep[4]+1, posSep[6]-1);
				posSep[7] = onceQsort(data, posSep[6]+1, posSep[8]-1);
			}
			else if(rank==2)
			{
				posSep[9] = onceQsort(data, posSep[8]+1, posSep[10]-1);
				posSep[11] = onceQsort(data, posSep[10]+1, posSep[12]-1);
			}
			else if(rank==3)
			{
				posSep[13] = onceQsort(data, posSep[12]+1, posSep[14]-1);
				posSep[15] = onceQsort(data, posSep[14]+1, lastPos);
			}
			#pragma omp barrier
		}
		/* Set first and last pos */
		posSep[16] = lastPos;
		posSep[0] = 0;
	}
	
}

void initThreadQsort(double *data, int lastPos, int posSep[])
{
	int rank = 0;
	/* First separation */
	posSep[2] = onceQsort(data, 0, lastPos);
	#pragma omp parallel shared(data) private(rank)
	{
		rank = omp_get_thread_num();
		/* Second & Third separation */
		if(rank == 0)
			posSep[1] = onceQsort(data, 0, posSep[2]-1);
		else if(rank==1)
			posSep[3] = onceQsort(data, posSep[2]+1, lastPos);
		#pragma omp barrier
	}
	/* Set first and last pos */
	posSep[4] = lastPos;
	posSep[0] = 0;
}

int main (int argc, char *argv[])
{
	int id = 0;								/* Rank id */
	int size = 0;							/* Size of MPI World, size-1 is number of c core */
	int numT = 0;							/* Get num thread */
	int posCSep[17] = {0};					/* Position separate after first initial for process */
	int posTSep[5] = {0};					/* Position separate after first initial for thread */
	int rank = 0;
	int numData = 0;
	int dataSize = 0;
	int sizeSep[16] = {0};
	double *data;
	MPI_Status status;						/* Status of MPI */
	MPI_Request request;

	/*MPI Initialize */
	MPI_Init(&argc, &argv); 
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
	omp_set_num_threads(atoi(argv[3])); 

	
	//printf("RUNNING\n");
	/*Only sort process & 1 thread*/
	//printf("RUNNING\n");
	FILE *file = NULL;					/* File read */
	if(id==0)
	{
		/* Read file section */
		file = fopen(argv[1],"r");
		if(file == NULL)
		{
			printf("ERROR: Cannot open input file!\n");
			exit(0);
		}
		fscanf(file, "%d", &dataSize);
		data = malloc(dataSize * sizeof(double));
		for(register int i = 0; i < dataSize;)
		{
			fscanf(file, "%lf", &data[i++]);
		}
		fclose(file);
		initCoreQsort(data, dataSize-1, posCSep, size);
		for(register int i = 1; i < size; i++)
			{
				sizeSep[i] = posCSep[i+1] - (posCSep[i]+1)+1;
				// printf("%d:%d - %d + 1 = %d\n",i,posCSep[i+1], posCSep[i]+1,sizeSep[i]);
				MPI_Send(&sizeSep[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
				MPI_Send(&data[posCSep[i]+1], sizeSep[i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
			}
		
		/* Start separated thread */
		numData = posCSep[1];
		initThreadQsort(data, numData, posTSep);
		#pragma omp parallel shared(data) private(rank)
		{
			rank = omp_get_thread_num();
			if(rank==0)
				quickSort(&data[posTSep[rank]] , posTSep[1] + 1);
			else
				quickSort(&data[posTSep[rank]+1] , posTSep[rank+1] - (posTSep[rank]+1) +1);
			#pragma omp barrier
		}

		for(register int i = 1; i < size; i++)
			MPI_Recv(&data[posCSep[i]+1], sizeSep[i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
		
		file = fopen(argv[2], "w");
		fprintf(file, "%d\n", dataSize);
		for(register int i = 0 ; i < dataSize ;i++)
			fprintf(file, "%.4lf\n", data[i]);
		fclose(file);
		free(data);
	}
	else
	{
		MPI_Recv(&numData, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		data = malloc(numData * sizeof(double));
		MPI_Recv(data, numData, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);

		initThreadQsort(data, numData-1, posTSep);
		#pragma omp parallel shared(data) private(rank)
		{
			rank = omp_get_thread_num();
			if(rank==0)
				quickSort(&data[posTSep[rank]] , posTSep[1] + 1);
			else
				quickSort(&data[posTSep[rank]+1] , posTSep[rank+1] - (posTSep[rank]+1) +1);
			#pragma omp barrier
		}
		MPI_Send(data, numData, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		free(data);
	}
	MPI_Finalize();
	return 0;
}
