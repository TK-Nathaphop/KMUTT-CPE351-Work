#include <omp.h> 
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


int main (int argc, char *argv[]) { 
	FILE *file = NULL;
	int rowA,colA,rowB,colB;
	double stime,etime;
	int nt;
	omp_set_num_threads(atoi(argv[4])); 
	
	stime = clock();
	/*open and get col row file A*/
	file = fopen(argv[1],"r"); 
	fscanf(file, "%d %d", &rowA, &colA);
	double *matA;
	matA = malloc(rowA * colA * sizeof(double));
	/*read file A*/
	int i;
	i = 0;
	while (fscanf(file, "%lf", &matA[i++]) == 1);
	fclose(file);

	/*open and get col row file B*/
	file = fopen(argv[2],"r"); 
	fscanf(file, "%d %d", &rowB, &colB);
	double *matB;
	matB = malloc(rowB * colB * sizeof(double));
	/*read file B*/

	i = 0;
	while (fscanf(file, "%lf", &matB[i++]) == 1);
	fclose(file);

	/*check size of matrix*/
	if( (colA != colB) || (rowA != rowB) )
		{
		printf("ERROR : row col not match\n");
		exit(0);
		}


	/*Add matrix by parallel*/
	int rowCol = rowA * colA;
	#pragma omp parallel for
	for(int i = 0; i < rowCol ; i++)
		{
		matA[i] = matA[i] + matB[i];
		}
	free(matB);

	/*write file*/
	file = fopen(argv[3], "w");
	fprintf(file, "%d %d", rowA, colA);
	for(int i = 0; i < rowCol ; i++)
		{
		if( i % colA == 0)
		fprintf(file, "\n");
		fprintf(file, "%.1f ",  matA[i]);
		
		}
	fprintf(file, "\n");

	fclose(file);
	free(matA);

	etime = clock();
	printf("RESULT : time = %lf\n",(double)(etime - stime)/(CLOCKS_PER_SEC));
} 
