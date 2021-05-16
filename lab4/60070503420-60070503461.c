#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[])
{
	int id = 0;								/* Rank id */
	int size = 0;							/* Size of MPI World, size-1 is number of c core */	
	int i,j,k;								/* Counter */
	FILE *file = NULL;						/* File read */
	MPI_Status status;						/* Status of MPI */

	/* Initialize */
	MPI_Init(&argc, &argv); 
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    if(size == 1)
    {
    	double *matrixA, *matrixB ,*product;	/* Matrix A data */
		int rowA, colA ,rowB, colB, temp1, temp2, temp3;

		file = fopen(argv[1],"r"); // Read matrix A
		fscanf(file, "%d %d", &rowA, &colA);
		matrixA = malloc(rowA * colA * sizeof(double));
		i = 0;
		while (fscanf(file, "%lf", &matrixA[i++]) == 1);
		fclose(file);
		// printf("MatrixA:\n");
		// printAll(matrixA, rowA, colA);

		file = fopen(argv[2],"r"); // Read matrix B
		fscanf(file, "%d %d", &rowB, &colB);
		matrixB = malloc(rowB * colB * sizeof(double));
		i = 0; //Row
		j = 0; //Col
		while (fscanf(file, "%lf", &matrixB[(j*rowB)+i]) == 1)
		{
            j++;
			if (j == colB)
            {
                i++;
                j=0;
            }
		}
		fclose(file);
		// printf("\nMatrixB:\n");
		// printAll(matrixB, colB, rowB);

		product = calloc(rowA * colB, sizeof(double));
		for(i = rowA; i--;)
		{
			temp1 = i*colA;
			temp2 = i*colB;
			for(j = colB; j--;)
			{
				temp3 = j*colA;
				for(k = rowB; k--;)
					product[temp2+j] += matrixA[temp1+k] * matrixB[temp3+k];
			}
		}
		// for(i = rowA; i--;)
		// {
		// 	int iColA = i*colA;
		// 	int iColB = i*colB;
		// 	for(j = colB; j--;)
		// 	{
		// 		int jColA = j*colA;
		// 		for(k = rowB; k--;)
		// 			product[iColB+j] += matrixA[iColA+k] * matrixB[jColA+k];
		// 	}
		// }

		file = fopen(argv[3], "w");
		fprintf(file, "%d %d\n", rowA, colB);
		int i = 0;
		int rowAcolB = rowA*colB;
		while(i < rowAcolB)
			{
			fprintf(file, "%.10lf ",product[i]);
			i++;
			if( i % colB == 0)
				fprintf(file, "\n");
			}
		// for(int i = 0; i < rowA*colB; i++)
		// {
		// 	fprintf(file, "%.10lf ",product[i]);
		// 	if((i+1) % colB == 0)
		// 		fprintf(file, "\n");
		// }
		fclose(file);
		free(matrixA);
		free(matrixB);
		free(product);

    }
    else if(id == 0)
    {
		double *matrixA, *matrixB ,*product;	/* Matrix A data */
		int rowA, colA, rowB, colB ,temp1, temp2, temp3;
		int *rowDiv;
		file = fopen(argv[1],"r"); // Read matrix A
		fscanf(file, "%d %d", &rowA, &colA);
		matrixA = malloc(rowA * colA * sizeof(double));
		i = 0;
		while (fscanf(file, "%lf", &matrixA[i++]) == 1);
		fclose(file);

		/* Divide Row */
		int numerator = 1, rowNum = rowA;
		rowDiv = malloc (size * sizeof(int));
		numerator = rowNum % size;
		rowNum = (rowNum - numerator)/size;
		for(i = 0; i < size; i++)
			if(numerator > 0)
			{
				rowDiv[i] = rowNum+1;
				numerator--;
			}
			else
				rowDiv[i] = rowNum;

		/* Send from A to another core */
		int sizeA;
		MPI_Scatter(rowDiv, 1, MPI_INT, &sizeA, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&colA, 1, MPI_INT, 0, MPI_COMM_WORLD);
		int *countSend = malloc (size * sizeof(int));
		int *displs = malloc (size * sizeof(int));
		countSend[0] = rowDiv[0] * colA;
		displs[0] = 0;
		for(i = 1; i < size; i++)
		{
			countSend[i] = rowDiv[i] * colA;
			displs[i] = displs[i-1] + countSend[i-1];
		}
		double *tempMat = malloc(rowDiv[0] * colA * sizeof(double));
		MPI_Scatterv(&matrixA[0], countSend, displs, MPI_DOUBLE, tempMat, countSend[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);
		free(tempMat);

		MPI_Bcast(&rowB, 1, MPI_INT, 1, MPI_COMM_WORLD);
		MPI_Bcast(&colB, 1, MPI_INT, 1, MPI_COMM_WORLD);
		matrixB = malloc(rowB * colB * sizeof(double));
		MPI_Bcast(matrixB, rowB*colB, MPI_DOUBLE, 1, MPI_COMM_WORLD);

		product = calloc(rowA * colB, sizeof(double));
		for(i = sizeA; i--;)
		{
			temp1 = i*colA;
			temp2 = i*colB;
			for(j = colB; j--;)
			{
				temp3 = j*colA;
				for(k = rowB; k--;)
					product[temp2+j] += matrixA[temp1+k] * matrixB[temp3+k];
			}
		}

		countSend[0] = sizeA * colB;
		displs[0] = 0;
		for(i = 1; i < size; i++)
		{
			countSend[i] = rowDiv[i] * colB;
			displs[i] = displs[i-1] + countSend[i-1];
		}
		MPI_Gatherv(NULL, 0, MPI_DOUBLE, product, countSend, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		file = fopen(argv[3], "w");
		fprintf(file, "%d %d\n", rowA, colB);

		int i = 0;
		int rowAcolB = rowA*colB;
		while(i < rowAcolB)
			{
			fprintf(file, "%.10lf ",product[i]);
			i++;
			if( i % colB == 0)
				fprintf(file, "\n");
			}
		// for(int i = 0; i < rowA*colB; i++)
		// {
		// 	fprintf(file, "%.10lf ",product[i]);
		// 	if((i+1) % colB == 0)
		// 		fprintf(file, "\n");
		// }
		fclose(file);
	}
	else if(id == 1)
	{
		double *matrixA, *matrixB ,*product;	
		int rowA, colA, rowB, colB ,temp1, temp2, temp3;

		file = fopen(argv[2],"r"); // Read matrix B
		fscanf(file, "%d %d", &rowB, &colB);
		matrixB = malloc(rowB * colB * sizeof(double));
		i = 0; //Row
		j = 0; //Col
		while (fscanf(file, "%lf", &matrixB[(j*rowB)+i]) == 1)
		{
            j++;
			if (j == colB)
            {
                i++;
                j=0;
            }
		}
		fclose(file);

		MPI_Scatter(NULL, 1, MPI_INT, &rowA, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&colA, 1, MPI_INT, 0, MPI_COMM_WORLD);
		matrixA = malloc(rowA * colA * sizeof(double));
		MPI_Scatterv(NULL, NULL, NULL, MPI_DOUBLE, matrixA, rowA*colA, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		
		MPI_Bcast(&rowB, 1, MPI_INT, 1, MPI_COMM_WORLD);
		MPI_Bcast(&colB, 1, MPI_INT, 1, MPI_COMM_WORLD);
		MPI_Bcast(matrixB, rowB * colB, MPI_DOUBLE, 1, MPI_COMM_WORLD);

		product = calloc(rowA * colB, sizeof(double));
		for(i = rowA; i--;)
		{
			temp1 = i*colA;
			temp2 = i*colB;
			for(j = colB; j--;)
			{
				temp3 = j*colA;
				for(k = rowB; k--;)
					product[temp2+j] += matrixA[temp1+k] * matrixB[temp3+k];
			}
		}
		MPI_Gatherv(product, rowA*colB, MPI_DOUBLE, NULL, NULL, NULL, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
	else
	{
		double *matrixA, *matrixB, *product;	/* Matrix data */
		int rowA, colA,rowB, colB ,temp1, temp2, temp3;

		MPI_Scatter(NULL, 1, MPI_INT, &rowA, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&colA, 1, MPI_INT, 0, MPI_COMM_WORLD);
		matrixA = malloc(rowA * colA * sizeof(double));
		MPI_Scatterv(NULL, NULL, NULL, MPI_DOUBLE, matrixA, rowA*colA, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		
		MPI_Bcast(&rowB, 1, MPI_INT, 1, MPI_COMM_WORLD);
		MPI_Bcast(&colB, 1, MPI_INT, 1, MPI_COMM_WORLD);
		matrixB = malloc(rowB * colB * sizeof(double));
		MPI_Bcast(matrixB, rowB * colB, MPI_DOUBLE, 1, MPI_COMM_WORLD);

		product = calloc(rowA * colB, sizeof(double));
		for(i = rowA; i--;)
		{
			temp1 = i*colA;
			temp2 = i*colB;
			for(j = colB; j--;)
			{
				temp3 = j*colA;
				for(k = rowB; k--;)
					product[temp2+j] += matrixA[temp1+k] * matrixB[temp3+k];
			}
		}
		MPI_Gatherv(product, rowA*colB, MPI_DOUBLE, NULL, NULL, NULL, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
	MPI_Finalize();
	return 0;
}