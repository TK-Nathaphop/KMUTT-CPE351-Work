#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define LEN 50000
#define FILE_1 "matAlarge.txt"
#define FILE_2 "matBlarge.txt"
//#define FILE_1 "testMat2.txt"
//#define FILE_2 "TestMat2.txt"
#define RESULT_FILE "result.txt"
#define TRUE 1
#define FALSE 0

/* Open file and check that can open or not */
FILE *openFile(char name[], char type[])
	{
	FILE* file = fopen(name,type);
	if(file == NULL)
		{
		printf("Error, can't read file\n");
		exit(0);
		}
	return file;
	}


/* Close the file */
void closeFile(FILE *file)
	{
	fclose(file);
	}

/* Print all data in 2D array */
void printAll(double **data, int row, int col)
	{
	for(int i = 0; i < row; i++)
		{
		for(int j = 0; j < col; j++)
			printf("%.01lf ",data[i][j]);
		printf("\n");
		}
	}

/* Allocate matrix memory */
double **allocMatrix(int row, int col)
	{
	double *data = (double *)calloc(row*col, sizeof(double));
	if(data == NULL)
		{
		printf("Error, can't allocate the 1d array\n");
		exit(0);
		}
    double **array= (double **)calloc(row, sizeof(double*));
    if(array == NULL)
			{
			printf("Error, can't allocate the 2d array\n");
			exit(0);
			}
    for (int i=0; i<row; i++)
        array[i] = &(data[col*i]);
    return array;
	}

void freeMatrix(double **array)
{
	free(array[0]);
	free(array);
}

/* Divide number of row for each core, return array for number of each row back*/
int *calcRow(int rowSize, int coreSize)
{
	int numerator = 1, size = rowSize;
	int *row;
	row = calloc(coreSize,sizeof(int));
	if(row == NULL)
		{	
		printf("Error! Cannot allocate memory\n");
		exit(0);
		}
	numerator = size%coreSize;
	size = (size-numerator)/coreSize;
	for(int i = 0; i < coreSize; i++)
		if(numerator > 0)
		{
			row[i] = size+1;
			numerator--;
		}
		else
			row[i] = size;
	return row;
}

/* Divide number of col for each core, return array for number of each col back*/
int *calcCol(int colSize, int coreSize)
{
	int *col;
	col = calloc(coreSize,sizeof(int));
	if(col == NULL)
		{	
		printf("Error! Cannot allocate memory\n");
		exit(0);
		}
	for(int i = 0; i < coreSize; i++)
		col[i] = colSize;
	return col;
}

double **sumMatrix(double **matrixA, double **matrixB, int row, int col)
{
	double **res = allocMatrix(row,col);
	for(int i = 0; i < row; i++)
		for(int j = 0; j < col; j++)
			res[i][j] = matrixA[i][j] + matrixB[i][j];
	return res;
}


/* Calculate amount of row and column in each core. Then allocate and prepare data for each core.
 * In the end send divided data
 */
double ***readData(FILE *file, int *rowData[], int *colData[], int size)
	{
	char input[LEN] = {0}; 			/* Input from reading file */
	char line[LEN] = {0};			/* Read line from text */
	char *token = NULL;
	double ***data = NULL;			/* Data that read from file */
	int row, col;					/* All row and column */
	int i, j, k;

	data = (double***)calloc(size, sizeof(double**));
	if(data == NULL)
		{
		printf("Error, can't allocate the 3d array\n");
		exit(0);
		}
	fgets(line, sizeof(input),file);
	sscanf(line,"%d %d", &row, &col);
	(*rowData) = calcRow(row, size);
	(*colData) = calcCol(col, size);
	/* Debugging */
	// for(i = 0; i < size; i++)
	// 	printf("=%d Row:%d, Col:%d\n",i,(*rowData)[i],(*colData)[i]);

	/* Allocate 3d array and 2d array to send to each core */
	for(i = 0; i < size; i++)
		data[i] = allocMatrix((*rowData)[i], (*colData)[i]);


	/* Initial value before read from line */
	i = 0;
	j = 0;

	/* Loop read and divide to keep in each matrix to send for computation
	 * i - matrix element to send to each core
	 * j - column
	 * k - row
	 */
	while(fgets(line, sizeof(input), file) != NULL)
		{
		token = strtok(line, " ");
		for (k = 0; k < (*colData)[i]; k++)
			{
			sscanf(token, "%lf", &data[i][j][k]);
			// printf("token: %s \n",token);	/* Debugging */
			token = strtok(NULL, " ");
			}		

		/* If it's last of row in current element, set to next element with new row */
		if(j == ((*rowData)[i] - 1))
			{
			i++;
			j = 0;
			}
		else
			j++;
		}
		/* For debugging */
		// for(int i = 0; i < size; i++)
		// 	printAll(data[i],(*rowData)[i],(*colData)[i]);
	return data;
	}

void writeFile(FILE *file, double **matrix, int row, int col)
{
	fprintf(file, "%d %d\r\n", row, col);
	for(int i = 0; i < row; i++)
		{
		for(int j = 0; j < col; j++)
			{
			fprintf(file, "%.01lf",matrix[i][j]);
			if(j == col-1)
				fprintf(file, "\r\n");
			else
				fprintf(file, " ");
			}
		}
	printf("Output file: %s\n",RESULT_FILE);
}

int setMatrixData(double **array, double **data, int posRow, int row, int col)
{
	int i,j;
	for(i = posRow; i < row+posRow; i++)
		{
		for(j = 0; j < col; j++)
			{
			array[i][j] = data[i-posRow][j];
			}

		}
	return i;
}

int main(int argc, char *argv[])
	{
	int id = 0;								/* Rank id */
	int size = 0;							/* Size of MPI World, size-1 is number of c core */	
	double ***matrixDiv[2];					/* Matrix divide into array for sending to each core */
	int *rowDiv[2], *colDiv[2];
	double **matrixA;	/* Matrix A data */
	int rowA, colA;
	double **matrixB;	/* Matrix B data */
	int rowB, colB;
	double **res;							/* Result matrix */
	int rowRes, colRes;
	double **sum;							/* Result matrix */
	int rowSum, colSum;
	int resultRow[16],resultCol[16];
	double StartTime,EndTime;				/* Execute time run */
	double StartSend, EndSend;
	double StartRecv,EndRecv;
	double StartCal, EndCal;
	int i,j,k;								/* Counter */
	int requestTag = 0;
	int countReceiveSum = 0;
	int divideRow[16] = {0};
	int insertPosition = 0;
	int readyReceive[16] = {0};

	MPI_Status status;						/* Status of MPI */
	MPI_Request request[24];				/* Request of MPI (non-blocking)*/
	MPI_Request requestSum[16];
	int flagRecv[64] = {0};
	int flagRecvSum[16] = {0};
	FILE *file = NULL;						/* File read */
	
	/* Initialize */
	MPI_Init(&argc, &argv); 
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
	
	if(id == 0)
		{

		/***************************
		 * Read first file and send to each core
		 ***************************/
		file = openFile(FILE_1,"r");
		matrixDiv[0] = readData(file, &rowDiv[0], &colDiv[0], size);
		closeFile(file);

		/* Sending divided matrix A to each core */

		/*Prepare the row and column of result */
		for(i = 0, rowRes = 0; i < size; i++)
			rowRes = rowRes + rowDiv[0][i];
		colRes = colDiv[0][0];
		// printf("Row Result: %d, Col Result: %d\n",rowRes,colRes); /* Debugging */

		/***************************
		 * Read second file and send to each core
		 ***************************/
		file = openFile(FILE_2,"r");
		matrixDiv[1] = readData(file, &rowDiv[1], &colDiv[1], size);
		closeFile(file);

		StartTime = MPI_Wtime(); /* Start Execution Time */
		StartSend = MPI_Wtime();
		/* Sending divided matrix A to each core */
		for(i = 1; i < size; i++)
			{
			//printf("RANK0 : sending Matrix A to rank %d / Dimension %d %d\n",i,rowDiv[0][i],colDiv[0][i]);
			divideRow[i] = rowDiv[0][i];
			MPI_Isend(&rowDiv[0][i], 1, MPI_INT, i, 11, MPI_COMM_WORLD ,&request[0]);
			MPI_Isend(&colDiv[0][i], 1, MPI_INT, i, 22, MPI_COMM_WORLD ,&request[1]);
			MPI_Isend(matrixDiv[0][i][0], colDiv[0][i] * rowDiv[0][i], MPI_DOUBLE, i, 33, MPI_COMM_WORLD ,&request[2]);
			}
		EndSend = MPI_Wtime();
		printf("Sending Time MatrixA to another core: %lf\n",EndSend - StartSend);

		StartSend = MPI_Wtime();
		/* Sending divided matrix B to each core */
		for(i = 1; i < size; i++)
			{
			//printf("RANK0 : sending Matrix B to rank %d / Dimension %d %d\n",i,rowDiv[1][i],colDiv[1][i]);
			MPI_Isend(&rowDiv[1][i], 1, MPI_INT, i, 44, MPI_COMM_WORLD ,&request[3]);
			MPI_Isend(&colDiv[1][i], 1, MPI_INT, i, 55, MPI_COMM_WORLD ,&request[4]);
			MPI_Isend(matrixDiv[1][i][0], colDiv[1][i] * rowDiv[1][i], MPI_DOUBLE, i, 66, MPI_COMM_WORLD ,&request[5]);
			}
		EndSend = MPI_Wtime();

		printf("Sending Time MatrixB to another core: %lf\n",EndSend - StartSend);


		/***************************
		 * Calculate
		 ***************************/
		StartCal = MPI_Wtime();
		matrixA = matrixDiv[0][0];
		rowA = rowDiv[0][0];
		colA = colDiv[0][0];

		matrixB = matrixDiv[1][0];
		rowB = rowDiv[1][0];
		colB = colDiv[1][0];
		// printf("ID:%d, Row:%d, Col:%d\n",id,rowSum,colSum);
		rowSum = rowA;
		colSum = colA;
		sum = sumMatrix(matrixA, matrixB, rowSum, colSum);
		EndCal = MPI_Wtime();
		printf("Calculation Time at rank0 : %lf\n", EndCal - StartCal);

		/***************************
		 * Get the result
		 ***************************/
		int rowTemp;		/* For row result from each core */
		int count = 0;		/* Count data sending */

		/* Allocate result array */
		res = allocMatrix(rowRes,colRes);

		/* Get data from current core */
		count = setMatrixData(res, sum, 0, rowSum, colSum);
		// printAll(res, rowRes, colRes);	/* Debugging */
		freeMatrix(sum);

		StartRecv = MPI_Wtime();
		/* Get data from another core */

		/*Get row and col from another core */
		for(i = 1; i < size; i++)
			{
			requestTag = (i - 1) * 2;
			MPI_Irecv(&resultRow[i], 1, MPI_INT, i, 111, MPI_COMM_WORLD, &request[requestTag]);
			requestTag++;
			MPI_Irecv(&resultCol[i], 1, MPI_INT, i, 122, MPI_COMM_WORLD, &request[requestTag]);
			}

		/*Check until receive all matrix from all core*/
		while(countReceiveSum < size - 1 )
			{
			for(i = 1; i < size; i++)
				{
				requestTag = (i - 1) * 2;
				MPI_Test(&request[requestTag],&flagRecv[requestTag],&status);
				requestTag++;
				MPI_Test(&request[requestTag],&flagRecv[requestTag],&status);
				}
			
			countReceiveSum = 0;
			for(i = 1; i < size; i++)
				{
				requestTag = (i - 1) * 2;
				if(flagRecv[requestTag] && flagRecv[requestTag + 1])
					{
					/*First time that row and col of each rank ready, will start receive matrix data*/
					if(readyReceive[i] == 0)
						{
						//fine insert position to result matrix
						insertPosition = count;
						for (j = 1 ; j < i ; j++)
							insertPosition += divideRow[j];
						//printf("insertPosition of rank %d is %d\n",i,insertPosition);

						//printf("Ready to receive from rank %d\n",i);
						MPI_Irecv(&res[insertPosition][0], (resultRow[i]) * (resultCol[i]), MPI_DOUBLE, i, 133, MPI_COMM_WORLD, &requestSum[i]);
						readyReceive[i] = 1;
						}
			
					/*Check and count, how many rank have done transfer matrix*/
					MPI_Test(&requestSum[i],&flagRecvSum[i],&status);
					countReceiveSum += flagRecvSum[i]; 

					//printf("countReceiveSum = %d\n",countReceiveSum);
					}
				}
			}

		EndRecv = MPI_Wtime();
		printf("Recieve Message fome all core: %lf\n",EndRecv-StartRecv);
		EndTime = MPI_Wtime();
		printf("Timings: %lf Sec\n", EndTime-StartTime);
		/* Print result here, comment */
		// printAll(res,rowRes,colRes);

		/* Write result to file (Optional). Just change btw TRUE and FALSE */
		if(TRUE)
			{
			file = openFile(RESULT_FILE, "w");
			writeFile(file, res, rowRes, colRes);
			closeFile(file);
			}
		
		/* Free divided matrix A*/
		for(int i = 1; i < size; i++)
			freeMatrix(matrixDiv[0][i]);
		free(matrixDiv[0]);
		free(rowDiv[0]);
		free(colDiv[0]);

		/* Free divided matrix B */
		for(i = 1; i < size; i++)
			freeMatrix(matrixDiv[1][i]);
		free(matrixDiv[1]);
		free(rowDiv[1]);
		free(colDiv[1]);

		// /* Free matrix */
		freeMatrix(res);

		}
	else
		{

		/***************************
		 * Recieve matrix A and B
		 ***************************/
		MPI_Irecv(&rowA, 1, MPI_INT, 0, 11, MPI_COMM_WORLD, &request[0]);
		MPI_Irecv(&colA, 1, MPI_INT, 0, 22, MPI_COMM_WORLD, &request[1]);
		MPI_Irecv(&rowB, 1, MPI_INT, 0, 44, MPI_COMM_WORLD, &request[2]);
		MPI_Irecv(&colB, 1, MPI_INT, 0, 55, MPI_COMM_WORLD, &request[3]);
		
		flagRecv[0] = 0;
		flagRecv[1] = 0;
		flagRecv[2] = 0;
		flagRecv[3] = 0;
		flagRecv[4] = 0;
		flagRecv[5] = 0;
		matrixA = NULL;
		matrixB = NULL;

		//printf("RANK%d : A == %p and B == %p\n",id,matrixA,matrixB);
		while(flagRecv[4] == 0 || flagRecv[5] == 0)
			{
			MPI_Test(&request[0],&flagRecv[0],&status);
			MPI_Test(&request[1],&flagRecv[1],&status);
			MPI_Test(&request[2],&flagRecv[2],&status);
			MPI_Test(&request[3],&flagRecv[3],&status);

			if(flagRecv[0] && flagRecv[1] && matrixA == NULL)
				{
				//printf("RANK%d : allocMatrix for A by %d %d\n", id, rowA, colA);
				matrixA = allocMatrix(rowA, colA);
				MPI_Irecv(matrixA[0], rowA*colA, MPI_DOUBLE, 0, 33, MPI_COMM_WORLD, &request[4]);
				MPI_Test(&request[4],&flagRecv[4],&status);
				}

			if(flagRecv[2] && flagRecv[3] && matrixB == NULL)
				{
				//printf("RANK%d : allocMatrix for B by %d %d\n",id, rowB, colB);
				matrixB = allocMatrix(rowB, colB);
				MPI_Irecv(matrixB[0], rowB*colB, MPI_DOUBLE, 0, 66, MPI_COMM_WORLD, &request[5]);
				MPI_Test(&request[5],&flagRecv[5],&status);
				}
			
			if(matrixA != NULL)
				MPI_Test(&request[4],&flagRecv[4],&status);
			if(matrixB != NULL)
				MPI_Test(&request[5],&flagRecv[5],&status);

			//printf("RANK%d : flagRecv[4] == %d and flagRecv[5] == %d\n",id,flagRecv[4] ,flagRecv[5]);
			}
		//printf("RANK%d : Get two matrix A == %p and B == %p\n",id,matrixA,matrixB);

		// printf("ID:%d, MatrixB: \n",id);
		// printAll(matrixB, rowB, colB); /* For debugging */

		/***************************
		 * Sum Matrix
		 ***************************/
		rowSum = rowA;
		colSum = colA;
		sum = sumMatrix(matrixA, matrixB, rowSum, colSum);
		// printAll(sum, rowSum, colSum); /* For debugging */

		//printf("RANK%d : I have my result, sending to rank0\n",id);
		MPI_Isend(&rowSum, 1, MPI_INT, 0, 111, MPI_COMM_WORLD, &request[6]);
		MPI_Isend(&colSum, 1, MPI_INT, 0, 122, MPI_COMM_WORLD, &request[7]);	
		MPI_Isend(sum[0], rowSum*colSum, MPI_DOUBLE, 0, 133, MPI_COMM_WORLD, &request[8]);

		MPI_Wait(&request[8], &status);

		freeMatrix(matrixA);
		freeMatrix(matrixB);
		freeMatrix(sum);

		//printf("RANK%d : I'm done my job, Bye!\n",id);
		}
	MPI_Finalize();
	}
