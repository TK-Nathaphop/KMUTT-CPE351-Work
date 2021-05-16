#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define LEN 50000
#define FILE_1 "matAlarge.txt"
#define FILE_2 "matBlarge.txt"
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
       	memcpy(array[i], data[i], col*sizeof(double));
    // printAll(array,row,col);
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
	double StartTime,EndTime;				/* Execute time run */
	double StartSend, EndSend;
	double StartRecv,EndRecv;
	double StartCal, EndCal;
	int i,j,k;								/* Counter */

	MPI_Status status;						/* Status of MPI */
	FILE *file = NULL;						/* File read */
	
	/* Initialize */
	MPI_Init(&argc, &argv); 
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
	
	if(id == 0)
		{

		/***************************
		 * Read first file
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
		 * Read second file
		 ***************************/
		file = openFile(FILE_2,"r");
		matrixDiv[1] = readData(file, &rowDiv[1], &colDiv[1], size);
		closeFile(file);

		StartTime = MPI_Wtime(); /* Start Execution Time */

		/***************************
		 * RSending divided matrix A to each core
		 ***************************/
		StartSend = MPI_Wtime();
		for(i = 1; i < size; i++)
			{
			MPI_Send(&rowDiv[0][i], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
			MPI_Send(&colDiv[0][i], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
			MPI_Send(matrixDiv[0][i][0], colDiv[0][i] * rowDiv[0][i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
			}
		EndSend = MPI_Wtime();
		printf("Sending Time MatrixA to another core: %lf\n",EndSend - StartSend);

		/***************************
		 * Sending divided matrix B to each core
		 ***************************/
		StartSend = MPI_Wtime();
		for(i = 1; i < size; i++)
			{
			MPI_Send(&rowDiv[1][i], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
			MPI_Send(&colDiv[1][i], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
			MPI_Send(matrixDiv[1][i][0], colDiv[1][i] * rowDiv[1][i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
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
		printf("Calculation Time: %lf\n", EndCal - StartCal);

		/***************************
		 * Get the result
		 ***************************/
		int rowTemp;		/* For row result from each core */
		int count = 0;		/* Count data sending */
		int pos = 0;		/* Position to free data after finished */

		/* Allocate result array */
		res = allocMatrix(rowRes,colRes);

		/* Get data from current core */
		StartCal = MPI_Wtime();
		pos = setMatrixData(res, sum, 0, rowSum, colSum);
		count = pos;
		EndCal = MPI_Wtime();
		printf("Set Matrix Time: %lf\n", EndCal - StartCal);
		// printAll(res, rowRes, colRes);	/* Debugging */

		/* Get data from another core */
		StartRecv = MPI_Wtime();
		for(i = 1; i < size; i++)
			{
			MPI_Recv(&rowSum, 1, MPI_INT, i, 1, MPI_COMM_WORLD, &status);
			MPI_Recv(&colSum, 1, MPI_INT, i, 1, MPI_COMM_WORLD, &status);
			MPI_Recv(&res[count][0], rowSum*colSum, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, &status);
			count += rowSum;
			// count = setMatrixData(res, sum, count, rowSum, colSum);
			// printf("count: %d, rowSum:%d, colSum:%d\n",count,rowSum,colSum); /* Debugging */
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

		freeMatrix(sum);

		// /* Free matrix */
		freeMatrix(res);

		}
	else
		{

		/***************************
		 * Recieve matrix A 
		 ***************************/
		MPI_Recv(&rowA, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&colA, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		matrixA = allocMatrix(rowA, colA);
		MPI_Recv(matrixA[0], rowA*colA, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
		
		// printf("ID:%d, MatrixA: \n",id);
		// printAll(matrixA, rowA, colA); /* For debugging */

		/***************************
		 * Recieve matrix B
		 ***************************/
		MPI_Recv(&rowB, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&colB, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		matrixB = allocMatrix(rowB, colB);
		MPI_Recv(matrixB[0], rowB*colB, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);

		// printf("ID:%d, MatrixB: \n",id);
		// printAll(matrixB, rowB, colB); /* For debugging */

		/***************************
		 * Sum Matrix
		 ***************************/
		rowSum = rowA;
		colSum = colA;
		sum = sumMatrix(matrixA, matrixB, rowSum, colSum);
		// printAll(sum, rowSum, colSum); /* For debugging */

		MPI_Send(&rowSum, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
		MPI_Send(&colSum, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);	
		MPI_Send(sum[0], rowSum*colSum, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);


		freeMatrix(matrixA);
		freeMatrix(matrixB);
		freeMatrix(sum);
		}
	MPI_Finalize();
	}
