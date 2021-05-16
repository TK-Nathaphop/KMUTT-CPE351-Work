#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define LEN 100000
#define SIZE 16

void printAll(float *data, int row, int col)
	{
	printf("Row:%d Col:%d\n",row ,col);
	for(int i = 0; i < row*col; i++)
		{
		printf("%4.0lf ",data[i]);
		if((i+1)%col == 0)
			printf("\n");
		}
	}

int main(int argc, char *argv[])
{
    char line[LEN], *token;
	int id = 0;								/* Rank id */
	int size = 0;							/* Size of MPI World, size-1 is number of c core */	
	int i,j,k;								/* Counter */
	FILE *file = NULL;						/* File read */
	MPI_Status status;						/* Status of MPI */
	MPI_Request reqRow, reqCol, reqPanel, request;

	/* Initialize */
	MPI_Init(&argc, &argv); 
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    if(size == 1)
    {
    	float *panel, *ans;
		int rowP, colP, rowAns, colAns;
		int iteration;

		/* Read Panel in heat transfer */
		file = fopen(argv[1],"r");
		fscanf(file, "%d %d", &rowP, &colP);
		panel = malloc(rowP * colP * sizeof(float));

		/* Get data in panel */
		int pos = 0;
		while (fscanf(file, "%f", &panel[i++]) == 1);
		fclose(file);
		/** Debugging */
		// printf("Matrix Panel:\n");
		// printAll(panel, rowP, colP);
		// printf("\n");

		/* Prepare memory for result */
		rowAns = rowP-2;
		colAns = colP-2;
		ans = malloc(rowAns * colAns * sizeof(float));

		sscanf(argv[3], "%d", &iteration);
		float StartCal = MPI_Wtime();
		for(i = 0; i < iteration; i++)
		{
			int rowJ, rowJ1,rowJ2,rowJ3;
			/* Sum first row */
			for(j = 2, pos = 0; j < rowP; j++)
			{
				rowJ1 = j * colP;
				rowJ2 = (j-1) * colP;
				rowJ3 = (j-2) * colP;
				for(k = 2; k < colP; k++)
				{
					ans[pos] = panel[rowJ1+k] + panel[rowJ1+k-1] + panel[rowJ1+k-2]
								+ panel[rowJ2+k] + panel[rowJ2+k-1] + panel[rowJ2+k-2]
								+ panel[rowJ3+k] + panel[rowJ3+k-1] + panel[rowJ3+k-2];
					pos++;

				}
			}

			/* Divide 9 and set answer to panel */
			int rowPSub = rowP-1;
			int colPSub = colP-1;
			
			for(j = 1, pos = 0; j < rowPSub; j++)
			{
				rowJ = j*colP;
				for(k = 1; k < colPSub; k++)
				{
					panel[rowJ+k] = ans[pos]/9;
					pos++;
				}
			}
		}
		// printf("\n");
		// printAll(panel,rowP,colP);

		file = fopen(argv[2], "w");
		fprintf(file, "%d %d\n", rowP, colP);
		i = 0;
		int rowPcolP = rowP*colP;
		while(i < rowPcolP)
			{
			fprintf(file, "%.0f ", panel[i]);
			i++;
			if(i % colP == 0)
				fprintf(file, "\n");
			}
		fclose(file);
		free(panel);
		free(ans);

		float End = MPI_Wtime();
		printf("Time All: %f\n",End-Start);

    }
    else if(id == 0)
    {
		float *panel ,*ans;
		int rowP, colP, rowCal, rowAns, colAns, pos;
		int iteration;
		int rowDiv[SIZE], countSend[SIZE], displs[SIZE] = {0};

		/* Read Panel in heat transfer */
		file = fopen(argv[1],"r");
		fscanf(file, "%d %d", &rowP, &colP);
		panel = malloc(rowP * colP * sizeof(float));

		/* Get data in panel */
		pos = 0;
		while (fscanf(file, "%f", &panel[i++]) == 1);
		fclose(file);

		/* Debugging */
		// printf("Matrix Panel:\n");
		// printAll(panel, rowP, colP);
		// printf("\n");

		/* Find amount of row to send in each core */
		int mod, rowNum = rowP;
		pos = 0;
		mod = rowNum % size;
		rowNum = (rowNum - mod)/size;

		for(i = 0; i < size; i++)
		{
			if(mod > 0)
			{
				rowDiv[i] = rowNum+1;
				mod--;
			}
			else
				rowDiv[i] = rowNum;

			/* Calculate for data that need to send each core and set displacement */
			displs[i] = pos - colP;
			pos += rowDiv[i] * colP;
			countSend[i] = (rowDiv[i]+2) * colP;
		}
		displs[0] = 0;
		/* Decrease amount of row for sending to first core and last core */
		countSend[0] -= colP;
		countSend[size-1] -= colP;

		// for(i = 0; i < size; i++)
		// 	printf("i:%d rowDiv:%d countSend:%d displs:%d\n", i, rowDiv[i], countSend[i], displs[i]);

		int temp;
		rowCal = rowDiv[0] + 1;
		float *pTemp = malloc(rowCal * colP * sizeof(float));
		MPI_Iscatter(rowDiv, 1, MPI_INT, &temp, 1, MPI_INT, 0, MPI_COMM_WORLD, &reqRow);
		MPI_Ibcast(&colP, 1, MPI_INT, 0, MPI_COMM_WORLD, &reqCol);
		MPI_Iscatterv(panel, countSend, displs, MPI_FLOAT, pTemp, rowCal * colP, MPI_FLOAT, 0, MPI_COMM_WORLD, &reqPanel);
		// printAll(panel,rowCal,colP);

		/* Calculation */
		rowAns = rowCal-2;
		colAns = colP-2;
		ans = malloc(rowAns * colAns * sizeof(float));

		sscanf(argv[3], "%d", &iteration);

		float StartCal = MPI_Wtime();
		for(i = 0; i < iteration; i++)
		{
			int rowJ, rowJ1,rowJ2,rowJ3;
			/* Sum first row */
			for(j = 2, pos = 0; j < rowCal; j++)
			{
				rowJ1 = j * colP;
				rowJ2 = (j-1) * colP;
				rowJ3 = (j-2) * colP;
				for(k = 2; k < colP; k++)
				{
					ans[pos] = panel[rowJ1+k] + panel[rowJ1+k-1] + panel[rowJ1+k-2]
								+ panel[rowJ2+k] + panel[rowJ2+k-1] + panel[rowJ2+k-2]
								+ panel[rowJ3+k] + panel[rowJ3+k-1] + panel[rowJ3+k-2];
					pos++;

				}
			}

			/* Divide 9 and set answer to panel */
			int rowCalSub = rowCal-1;
			int colPSub = colP-1;
			for(j = 1, pos = 0; j < rowCalSub; j++)
			{
				rowJ = j*colP;
				for(k = 1; k < colPSub; k++)
				{
					panel[rowJ+k] = ans[pos]/9;
					pos++;
				}
			}

			/* Send last answer to another core */
			if(i+1 < iteration)
			{
				MPI_Send(&panel[(rowCal-2) * colP], colP, MPI_FLOAT, 1, 0, MPI_COMM_WORLD);
				MPI_Recv(&panel[(rowCal-1) * colP], colP, MPI_FLOAT, 1, 0, MPI_COMM_WORLD, &status);
			}
		}

		countSend[0] = rowDiv[0] * colP;
		displs[0] = 0;
		for(i = 1; i < size; i++)
		{
			countSend[i] = rowDiv[i] * colP;
			displs[i] = displs[i-1] + countSend[i-1];
			// printf("i:%d rowDiv:%d countsend:%d displs:%d\n",i, rowDiv[i] ,countSend[i],displs[i]);
		}

		MPI_Gatherv(NULL, 0, MPI_FLOAT, panel, countSend, displs, MPI_FLOAT, 0, MPI_COMM_WORLD);
		// printf("ID:%d\n",id);
		// printAll(panel, rowP, colP);
		
		file = fopen(argv[2], "w");
		fprintf(file, "%d %d\n", rowP, colP);
		i = 0;
		int rowPcolP = rowP*colP;
		while(i < rowPcolP)
			{
			fprintf(file, "%.0f ", panel[i]);
			i++;
			if(i % colP == 0)
				fprintf(file, "\n");
			}
		fclose(file);

		free(panel);
		free(pTemp);

	}
	else
	{
		float *panel, *ans;
		int rowP, colP, rowAns, colAns, pos;
		int iteration;

		MPI_Iscatter(NULL, 1, MPI_INT, &rowP, 1, MPI_INT, 0, MPI_COMM_WORLD, &reqRow);
		MPI_Ibcast(&colP, 1, MPI_INT, 0, MPI_COMM_WORLD, &reqCol);
		MPI_Wait(&reqRow, &status);
		MPI_Wait(&reqCol, &status);

		if(id == size-1)
			rowP += 1;
		else
			rowP += 2;

		panel = malloc(rowP *  colP * sizeof(float));
		MPI_Iscatterv(NULL, NULL, NULL, MPI_FLOAT, panel, rowP * colP, MPI_FLOAT, 0, MPI_COMM_WORLD, &reqPanel);
		MPI_Wait(&reqPanel, &status);

		rowAns = rowP-2;
		colAns = colP-2;
		ans = malloc(rowAns * colAns * sizeof(float));

		sscanf(argv[3], "%d", &iteration);
		for(i = 0; i < iteration; i++)
		{
			int rowJ, rowJ1,rowJ2,rowJ3;
			for(j = 2, pos = 0; j < rowP; j++)
			{
				rowJ1 = j * colP;
				rowJ2 = (j-1) * colP;
				rowJ3 = (j-2) * colP;
				for(k = 2; k < colP; k++)
				{
					ans[pos] = panel[rowJ1+k] + panel[rowJ1+k-1] + panel[rowJ1+k-2]
								+ panel[rowJ2+k] + panel[rowJ2+k-1] + panel[rowJ2+k-2]
								+ panel[rowJ3+k] + panel[rowJ3+k-1] + panel[rowJ3+k-2];
					pos++;

				}
			}

			/* Divide 9 and set answer to panel */
			int rowPSub = rowP-1;
			int colPSub = colP-1;
			
			for(j = 1, pos = 0; j < rowPSub; j++)
			{
				rowJ = j*colP;
				for(k = 1; k < colPSub; k++)
				{
					panel[rowJ+k] = ans[pos]/9;
					pos++;
				}
			}

			if(i+1 < iteration)
			{
				MPI_Recv(panel, colP, MPI_FLOAT, id-1, 0, MPI_COMM_WORLD, &status);
				MPI_Send(&panel[colP], colP, MPI_FLOAT, id-1, 0, MPI_COMM_WORLD);
				if(id != size-1)
				{
					MPI_Send(&panel[(rowP-2) * colP], colP, MPI_FLOAT, id+1, 0, MPI_COMM_WORLD);
					MPI_Recv(&panel[(rowP-1) * colP], colP, MPI_FLOAT, id+1, 0, MPI_COMM_WORLD, &status);
				}
			}
		}
		if(id == size-1)
			rowP -= 1;
		else
			rowP -= 2;
		MPI_Gatherv(&panel[colP], rowP * colP, MPI_FLOAT, NULL, NULL, NULL, MPI_FLOAT, 0, MPI_COMM_WORLD);
		free(panel);
		free(ans);
	}
	MPI_Finalize();
	return 0;
}