#include <mpi.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
    int id = 0;
    int size = 0;
    MPI_Status status;
    char msg[32];
    MPI_Init(&argc, &argv); 
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    if(id == 0)
    {
        for(int i = 1; i < size; i++)
        {
                MPI_Recv(msg,sizeof(msg), MPI_BYTE, i, 0, MPI_COMM_WORLD, &status);
                printf("%s\n",msg);
        }

    }
    else
    {
            sprintf(msg,"Hello rank 0. I'm rank %d.",id);
            MPI_Send(msg, sizeof(msg), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
    }
    MPI_Finalize();
}

