#include <mpi.h>
#include <iostream>
#include <cstdlib>

using namespace std;

int main (int argc, char* argv[])
{
    MPI_Init(&argc,&argv);

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    int appnum;
#ifdef SLURM
    if (argc<3) MPI_Abort(MPI_COMM_WORLD,0);
    int task = atoi(argv[1]);
    int offset = atoi(argv[2]);
    appnum = task-offset;
#else
    void* v;
    int flag;
    MPI_Comm_get_attr(MPI_COMM_WORLD,MPI_APPNUM,&v,&flag);
    appnum = *(int*)v;
#endif

    MPI_Comm domain;
    MPI_Comm_split(MPI_COMM_WORLD,appnum,rank,&domain);

    cout << rank << "/" << size << ": appnum " << appnum << " domain " << domain << "\n" << flush;

    MPI_Finalize();

    return 0;
}
