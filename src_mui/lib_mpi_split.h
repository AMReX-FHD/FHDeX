/*
 * lib_mpi_split.h
 *
 *  Created on: Mar 14, 2014
 *      Author: ytang
 */

#ifndef LIB_MPI_SPLIT_H_
#define LIB_MPI_SPLIT_H_

using namespace std;

namespace mui {

inline void mpi_finalize_after_split() {
    int flag;
    MPI_Finalized(&flag);
    if (!flag) MPI_Finalize();
}

inline MPI_Comm mpi_split_by_app( int argc=0, char **argv=NULL )
{
    {
        int flag;
        MPI_Initialized(&flag);
        if( !flag ) {
            MPI_Init( &argc, &argv );
            atexit( mpi_finalize_after_split );
        }
    }

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    int appnum;
#ifdef SLURM
    // expected args: exec %t %o ...
    // consecutive task numbers are assumed for each app
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
    return domain;
}

}

#endif /* LIB_MPI_SPLIT_H_ */