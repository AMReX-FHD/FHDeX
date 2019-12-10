#include <cstdio>
#include <cstdlib>
#include <mpi.h>
#include <lib_mpi_split.h>
#include <mui.h>
#include <cstring>

int main(int argc, char **argv)
{
  MPI_Comm comm = mui::mpi_split_by_app(argc,argv);
  int myrank;
  MPI_Comm_rank(comm,&myrank);

  int nt = 5;

  // FHD discretization
  int nx1 = 4;
  int ny1 = 4;
  double dx1 = 8e-6;
  double dy1 = 8e-6;

  // finer descritization for KMC
  int nx2 = 3;
  int ny2 = 3;
  double dx2 = dx1/nx2;
  double dy2 = dy1/ny2;
  int nx = nx1*nx2;
  int ny = ny1*ny2;

  mui::uniface2d uniface( "mpi://KMC-side/FHD-KMC-coupling" );
  mui::sampler_kmc_fhd2d<double> u({dx1,dy1});
  mui::chrono_sampler_exact2d t;

  char rCH1[100] = "CH_density1";
  char rCH2[100] = "CH_density2";
  char rCH3[100] = "CH_density3";
  char rCH4[100] = "CH_density4";

  char sCH1[100] = "CH_dn1";
  char sCH2[100] = "CH_dn2";
  char sCH3[100] = "CH_dn3";
  char sCH4[100] = "CH_dn4";

  for(int step=1;step<=nt;step++)
  {
    // mui push
    for (int i=0;i<nx;i++)
    {
      for (int j=0;j<ny;j++)
      {
	     double x = (i+0.5)*dx2;
	     double y = (j+0.5)*dy2;

         int dn1 = 1;
         int dn2 = 2;
         int dn3 = 3;
         int dn4 = 4;

         uniface.push(sCH1,{x,y},dn1);
         uniface.push(sCH2,{x,y},dn2);
         uniface.push(sCH3,{x,y},dn3);
         uniface.push(sCH4,{x,y},dn4);
      }
    }

    uniface.commit(step);

    // mui fetch
    for (int i=0;i<nx;i++)
    {
      for (int j=0;j<ny;j++)
      {
	     double x = (i+0.5)*dx2;
	     double y = (j+0.5)*dy2;

         printf("receiver=%d, step=%d: (%e,%e,%e,%e) at (%e,%e)\n",myrank,step,
           uniface.fetch(rCH1,{x,y},step,u,t),
           uniface.fetch(rCH2,{x,y},step,u,t),
           uniface.fetch(rCH3,{x,y},step,u,t),
           uniface.fetch(rCH4,{x,y},step,u,t),
           x,y); 
      }
    }

    uniface.forget(step);
  }

  return 0;
}
