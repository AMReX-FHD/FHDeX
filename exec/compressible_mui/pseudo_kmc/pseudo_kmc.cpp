#include <cstdio>
#include <cstdlib>
#include <random>
#include <mpi.h>
#include <lib_mpi_split.h>
#include <mui.h>
#include <cstring>

#define MAXCHLEN 100

int main(int argc, char **argv)
{
  MPI_Comm comm = mui::mpi_split_by_app(argc,argv);
  int myrank;
  MPI_Comm_rank(comm,&myrank);

  int nt = 2;

  // FHD discretization
  int nx1 = 4;
  int ny1 = 4;
  double dx1 = 1.e-5;
  double dy1 = 1.e-5;

  // finer descritization for KMC
  int nx2 = 2;
  int ny2 = 2;
  double dx2 = dx1/nx2;
  double dy2 = dy1/ny2;
  int nx = nx1*nx2;
  int ny = ny1*ny2;

  mui::uniface2d uniface( "mpi://KMC-side/FHD-KMC-coupling" );
  mui::sampler_kmc_fhd2d<double> u({dx1,dy1});
  mui::chrono_sampler_exact2d t;

  char rCH1[MAXCHLEN] = "CH_density1";
  char rCH2[MAXCHLEN] = "CH_density2";
  char rCH3[MAXCHLEN] = "CH_density3";
  char rCH4[MAXCHLEN] = "CH_density4";
  char rCH5[MAXCHLEN] = "CH_temp";

  char sCH1[MAXCHLEN] = "CH_ac1";
  char sCH2[MAXCHLEN] = "CH_ac2";
  char sCH3[MAXCHLEN] = "CH_ac3";
  char sCH4[MAXCHLEN] = "CH_ac4";
  char sCH5[MAXCHLEN] = "CH_dc1";
  char sCH6[MAXCHLEN] = "CH_dc2";
  char sCH7[MAXCHLEN] = "CH_dc3";
  char sCH8[MAXCHLEN] = "CH_dc4";

  FILE *fp = fopen("log.kmc","w");

  std::default_random_engine rng_eng;
  std::poisson_distribution<int> poiss_dist(2.);

  for(int step=1;step<=nt;step++)
  {
    // mui fetch

    fprintf(fp,"## step %d - fetch (x y dens1 dens2 dens3 dens4 temp)\n",step);

    for (int i=0;i<nx;i++)
    {
      for (int j=0;j<ny;j++)
      {
        double x = (i+0.5)*dx2;
        double y = (j+0.5)*dy2;

        fprintf(fp,"%e\t%e\t",x,y);
        fprintf(fp,"%e\t",uniface.fetch(rCH1,{x,y},step,u,t));
        fprintf(fp,"%e\t",uniface.fetch(rCH2,{x,y},step,u,t));
        fprintf(fp,"%e\t",uniface.fetch(rCH3,{x,y},step,u,t));
        fprintf(fp,"%e\t",uniface.fetch(rCH4,{x,y},step,u,t));
        fprintf(fp,"%e\n",uniface.fetch(rCH5,{x,y},step,u,t));
      }
    }

    uniface.forget(step);

    // mui push

    fprintf(fp,"## step %d - push (x y ac1 dc1 ac2 dc2 ac3 dc3 ac4 dc4\n",step);

    for (int i=0;i<nx;i++)
    {
      for (int j=0;j<ny;j++)
      {
        double x = (i+0.5)*dx2;
        double y = (j+0.5)*dy2;

        int ac1 = poiss_dist(rng_eng);
        int ac2 = poiss_dist(rng_eng);
        int ac3 = poiss_dist(rng_eng);
        int ac4 = poiss_dist(rng_eng);
        int dc1 = poiss_dist(rng_eng);
        int dc2 = poiss_dist(rng_eng);
        int dc3 = poiss_dist(rng_eng);
        int dc4 = poiss_dist(rng_eng);

        uniface.push(sCH1,{x,y},ac1);
        uniface.push(sCH2,{x,y},ac2);
        uniface.push(sCH3,{x,y},ac3);
        uniface.push(sCH4,{x,y},ac4);
        uniface.push(sCH5,{x,y},dc1);
        uniface.push(sCH6,{x,y},dc2);
        uniface.push(sCH7,{x,y},dc3);
        uniface.push(sCH8,{x,y},dc4);

        fprintf(fp,"%e\t%e\t",x,y);
        fprintf(fp,"%d\t%d\t",ac1,dc1);
        fprintf(fp,"%d\t%d\t",ac2,dc2);
        fprintf(fp,"%d\t%d\t",ac3,dc3);
        fprintf(fp,"%d\t%d\n",ac4,dc4);
      }
    }

    uniface.commit(step);
  }

  fclose(fp);

  return 0;
}