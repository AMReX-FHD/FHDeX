#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mpi.h>
#include <lib_mpi_split.h>
#include <mui.h>

int main(int argc, char **argv)
{
  MPI_Comm comm = mui::mpi_split_by_app(argc,argv);
  int myrank;
  MPI_Comm_rank(comm,&myrank);

  int ntimestamp = 1000;

  double dx = 8.e-6;
  double offset = 0.5*dx;
  int nx = 4;
  int ny = 4;

  mui::uniface2d uniface("mpi://FHD-side/FHD-KMC-coupling");
  mui::sampler_kmc_fhd2d<int> u({dx,dx});
  mui::chrono_sampler_exact2d t;

  char channelp1[100] = "CH_density1";
  char channelp2[100] = "CH_density2";
  char channelp3[100] = "CH_density3";
  char channelp4[100] = "CH_density4";
  char channelp5[100] = "CH_temp";

  char channelf1[100] = "CH_ac1";
  char channelf2[100] = "CH_dc1";
  char channelf3[100] = "CH_ac2";
  char channelf4[100] = "CH_dc2";
  char channelf5[100] = "CH_ac3";
  char channelf6[100] = "CH_dc3";
  char channelf7[100] = "CH_ac4";
  char channelf8[100] = "CH_dc4";

  for (int timestamp=1;timestamp<=ntimestamp;timestamp++)
  {
    // push

    for (int i=0;i<nx;i++) {
      for (int j=0;j<ny;j++) {
        double x = i*dx+offset;
        double y = j*dx+offset;
        double density1 = 1.817055e+19;
        double density2 = 3.604090e+18;
        double density3 = 1.820603e+18;
        double density4 = 8.678933e+17;
        double temp = 300.;
        uniface.push(channelp1,{x,y},density1);
        uniface.push(channelp2,{x,y},density2);
        uniface.push(channelp3,{x,y},density3);
        uniface.push(channelp4,{x,y},density4);
        uniface.push(channelp5,{x,y},temp);
      }
    }

    uniface.commit(timestamp);

    // fetch

    bool file_output = false;
    FILE *fp;

    if (file_output)
    {
      char filename[100];
      sprintf(filename,"res.acdc_ts%d",timestamp);
      fp = fopen(filename,"w");

      fprintf(fp,"# timestamp %d\n",timestamp);
      fprintf(fp,"# x y ac1 dc1 ac2 dc2 ac3 dc3 ac4 dc4\n");
    }

    int sum1 = 0;
    int sum2 = 0;
    int sum3 = 0;
    int sum4 = 0;
    int sum5 = 0;
    int sum6 = 0;
    int sum7 = 0;
    int sum8 = 0;

    for (int i=0;i<nx;i++)
    {
      for (int j=0;j<ny;j++)
      {
        double x = i*dx+offset;
        double y = j*dx+offset;

        int ac1 = uniface.fetch(channelf1,{x,y},timestamp,u,t);
        int dc1 = uniface.fetch(channelf2,{x,y},timestamp,u,t);
        int ac2 = uniface.fetch(channelf3,{x,y},timestamp,u,t);
        int dc2 = uniface.fetch(channelf4,{x,y},timestamp,u,t);
        int ac3 = uniface.fetch(channelf5,{x,y},timestamp,u,t);
        int dc3 = uniface.fetch(channelf6,{x,y},timestamp,u,t);
        int ac4 = uniface.fetch(channelf7,{x,y},timestamp,u,t);
        int dc4 = uniface.fetch(channelf8,{x,y},timestamp,u,t);

        sum1 += ac1;
        sum2 += dc1;
        sum3 += ac2;
        sum4 += dc2;
        sum5 += ac3;
        sum6 += dc3;
        sum7 += ac4;
        sum8 += dc4;

        if (file_output)
        {
          fprintf(fp,"%e\t%e\t",x,y);
          fprintf(fp,"%d\t%d\t",ac1,dc1);
          fprintf(fp,"%d\t%d\t",ac2,dc2);
          fprintf(fp,"%d\t%d\t",ac3,dc3);
          fprintf(fp,"%d\t%d\n",ac4,dc4);
        }
      }
    }

    if (file_output) fclose(fp);

    printf("timestamp=%d:\t",timestamp);
    printf("sum_ac1= %d\tsum_dc1= %d\t",sum1,sum2);
    printf("sum_ac2= %d\tsum_dc2= %d\t",sum3,sum4);
    printf("sum_ac3= %d\tsum_dc3= %d\t",sum5,sum6);
    printf("sum_ac4= %d\tsum_dc4= %d\n",sum7,sum8);

    uniface.forget(timestamp);
  }

  return 0;
}
