#include <cstdio>
#include <iostream>
#include <cmath>
#include <random>
#include <cassert>

using namespace std;

default_random_engine rng;
normal_distribution<double> normal(0.,1.);
uniform_real_distribution<double> unif(0.,1.);

double sample_Maxwell_inflow_normal(double a)
{
    double z;

    if (a<=0)
    {
        while (true)
        {
            z = -sqrt(a*a-log(unif(rng)));
            if ((a-z)/(-z)>unif(rng)) break;
        }
    }
    else
    {
        double arpi = a*sqrt(M_PI);

        while (true)
        {
            double u = unif(rng);

            if (arpi/(arpi+1+a*a) > u)
            {
                z = -fabs(normal(rng))/sqrt(2);
                break;
            }
            else if ((arpi+1)/(arpi+1+a*a) > u)
            {
                z = -sqrt(-log(unif(rng)));
                break;
            }
            else
            {
                z = (1-sqrt(unif(rng)))*a;

                if (exp(-z*z)>unif(rng)) break;
            }
        }
    }

    return z;
}

int main(int argc,char **argv)
{
    assert(argc==2);

    double kB = 1.38064852e-16;
    double Navo = 6.0221409e+23;

    double Temp = 300;

    double M = 83.8000;
    double m = M/Navo;

    double vT = sqrt(2*kB*Temp/m);
    double vT2 = sqrt(kB*Temp/m);

    double Vx = 0.;
    double Vy = 0.;

    double a = -atof(argv[1]);
    double Vz = a*vT;

    int Nsample = 1e6;

    for (int i=0;i<Nsample;i++)
    {
        double vx = vT2*normal(rng);
        double vy = vT2*normal(rng);

        double z = sample_Maxwell_inflow_normal(a);
        double vz = -(a-z)*vT;

        cout <<  vx << '\t';
        cout <<  vy << '\t';
        cout <<  vz << endl;
    }

    return 0;
}
