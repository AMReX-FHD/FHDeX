#include <cstdio>
#include <iostream>
#include <cmath>
#include <random>
#include <cassert>

using namespace std;

int main()
{
    double kB = 1.38064852e-16;
    double Navo = 6.0221409e+23;

    double Temp = 300;

    double M = 83.8000;
    double m = M/Navo;

    double vT = sqrt(2*kB*Temp/m);
    double vT2 = sqrt(kB*Temp/m);

    double Vx = 0.;
    double Vy = 0.;
    double Vz = 0.1*vT;

    double a = Vz/vT;
    assert(a>0.);

    default_random_engine rng;
    normal_distribution<double> normal(0.,1.);
    uniform_real_distribution<double> unif(0.,1.);

    int Nsample = 1e6;

    double arpi = a*sqrt(M_PI);

    for (int i=0;i<Nsample;i++)
    {
        double vx = vT2*normal(rng);
        double vy = vT2*normal(rng);

        double z;

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

        double vz = (a-z)*vT;

        cout <<  vx << '\t';
        cout <<  vy << '\t';
        cout <<  vz << endl;
    }

    return 0;
}
