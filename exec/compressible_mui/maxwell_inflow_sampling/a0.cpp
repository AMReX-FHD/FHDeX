#include <cstdio>
#include <iostream>
#include <cmath>
#include <random>

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

    default_random_engine rng;
    normal_distribution<double> normal(0.,1.);
    uniform_real_distribution<double> unif(0.,1.);

    int Nsample = 1e6;

    for (int i=0;i<Nsample;i++)
    {
        double vx = vT2*normal(rng);
        double vy = vT2*normal(rng);
        double vz = vT*sqrt(-log(unif(rng)));

        cout <<  vx << '\t';
        cout <<  vy << '\t';
        cout <<  vz << endl;
    }

    return 0;
}
