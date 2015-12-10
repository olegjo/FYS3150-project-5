#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>

using namespace std;

// function for gaussian random numbers
double gaussian_deviate(long *);
// ran2 for uniform deviates, initialize with negative seed.
double ran2(long *);



int main(int argc, char* argv[])
{

    double ZERO = 1e-4;
    long idum = -1;
    double dt = 0.0005;
    double tStop = atof(argv[1]);
    double maxPosition = 1;
    double l0 = sqrt(2*dt); // Defining D = 1

    int nWalkers = 10000;

    vector<double> walkers(nWalkers, 0);

    for (double T = 0; T < tStop; T += dt) {
        int N = walkers.size();
        for (int i = 0; i < N; i++) {

            double dx = l0*gaussian_deviate(&idum);
            double oldPosition = walkers[i];
            double newPosition = walkers[i] + dx;

            if (newPosition > 0) {
                walkers[i] = newPosition;
            }
            if (newPosition > 0 && oldPosition == 0) {
                walkers.push_back(0);
            }
            if (newPosition >= maxPosition) {
                walkers.erase(walkers.begin()+i);
                if (i != N-1) {
                    i--;
                    N--;
                }
            }

            if (newPosition > 0 && newPosition < ZERO && oldPosition > ZERO) {
                walkers.erase(walkers.begin()+i);
                if (i != N-1) {
                    i--;
                    N--;
                }
            }
        }

       cout << "Time = " << T << endl;
    }

    ofstream outfile;
    string filename = argv[2];
    outfile.open(filename, ios::out);
    outfile << "# " << l0 << endl;
    for (int i = 0; i < walkers.size(); i++) {
        outfile << setw(15) << setprecision(8) << walkers[i] << endl;
    }
    outfile.close();





    return 0;
}



// random numbers with gaussian distribution
double gaussian_deviate(long* idum)
{
    static int iset = 0;
    static double gset;
    double fac, rsq, v1, v2;
    if ( idum < 0) iset =0;
    if (iset == 0) {
        do {
            v1 = 2.*ran2(idum) -1.0;
            v2 = 2.*ran2(idum) -1.0;
            rsq = v1*v1+v2*v2;
        } while (rsq >= 1.0 || rsq == 0.);
        fac = sqrt(-2.*log(rsq)/rsq);
        gset = v1*fac;
        iset = 1;
        return v2*fac;
    } else {
        iset =0;
        return gset;
    }
} // end function for gaussian deviate

/*
** The function
**         ran2()
** is a long periode (> 2 x 10^18) random number generator of
** L’Ecuyer and Bays-Durham shuffle and added safeguards.
** Call with idum a negative integer to initialize; thereafter,
** do not alter idum between sucessive deviates in a
** sequence. RNMX should approximate the largest floating point value
** that is less than 1.
** The function returns a uniform deviate between 0.0 and 1.0
** (exclusive of end-point values).
*/
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran2(long *idum)
{
    int             j;
    long            k;
    static long     idum2 = 123456789;
    static long     iy = 0;
    static long     iv[NTAB];
    double          temp;
    if(*idum <= 0) {
        if(-(*idum) < 1) *idum = 1;
        else             *idum = -(*idum);
        idum2 = (*idum);
        for(j = NTAB + 7; j >= 0; j--) {
            k     = (*idum)/IQ1;
            *idum = IA1*(*idum - k*IQ1) - k*IR1;
            if(*idum < 0) *idum +=  IM1;
            if(j < NTAB)  iv[j]  = *idum;
        }
        iy=iv[0];
    }
    k   = (*idum)/IQ1;
    *idum = IA1*(*idum - k*IQ1) - k*IR1;
    if(*idum < 0) *idum += IM1;
    k     = idum2/IQ2;
    idum2 = IA2*(idum2 - k*IQ2) - k*IR2;
    if(idum2 < 0) idum2 += IM2;
    j     = iy/NDIV;
    iy    = iv[j] - idum2;
    iv[j] = *idum;
    if(iy < 1) iy += IMM1;
    if((temp = AM*iy) > RNMX) return RNMX;
    else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
// End: function ran2()
