#include <iostream>
#include "../partial_ODE.h"
using namespace std;

int main(int argc, char* argv[])
{

    // initialization
    double dx = stod(argv[1]);
    double t_stop = stod(argv[2]);    // stop computing after this much time
    int N_x = 1.0/dx + 1;   // Number of steps in the x-dimension
    double* v_initial = new double[N_x];
    double* v_fEuler  = new double[N_x];
    double* v_bEuler  = new double[N_x];
    double* v_CN      = new double[N_x];

    // inserting initial and boundary conditions:
    v_initial[0] = v_initial[N_x-1] = 0; // boundary
    for (int i = 1; i < N_x-1; i++) {
        v_initial[i] = dx*i - 1; // initial conditions
    }

    // The explicit, forward Euler, scheme
    forwardEuler(v_initial, v_fEuler, N_x, dx, t_stop);

    // The implicit, backward Euler, scheme
    double dt = 5e-9;
    backwardEuler(v_initial, v_bEuler, N_x, dx, dt, t_stop);

    // The Crank-Nickolson scheme
    Crank_Nicolson(v_initial, v_CN, N_x, dx, dt, t_stop);

    // remember, v = u - u_s -> u = v + u_s = v + 1 - x. Changing to this value
    for (int i = 0; i < N_x; i++) {
        v_fEuler[i] += 1 - dx*i;
        v_bEuler[i] += 1 - dx*i;
        v_CN[i]     += 1 - dx*i;
    }

    // print to files
    string filename;

    string extension = argv[3];
    filename = "ForwardEuler_"+extension+".txt";
    outputToFile(filename, v_fEuler, dx, N_x);

    filename = "BackwardEuler_"+extension+".txt";
    outputToFile(filename, v_bEuler, dx, N_x);

    filename = "CrankNickolson_"+extension+".txt";
    outputToFile(filename, v_CN, dx, N_x);


    delete[] v_initial;
    delete[] v_fEuler;
    delete[] v_bEuler;
    delete[] v_CN;
    return 0;
}











