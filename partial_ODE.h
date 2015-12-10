#include <iostream>
using namespace std;



void forwardEuler(double* u, double* unew, int N_x, double dx, double t_stop);


void outputToFile(string filename, double* u, double dx, int N_x);

void backwardEuler(double* u, double* unew, int N_x, double dx, double dt, double t_stop);

void Crank_Nicolson(double* u, double* unew, int N_x, double dx, double dt, double t_stop);
