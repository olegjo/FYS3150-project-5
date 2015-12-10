#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;


void forward_substitution(double* &a, double* &b, double* &c, double* &b_tilde, int n);
void backward_substitution(double* &b, double* &c, double* &b_tilde, double* &v, int n);
void copy(double* uold, double *unew, int N);

void forwardEuler(double* u, double* unew, int N_x, double dx, double t_stop)
{
    double T = 0;
    // dont want to change the original u, because it needs to be used later as is. Therefore, make a copy.
    double* u_internal = new double[N_x];
    copy(u, u_internal, N_x);

    // boundary conditions
    unew[0] = u[0]; unew[N_x - 1] = u[N_x - 1];

    double alpha = 0.5; // alpha = dt/dx^2, dt = 0.5*dx^2 -> alpha = 0.
    double dt = 0.5*dx*dx;
    for (int j = 0; T <= t_stop; j++) {
        for (int i = 1; i < N_x - 1; i++) {
            unew[i] = alpha*(u_internal[i-1] + u_internal[i+1]) + (1 - 2*alpha)*u_internal[i];
        }

        copy(unew, u_internal, N_x);

        T += dt;
    }
    delete[] u_internal;
} // End: function forwardEuler


void backwardEuler(double* u, double* unew, int N_x, double dx, double dt, double t_stop)
{
    // dont want to change the original u, because it needs to be used later as is
    double* u_internal = new double[N_x];
    copy(u, u_internal, N_x);

    double* a = new double[N_x];
    double* b = new double[N_x];
    double* c = new double[N_x];
    double alpha = dt/dx/dx;
    double T = 0;


    // iterate over time
    // for each iteration, we wish to solve the set of linear equations
    // A*unew = u_internal
    // the input for next time iteration is then unew, so need to change u_internal = u_new
    for (int i = 0; T <= t_stop; i++) {
        for (int j = 0; j < N_x; j++) {
            a[j] = -alpha;
            b[j] = 1 + 2*alpha;
            c[j] = -alpha;
        }
        forward_substitution(a, b, c, u_internal, N_x);
        backward_substitution(b, c, u_internal, unew, N_x);
        // the boundary conditions may be changed in the tridiag solver. Need to change them back
        unew[0] = u[0];
        unew[N_x-1] = u[N_x-1];

        // update the known vector
        copy(unew, u_internal, N_x);
        T += dt;
    }

    delete[] u_internal;
    delete[] a;
    delete[] b;
    delete[] c;

} // End: function backwardEuler



void Crank_Nicolson(double* u, double* unew, int N_x, double dx, double dt, double t_stop)
{
    double* RHS = new double[N_x];
    double* a = new double[N_x];
    double* b = new double[N_x];
    double* c = new double[N_x];
    double* u_internal = new double[N_x];

    copy(u, u_internal, N_x);

    // boundary conditions:
    unew[0] = u[0]; unew[N_x - 1] = u[N_x - 1];
    RHS[0] = u[0]; RHS[N_x - 1] = u[N_x - 1];
    double alpha = dt/dx/dx;
    double T = 0;
    for (int i = 0; T <= t_stop; i++) {
        // calculating RHS
        for (int j = 1; j < N_x - 1; j++) {
            RHS[j] = alpha*(u_internal[j-1] + u_internal[j+1]) + (2 - 2*alpha)*u_internal[j];
        }



        // now calculate the tridiagonal linear eq.
        // first, set up the matrix
        for (int j = 0; j < N_x; j++) {
            a[j] = -alpha;
            b[j] = 2 + 2*alpha;
            c[j] = -alpha;
        }
        forward_substitution(a, b, c, RHS, N_x);
        backward_substitution(b, c, RHS, unew, N_x);
        // The boundary conditions are changed by forward/backword substitution. Need to change them back.
        unew[0] = u[0]; unew[N_x - 1] = u[N_x - 1];
        RHS[0] = u[0]; RHS[N_x - 1] = u[N_x - 1];

        copy(unew, u_internal, N_x);

        T += dt;
    }

    delete[] u_internal;
    delete[] RHS;
    delete[] a;
    delete[] b;
    delete[] c;

} // End: function Crank_Nickolson











void outputToFile(string filename, double *u, double dx, int N_x)
{
    ofstream outfile;
    outfile.open(filename, ios::out);
    for (int i = 0; i < N_x; i++) {
        outfile << setw(15) << setprecision(8) << i*dx;
        outfile << setw(15) << setprecision(8) << u[i] << endl;
    }
    outfile.close();
} // End: function otputToFile





void forward_substitution(double* &a, double* &b, double* &c, double* &b_tilde, int n){
    for (int i = 1; i < n; i++){
        b[i] = b[i-1]*b[i] - a[i]*c[i-1];
        c[i] = b[i-1]*c[i];
        b_tilde[i] = b[i-1]*b_tilde[i] - a[i]*b_tilde[i-1];
    }
} // End: function forward_substitution

void backward_substitution(double* &b, double* &c, double* &b_tilde, double* &v, int n){
    v[n-1] = b_tilde[n-1]/b[n-1];
    for (int i = n-2; i >= 0; i--){
        v[i] = (b_tilde[i] - c[i]*v[i+1])/b[i];
    }
} // End: function backward_substitution



void copy(double *uold, double* unew, int N)
{
    for (int i = 0; i < N; i++) {
        unew[i] = uold[i];
    }
} // End: function copy
