/*
 - AE3-422 High Performance Computing
 - Author: Kieran Downie (00964704)
 
 - Time Integration Schemes in Serial.
 
 */

#include <iostream>
#include <fstream>

#include "integration.h"
#include "assembly.h"

using namespace std;

#define F77NAME(x) x##_
extern "C" {
    void F77NAME(dpbtrf) (const char& UPLO, const int& N, const int& KD, double* AB, const int& ldab, int& INFO);
    void F77NAME(dpbtrs) (const char& UPLO, const int& N, const int& KD, const int& nrhs, const double* AB,
                          const int& lda, double* B, const int& ldb, int& INFO);
    void F77NAME(daxpy) (const int& N, const double& alpha, double* A, const int& INCX, double* B, const int& INCY);
    void F77NAME(dsbmv) (const char& UPLO, const int& N, const int& K, const double& alpha, double* A,
                         const int& lda, double* x, const int& INCX, const double& beta, double* y, const int& INCY);
    void F77NAME(dcopy) (const int& N, double* X, const int& INCX, double* Y, const int& INCY);
}

int timeIntegrationExplicit (double* M, double* K, double L, double l, double Fy, int pos, double qy, double qx, const int dof, const double dt, const double Tl, const double maxTime) {
    
    double t = 0;
    
    // Intialise the current and past displacement vectors as 0.
    double   u [dof] = {0};
    double  _u [dof] = {0};
    
    double  f [dof];
    
    // Define index of the tranverse degree of freedom of the central node.
    int index = dof / 2;
    
    for (int i = 0; i < dof; ++i) {
        K[6 * i + 5] -= 2 * M[i];
    }
    
    ofstream fOutput1;
    fOutput1.open("./output/dynamic.txt");
    
    if (!fOutput1.good()) {
        cout << "Error: Unable to write to output file: dynamic.txt" << endl;
        return 2;
    }
    
    double p;
    
    while (t <= maxTime) {
        
        if (t < Tl) {
            p = t / Tl;
            assembleForceVector(f, dof, l, pos, p * Fy, p * qy, qx);
        } else {
            assembleForceVector(f, dof, l, pos, Fy, qy, qx);
        }
        
        F77NAME(dsbmv) ('U', dof, 5, -1, K, 6,  u, 1, 1, f, 1);
        F77NAME(dsbmv) ('U', dof, 0, -1, M, 1, _u, 1, 1, f, 1);
        
        // Copy the current displacement vector to the past.
        F77NAME(dcopy) (dof, u, 1, _u, 1);
        
        // Multiply the RHS of eq. 9 in handout by the inverse of the mass matrix.
        for (int i = 0; i < dof; ++i) {
            u[i] = f[i] / M[i];
        }
        
        fOutput1.precision(5);
        fOutput1.width(6);
        fOutput1 << fixed << t << " ";
        
        fOutput1.precision(10);
        fOutput1.width(11);
        fOutput1 << u[index] << "\n";
        
        t += dt;
    }
    
    fOutput1.close();
    
    return 0;
    
}

int timeIntegrationImplicit (double* M, double* K, double L, double l, double Fy, int pos, double qy, double qx, const int dof, const double dt, const double Tl, const double maxTime) {
    
    double beta = 0.25;
    double gamma = 0.5;
    
    // Declare displacement, velocity and acceleeration terms, and pre-alloacate as 0.
    double u    [dof] = {0};
    double u_d  [dof] = {0};
    double u_dd [dof] = {0};
    
    // Declare past acceleration term and pre-allocate 0.
    double _u_dd [dof] = {0};
    
    // Define certain coefficients for Newark scheme.
    double bdtdt = 1.0 / (beta * dt * dt);
    double bdt   = 1.0 / (beta * dt);
    double b21   = (1.0 / (2 * beta) - 1);
    
    double dt1g = dt * (1.0 - gamma);
    double dtg  = dt * gamma;
    
    // Declare intermediate variable to hold the parenthetical expression in eq. 12 (handout).
    double lambda [dof] = {0};
    
    // Declare time-varying force vector, and it's scaling factor.
    double f [dof];
    double p;
    
    // Define the position of the central node.
    int index = dof / 2;
    
    // Re-define K as the effective stiffness matrix.
    for (int i = 0; i < dof; ++i) {
        K[6 * i + 5] += bdtdt * M[i];
    }
    
    int INFO;
    
    // Pre-factorise the effective K matrix
    F77NAME(dpbtrf) ('U', dof, 5, K, 6, INFO);
    
    if (INFO) {
        cout << "Error: Unable to pre-factor for solution of the linear system: " << INFO << "." << endl;
        return 1;
    }
    
    // Open text file to output results.
    ofstream fOutput1;
    fOutput1.open("./output/dynamic.txt");
    
    if (!fOutput1.good()) {
        cout << "Error: Unable to write to the output file: dynamic.txt" << endl;
        return 2;
    }
    
    
    double t = 0;
    
    while (t < maxTime) {
        
        // Define the future force vector.
        if (t < Tl) {
            p = (t + dt) / Tl;
            assembleForceVector(f, dof, l, pos, p * Fy, p * qy, qx);
        } else {
            assembleForceVector(f, dof, l, pos, Fy, qy, qx);
        }
        
        // Form the rhs of equation 12 (handout), stored in f.
        for (int i = 0; i < dof; ++i) {
            f[i] += M[i] * lambda[i];
        }
        
        // Solve equation 12, using pre-factored matrix in K, storing the solution (u) in f:
        F77NAME(dpbtrs) ('U', dof, 5, 1, K, 6, f, dof, INFO);
        
        if (INFO) {
            cout << "Error: Unable to compute the solution of the linear system: " << INFO << "." << endl;
            return 1;
        }
        
        // Define the future acceleration and velocity terms.
        for (int i = 0; i < dof; ++i) {
            u_dd[i] = bdtdt * (f[i] - u[i])  - bdt * u_d[i] - b21 * _u_dd[i];
            u_d [i] = u_d[i] + dt1g * _u_dd[i] + dtg * u_dd[i];
        }
        
        // Copy the future values into the past.
        F77NAME(dcopy) (dof, u_dd, 1, _u_dd, 1);
        F77NAME(dcopy) (dof, f, 1, u, 1);
        
        // Solve for the parenthetical expression stored in lambda.
        for (int i = 0; i < dof; ++i) {
            lambda[i] = bdtdt * u[i] + bdt * u_d[i] + b21 * _u_dd[i];
        }
        
        // Write to output file.
        fOutput1.precision(5);
        fOutput1.width(6);
        fOutput1 << fixed << t << " ";
        
        fOutput1.precision(10);
        fOutput1.width(11);
        fOutput1 << u[index] << "\n";
        
        t += dt;
    }
    
    fOutput1.close();
    
    return 0;
}
