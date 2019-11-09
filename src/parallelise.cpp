/*
 - AE3-422 High Performance Computing
 - Author: Kieran Downie (00964704)
 
 - Time Integration Schemes in Parallel.
 
 */

#include <mpi.h>
#include <fstream>

#include "parallelise.h"
#include "assembly.h"

using namespace std;

#define F77NAME(x) x##_
extern "C" {
    void F77NAME(pdpbtrf) (const char& UPLO, const int& N, const int& K, double* a, int& ja, int* desca,
                          double* af, int& ldaf, double* work, int& lwork, int& INFO);
    void F77NAME(pdpbtrs)(const char& UPLO, const int& N, int& K, int& NRHS, double* a, int& ja, int* desca, double* b,
                          int& ib, int* descb, double* af, int& ldaf, double* work, int& lwork, int& INFO);

    void F77NAME(dsbmv) (const char& UPLO, const int& N, const int& K, const double& alpha, double* A,
                         const int& lda, double* x, const int& INCX, const double& beta, double* y,
                         const int& INCY);
    void F77NAME(dgemv) (const char& TRANS, const int& M, const int& N, const double& alpha, double* A, const int& lda,
                         double* X, const int& INCX, const double& beta, double* Y, const int& INCY);
    void F77NAME(dcopy) (const int& N, double* X, const int& INCX, double* Y, const int& INCY);
    
    void Cblacs_pinfo       (int&, int&);
    void Cblacs_get         (int, int, int&);
    void Cblacs_gridinit    (int&, char&, int, int);
    void Cblacs_gridinfo    (int, int&, int&, int&, int&);
    void Cblacs_gridexit    (int);
    void Cblacs_exit        (int);
    
}


int parallelisedExplicitIntegration (double E, double A, double I, double rho, double L, double l, double Fy, int pos, double qy, double qx, const int Nx, const double dt, const double Tl, const double maxTime) {
    
    int err = 0;
    // Extend for parallelisation over for more processes?
    
    // Initialise the MPI environment.
    err = MPI_Init(NULL, NULL);
    
    if (err != MPI_SUCCESS) {
        cout << "Failed to initialise MPI environment." << endl;
        return 3;
    }
    
    // Number of processes.
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    
    // Rank of local process.
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
    ofstream fOutput2;
    
    // Position of central node in each process.
    int index;
    
    // Define number of degrees of freedom in each partition.
    const int dof = ((Nx / world_size) + (1 - world_rank)) * 3;
    
    if (world_rank == 0) {
    
        fOutput2.open("./output/dynamic.txt");
        
        if (!fOutput2.good()) {
            cout << "Error: Unable to write to output file: dynamic.txt" << endl;
            return 2;
        }
        
        index = dof - 5;
        
    } else if (world_rank == 1) {
        
        index = 1;
        
    }
    
    double K [6 * dof];
    double M [dof];
    double f [dof];
    
    double  u [dof] = {0};
    double _u [dof] = {0};
    
    double multiplier = rho * A * l / (dt * dt);
    
    assembleStiffnessMatrix(K, dof, E, I, A, l);
    assembleMassMatrix(M, dof, l, multiplier);
    
    for (int i = 0; i < dof; ++i) {
        K[6 * i + 5] -= 2 * M[i];
    }
    
    double t = 0;
    double p;
    
    while (t <= maxTime) {
        
        if (t < Tl) {
            p = t / Tl;
            assembleForceVector(f, dof, l, pos, p * Fy, p * qy, qx);
        } else {
            assembleForceVector(f, dof, l, pos, Fy, qy, qx);
        }

        F77NAME(dsbmv) ('U', dof, 5, -1, K, 6, u, 1, 1, f, 1);
        
        if (world_rank == 0) {
            
            MPI_Ssend(&f[dof - 6], 3, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
            MPI_Recv(&f[dof - 3], 3, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
        } else if (world_rank == 1) {
            
            MPI_Recv(&f[0], 3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Ssend(&f[3], 3, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
            
        }
        
        F77NAME(dsbmv) ('U', dof, 0, -1, M, 1, _u, 1, 1, f, 1);
        F77NAME(dcopy) (dof, u, 1, _u, 1);
        
        for (int i = 0; i < dof; ++i) {
            u[i] = f[i] / M[i];
        }
        
        if (world_rank == 0) {
            
            fOutput2.precision(5);
            fOutput2.width(6);
            fOutput2 << fixed << t << " ";
            
            fOutput2.precision(10);
            fOutput2.width(11);
            fOutput2 << u[index] << "\n";

        }
        
        t += dt;
    }
    
    if (world_rank == 0) {
        fOutput2.close();
    }
    
    MPI_Finalize();
    
    return 0;
}

int parallelisedImplicitIntegration (double E, double A, double I, double rho, double L, double l, double Fy, int pos, double qy, double qx, const int Nx, const double dt, const double Tl, const double maxTime) {
    
    // Global number of degrees of freedom
    const int DOF = (Nx - 1) * 3;
    
    // Explicit call to MPI_Init required.
    int err = MPI_Init(NULL, NULL);
    
    if (err != MPI_SUCCESS) {
        cout << "Failed to initialise MPI environment." << endl;
        //return -1;
    }
    
    // Initialise BLACS environment.
    int mype, npe, ctx, nrow, ncol, myrow, mycol;
    char order = 'C';
    
    // Initialises the BLACS world communicator.
    Cblacs_pinfo(mype, npe);
    
    // Get the default system context.
    Cblacs_get(0, 0, ctx);
    
    // Initialise a process grid of 1 rows and npe columns
    Cblacs_gridinit(ctx, order, 1, npe);
    
    // Get info about the grid to verify it is set up
    Cblacs_gridinfo(ctx, nrow, ncol, myrow, mycol);

    // Local number of degrees of freedom.
    int dof = 3 * (Nx / npe);
    
    // End process has one less node to satisfy boundary conditions.
    if (mype == npe - 1) {
        dof -= 3;
    }
    
    // Define variables associated with Newmark Scheme.
    double beta = 0.25;
    double gamma = 0.5;
    
    double K [6 * dof];
    double M [dof];
    
    double multiplier = rho * A * l;
    
    // Assemble local stiffness & mass matrices.
    assembleStiffnessMatrix(K, dof, E, I, A, l);
    assembleMassMatrix(M, dof, l, multiplier);
    
    // Declare displacement, velocity and acceleeration terms, and pre-alloacate as 0.
    double u    [dof] = {0};
    double u_d  [dof] = {0};
    double u_dd [dof] = {0};
    
    // Declare past acceleration term and pre-allocate 0.
    double _u_dd [dof] = {0};
    
    // Define certain coefficients for Newark scheme.
    double bdtdt = 1.0 / (beta * dt * dt);
    double bdt   = 1.0 / (beta * dt);
    double b21   = 1.0 / (2 * beta) - 1;
    
    double dt1g = dt * (1.0 - gamma);
    double dtg  = dt * gamma;
    
    // Declare intermediate variable to hold the parenthetical expression in eq. 12 (handout).
    double lambda [dof] = {0};
    
    // Declare time-varying force vector, and it's scaling factor.
    double f [dof];
    double p;
    
    // Re-define K as the effective stiffness matrix.
    for (int i = 0; i < dof; ++i) {
        K[6 * i + 5] += bdtdt * M[i];
    }
    
    // Define the descriptor for implementing SCALAPACK routine.
    char UPLO = 'U';
    int     n = DOF;                        // Number of columns in the K matrix to be pre-factored.
    int  NB_A = 3 * (Nx / npe);
    int    bw =   5;                        // Half- bandwidth of the K matrix (number of superdiagonals).
    int   laf = (NB_A + (2 * bw)) * bw;     // length of space reserved for use by the function.
    int   col = 1; //
    int lwork = bw * bw; // = 1000; Arbitrarily large workspace for process above np = 4
    int  nRHS = 1;
    int   row = 1;
    
    int desc_K [7];
    
    desc_K [0] =  501;
    desc_K [1] =  ctx;
    desc_K [2] =  DOF;
    desc_K [3] = NB_A;
    desc_K [4] =    0;
    desc_K [5] =    6;
    desc_K [6] =    0;
    
    int desc_RHS [7];
    
    desc_RHS [0] =  502;
    desc_RHS [1] =  ctx;
    desc_RHS [2] =  DOF;
    desc_RHS [3] = NB_A;
    desc_RHS [4] =    0;
    desc_RHS [5] = NB_A;
    desc_RHS [6] =    0;
    
    int INFO;
    
    double af   [laf]  ;
    double work [lwork];
    
    // Pre-factorise the effective K matrix.
    F77NAME(pdpbtrf) (UPLO, n, bw, K, col, desc_K, af, laf, work, lwork, INFO);
    
    if (INFO) {
        cout << "Error: Unable to pre-factor for solution of the linear system: " << INFO << "." << endl;
        return 1;
    }
    
    ofstream fOutput3;
    
    if (pos >= mype * (Nx * 3 / npe) && pos <= (mype + 1) * (Nx * 3 / npe)) {
        
        // Convert position from global to local co-ordinates in the sub-domain.
        pos -= mype * (Nx * 3 / npe);
        
        // Open text file to output results.
        fOutput3.open("./output/dynamic.txt");
        
        if (!fOutput3.good()) {
            cout << "Error: Unable to write to the output file: dynamic.txt" << endl;
            return 2;
        }

    } else {
        Fy  = 0;
        pos = 0;
    }
    
    double t = 0;
    
    while (t <= maxTime) {
        
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
        F77NAME(pdpbtrs) (UPLO, n, bw, nRHS, K, col, desc_K, f, row, desc_RHS, af, laf, work, lwork, INFO);
        
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
        
        if (fOutput3.good()) {
            
            // Write to output file.
            fOutput3.precision(5);
            fOutput3.width(6);
            fOutput3 << fixed << t << " ";
            
            fOutput3.precision(10);
            fOutput3.width(11);
            fOutput3 << u[pos] << "\n";

        }
        
        
        t += dt;
    }
    
    // Finalise BLACS grid context.
    Cblacs_gridexit(ctx);
    
    // Finalise MPI.
    Cblacs_exit(0);
    
    if (fOutput3.good()) { fOutput3.close(); }
    
    return 0;
}
