/*
 - AE3-422 High Performance Computing
 - Author: Kieran Downie (00964704)
 
 - Functions defined to assemble system matrices.
 
 */


using namespace std;

#include "assembly.h"

void assembleStiffnessMatrix(double *K, const int _dof, double E, double I, double A, double l) {
    
    // Define material properties as the product of E, I & A to save repeating computation later.
    const double EI = E * I;
    const double EA = E * A;
    
    // Split the elemental K matrix into 4 sub-matrices, and store these to easily assign values in global K matrix below.
    
    //      [[1]|[2]]       [4] = [1] + [3]
    //      ---------
    //      [[2]|[3]]
    
    // Note that only [2] and [4] are required to assemble the reduced stiffness matrix.
    
    // Submatrix 2, stored column major
    double vMatrix2[9] = {-1 * A * E / l, 0, 0, 0, -12 * E * I / (l * l * l), -6 * E * I / (l * l), 0, 6 * E * I / (l * l), 2 * E * I / l};
    // Submatrix 4 stored as the sum of submatrices 1 & 3
    double vMatrix4[6] = {2 * EA / l, 0, 24 * EI / (l * l * l), 0, 0, 8 * EI / l};
    
    // Declare counters for each sub-matrix, pre-allocate count2 & count22 to prevent extra conditional statements later;
    int count1, count2;
    count1  = 0;
    count2  = 9;
    
    // Global K matrix is assembled by first identifying regions of the global matrix which are filled by each sub-matrix, then populating them with an independent counter for each such that the matrix is populated uniformly. Column-major for blas compatibility.
    for (int i = 0; i < _dof; ++i) {
        for (int j = 0; j < 6; ++j) {
            if (j >= 5 - (i % 3)) {
                // Submatrix 4
                if (i % 3 == 0) {count1 = 0;}
                K[i * 6 + j] = vMatrix4[count1];
                ++count1;
            } else if (j >= 2 - (i % 3)) {
                // Submatrix 2
                if (i % 3 == 0 && count2 == 9) {count2 = 0;}
                K[i * 6 + j] = vMatrix2[count2];
                ++count2;
            } else {
                K[i * 6 + j] = 0;
            }
        }
    }

}

void assembleForceVector(double* F, const int dof, double l, int position, double Fy, double qy, double qx) {

    for (int i = 0; i < dof; ++i) {
        if (i % 3 == 0) {
            F[i] = qx * l;
        } else if (i % 3 == 1) {
            F[i] = qy * l;
        } else {
            F[i] = 0;
        }
    }

    F[position] += Fy;
    
}

void assembleMassMatrix (double* M, const int dof, double l, double multiplier) {
    
    double alpha = 1.0 / 24.0;
    
    // Mass matrix is assembled as a single line for efficiency of storage.
    for (int i = 0; i < dof; ++i) {
        if (i % 3 < 2) {
            M[i] = multiplier;
        } else if (i % 3 == 2) {
            M[i] = 2 * multiplier * alpha * l * l;
        }
    }
}
