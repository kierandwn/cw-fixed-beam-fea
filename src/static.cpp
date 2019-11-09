/*
 - AE3-422 High Performance Computing
 - Author: Kieran Downie (00964704)
 
 - Solves the Static Problem in Serial.
 
 */

#include <iostream>
#include <fstream>

using namespace std;

#include "static.h"
#include "assembly.h"

#define F77NAME(x) x##_
extern "C" {
    void F77NAME(dpbsv) (const char& UPLO, const int& N, const int& KD, const int& NRHS, double *AB,
                         const int& ldab, const double *B, const int& LDB, int& INFO);
}


int staticSolve (int dof, double E, double I, double A, double l, double Fy, int pos, double qy, double qx) {
    
    // Declare (reduced) K matrix (symmetric banded) and force vector in memory using (reduced) number of degrees of freedom.
    double K [6 * dof];
    double F [dof];
    
    assembleStiffnessMatrix(K, dof, E, I, A, l);
    assembleForceVector(F, dof, l, pos, Fy, qy, qx);

    int INFO;
    F77NAME(dpbsv) ('U', dof, 5, 1, K, 6, F, dof, INFO);

    if (INFO) {
        cout << "An error occurred during the solution of the linear system: " << INFO << "." << endl;
        return 1;
    }

    ofstream fOutput;
    fOutput.open("./output/static.txt");

    // Perform check if file is writable.
    if (!fOutput.good()) {
        cout << "Error: Unable to write to output file: output.txt" << endl;
        return 2;
    }

    double x = l;

    for (int i = 0; i < dof; i += 3) {
        
        fOutput.precision(3);
        fOutput.width(4);
        
        fOutput << fixed << x << " ";
        
        fOutput.precision(10);
        fOutput.width(11);
        fOutput << fixed << F[i] << " " << F[i + 1] << " " << F[i + 2] << "\n";
        
        x += l;
    }

    fOutput.close();

    cout << "Static solve: Complete" << endl;
    return 0;

}
