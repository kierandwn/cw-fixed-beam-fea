void assembleStiffnessMatrix(double *K, const int dof, double E, double I, double A, double l);
void assembleForceVector(double* F, const int dof, double l, int position, double Fy, double qy, double qx);
void assembleMassMatrix (double* M, const int dof, double l, double multiplier);
