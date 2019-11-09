int timeIntegrationExplicit (double* M, double* K, double L, double l, double Fy, int pos, double qy, double qx, const int _dof, const double dt, const double Tl, const double maxTime);
int timeIntegrationImplicit (double* M, double* K, double L, double l, double Fy, int pos,  double qy, double qx, const int dof, const double dt, const double Tl, const double maxTime);
