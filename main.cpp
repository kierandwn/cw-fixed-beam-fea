/*
 - AE3-422 High Performance Computing
 - Author: Kieran Downie (00964704)
 
 - Main calling function. Handles all input.
 
 */

#include <iostream>
#include <cstdlib>

using namespace std;

#include <boost/program_options.hpp>

// Alias for ease of typing.
namespace po = boost::program_options;

// Contains declarations for the functions defined in other source files.
#include "assembly.h"
#include "static.h"
#include "integration.h"
#include "parallelise.h"

int main (int argc, char** argv) {
    
    // Parse inputs using boost, if invalid these are assigned to be those given in the handout.
    po::options_description desc("This program computes a Finite Element Analysis of a beam fixed at both ends. \n Options: ");
    desc.add_options()
    // Problem geometry.
    ("length"  , po::value<double>(), "Length of Beam. (>= 0)")
    ("force"   , po::value<double>(), "Position of Concentrated Force along the beam. (0 <= x <= 1)")
    // Material properties.
    ("area"    , po::value<double>(), "Area of Beam Cross-section. (> 0)")
    ("I"       , po::value<double>(), "Beam Cross-sectional Second Moment of Area. (> 0)")
    ("E"       , po::value<double>(), "Material Modulus of Elasticity. (> 0)")
    ("rho"     , po::value<double>(), "Material Density. (> 0)")
    // Mesh definition.
    ("Nx"      , po::value<int>()   , "Number of Elements. (> 0)")
    ("Nt"      , po::value<int>()   , "Number of Time Increments. (> 0)")
    // Time domain.
    ("T"       , po::value<double>(), "Maximmum Time in Simulation. (>= 0)")
    ("Tl"      , po::value<double>(), "Maximum Time in Transient Force Phase. (>= 0)")
    // Solver options.
    ("static"  ,                      "Solve Static Equation only. (default)")
    ("dynamic" ,                      "Solve Dynamic Eqaution using specified integration scheme.")
    ("explicit",                      "Solve Dynamic Eqaution using explicit time integration scheme.")
    ("implicit",                      "Solve Dynamic Eqaution using implicit time integration scheme.")
    ("parallel",                      "Solve the specified equation using the specified integration scheme in parallel.")
    
    ("help"    ,                      "Print help message.");
    
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    
    if (vm.count("help")) {
        cout << desc << endl;
        return 1;
    }

    // Static Inputs.
    const double L   = vm.count("length") ? vm["length"].as<double>() :    10000.0   ;
    const double pos = vm.count("force" ) ? vm["force" ].as<double>() :        0.5   ;
    const int    Nx  = vm.count("Nx"    ) ? vm["Nx"    ].as<int>()    :       24     ;
    const double A   = vm.count("area"  ) ? vm["area"  ].as<double>() :    12000.0   ;
    const double I   = vm.count("I"     ) ? vm["I"     ].as<double>() :     14.4e+06 ;
    const double E   = vm.count("E"     ) ? vm["E"     ].as<double>() :      210e+03 ; // Pa
    
    // Dynamic Inputs.
    const double rho = vm.count("rho"   ) ? vm["rho"   ].as<double>() : 7850e-12      ;
    const double T   = vm.count("T"     ) ? vm["T"     ].as<double>() :    10.0       ;
    const double Tl  = vm.count("Tl"    ) ? vm["Tl"    ].as<double>() :     1.0       ;
    const int    Nt  = vm.count("Nt"    ) ? vm["Nt"    ].as<int>()    : 10000         ;
    
    const double l = L / Nx;
    
    if (Nx % 2 != 0) {
        cout << "There is not a node under the central concentrated force: solution will not be exact." << endl;
    }
    
    // Boost library will throw an error if inputs are not of correct data type, further validation required to make sure inputs are physical.
    if (L <= 0 || Nx <= 0 || A <= 0 || I <= 0 || E <= 0 || T < 0 || Tl < 0 || Nt <= 0 || rho <= 0) {
        cout << "Invalid input, see option constraints below." << endl;
        cout << desc << endl;
    }

    // Define loading conditions.
    const double qx =       0; // N / mm (Trivial for this exercise, but included to keep the solution generalised.)
    const double qy = -   1.0; // N / mm
    const double Fy = -1000.0; // N
    
    // Reduced degrees of freedom: 3 dof per node, Nx - 1 nodes after removal of boundary nodes.
    const int dof = (Nx - 1) * 3;
    
    // Global position of force in dof domain (rounded down to integer value).
    int POS = 3.0 * (pos * Nx) - 2;
    
    int err;

    if (vm.count("dynamic")) {
        
        const double dt = Tl / Nt;
        
        double M[dof];
        
        if (vm.count("explicit")) {
            
            if (vm.count("parallel")) {
                
                // Check inputs properly.
                err = parallelisedExplicitIntegration(E, A, I, rho, L, l, Fy, POS, qy, qx, Nx, dt, Tl, T);
                
                if (err) {
                    return err;
                }
                
            } else {
                
                // Include the dt^2 denominator in the mass matrix since is a common factor in each mass term of the explicit solver.
                double multiplier = rho * A * l / (dt * dt);
                double K [6 * dof];
                
                assembleMassMatrix(M, dof, l, multiplier);
                assembleStiffnessMatrix(K, dof, E, I, A, l);
                
                err = timeIntegrationExplicit(M, K, L, l, Fy, POS, qy, qx, dof, dt, Tl, T);
                
                if (err) {
                    return err;
                }
            }
            
            
        } else if (vm.count("implicit")) {
            
            const double dt = Tl / Nt;
            
            if (vm.count("parallel")) {
                
                err = parallelisedImplicitIntegration(E, A, I, rho, L, l, Fy, POS, qy, qx, Nx, dt, Tl, T);
                
                if (err) {
                    return err;
                }
                
            } else {
                
                double multiplier = rho * A * l;
                double K [6 * dof];
                
                assembleMassMatrix(M, dof, l, multiplier);
                assembleStiffnessMatrix(K, dof, E, I, A, l);

                err = timeIntegrationImplicit(M, K, L, l, Fy, POS, qy, qx, dof, dt, Tl, T);
                
                if (err) {
                    return err;
                }
            }
            
        }
        
    }  else {
        
        // Run the Static solver as default if no scheme identified.
        err = staticSolve (dof, E, I, A, l, Fy, POS, qy, qx);
        
        if (err) {
            return err;
        }
        
    }

    
    return 0;
}
