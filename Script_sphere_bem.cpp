// Standard C libraries
#include <cmath> // for functions like pow, csqrt, sqrt, complex math, etc. from C library
#include <cstdlib> // for functions like atoi etc. from C library
#include <ctime> // for time from C library
#include <cstdio> // for functions like printf, sprintf, and file I/O, etc from C library
#include <cstring> // for functions like strcat, strcpy, etc from C library
#include <cassert> // for assert function

// Standard C++ libraries
#include <iostream> // for cin, cout
#include <fstream> // for file  I/O
#include <sstream> // for string I/O ( good for num2str conversion)
#include <complex>
#include <string>
#include <iomanip>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>


// External C++ libraries
#include <libscuff.h>


using namespace scuff;
using namespace std;
typedef complex<double> cdouble;

#define pi 3.14159265358979323
#define MASTER 0


// Global variables

// *********
// SI units
// ********


// Constants
extern const double epsilon_0 = 8.85e-12;
extern const double mu_0 = 4.0*pi*1.0e-7;
extern const double c = 3e8;
extern const double hbar = 1.05*1e-34 ;
extern const double me = 9.109*1e-31 ;
extern const double kb = 1.3806503e-23 ;
extern const double ec = 1.6e-19 ; 
extern const double aB = 5.3e-11 ;  
extern const double a = 1e-6;


// Source parameters
extern const double lambda_SI = 2.0*pi*a;
extern const double freq_SI = c/lambda_SI;
extern const double omega_SI = 2.0*pi*freq_SI;


// Material parameters
extern const double xc_SI = 0.0; 
extern const double yc_SI = 0.0;
extern const double zc_SI = 0.0;
extern const double rad_SI = 10.0*a;

// Parameters (for Aluminum)
extern const double nu_SI = 7.596e13 ; // in rad/s  // For Copper: 4*pi*c*1e4 ; 
extern const double freq_plasma_SI = 1.747e16/(2*pi); //  in rad/s  // \omega_p = sqrt{\sigma/(\tau*\epsilon_0)}  // For Copper: 1.6/(2*pi)*2.42e14; 


// ************
// Scuffem units
// ************


// Source parameters
extern const double lambda = lambda_SI/a;
extern const double freq = freq_SI*a/c;
extern const double omega = omega_SI*a/c;


// Material parameters
extern const double xc = xc_SI/a; 
extern const double yc = yc_SI/a;
extern const double zc = zc_SI/a;
extern const double rad = rad_SI/a;


// Source/Observation location (angle in radians)
extern const double THETA_s = pi/4.0;
extern const double PHI_s = 0.0;
extern const double THETA_o = THETA_s;
extern const double PHI_o = PHI_s;




// Main
int main(int argc, char **argv) {


  // Start time
	time_t current_time = time(NULL) ;
	char* pcurrent_time_string ;
	
	// Convert to local time format and print the string
	pcurrent_time_string = ctime(&current_time) ;
	cout << "START TIME: " << pcurrent_time_string << endl ;

    
    
	  // Radial distances from the surface of the sphere
    int n;
    int dlen = 10; 
    double dstep = 1.0*a/a; 
    double dstart = 1.0*a/a; 

    double pDspan[dlen]; 
    for (n=0;n<dlen;n++){
        pDspan[n] = dstart +  n*dstep;
    }
    
    // Source and field locations 
    double sx, sy, sz, ox, oy, oz, R_s, R_o;
    R_s = rad + pDspan[0];
    sx = xc + R_s*cos(PHI_s)*sin(THETA_s);
    sy = yc + R_s*sin(PHI_s)*sin(THETA_s);
    sz = zc + R_s*cos(THETA_s);
    
    R_o = R_s;
    ox = xc + R_o*cos(PHI_o)*sin(THETA_o);
    oy = yc + R_o*sin(PHI_o)*sin(THETA_o);
    oz = zc + R_o*cos(THETA_o);

    double SLoc[3] = {sx, sy, sz};  // point source location 
    double OLoc[3] = {ox, oy, oz}; // observation location



    // Store output E-fields in a matrix
	HMatrix FxMatrix(dlen, 1, LHM_COMPLEX, LHM_NORMAL); //  complex, non-symmetric matrix
	HMatrix FyMatrix(dlen, 1, LHM_COMPLEX, LHM_NORMAL); //  complex, non-symmetric matrix
	HMatrix FzMatrix(dlen, 1, LHM_COMPLEX, LHM_NORMAL); //  complex, non-symmetric matrix


	HMatrix Fx0Matrix(dlen, 1, LHM_COMPLEX, LHM_NORMAL); //  complex, non-symmetric matrix
	HMatrix Fy0Matrix(dlen, 1, LHM_COMPLEX, LHM_NORMAL); //  complex, non-symmetric matrix
	HMatrix Fz0Matrix(dlen, 1, LHM_COMPLEX, LHM_NORMAL); //  complex, non-symmetric matrix
	
	
	// Store temporary E- and H-fields in array
	cdouble EH[6];
	cdouble E0H0[6];
	


    // Define geometry
	RWGGeometry *Geom=new RWGGeometry("Sphere.scuffgeo", a); // <--- I modified the source code so it takes lengthscale "a" as input
	
	
	SetLogFileName("Sphere.log");
	Geom->SetLogLevel(SCUFF_VERBOSELOGGING);

	// Preallocate BEM matrix and RHS vector
	HMatrix *M  = Geom->AllocateBEMMatrix();
	HVector *KN = Geom->AllocateRHSVector();
	

	// Point source
	cdouble SPol[3] = {0.0, 0.0, 1.0};  // point source polarization and strength
	PointSource *PS = new PointSource(SLoc, SPol, 0); // use 0 for electric dipole and 1 for magnetic dipole source


	// Material property
	MatProp* HSMP = new MatProp("ALUMINUM"); // sphere material property



	// Assemble BEM matrix
	Geom->AssembleBEMMatrix(static_cast<cdouble>(omega), M);
	

	// Factorize BEM matrix
	M->LUFactorize();
	


	// Loop through distances
	for (n=0;n<dlen;n++){

		// Assemble RHS vector
		Geom->AssembleRHSVector(static_cast<cdouble>(omega), PS, KN);
		

		// Solve BEM equation
		M->LUSolve(KN);
		
		
		// Obtain E- and H-fields. These fields will be used to compute dyadic greens function
		Geom->GetFields(0, KN, static_cast<cdouble>(omega), OLoc, E0H0);
		Geom->GetFields(PS, KN, static_cast<cdouble>(omega), OLoc, EH);
		// Using PS gives total field (scattered + incident). Replace PS by 0 to get scattered field only

		

		// Store E-fields in HMatrix (H-fields are not stored)
		Fx0Matrix.SetEntry(n, 1, E0H0[0]);
		Fy0Matrix.SetEntry(n, 1, E0H0[1]);
		Fz0Matrix.SetEntry(n, 1, E0H0[2]);
		FxMatrix.SetEntry(n, 1, EH[0]);
		FyMatrix.SetEntry(n, 1, EH[1]);
		FzMatrix.SetEntry(n, 1, EH[2]);


		// Update source/observation location
        R_s = rad + pDspan[n+1];
        sx = xc + R_s*cos(PHI_s)*sin(THETA_s);
        sy = yc + R_s*sin(PHI_s)*sin(THETA_s);
        sz = zc + R_s*cos(THETA_s);
    
        R_o = R_s;
        ox = xc + R_o*cos(PHI_o)*sin(THETA_o);
        oy = yc + R_o*sin(PHI_o)*sin(THETA_o);
        oz = zc + R_o*cos(THETA_o);

		SLoc[0] = sx;
		SLoc[1] = sy;
		SLoc[2] = sz;
		OLoc[0] = ox;
		OLoc[1] = oy;
		OLoc[2] = oz;

		// Update point source location
		PS->SetX0(SLoc);
	}


	// Delete mat prop
	delete HSMP;
	

  // End time
	current_time = time(NULL) ;
	
	// Convert to local time format and print the string
	pcurrent_time_string = ctime(&current_time) ;
  printf("END TIME: %s", pcurrent_time_string) ;


	return 0;

}



