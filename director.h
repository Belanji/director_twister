
struct lc_cell
{

  double k11,k22,k33;
  double n0, ne;
  double viscosity, surf_viscosity[2];
  double pretwist[2], pretilt[2];
  double wa[2], omega_d[2];
  double cell_length;
  
};

int frank_energy (double t, const double phi[], double f[], struct lc_cell *params);

int jacobian(double t, const double phi[], double * dfdphi, double dfdt[], void * params);


int print_phi_time( const double *,
		    const double  ,
		    const double  );

int print_snapshot_to_file(const double *,
			   const double  ,
			   const double  ,
			         FILE   *);
