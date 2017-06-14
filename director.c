#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <math.h>
#include "./director.h"

const int nz=50;
const int time_relvant_digits=2;





void print_snapshot_with_time( const double * phi,
			       const double time)
			       
{
  FILE * snapshot_file;
  
  //open()

}


int main (int argc, char * argv[]) {

  double * theta, * phi;
  struct lc_cell lc_environment;
  double ti=0.0, tf=200.0;
  double time=ti;
  double timeprint=0.2;
  FILE * time_file;
  const char * time_file_name="middle_sin.dat";    
  
  lc_environment.k11=16.7;
  lc_environment.k22=7.0;
  lc_environment.k33=18.1;
  lc_environment.n0=1.4774;
  lc_environment.ne=1.5578;
  lc_environment.viscosity=186.0;
  lc_environment.cell_length=5.0;
  
  lc_environment.surf_viscosity[0]=0.1*186.0*5.0;
  lc_environment.surf_viscosity[1]=0.1*186.0*5.0;

  lc_environment.pretilt[0]=0.0;
  lc_environment.pretilt[1]=0.0;

  lc_environment.pretwist[0]=0.0;
  lc_environment.pretwist[1]=0.0;
  

  lc_environment.wa[0]=atof(argv[1]);
  lc_environment.wa[1]=atof(argv[2]);

  lc_environment.omega_d[0]=0.2;
  lc_environment.omega_d[1]=0.0;

  
  gsl_odeiv2_system sys = {frank_energy, jacobian, nz+1, &lc_environment};


  //gsl_odeiv2_driver * pde_driver =gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, 1e-9, 1e-6, 0.0);
  gsl_odeiv2_driver * pde_driver =gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_bsimp, 1e-6, 1e-6, 0.0);

  time_file=fopen(time_file_name,"w");
  
  phi= (double *) malloc( (nz+1)*sizeof(double) );

  for (int ii=0; ii<=nz;ii++)
    {
      phi[ii]=0.0;

    };

  fprintf(time_file,"%f  %f \n",time, sin(phi[nz/2]) );

  while(time <tf)
    {
      int status = gsl_odeiv2_driver_apply (pde_driver, &time, time+timeprint, phi);      
	
      if (status != GSL_SUCCESS)
	{
	  printf ("error, return value=%d\n", status);
   
	};

  
        fprintf(time_file,"%f  %f \n",time, sin(phi[nz/2]) );

    };
  

  gsl_odeiv2_driver_free (pde_driver);
  free(phi);
  fclose(time_file);
  return 0;


};
     

    
int frank_energy (double t, const double phi[], double f[], struct lc_cell *params)
{ 
  struct lc_cell mu = *(struct lc_cell *)params;
  double dz = mu.cell_length/nz;
  double k22= mu.k22;
  double omega_d[2],wa[2], pretwist[2];
  double viscosity, surf_viscosity[2];
  double dphi, d2phi;
  omega_d[0]=mu.omega_d[0];
  omega_d[1]=mu.omega_d[1];
  wa[0]=mu.wa[0];
  wa[1]=mu.wa[1];
  pretwist[0]=mu.pretwist[0];
  pretwist[1]=mu.pretwist[1];
  viscosity=mu.viscosity;
  surf_viscosity[0]=mu.surf_viscosity[0];
  surf_viscosity[1]=mu.surf_viscosity[1];
  

  /*Bulk equations */  
  for(int ii=1; ii<nz; ii++)
    {

      d2phi=(phi[ii+1]+phi[ii-1]-2.0*phi[ii])/(dz*dz);

      f[ii]= d2phi*k22/viscosity;
          

    };

  
  f[0]=(((phi[1]-phi[0])/dz)*k22 + wa[0]*cos( omega_d[0]*t - phi[0] )*sin( omega_d[0]*t - phi[0] ))/surf_viscosity[0];

  f[nz]=(-(((phi[nz]-phi[nz-1])/dz)*k22) + wa[1]*cos( omega_d[1]*t - phi[nz] )*sin( omega_d[1]*t - phi[nz] ))/surf_viscosity[1];
  
  

  //Boundary conditions n=1

  //f[0]=( (phi[1]-phi[0])*k22/dz + wa[0]*(t*omega_d[0]+pretwist[0] - phi[0]) )/surf_viscosity[0]; 
  //f[nz]=( -(phi[nz]-phi[nz-1])*k22/dz + wa[1]*(t*omega_d[1] + pretwist[1] - phi[nz]) )/surf_viscosity[1];
      

  return GSL_SUCCESS;

};


int jacobian(double t, const double phi[], double * dfdphi, double dfdt[], void * params)
{
  struct lc_cell mu = *(struct lc_cell *)params;
  double dz = mu.cell_length/nz;
  double k22= mu.k22;
  double omega_d[2],wa[2], pretwist[2];
  double viscosity, surf_viscosity[2];
  double dphi, d2phi;
  gsl_matrix_view dfdphi_mat= gsl_matrix_view_array (dfdphi, nz+1, nz+1);
  
  omega_d[0]=mu.omega_d[0];
  omega_d[1]=mu.omega_d[1];
  wa[0]=mu.wa[0];
  wa[1]=mu.wa[1];
  pretwist[0]=mu.pretwist[0];
  pretwist[1]=mu.pretwist[1];
  viscosity=mu.viscosity;
  surf_viscosity[0]=mu.surf_viscosity[0];
  surf_viscosity[1]=mu.surf_viscosity[1];
  
  gsl_matrix_set_zero( &dfdphi_mat.matrix );
  
  for(int i=0; i<=nz;i++)
    {

      dfdt[i]=0;
      
    };
  
  dfdt[0]=(omega_d[0]*wa[0]*pow(cos(pretwist[0] + omega_d[0]*t - phi[0]),2) - omega_d[0]*wa[0]*pow( sin(pretwist[0] + omega_d[0]* t - phi[0]), 2))/surf_viscosity[0];
  dfdt[nz]=(omega_d[1]*wa[1]*pow(cos(pretwist[1] + omega_d[1]*t - phi[nz]),2) - omega_d[1]*wa[1]*pow( sin(pretwist[1] + omega_d[1]* t - phi[nz]), 2))/surf_viscosity[1];

  
  for(int i=1;i<nz;i++){

    gsl_matrix_set ( &dfdphi_mat.matrix,i,i-1,k22/(dz*dz*viscosity) );
    gsl_matrix_set ( &dfdphi_mat.matrix,i,i  ,-2.0*k22/(dz*dz*viscosity));
    gsl_matrix_set ( &dfdphi_mat.matrix,i,i+1,k22/(dz*dz*viscosity) );

  };


  gsl_matrix_set ( &dfdphi_mat.matrix,0,0,( -(k22/dz) - wa[0]*pow(cos(pretwist[0] + omega_d[0]*t - phi[0]),2) + wa[0]*pow(sin(pretwist[0] + omega_d[0]*t - phi[0]),2) )/surf_viscosity[0] );
  gsl_matrix_set ( &dfdphi_mat.matrix,0,1,k22/(dz*surf_viscosity[0]) );

  gsl_matrix_set( &dfdphi_mat.matrix,nz,nz-1,k22/(dz*surf_viscosity[1]) );		  
  gsl_matrix_set( &dfdphi_mat.matrix,nz,nz,(-(k22/dz) - wa[1]*pow(cos(pretwist[1] + omega_d[1]*t - phi[nz]),2) + wa[1]*pow(sin(pretwist[1] + omega_d[1]*t - phi[nz]),2))/surf_viscosity[1] );
    
  return GSL_SUCCESS;
  
};

