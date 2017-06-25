#include <stdio.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <math.h>
#include <stdlib.h>
#include "./director.h"
#include "./parser.h"
const int nz=60;


int main (int argc, char * argv[]) {

  const double pi=3.141592653589793;
  double * theta, * phi;
  struct lc_cell lc_environment;
  double ti=0.0, tf=50.0;
  double time=ti, dz,dt=1e-6;
  double timeprint=0.2;
  FILE * time_file, * snapshot_file;
  const char * time_file_name="middle_sin.dat";    
  int timesteper_kind_flag=0;

  //Standard values:
  lc_environment.k11=1.0;
  lc_environment.k22=1.0;
  lc_environment.k33=1.0;
  lc_environment.q=0.0;
  lc_environment.n0=1.0;
  lc_environment.ne=1.0;
  lc_environment.viscosity=1.0;
  lc_environment.cell_length=1.0;

  lc_environment.surf_viscosity[0]=0.1*lc_environment.viscosity*lc_environment.cell_length;
  lc_environment.surf_viscosity[1]=0.1*lc_environment.viscosity*lc_environment.cell_length;

  lc_environment.pretilt[0]=0.0;
  lc_environment.pretilt[1]=0.0;

  lc_environment.pretwist[0]=0.0;
  lc_environment.pretwist[1]=0.0;
  

  lc_environment.wa[0]=1.0;
  lc_environment.wa[1]=1.0;

  lc_environment.omega_d[0]=0.0;
  lc_environment.omega_d[1]=0.0;


  //Read the parameter values form the input file:
  parse_input_file(  & lc_environment, & tf, & timeprint , & dt );
  print_log_file( lc_environment, tf, dt, "log.file");
  
  dz=lc_environment.cell_length/nz;


  //Starting the PDE solver:
  gsl_odeiv2_system sys = {frank_energy, jacobian, nz+1, &lc_environment};


  gsl_odeiv2_driver * pde_driver =gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-9, 0.0);
  //gsl_odeiv2_driver * pde_driver =gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_msbdf, 1e-8, 1e-8, 0.0);


  time_file=fopen(time_file_name,"w");
  snapshot_file=fopen("phi_time.dat","w");

  phi= (double *) malloc( (nz+1)*sizeof(double) );

  for (int ii=0; ii<=nz;ii++)
    {

      phi[ii]=(lc_environment.q*ii)/nz;
      
    };

 
  print_snapshot_to_file(phi,time,dz,snapshot_file);
  fprintf(time_file,"%f  %f \n",time, sin(phi[nz/2]) );

  while(time <tf)
    {
      int status = gsl_odeiv2_driver_apply (pde_driver, &time, time+timeprint, phi);      

      
      if (status != GSL_SUCCESS)
	{
	  printf ("error, return value=%d\n", status);
   
	};

      
      print_snapshot_to_file(phi,time,dz,snapshot_file);
      fprintf(time_file,"%f  %f \n",time, sin(phi[nz/2]) );

    };
  

  gsl_odeiv2_driver_free (pde_driver);
  free(phi);
  fclose(time_file);
  fclose(snapshot_file);
  return 0;


};
     

    
int frank_energy (double t, const double phi[], double f[], void * params)
{ 
  struct lc_cell mu = *(struct lc_cell *)params;
  double dz = mu.cell_length/nz;
  double k22= mu.k22;
  double omega_d[2],wa[2], pretwist[2];
  double viscosity, surf_viscosity[2];
  double dphi, d2phi, q;
  omega_d[0]=mu.omega_d[0];
  omega_d[1]=mu.omega_d[1];
  wa[0]=mu.wa[0];
  wa[1]=mu.wa[1];
  pretwist[0]=mu.pretwist[0];
  pretwist[1]=mu.pretwist[1];
  viscosity=mu.viscosity;
  surf_viscosity[0]=mu.surf_viscosity[0];
  surf_viscosity[1]=mu.surf_viscosity[1];
  q=mu.q;

  /*Bulk equations */  
  for(int ii=1; ii<nz; ii++)
    {

      d2phi=(phi[ii+1]+phi[ii-1]-2.0*phi[ii])/(dz*dz);

      f[ii]= d2phi*k22/viscosity;
          

    };

  f[0]=(-(k22*(-((phi[1]-phi[0])/dz) + q)) + wa[0]*cos(pretwist[0] + omega_d[0]*t - phi[0])*sin(pretwist[0] + omega_d[0]*t - phi[0]))/surf_viscosity[0];


f[nz]=(k22*(-((phi[nz]-phi[nz-1])/dz) + q) + wa[1]*cos(pretwist[1] + omega_d[1]*t - phi[nz])*sin(pretwist[1] + omega_d[1]*t - phi[nz]))/surf_viscosity[1];

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
  
  for(int i=1; i<nz;i++)
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


int print_snapshot_to_file(const double * phi,
			   const double time,
			   const double dz,
			   FILE * snapshot_file)
{

      fprintf(snapshot_file,"#time=%f\n",time);

      for(int ii=0;ii<=nz;ii++)
	{
	  
	  fprintf(snapshot_file,"%f  %f\n",ii*dz,phi[ii]);
      

	};
      fprintf(snapshot_file,"\n\n");

};

int print_phi_time( const double * phi,
		    const double time,
		    const double dz )
			       
{
  FILE * snapshot_file;
  char snap_name[50];

  sprintf(snap_name,"time=%f.dat",time);
  snapshot_file=fopen(snap_name,"w");

  fprintf(snapshot_file,"#i  nx        ny         nz         phi   time=%f, dz=%f\n",time,dz);

  for(int ii=0;ii<=nz;ii++)
    {

      fprintf(snapshot_file,"%i  %f  %f  %f  %f\n",ii,cos(phi[ii]),sin(phi[ii]), 0.0, phi[ii]);
      

    };

  fclose(snapshot_file);
  
  return 0;
};


void print_log_file(const struct lc_cell lc,
		    const double  tf,
		    const double  dt,
		    const char something[])
{

  printf("\n\nValues used for the parameters:\n\n");
  printf( "Kii:                        %lf  %lf  %lf\n",lc.k11,lc.k22,lc.k33 );
  printf( "timestep(initial dt):       %lf \n",dt);
  printf( "cell length:                %lf \n",lc.cell_length);
  printf( "Chiral power (q):           %lf \n",lc.q);
  printf( "bulk viscosity:             %lf \n",lc.viscosity);
  printf( "surface viscosity:          %lf  %lf \n",lc.surf_viscosity[0], lc.surf_viscosity[1]);
  printf( "anchoring energy(wa):       %lf  %lf \n",lc.wa[0], lc.wa[1]);
  printf( "twsiting velocity(omega_d): %lf  %lf \n\n",lc.omega_d[0], lc.omega_d[1]);
    
    
};

