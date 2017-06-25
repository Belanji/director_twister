#include "director.h"
#include "parser.h"
#include <string.h>
#include <stdlib.h>
const double pi=3.141592653589793;
 
void parse_input_file(struct lc_cell  * lc,
		      double * tf,
		      double * timeprint,
		      double * dt )
{

  char parser[80];
  char garbage[400];
  int error_handler;

  while (   scanf("%79s",parser) != EOF )
    {


      if ( strcmp(parser,"k11") == 0 || strcmp(parser,"K11") == 0 )
	{

	  error_handler=scanf("%lf",&lc->k11);

	  if (error_handler <= 0 )
	    {

	      printf("You placed a comment or a non numeric value after %s in your input file.\n",parser);
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);
	    };
		
	  
	  fgets(garbage,400,stdin);


	}
      else if ( strcmp(parser,"k22")== 0 || strcmp(parser,"K22")== 0 )
	{
	

	  error_handler=scanf("%lf",&lc->k22);


	  if (error_handler <= 0 )
	    {

	      printf("You placed a comment or a non numeric value after %s in your input file.\n",parser);
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);
	    };
		  
	  
	  fgets(garbage,400,stdin);
	  
	  
      
	}
      else if ( strcmp(parser,"k33") == 0 || strcmp(parser,"K33")== 0 )
	{

	  error_handler=scanf("%lf",&lc->k33);

	  if (error_handler <= 0 )
	    {

	      printf("You placed a comment or a non numeric value after %s in your input file.\n",parser);
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);
	    };
	  
	  fgets(garbage,400,stdin);


	}
      else if (strcmp(parser,"tf") == 0 || strcmp(parser,"run_time") ==0 )
	{


	  error_handler=scanf("%lf",&tf);


	  if (error_handler <= 0 )
	    {

	      printf("You placed a comment or a non numeric value after %s in your input file.\n",parser);
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);

	    };


	  fgets(garbage,400,stdin);


	}
      else if (strcmp(parser,"dt") == 0 || strcmp(parser,"maximum_timestep") ==0 )
	{


	  error_handler=scanf("%lf",&dt);


	  if (error_handler <= 0 )
	    {

	      printf("You placed a comment or a non numeric value after %s in your input file.\n",parser);
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);

	    };


	  fgets(garbage,400,stdin);


	}
      else if ( strcmp(parser,"timeprint") == 0 || strcmp(parser,"print_every") ==0 )
	{


	  error_handler=scanf("%lf",&timeprint);


	  if (error_handler <= 0 )
	    {

	      printf("You placed a comment or a non numeric value after %s in your input file.\n",parser);
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);
	    };


	  fgets(garbage,400,stdin);


	}
      else if (strcmp(parser,"bulk_viscosity") == 0  || strcmp(parser,"bviscosity") == 0 || strcmp(parser,"bulk_visc") == 0 || strcmp(parser,"gamma1") == 0 || strcmp(parser,"gamma_1") == 0 )
	{


	  error_handler=scanf("%lf",&lc->viscosity);

	  if (error_handler <= 0 )
	    {

	      printf("You placed a comment or a non numeric value after %s in your input file.\n",parser);
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);
	    };

	  
	  fgets(garbage,400,stdin);
	  
	}
      else if (strcmp(parser,"bottom_surface_viscosity") == 0  || strcmp(parser,"bottom_sviscosity") == 0 ||  strcmp(parser,"bsurf_visc") == 0 || strcmp(parser,"b_gammas") == 0 )
	{

	  
	  error_handler=scanf("%lf",&lc->surf_viscosity[0]);

	  if (error_handler <= 0 )
	    {

	      printf("You placed a comment or a non numeric value after %s in your input file.\n",parser);
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);
	    };
	  
	  fgets(garbage,400,stdin);
	  

	}
      else if (strcmp(parser,"top_surface_viscosity") == 0 || strcmp(parser,"top_sviscosity") == 0 ||  strcmp(parser,"tsurf_visc") == 0 || strcmp(parser,"t_gammas") == 0 )
	{

	  
	  error_handler=scanf("%lf",&lc->surf_viscosity[1]);

	  if (error_handler <= 0 )
	    {

	      printf("You placed a comment or a non numeric value after %s in your input file.\n",parser);
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);
	    };
	  
	  fgets(garbage,400,stdin);
	  

	}
      else if (strcmp(parser,"surface_viscosity") == 0 || strcmp(parser,"sviscosity") == 0 ||  strcmp(parser,"surf_visc") == 0 || strcmp(parser,"gammas") == 0 || strcmp(parser,"gamma_s") == 0 )
	{

	  
	  error_handler=scanf("%lf",&lc->surf_viscosity[0]);
	  if (error_handler <= 0 )
	    {

	      printf("You placed a comment or a non numeric value after %s in your input file.\n",parser);
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);

	    };
	  


	  error_handler=scanf("%lf",&lc->surf_viscosity[1]);
	  if (error_handler <= 0 )
	    {

	      printf("The keyword '%s' must be followed by 2 numeric arguments.\n",parser);
	      printf("The first is the bottom %s, while the second is the top %s.\n",parser,parser);	      
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);

	    };
	  
	  fgets(garbage,400,stdin);
	  

	}
      else if (strcmp(parser,"cell_length") == 0  || strcmp(parser,"cell_height") == 0 ||  strcmp(parser,"lenght_z") == 0 || strcmp(parser,"d")==0  )
	{

	  
	  error_handler=scanf("%lf",&lc->cell_length);
	  
	  if (error_handler <= 0 )
	    {

	      printf("You placed a comment or a non numeric value after %s in your input file.\n",parser);
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);
	    };
	  
	  fgets(garbage,400,stdin);
	  

	}
            else if (strcmp(parser,"chiral_power") == 0  || strcmp(parser,"q") == 0 )
	{

	  
	  error_handler=scanf("%lf",&lc->q);
	  
	  if (error_handler <= 0 )
	    {

	      printf("You placed a comment or a non numeric value after %s in your input file.\n",parser);
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);
	    };
	  
	  fgets(garbage,400,stdin);
	  

	}
      else if (strcmp(parser,"pitch") == 0 || strcmp(parser,"Pitch") == 0  || strcmp(parser,"P") == 0 )
	{

	  double lc_pitch;
	  error_handler=scanf("%lf",&lc_pitch);
	  
	  if (error_handler <= 0 )
	    {

	      printf("You placed a comment or a non numeric value after %s in your input file.\n",parser);
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);
	    };

	  lc->q=2.0*pi/lc_pitch;
	  
	  fgets(garbage,400,stdin);
	  

	}
      else if (strcmp(parser,"azimutal_ancoring") == 0  || strcmp(parser,"wa" ) == 0  )
	{

	  
	  error_handler=scanf("%lf",&lc->wa[0]);
	  
	  if (error_handler <= 0 )
	    {

	      printf("You placed a comment or a non numeric value after %s in your input file.\n",parser);
	      
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);
	    };

	  error_handler=scanf("%lf",&lc->wa[1]);


	  if (error_handler <= 0 )
	    {

	      printf("The keyword '%s' must be followed by 2 numeric arguments.\n",parser);
	      printf("The first is the bottom %s, while the second is the top %s.\n",parser,parser);	      
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);
	    };
	  
	  fgets(garbage,400,stdin);
	  
	}
      else if (strcmp(parser,"omega_d") == 0  || strcmp(parser,"omega" ) == 0 || strcmp(parser,"omega" ) == 0 || strcmp(parser,"omegad" ) == 0  )
	{
	  
	  error_handler=scanf("%lf",&lc->omega_d[0]);
	  
	  if (error_handler <= 0 )
	    {

	      printf("You placed a comment or a non numeric value after %s in your input file.\n",parser);
	      
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);
	    };

	  error_handler=scanf("%lf",&lc->omega_d[1]);

	  if (error_handler <= 0 )
	    {

	      printf("The keyword '%s' must be followed by 2 numeric arguments.\n",parser);
	      printf("The first is the bottom %s, while the second is the top %s.\n",parser,parser);	      
	      printf("Please review your input file.\n Aborting the program\n");
	      exit(0);
	    };
	  
	  fgets(garbage,400,stdin);
	  
	}
      else if (parser[0]=='#')
	{

	  fgets(garbage,400,stdin);

	}	  
      else
	{

	  printf("The parser did not recognize the option '%s'. Please review your input file\n", parser);
	  printf("Aborting the program\n");
	  exit(0);
	};
      
    };
};
