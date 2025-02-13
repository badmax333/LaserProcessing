/*

Contact mechanics simulation by Green's function molecular dynamics (GFMD) method in continuum formulation (reference implementation).

Copyright 2025 Leonid Dorogin.

*/

// ###############################################################
// Programming language features common for many scientific tasks.
// ###############################################################


#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <math.h>

#include <time.h>

typedef double _Complex Complex;

// "Pi" constant might already be defined in the C/C++ library, or use the following line to define it.
// #define M_PI 3.14159265358979323846

// Include the programming definitions to use the external library "FFTW" verion 3.* for the Fourier transformation procedures.
#include <fftw3.h>




// ##################################
// Definition of the physical system.
// ##################################

/*

Consider contact of two elastic bodies with defined elastic moduli.
Then we will redefine the system as rigid punch against elastic substrate.

*/


double punch_Young_modulus = 73.1e+9; // INPUT; Young modulus of the punch; measures in Newton/(meter*meter); 73.1 GPa for silica glass
double punch_Poisson_ratio = 0.17; // INPUT; Poisson's ratio of the punch; 0.17 for silica glass.
double substrate_Young_modulus = 1.6e+6; // INPUT; Young modulus of the substrate; measures in Newton/(meter*meter); ~1.6 MPa static limit for PDMS Sylgard 184.
double substrate_Poisson_ratio = 0.5; // INPUT; Poisson's ratio of the substrate; ~0.5 typical for PDMS Sylgard 184.
double X_length = 0.1; // INPUT; size of the system, the system is periodic, so same as period; measures in meters.
double punch_pressure = 1.0e+5; // INPUT; magnitude of applied load as pressure; measures in Newton/(meter*meter).

// INPUT: two geometries are possible:
//     (1) cylindrical punch on a flat base surface;
//     (2) sinusoidal wave punch.
// Comment out the following line to solve the (2)nd problem.
//#define cylindrical_punch

#ifdef cylindrical_punch
double punch_radius = 0.020; // INPUT; radius of cylindrical rigid punch; measures in meters.
#endif

double punch_amplitude = 0.010; // INPUT; height of the rigid punch or wave amplitude; measures in meters.



// INPUT: two out-of-plane boundary conditions are possible:
//    (1) zero-displacement, i.e. zero-strain, i.e. fixed ("pipe") boundary, also called "plane strain" problem;
//    (2) zero-stress, i.e. free ("film") boundary, also called "plane stress" problem.
// Comment out the following line to solve the (2)nd problem.
//#define outplane_is_fixed


/*

Below we use the redefined system as rigid punch against elastic substrate.
It is assumed that the contact mechanics properties such as contact area are preserved.
For more details we refer to the book: K. L. Johnson, "Contact Mechanics", Cambridge University Press, New York (1985).

*/

#ifndef outplane_is_fixed
// Plane stress solution can be obtained from plane strain solution by the following substitution.
// For details, see, e.g. the book of A. I. Lurie "Theory of Elasticity", Springer (2005).
double effective_punch_Poisson_ratio = punch_Poisson_ratio/(1+punch_Poisson_ratio);
double effective_substrate_Poisson_ratio = substrate_Poisson_ratio/(1+substrate_Poisson_ratio);
#endif
#ifdef outplane_is_fixed
// If the plain strain problem is desired, keep them intact.
double effective_punch_Poisson_ratio = punch_Poisson_ratio;
double effective_substrate_Poisson_ratio = substrate_Poisson_ratio;
#endif

/*
Hereby the reformulated problem assumes incompressible substrate.
The closer "effective_punch_Poisson" and "effective_substrate_Poisson_ratio" to 0.5, the more accuracy is achieved.
*/

double effective_elastic_modulus=1/((1-effective_punch_Poisson_ratio*effective_punch_Poisson_ratio)/punch_Young_modulus+(1-effective_substrate_Poisson_ratio*effective_substrate_Poisson_ratio)/substrate_Young_modulus); // OUTPUT; measures in Newton/(meter*meter).
double punch_force = punch_pressure*X_length; // OUTPUT; magnitude of applied load; measures in Newton/meter.
double contact_width; // OUTPUT; measures in meters.
double indentation_depth; // OUTPUT; measures in meters.
double average_contact_gap; // OUTPUT; measures in meters.
double Hertzian_contact_width; // OUTPUT; measures in meters.

#ifndef cylindrical_punch
double punch_radius; // OUTPUT; radius of wave punch asperity; measures in meters.
#endif



// ###############################
// Discretization and fine tuning.
// ###############################

//
// These are the parameters you need to adjust to get an accurate solution.
// If the result is noisy, crazy orders of magnitude, it often means the mathematical time increment is not fine (small) enough.
//
long X_points=128; // INPUT; number of spacial discretization points; it should be a power of 2 (e.g. 128, 256, 512 ...).
double mathematical_time_increment_prefactor = 0.001; // INPUT; adjusts energy-minimization virtual time increment (step).
long mathematical_time_steps_prefactor = 4; // INPUT; adjusts the number of mathematical time steps to be taken during the simulation.
double damping_prefactor = 1; // INPUT.


long half_X_points = X_points/2;
long Fourier_X_points = X_points/2+1; // number of degrees of freedom in Fourier space. See: http://www.fftw.org/fftw3_doc/The-1d-Real_002ddata-DFT.html#The-1d-Real_002ddata-DFT
double mathematical_time_increment; // energy-minimization virtual time increment (step).
long mathematical_time_steps; // how many mathematical time steps will be taken during the simulation.

//
// Dynamically allocated arrays.
//
double* punch_surface_profile;
double* surface_displacement;
double* surface_force;
Complex* Fourier_surface_displacement;
Complex* Fourier_surface_displacement_old;
Complex* Fourier_surface_force;
double* damping_factors;


//
// Output the state of the surface of the substrate to files.
//
void write_surface_state_output(const char* filename1, const char* filename2)
{

	//
	// Open the files.
	//
	FILE * displacement_output_file = fopen (filename1,"a");
	assert (displacement_output_file!=NULL);
	FILE * force_output_file = fopen (filename2,"a");
	assert (force_output_file!=NULL);


	double X;

	for (long ix=0; ix<X_points; ix++)
	{
			X=(double)ix*X_length/(double)X_points;

			fprintf (displacement_output_file, "%E\t%E\t%E\n", X, *(punch_surface_profile+ix), *(surface_displacement+ix));
			fprintf (force_output_file, "%E\t%E\n", X, *(surface_force+ix));

	}


	//
	// Close the files.
	//
	fclose (displacement_output_file);
	fclose (force_output_file);
}

double calc_time_estimate(clock_t steps_time, clock_t init_duration,
                          int current_step_count, int total_step_count)
{
	double avg_step_time = steps_time / (double) current_step_count;
	double estimate = init_duration + avg_step_time * total_step_count;
	estimate /= CLOCKS_PER_SEC;
	return estimate;
}

// File handler for the report file.
FILE * report_file;



//
// The entry point of the program.
//
int main(int argc, char** argv)
{
	// ***************************
	// Declare internal variables.
	// ***************************
	double X, Q;
	long ix;
	double sphere_boundary;
	Complex Fourier_surface_displacement_new;
	long contact_points;
	double surface_height_max, surface_height_min, gap_sum;


	clock_t computation_start_time = clock();


	// **************************************************************
	// Allocate memory for the spacial arrays
	//  and initialize all necessary variables, including the arrays.
	// **************************************************************


	//
	// Define the rigid punch.
	//
	punch_surface_profile = (double*)fftw_malloc(sizeof(double)*X_points); // Shape of the rigid punch as a height function.

#ifdef cylindrical_punch
	sphere_boundary = sqrt(punch_radius*punch_radius-(punch_radius-punch_amplitude)*(punch_radius-punch_amplitude));

    // Tip: The following "for" loop can be parallelized for private X variable.
	for (ix = 0; ix < X_points; ix++)
	{
		X = X_length*(double)ix/(double)X_points;

		if (X<(X_length/2-sphere_boundary))
		{
			*(punch_surface_profile+ix) = punch_amplitude;
			continue;
		}

		if (X>(X_length/2+sphere_boundary))
		{
			*(punch_surface_profile+ix) = punch_amplitude;
			continue;
		}

		*(punch_surface_profile+ix) = punch_radius - sqrt(punch_radius*punch_radius - (X-X_length/2)*(X-X_length/2));
	}
#endif
#ifndef cylindrical_punch
    // Tip: The following "for" loop can be parallelized for private X variable.
	for (ix = 0; ix < X_points; ix++)
	{
		X = X_length*(double)ix/(double)X_points;
		*(punch_surface_profile+ix) = punch_amplitude*(cos(2*M_PI*X/X_length)+1)/2;
	}
	punch_radius = X_length*X_length/(2*M_PI*M_PI*punch_amplitude);
#endif

	// Hertzian analytical solution.
	Hertzian_contact_width = 4*sqrt(punch_radius*punch_force/(M_PI*effective_elastic_modulus));


	surface_displacement = (double*)fftw_malloc(sizeof(double)*X_points); // Normal displacement (height profile) at the surface of the substrate.
	memset(surface_displacement,0,sizeof(double)*X_points);



	surface_force = (double*)fftw_malloc(sizeof(double)*X_points); // Normal force at the surface of the substrate.
	memset(surface_force,0,sizeof(double)*X_points);


	Fourier_surface_displacement = (Complex*)fftw_malloc(sizeof(Complex)*X_points); // Normal displacement on the surface of the substrate as Fourier image.
	memset(Fourier_surface_displacement,0,sizeof(Complex)*X_points);



	Fourier_surface_displacement_old = (Complex*)fftw_malloc(sizeof(Complex)*X_points); // Temporary data storage for the "Verlet" loop.
	memset(Fourier_surface_displacement_old,0,sizeof(Complex)*X_points);


	Fourier_surface_force = (Complex*)fftw_malloc(sizeof(Complex)*X_points); // Force acting on the body from the counterbody converted to Fourier image and added with the non-physical damping force.
	memset(Fourier_surface_force,0,sizeof(Complex)*X_points);



	//
	// Find the mathematical time increment and the number of mathemtical time steps to take for the simulation.
	//
	mathematical_time_increment = mathematical_time_increment_prefactor/sqrt(effective_elastic_modulus/((double)X_points*X_length));
	mathematical_time_steps = 256*(double)mathematical_time_steps_prefactor*sqrt((double)X_points);


	damping_factors = (double*)fftw_malloc(sizeof(double)*Fourier_X_points); // Define non-physical damping factors.
	*(damping_factors+0) = 0.75*sqrt(1.0/((double)X_points*X_length));

    // Tip: The following "for" loop can be parallelized.
	for (ix = 1; ix < Fourier_X_points; ix++)
	{
		*(damping_factors+ix) = damping_prefactor*(2*sqrt(M_PI*(double)ix*effective_elastic_modulus/((double)X_points*X_length))-M_PI*(double)ix*effective_elastic_modulus*mathematical_time_increment/((double)X_points*X_length));
	}


  // Plan FFTW operations.

	fftw_plan surface_displacement_fftw_forward = fftw_plan_dft_r2c_1d(X_points, surface_displacement, (fftw_complex*) Fourier_surface_displacement, FFTW_MEASURE);
	fftw_plan surface_displacement_fftw_backward = fftw_plan_dft_c2r_1d(X_points, (fftw_complex*) Fourier_surface_displacement, surface_displacement, FFTW_MEASURE);
	fftw_plan surface_force_fftw_backward = fftw_plan_dft_c2r_1d(X_points, (fftw_complex*) Fourier_surface_force, surface_force, FFTW_MEASURE);



	clock_t initialization_time = clock() - computation_start_time;
	int estimate_iters = 100;
	int only_time_estimate = (argc == 2 && !strcmp("--benchmark", argv[1]));
	clock_t loop_start_time = clock();

	// *************************************************************************
	// Iterate for the energy minimization of the body. It is not physical time.
	// *************************************************************************

	for (int current_mathematical_time_step=0; current_mathematical_time_step<mathematical_time_steps; current_mathematical_time_step++)
	{
		if (only_time_estimate)
		{
			if (current_mathematical_time_step == estimate_iters)
			{
				printf("estimated computation time: %f s\n",
				       calc_time_estimate(clock() - loop_start_time,
				                          initialization_time,
				                          estimate_iters,
				                          mathematical_time_steps));
				return 0;
			}
		}

		// Obtain "Fourier_surface_displacement" from "surface_displacement".
		fftw_execute(surface_displacement_fftw_forward);


		*(Fourier_surface_force+0) = punch_force/X_length;

    // Tip: The following "for" loop can be parallelized for private Q variable.
		for (ix=1; ix<Fourier_X_points; ix++)
		{
			Q = 2*M_PI*ix/X_length;
			*(Fourier_surface_force+ix) = -Q*effective_elastic_modulus*(*(Fourier_surface_displacement+ix))/(2*(double)X_points);
		}


		//
		// "Verlet" process to the "Fourier_surface_displacement" and "Fourier_surface_displacement_old" arrays.
		//
		Complex force;

      // Tip: The following "for" loop can be parallelized
      //  for private force and Fourier_surface_displacement_new variables.
		for (ix=0; ix<Fourier_X_points;ix++)
		{
			force = *(Fourier_surface_force+ix) - (*(damping_factors+ix))*(*(Fourier_surface_displacement+ix)-*(Fourier_surface_displacement_old+ix))/mathematical_time_increment;
			Fourier_surface_displacement_new = 2.0*(*(Fourier_surface_displacement+ix)) - *(Fourier_surface_displacement_old+ix) + force*mathematical_time_increment*mathematical_time_increment;
			*(Fourier_surface_displacement_old+ix) = *(Fourier_surface_displacement+ix);
			*(Fourier_surface_displacement+ix) = Fourier_surface_displacement_new;
		}


		// Obtain displacements in the real coordinates from their Fourier image.
		fftw_execute(surface_displacement_fftw_backward);


		//
		// Apply the punch interaction as a hardwall.
		//
		contact_points = 0;

        // Tip: The following "for" loop can be parallelized
        // for shared (+reduction) contact_points variable.
		for (ix=0; ix<X_points;ix++)
		{
			// First normalize the value after the FFTW transformation.
			*(surface_displacement+ix) = *(surface_displacement+ix) / (double)X_points;

			if ((*(punch_surface_profile+ix)) <= (*(surface_displacement+ix)))
			{
				*(surface_displacement+ix) = *(punch_surface_profile+ix);
				contact_points = contact_points + 1;
			}
		}
	}
	// End of the iterations.




	// *********************
	// Finalize the results.
	// *********************


	fftw_execute(surface_force_fftw_backward);
	// Note, no need to normalize this array, because it is defined directly as a Fourier image according to physical units.


	contact_width = (double)contact_points*X_length/(double)X_points;


	//
	// Find the indentation depth as a distance between the maximum and minimum heights, and the mean gap between the bodies.
	//
	surface_height_max = *(surface_displacement+0);
	surface_height_min = *(surface_displacement+0);
	gap_sum = 0;

    // Tip: The following "for" loop can be parallelized
    //  with few properly shared variables.	
    for (ix=0; ix<X_points; ix++)
	{
		if (*(surface_displacement+ix)>surface_height_max)
		{
			surface_height_max = *(surface_displacement+ix);
		}
		if (*(surface_displacement+ix)<surface_height_min)
		{
			surface_height_min = *(surface_displacement+ix);
		}
		gap_sum = gap_sum + (*(punch_surface_profile+ix) - *(surface_displacement+ix));
	}
	indentation_depth = surface_height_max-surface_height_min;
	average_contact_gap = gap_sum/(double)X_points;




	// **************************************************************
	// Display the results and save the state of the system to files.
	// **************************************************************


	// Open a file for the textual report as HTML.
	report_file = fopen("output/textual_report.htm","w");
	assert(report_file!=NULL);


	//
	// Print the summary of input parameters and output results.
	//
	printf("\nINPUT for the GFMD computer simulation:\n"); // for the display output
	fprintf(report_file, "<P><B><U>Input</U></B> for the GFMD computer simulation:<BR>\n"); // for the report file and so on.
#ifdef outplane_is_fixed
	printf("zero out-of-plane displacements (\"pipe\" mode).\n");
	fprintf(report_file, "zero out-of-plane displacements (\"pipe\" mode).<BR>\n");
#endif
#ifndef outplane_is_fixed
	printf("zero out-of-plane stress (\"free film\" mode).\n");
	fprintf(report_file, "zero out-of-plane stress (\"free film\" mode).<BR>\n");
#endif
#ifdef cylindrical_punch
	printf("cylindrical punch of radius [m] = %E\n", punch_radius);
	fprintf(report_file, "cylindrical punch of radius [m] = %E<BR>\n", punch_radius);
	printf("punch height [m] = %E\n", punch_amplitude);
	fprintf(report_file, "punch height [m] = %E<BR>\n", punch_amplitude);
#endif
#ifndef cylindrical_punch
	printf("wave punch asperity radius [m] = %E\n", punch_radius);
	fprintf(report_file, "wave punch asperity radius [m] = %E<BR>\n", punch_radius);
	printf("wave punch amplitude [m] = %E\n", punch_amplitude);
	fprintf(report_file, "wave punch amplitude [m] = %E<BR>\n", punch_amplitude);
#endif
	printf("Young modulus of the punch [Pa] = %E\n", punch_Young_modulus);
	fprintf(report_file, "Young modulus of the punch [Pa] = %E<BR>\n", punch_Young_modulus);
	printf("Poisson's ratio of the punch = %E\n", punch_Poisson_ratio);
	fprintf(report_file, "Poisson's ratio of the punch = %E<BR>\n", punch_Poisson_ratio);
	printf("Young modulus of the substrate [Pa] = %E\n", substrate_Young_modulus);
	fprintf(report_file, "Young modulus of the substrate [Pa] = %E<BR>\n", substrate_Young_modulus);
	printf("Poisson's ratio of the substrate = %E\n", substrate_Poisson_ratio);
	fprintf(report_file, "Poisson's ratio of the substrate = %E<BR>\n", substrate_Poisson_ratio);
	printf("punch force per unit area [Pa] = %E\n", punch_pressure);
	fprintf(report_file, "punch force per unit area [Pa] = %E<BR>\n", punch_pressure);
#ifdef cylindrical_punch
	printf("period of the system in X direction [m] = %E\n", X_length);
	fprintf(report_file, "period of the system in X direction [m] = %E<BR>\n", X_length);
#endif
#ifndef cylindrical_punch
	printf("period of the system in X direction [m] = wave punch period [m] = %E\n", X_length);
	fprintf(report_file, "period of the system in X direction [m] = wave punch period [m] = %E<BR>\n", X_length);
#endif
	printf("number of the discretization points in X direction [m] = %E\n", (double)X_points);
	fprintf(report_file, "number of the discretization points in X direction [m] = %E<BR></P>\n", (double)X_points);

	printf("\nOUTPUT of the computer simulation:\n");
	fprintf(report_file, "\n<P><B><U>Output</U></B> of the computer simulation:<BR>\n");
	printf("converted to the problem of a rigid punch on an elastic substrate.\n");
	fprintf(report_file, "converted to the problem of a rigid punch on an elastic substrate.<BR>\n");
	printf("effective elastic modulus of the substrate [Pa] = %E\n", effective_elastic_modulus);
	fprintf(report_file, "effective elastic modulus of the substrate [Pa] = %E<BR>\n", effective_elastic_modulus);
	printf("contact width  [m] = %E\n", contact_width);
	fprintf(report_file, "contact width  [m] = %E<BR>\n", contact_width);
	printf("contact width from Hertzian solution [m] = %E\n", Hertzian_contact_width);
	fprintf(report_file, "contact width from Hertzian solution [m] = %E<BR>\n", Hertzian_contact_width);
	printf("relative contact area = %E\n", contact_width/X_length);
	fprintf(report_file, "relative contact area = %E<BR>\n", contact_width/X_length);
	printf("indentation depth  [m] = %E\n", indentation_depth);
	fprintf(report_file, "indentation depth  [m] = %E<BR>\n", indentation_depth);
	printf("average contact gap  [m] = %E\n", average_contact_gap);
	fprintf(report_file, "average contact gap  [m] = %E<BR></P>\n", average_contact_gap);

	printf("\nINTERNALS of the computer simulation:\n");
	fprintf(report_file, "\n<P><B>Internals</B> of the computer simulation:<BR>\n");
	printf("mathematical time increment = %E\n", mathematical_time_increment);
	fprintf(report_file, "mathematical time increment = %E<BR>\n", mathematical_time_increment);
	printf("mathematical time steps = %E\n", (double)mathematical_time_steps);
	fprintf(report_file, "mathematical time steps = %E<BR></P>\n", (double)mathematical_time_steps);


	printf("\n");


	// Close the report file.
	fclose(report_file);


	write_surface_state_output("output/x_displacement_output.txt", "output/x_force_output.txt");




	// ****************************************************
	// Release the allocated memory for the spacial arrays.
	// ****************************************************


	//
	// Free the memory for FFTW metadata structures.
	//
	fftw_destroy_plan(surface_displacement_fftw_forward);
	fftw_destroy_plan(surface_displacement_fftw_backward);
	fftw_destroy_plan(surface_force_fftw_backward);


	//
	// Free the memory for the data arrays.
	//
	fftw_free(punch_surface_profile);
	fftw_free(surface_displacement);
	fftw_free(surface_force);
	fftw_free(Fourier_surface_displacement);
	fftw_free(Fourier_surface_displacement_old);
	fftw_free(Fourier_surface_force);
	fftw_free(damping_factors);

}
