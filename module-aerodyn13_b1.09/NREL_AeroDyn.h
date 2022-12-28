/* $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/modules/module-aerodyn/NREL_AeroDyn.h,v 1.8 2017/01/12 14:47:15 masarati Exp $ */
/* 
 * Copyright (C) 2003-2017
 *
 * Pierangelo Masarati	<masarati@aero.polimi.it>
 * Paolo Mantegazza	<mantegazza@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 *
 * This header file is free software; you can redistribute it at will,
 * under the same license conditions of the AeroDyn package.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 */

#ifndef AERO_DYN_H
#define AERO_DYN_H

/*
 * NOTE to gfortran users:
 *
 * - get AeroDyn 13
 * - run
	gfortran -O -c *.f90
	ar ru libAeroDyn.a *.o
 * - place libAeroDyn.a where the linker can find it,
 *   or tweak Makefile.inc as appropriate
 */

/*
 * NOTE to icc users:
 *
 * compile f90 files with

	ifc -r8

 * for double precision; defaults to single precision
 * link C++ executable with

	g++ -L /opt/intel/ia32/lib/ -lF90 -lCEPCF90 -lintrins

 */


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

//#if defined(USE_SINGLE_PRECISION)
//typedef float F_REAL;
//#elif defined(USE_DOUBLE_PRECISION)
//typedef double F_REAL;
//#else /* !USE_SINGLE_PRECISION && !USE_DOUBLE_PRECISION */
//#error "define either USE_SINGLE_PRECISION or USE_DOUBLE_PRECISION"
//#endif /* !USE_SINGLE_PRECISION && !USE_DOUBLE_PRECISION */
typedef float F_REAL;
typedef long int F_LOGICAL;
typedef char F_CHAR;
typedef long int F_INTEGER;

/*
 * Info from:

			USER'S GUIDE
	to the Wind Turbine Aerodynamics Computer Software
			AeroDyn
			
	David J. Laino
	A. Craig Hansen
	Windward Engineering, LC
	Salt Lake City, UT 84117
	www.windwardengineering.com

	Phone:	801-278-7852
	Fax:	801-272-4132
	email:	dlaino@windwardengineering.com
		chansen@windwardengineering.com
	
	Software date and version
	AeroDyn 12.43, 26-Apr-2002
	
	Prepared for the
	National Renewable Energy Laboratory
	under Subcontract No. TCX-9-29209-01

*/

/*
 * AeroDyn initialization; must be called as early as possible.
 */

extern int
__FC_DECL__(mbdyn_ad_inputgate)(F_REAL *b_length, F_CHAR *input_file_name, F_INTEGER *ifnamelen, F_CHAR *elem_file_name, F_INTEGER *efnamelen);

/*
 * Returns the force and moment for a given element.
 */
extern int
__FC_DECL__(call_ad_calculateloads)(F_REAL *curr_time,F_INTEGER *curr_blade,F_INTEGER *curr_elem,F_REAL * force,F_REAL *moment);

/*
 * interface to set data to aeordyn from mbdyn
 */

extern int
__FC_DECL__(set_aeromarkers)(F_INTEGER *curr_blade, F_INTEGER *curr_elem, F_REAL *position, F_REAL orientation[][3], F_REAL *translationvel, F_REAL *rotationvel);
extern int
__FC_DECL__(set_ad_interfacecomponent)( F_INTEGER *num_component,  F_INTEGER *curr_blede, F_REAL *position, F_REAL orientation[][3], F_REAL *transVel, F_REAL *rotVel);
/*
 * inteface to get data from aerodyn
 */
extern int
__FC_DECL__(get_horizonal_wind_speed)(F_REAL *ZTime, F_REAL *HorWindV);

/*
 * MBDyn stuff initialization
 */
extern int
__FC_DECL__(mbdyn_init)(F_CHAR *Version, F_INTEGER *nblades, F_INTEGER *nelems);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* AERO_DYN_H */

