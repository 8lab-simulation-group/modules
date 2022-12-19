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

struct Marker{
	F_REAL			Position[3];
	F_REAL 			Orientation[3][3];
	F_REAL			TranslationVel[3];
	F_REAL  		RotationVel[3];
};

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
__FC_DECL__(ad_inputgate)(F_CHAR *input_file);
extern int
__FC_DECL__(adinputgate)(void);
extern int
__FC_DECL__(mbdyn_ad_inputgate)(F_REAL *b_length, F_CHAR *input_file_name, F_INTEGER *ifnamelen, F_CHAR *elem_file_name, F_INTEGER *efnamelen);

// ADDED BY JENS VAN SCHELVE TO PROVIDE AERODYN ELEMENT DATA OUTPUT
extern int
__FC_DECL__(elemout)(void);

/*
 * Returns the force and moment for a given element.
 */
extern int
__FC_DECL__(call_ad_calculateloads)(F_REAL *c_time, F_REAL *b_length,  float DFN[], float DFT[], float PMA[]);

/*
 * interface to set data to aeordyn from mbdyn
 */
extern int
__FC_DECL__(set_nacelle_components)(F_REAL *nacelle_position, F_REAL *nacelle_orientation, F_REAL *nacelle_transvel, F_REAL *nacelle_rotvel);
extern int
__FC_DECL__(set_hub_components)(F_REAL *hub_position, F_REAL *hub_orientation, F_REAL *hub_transvel, F_REAL *hub_rotvel);
extern int
__FC_DECL__(set_rotorfurl_components)(F_REAL *rotorfurl_position, F_REAL *rotorfurl_orientation, F_REAL *rotorfurl_transvel, F_REAL *rotorfurl_rotvel);
extern int
__FC_DECL__(set_tower_components)(F_REAL *tower_position, F_REAL *tower_orientation, F_REAL *tower_transvel, F_REAL *tower_rotvel);
extern int
__FC_DECL__(set_blades_components)(F_INTEGER *c_blade, F_REAL *bladeroot_position, F_REAL *bladeroot_orientation, F_REAL *bladeroot_transvel, F_REAL *bladeroot_rotvel);
extern int
__FC_DECL__(set_bladeelem_component)(F_INTEGER *c_blade, F_INTEGER *c_elem, F_REAL * bladeelem_position, F_REAL * bladeelem_orientation, F_REAL *bladeelem_transvel, F_REAL *bladeelem_rotvel);

/*
 * inteface to get data from aerodyn
 */
extern int
__FC_DECL__(get_horizonal_wind_speed)(F_REAL *ZTime, F_REAL *HorWindV);

/*
 * Rotor parameters - called once per time step.
 */
/*
 * Write an error message to the appropriate stream
 *
 * FIXME: the "msg" and "level" arrays should be reset by the caller
 * before writing the message, otherwise they're not '\0' terminated
 */
extern int
__FC_DECL__(usrmes)(F_LOGICAL *Logical, F_CHAR msg[],
		F_INTEGER *code, F_CHAR level[]);

/*
 * self explanatory :)
 */
extern int
__FC_DECL__(ad_abort)(void);

/*
 * MBDyn stuff initialization
 */
extern int
__FC_DECL__(mbdyn_init)(F_CHAR *Version, F_INTEGER *nblades, F_INTEGER *nelems);

/*

/*
 * MBDyn logical true
 */
extern int
__FC_DECL__(mbdyn_true)(F_LOGICAL *val);

/*
 * MBDyn logical false
 */
extern int
__FC_DECL__(mbdyn_false)(F_LOGICAL *val);

/* 
 * This subroutine is to pass the current simulation time 
 * from MBDyn to AeroDyn!
 * c_time: current time
 * By Fanzhong MENG 19 June 2008
 */
extern int
__FC_DECL__(mbdyn_sim_time)(doublereal *c_time);

/* 
 * This subroutine is to pass the current simulation time step 
 * from MBDyn to AeroDyn!
 * dt: time step
 * By Fanzhong MENG 19 June 2008
 */
extern int
__FC_DECL__(mbdyn_time_step)(F_REAL *dt);



#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* AERO_DYN_H */

