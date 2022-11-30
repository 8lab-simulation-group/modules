/* $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/modules/module-aerodyn/module-aerodyn.cc,v 1.46 2017/01/12 14:47:15 masarati Exp $ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation (version 2 of the License).
 * 
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */
/* module-aerodyn
 * Authors: Fanzhong Meng, Pierangelo Masarati
 *
 * Copyright (C) 2008-2017
 *
 * Fanzhong Meng <f.meng@tudelft.nl>
 * 
 * Faculty of Aerospace Engineering - Delft University of Technology
 * Kluyverweg 1, 2629HS Delft, the Netherlands
 * http://www.tudelft.nl
 *
 */
/*
 * AeroDynModule: interface between MBDyn and NREL's AeroDyn.
 *
 * Author: Pierangelo Masarati
 *
 * Copyright (C) 2008-2017 All rights reserved
 *
 * Refactoring of module-aerodyn using the UserDefinedElem interface.
 * Based on module-aerodyn.cc by Fanzhong Meng and Pierangelo Masarati.
 */

#include "module-aerodyn13_b1.03.h"
#include "module-discon.h"

AeroDynModule::AeroDynModule(
	unsigned uLabel,
	const DofOwner *pDO,
	DataManager* pDM,
	MBDynParser& HP)
: Elem(uLabel, 0),
UserDefinedElem(uLabel, pDO),
bFirst(true)
{
	if (HP.IsKeyWord("help")) {
		/* NOTE: add help message */
		silent_cout(
		"Module: AeroDyn" << std::endl
		<< std::endl
		<< "Author: Fanzhong Meng, Pierangelo Masarati" << std::endl
		<< std::endl
		<< "This is the MBDyn interface to AeroDyn, the aerodynamic routines" << std::endl
		<< "developed by NREL <http://www.nrel.gov/> to model the aerodynamic" << std::endl
		<< "forces acting on wind turbines" << std::endl
		<< std::endl
		<< "usage:" << std::endl
		<< "#Nacelle node; requirements:" << std::endl
		<< "#\t\t- axis 3 is the shaft axis, pointing downstream the wind" << std::endl
		<< "\t<nacelle node label> ," << std::endl
		<< "\t<hub node label> ," << std::endl
		<< "\t<pylon top-hub xy distance> ," << std::endl
		<< "\t<hub radius> ," << std::endl
		<< "\t<rotor radius> ," << std::endl
		<< "\t<number of blades> ," << std::endl
		<< "\t<number of elements per blade> ," << std::endl
		<< "\t# for each blade..." << std::endl
		<< "\t#\t- axis 1 in the blade spanwise direction, towards blade tip" << std::endl
		<< "\t#\t- axis 2 in coord direction, towards leading edge" << std::endl
		<< "\t#\t- axis 3 in thickness direction, towards low pressure side" << std::endl
		<< "\t\t<i-th blade root orientation matrix> ," << std::endl
		<< "\t\t# for each blade element..." << std::endl
		<< "\t\t\t<i-th blade j-th node label>" << std::endl
		<< "\t\t\t[ , orientation , <i-th blade j-th node orientation> ]" << std::endl
		<< "\t[ , force scale factor , (DriveCaller)<factor> ]" << std::endl
		<< "\t[ , output file name , \"<file name>\" ]" << std::endl
		<< "\t[ , AeroDyn input file name , \"<file name>\" ]" << std::endl
		<< "\t[ , element file name , \"<file name>\" ]" << std::endl
		);
	}
	if (::module_aerodyn != 0) {
		silent_cerr("AeroDynModule::AeroDynModule(" << uLabel << "): "
			"AeroDyn already in use" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	/* read nacelle node */
	pNacelle = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));
	if (pNacelle == 0) {
		silent_cerr("AeroDynModule(" << GetLabel() << "): "
			"nacelle node not defined "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	/* read hub node */
	pHub = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));
	if (pHub == 0) {
		silent_cerr("AeroDynModule(" << GetLabel() << "): "
			"hub node not defined "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	
	/* read rotorfurl node */
	pRotorFurl = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));
	if (pRotorFurl == 0) {
		silent_cerr("AeroDynModule(" << GetLabel() << "): "
			"rotorfurl node not defined "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	
	/* read Tower node */
	pTower = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));
	if (pTower == 0) {
		silent_cerr("AeroDynModule(" << GetLabel() << "): "
			"tower node not defined "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	// read blade length
	b_length = HP.GetReal();
	/*
	 * Initialize AeroDyn package
	 *
	 * FIXME: does it return an error code?
	 */
	char Version[26 + 1];
	snprintf(Version, sizeof(Version), "(%s)", VERSION);
	for (unsigned i = strlen(Version); i < sizeof(Version); i++) {
		Version[i] = ' ';
	}
	Version[sizeof(Version) - 1] = '\0';

	// number of blades
	F_INTEGER NBlades = HP.GetInt();
	if (NBlades <= 0) {
		silent_cerr("AeroDynModule(" << GetLabel() << "): "
			"invalid number of blades " << NBlades
			<< " at line " << HP.GetLineData()
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	
	
	// number of elements per blade
	F_INTEGER NElems = HP.GetInt();
	if (NElems <= 0) {
		silent_cerr("AeroDynModule(" << GetLabel() << "): "
			"invalid number of elements per blade " << NElems
			<< " at line " << HP.GetLineData()
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	nblades = NBlades;
	nelems = NElems;

	// blade root node resize
	pBladeroot.resize(NBlades);
	// read blade root node lavel 
	for (blade = 0; blade < nblades; blade++){
		pBladeroot[blade].pNode = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));
		if (pBladeroot[blade].pNode == 0) {
			silent_cerr("AeroDynModule(" << GetLabel() << "): "
				"hub node not defined "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}
		
	
	// get blade root orientation(without pitch angle).

    nodes.resize(NBlades*NElems);
    bladeR.resize(NBlades);
	ReferenceFrame rf(pHub);
	for (elem = 0; elem < NBlades*NElems; elem++) {
		if ((elem % NElems) == 0) {
		    	/*
			 * Get the orientation matrix of blade root with respect to hub reference frame
			 */ 
			bladeR[elem/NElems] = HP.GetRotRel(rf);
#if 0
		std::cerr << "Root[" << elem/NElems << "] Rotation Matrix:" << bladeR[elem/NElems] << std::endl;
		std::cerr << "Hub Rotation Matrix:" << pHub->GetRCurr() << std::endl;
		std::cerr << "Nacelle Rotation Matrix:" << pNacelle->GetRCurr() << std::endl<< std::endl;
#endif
		}

		nodes[elem].pNode = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));

		// (optional) aerodynamics offset with respect to the node,
		// constant in the reference frame of the node
		if (HP.IsKeyWord("position")) {
			nodes[elem].f = HP.GetPosRel(ReferenceFrame(nodes[elem].pNode));

		} else {
			nodes[elem].f = Zero3;
		}

		// (optional) relative orientation between the aerodynamics
		// and the node
		if (HP.IsKeyWord("orientation")) {
			nodes[elem].Ra = HP.GetRotRel(ReferenceFrame(nodes[elem].pNode));

		} else {
			nodes[elem].Ra = Eye3;
			nodes[elem].dBuiltInTwist = 0.;
		}
	}

	if (HP.IsKeyWord("force" "scale" "factor")) {
		FSF.Set(HP.GetDriveCaller());

	} else {
		FSF.Set(new OneDriveCaller);
	}

	if (HP.IsKeyWord("output" "file" "name")) {
		const char *ofname = HP.GetFileName();
		if (ofname == 0) {
			silent_cerr("AeroDynModule(" << GetLabel() << "): "
				"unable to get file name "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		ofname = ofname;
		out.open(ofname);
		if (!out) {
			silent_cerr("AeroDynModule(" << GetLabel() << "): "
				"unable to open file \"" << ofname << "\" "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		} else {
			out << std::setw(16) << "label"
				<< std::setw(16) << "Time s"
			    << std::setw(16) << "Thrust kN"
			    << std::setw(16) << "Torque kNm"
			    << std::setw(16) << "R.Speed Rad/s"
			    << std::setw(16) << "Power kW"
				<< std::setw(16) << "Wind Speed m/s"
			    << std::endl;
		}
	}

	// Aerodynに引数を渡す
	__FC_DECL__(mbdyn_init)(Version, &NBlades, &NElems);
	
	std::string input_file_name;
	std::string elem_file_name;

	if (HP.IsKeyWord("input" "file" "name")) {
		const char *fname = HP.GetStringWithDelims();
		if (fname == 0) {
			silent_cerr("unable to get input file name "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		input_file_name = fname;
	}

	if (HP.IsKeyWord("element" "file" "name")) {
		const char *fname = HP.GetStringWithDelims();
		if (fname == 0) {
			silent_cerr("unable to get element file name "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		elem_file_name = fname;
	}

	F_INTEGER input_file_name_len = input_file_name.size();
	F_INTEGER elem_file_name_len = elem_file_name.size();


	// Aerodyn 初期化
	const Vec3& nacelle_p = pNacelle->GetXCurr();
	const Vec3& nacelle_t = pNacelle->GetVCurr();
	const Vec3& nacelle_r = pNacelle->GetWCurr();
	const Vec3& hub_p = pHub->GetXCurr();
	const Vec3& hub_t = pHub->GetVCurr();
	const Vec3& hub_r = pHub->GetWCurr();
	const Vec3& rotorfurl_p = pRotorFurl->GetXCurr();
	const Vec3& rotorfurl_t = pRotorFurl->GetVCurr();
	const Vec3& rotorfurl_r = pRotorFurl->GetWCurr();
	const Vec3& tower_p = pTower->GetXCurr();
	const Vec3& tower_t = pTower->GetVCurr();
	const Vec3& tower_r = pTower->GetWCurr();
	const Mat3x3& nacelle_o = pNacelle->GetRCurr();
	const Mat3x3& hub_o = pHub->GetRCurr();
	const Mat3x3& rotorfurl_o = pRotorFurl->GetRCurr();
	const Mat3x3& tower_o = pTower->GetRCurr();

	for (int i = 0; i < 3; i++){
		nacelle_position[i] = nacelle_p.dGet(i+1);
		nacelle_transvel[i] = nacelle_t.dGet(i+1);
		nacelle_rotvel[i] = nacelle_r.dGet(i+1);
		hub_position[i] = hub_p.dGet(i+1);
    	hub_transvel[i] = hub_t.dGet(i+1);
		hub_rotvel[i] = hub_r.dGet(i+1);
		rotorfurl_position[i] = rotorfurl_p.dGet(i+1);
    	rotorfurl_transvel[i] = rotorfurl_t.dGet(i+1);
		rotorfurl_rotvel[i] = rotorfurl_r.dGet(i+1);
		tower_position[i] = tower_p.dGet(i+1);
    	tower_transvel[i] = tower_t.dGet(i+1);
		tower_rotvel[i] = tower_r.dGet(i+1);
		for(int	j = 0; j<3; j++){
					nacelle_orientation[j+i*3] = nacelle_o.dGet(i+1,j+1);
					hub_orientation[j + i * 3] = hub_o.dGet(i + 1, j + 1);
					rotorfurl_orientation[j + i * 3] = rotorfurl_o.dGet(i + 1, j + 1);
					tower_orientation[j + i * 3] = tower_o.dGet(i + 1, j + 1);
		}
	}

	__FC_DECL__(set_nacelle_components)(nacelle_position, nacelle_orientation, nacelle_transvel, nacelle_rotvel);
	__FC_DECL__(set_hub_components)(hub_position, hub_orientation, hub_transvel, hub_rotvel);
	__FC_DECL__(set_rotorfurl_components)(rotorfurl_position, rotorfurl_orientation, rotorfurl_transvel, rotorfurl_rotvel);
	__FC_DECL__(set_tower_components)(tower_position, tower_orientation, tower_transvel, tower_rotvel);
	
	
	c_blade = 0;
	for (blade = 0; blade < NBlades; blade++){
		const Vec3& blade_p = pBladeroot[blade].pNode->GetXCurr();
		const Vec3& blade_t = pBladeroot[blade].pNode->GetVCurr();
		const Vec3& blade_r = pBladeroot[blade].pNode->GetWCurr();
		const Mat3x3& blade_o = pBladeroot[blade].pNode->GetRCurr();
		
		for (int i = 0; i < 3; i++){
			bladeroot_position[i] = blade_p.dGet(i+1);
			bladeroot_transvel[i] = blade_t.dGet(i+1);
			bladeroot_rotvel[i] = blade_r.dGet(i+1);
			for (int j = 0; j<3; j++){
				bladeroot_orientation[j + i * 3] = blade_o.dGet(i + 1, j + 1);
			}
		}

		c_blade = c_blade + 1;

		__FC_DECL__(set_blades_components)(&c_blade, bladeroot_position, bladeroot_orientation, bladeroot_transvel, bladeroot_rotvel);

	}
	
	
	int rc = __FC_DECL__(mbdyn_ad_inputgate)(&b_length, (F_CHAR *)input_file_name.c_str(), &input_file_name_len, (F_CHAR *)elem_file_name.c_str(), &elem_file_name_len);
	if (rc != 0) {
		silent_cerr("AeroDynModule(" << GetLabel() << "): "
			"initialization failed "
			"(err=" << rc << ")" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// FIXME: only needed when Beddoes and also Dynamic Inflow model are enabled 
	Time.Set(new TimeDriveCaller(pDM->pGetDrvHdl()));
	dCurTime = Time.dGet();
	(void)__FC_DECL__(mbdyn_sim_time)(&dCurTime);

	::module_aerodyn = this;


}

AeroDynModule::~AeroDynModule(void)
{
	NO_OP;
}

unsigned int
AeroDynModule::iGetNumDof(void) const
{
	return 0;
}

DofOrder::Order
AeroDynModule::GetDofType(unsigned int i) const
{
	return DofOrder::UNKNOWN;
}

void
AeroDynModule::Output(OutputHandler& OH) const
{
	if (out) {
		out << std::scientific
			<< std::setw(16) << GetLabel()
		    << std::setw(16) << dCurTime 
		    << std::setw(16) << Thrust/1000.0
		    << std::setw(16) << Torque/1000.0
		    << std::setw(16) << Rotor_speed
		    << std::setw(16) << Torque*Rotor_speed/1000.0
			<< std::setw(16) << Wind_speed
		    << std::endl; 
	}
}

std::ostream&
AeroDynModule::Restart(std::ostream& out) const
{
	return out << "not implemented yet;" << std::endl;
}

void
AeroDynModule::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 6*nblades*nelems;
	*piNumCols = 1;
}

VariableSubMatrixHandler& 
AeroDynModule::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	// should do something useful
	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler& 
AeroDynModule::AssRes(
	SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);

	WorkVec.ResizeReset(iNumRows);

	/*
	 * set sub-vector indices and coefs
	 */

	if (bFirst) {
		
		Mat3x3 BladeR;
		Vec3 BladeAxis;

		// force scale factor to "ramp up" loads
		doublereal dFSF = FSF.dGet();

		// sanity check
		ASSERT(dFSF >= 0.);
		ASSERT(dFSF <= 1.);

		// 空力計算に必要なnode情報の取得
		const Vec3& nacelle_p = pNacelle->GetXCurr();
		const Vec3& nacelle_t = pNacelle->GetVCurr();
		const Vec3& nacelle_r = pNacelle->GetWCurr();
		const Vec3& hub_p = pHub->GetXCurr();
		const Vec3& hub_t = pHub->GetVCurr();
		const Vec3& hub_r = pHub->GetWCurr();
		const Vec3& rotorfurl_p = pRotorFurl->GetXCurr();
		const Vec3& rotorfurl_t = pRotorFurl->GetVCurr();
		const Vec3& rotorfurl_r = pRotorFurl->GetWCurr();
		const Vec3& tower_p = pTower->GetXCurr();
		const Vec3& tower_t = pTower->GetVCurr();
		const Vec3& tower_r = pTower->GetWCurr();
		const Mat3x3& nacelle_o = pNacelle->GetRCurr();
		const Mat3x3& hub_o = pHub->GetRCurr();
		const Mat3x3& rotorfurl_o = pRotorFurl->GetRCurr();
		const Mat3x3& tower_o = pTower->GetRCurr();

		for (int i = 0; i < 3; i++){
			nacelle_position[i] = nacelle_p.dGet(i + 1);
			nacelle_transvel[i] = nacelle_t.dGet(i + 1);
			nacelle_rotvel[i] = nacelle_r.dGet(i + 1);
			hub_position[i] = hub_p.dGet(i + 1);
			hub_transvel[i] = hub_t.dGet(i + 1);
			hub_rotvel[i] = hub_r.dGet(i + 1);
			rotorfurl_position[i] = rotorfurl_p.dGet(i + 1);
			rotorfurl_transvel[i] = rotorfurl_t.dGet(i + 1);
			rotorfurl_rotvel[i] = rotorfurl_r.dGet(i + 1);
			tower_position[i] = tower_p.dGet(i + 1);
			tower_transvel[i] = tower_t.dGet(i + 1);
			tower_rotvel[i] = tower_r.dGet(i + 1);
			
			for (int j = 0; j<3; j++){
				nacelle_orientation[j + i * 3] = nacelle_o.dGet(i + 1, j + 1);
				hub_orientation[j + i * 3] = hub_o.dGet(i + 1, j + 1);
				rotorfurl_orientation[j + i * 3] = rotorfurl_o.dGet(i + 1, j + 1);
				tower_orientation[j + i * 3] = tower_o.dGet(i + 1, j + 1);
			}
		}

		__FC_DECL__(set_nacelle_components)(nacelle_position, nacelle_orientation, nacelle_transvel, nacelle_rotvel);
		__FC_DECL__(set_hub_components)(hub_position, hub_orientation, hub_transvel, hub_rotvel);
		__FC_DECL__(set_rotorfurl_components)(rotorfurl_position, rotorfurl_orientation, rotorfurl_transvel, rotorfurl_rotvel);
		__FC_DECL__(set_tower_components)(tower_position, tower_orientation, tower_transvel, tower_rotvel);


		c_blade = 0;
		for (blade = 0; blade < nblades; blade++){
			const Vec3& blade_p = pBladeroot[blade].pNode->GetXCurr();
			const Vec3& blade_t = pBladeroot[blade].pNode->GetVCurr();
			const Vec3& blade_r = pBladeroot[blade].pNode->GetWCurr();
			const Mat3x3& blade_o = pBladeroot[blade].pNode->GetRCurr();

			for (int i = 0; i < 3; i++){
				bladeroot_position[i] = blade_p.dGet(i + 1);
				bladeroot_transvel[i] = blade_t.dGet(i + 1);
				bladeroot_rotvel[i] = blade_r.dGet(i + 1);
				for (int j = 0; j<3; j++){
					bladeroot_orientation[j + i * 3] = blade_o.dGet(i + 1, j + 1);
					
				}
			}
			c_blade = c_blade + 1;

			__FC_DECL__(set_blades_components)(&c_blade, bladeroot_position, bladeroot_orientation, bladeroot_transvel, bladeroot_rotvel);

		}

		//AllAeroMaelersの取得に使用
		c_blade = 0;
		c_elem = 0;
		for (elem = 0; elem < nelems*nblades; elem++) {
			/*
			 * get the current blade number.
			 */
			if ((elem % nelems) == 0) {
				c_blade++;
				BladeR = pHub->GetRCurr()*bladeR[c_blade];
				BladeAxis = BladeR.GetVec(1);
			}

			/*
			 * get force/moment
			 */

			c_elem = elem%nelems + 1;

			const Mat3x3& bladeE_o = nodes[elem].pNode->GetRCurr()*nodes[elem].Ra;
			const Vec3& bladeE_p = nodes[elem].pNode->GetXCurr() + bladeE_o*nodes[elem].f;
			const Vec3& bladeE_t = nodes[elem].pNode->GetVCurr();
			const Vec3& bladeE_r = nodes[elem].pNode->GetWCurr();

			for (int i = 0; i < 3; i++){
				bladeelem_position[i] = bladeE_p.dGet(i + 1);
				bladeelem_transvel[i] = bladeE_t.dGet(i + 1);
				bladeelem_rotvel[i] = bladeE_r.dGet(i + 1);
				for (int j = 0; j < 3; j++){
					bladeelem_orientation[j + i * 3] = bladeE_o.dGet(i + 1,j + 1);
				}

			}
			__FC_DECL__(set_bladeelem_component)(&c_blade, &c_elem, bladeelem_position, bladeelem_orientation, bladeelem_transvel, bladeelem_rotvel);


		}

		F_REAL		DFN[nblades*nelems];
		F_REAL		DFT[nblades*nelems];
		F_REAL		PMA[nblades*nelems];

		c_time = Time.dGet();

		__FC_DECL__(call_ad_calculateloads)(&c_time, &b_length, DFN, DFT, PMA);

		
		c_blade = 0;
		int elem = 0;
		for (elem = 0; elem < nelems*nblades; elem++) {
			// ramp forces in the very first second up to ease convergence
			DFN[elem] *= dFSF;
			DFT[elem] *= dFSF;
			PMA[elem] *= dFSF;


			nodes[elem].FN = DFN[elem];
			nodes[elem].FT = DFT[elem];
			nodes[elem].AM = PMA[elem];


#ifdef MODULE_AERODYN_DEBUG
			silent_cerr("aerodyn[" << elem << "] in rotation plane: DFN=" << DFN << " DFT=" << DFT << " PMA=" << PMA << std::endl);
#endif // MODULE_AERODYN_DEBUG
			/*
			 * turn force/moment into the node frame (blade element frame used by Aerodyn to calculation forces)
			 * (passing thru local element frame),
			 * including offset
			 */

			nodes[elem].F = Vec3(DFN[elem], DFT[elem], 0.);
			nodes[elem].M = nodes[elem].Ra.GetVec(3)*PMA[elem] + nodes[elem].f.Cross(nodes[elem].F);
		
#ifdef MODULE_AERODYN_DEBUG

			silent_cerr("Moment on Node[" << elem << "]: " << nodes[elem].f.Cross(nodes[elem].F) << std::endl);

#endif // MODULE_AERODYN_DEBUG

#ifdef MODULE_AERODYN_DEBUG
			silent_cerr(std::setprecision(3) << std::fixed
				<< " aerodyn[" << elem << "]"
				<< " F = " << nodes[elem].F
				<< " M = " << nodes[elem].M
				<< " FN = " << nodes[elem].FN
				<< " FN = " << nodes[elem].FT
				<< " AM = " << nodes[elem].AM
				<< std::endl << std::endl);
#endif // MODULE_AERODYN_DEBUG

		}

#ifdef MODULE_AERODYN_DEBUG
		silent_cerr(std::endl);
#endif // MODULE_AERODYN_DEBUG

		bFirst = false;
	}


	TF = Vec3(Zero3);
	TM = Vec3(Zero3);
	for (elem = 0; elem < nblades*nelems; elem++) {
		/*
		 * set indices where force/moment need to be put
		 */
		integer iFirstIndex = nodes[elem].pNode->iGetFirstMomentumIndex();
		for (int i = 1; i <= 6; i++) {
			WorkVec.PutRowIndex(6*elem + i, iFirstIndex + i);
		}

		/*
		 * add force/moment to residual, after rotating them
		 * into the global frame
		 */
		/*
		 * Transfer the force/moment to global frame. 
	     */	 
  		WorkVec.Add(6*elem + 1, nodes[elem].pNode->GetRCurr()*(nodes[elem].Ra*nodes[elem].F));
  		WorkVec.Add(6*elem + 4, nodes[elem].pNode->GetRCurr()*(nodes[elem].Ra*nodes[elem].M));
		/*
		 * calculate the force and moment contributions on the rotor in absolute frame.
	     */
		const Vec3& Xp = nodes[elem].pNode->GetXCurr();
		const Vec3& Xh = pHub->GetXCurr();

		TF = TF + nodes[elem].pNode->GetRCurr()*nodes[elem].Ra*nodes[elem].F;
		TM = TM + nodes[elem].pNode->GetRCurr()*nodes[elem].Ra*nodes[elem].M
			+ (Xp-Xh).Cross(nodes[elem].pNode->GetRCurr()*nodes[elem].Ra*nodes[elem].F);
  	}

	/*
	 * Transfer the force/moment to hub reference frame.
	 */
	const Mat3x3& Rh = pHub->GetRCurr();
	TF_h = Rh.MulTV(TF);
	TM_h = Rh.MulTV(TM);

	/*
	 * The 1st component of TF_h and TM_h are the rotor Thrust and rotor Torque.
	 */
	Thrust = TF_h(1);
	Torque = TM_h(1);
	
	/*
	 * Calculating hub anguler velocity
	 */
	const Vec3& Wh = pHub->GetWCurr();
	const Vec3& Wn = pNacelle->GetWCurr();
	Rotor_speed =Rh.MulTV(Wh - Wn).dGet(1);
	Wind_speed = dGetWindSpeed();

	/* 
	 * make sure next time FirstLoop will be false
	 */
	if (FirstLoop) {
		(void)__FC_DECL__(mbdyn_false)(&FirstLoop);
	}
	
	return WorkVec;
}

#if 0
void
module_aerodyn_before_predict(const LoadableElem* pEl, 
	VectorHandler& X, VectorHandler& XP,
	VectorHandler& XPrev, VectorHandler& XPPrev)
{
	DEBUGCOUTFNAME("module_aerodyn_before_predict");
}
#endif

void
AeroDynModule::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
	bFirst = true;

	/*
	 * Get the current simulation time and time step!
	 * MENG: 19 June 2008.
	 */ 
	dDT = Time.dGet() - dOldTime;
	dCurTime = Time.dGet();
	(void)__FC_DECL__(mbdyn_sim_time)(&dCurTime);
	(void)__FC_DECL__(mbdyn_time_step)(&dDT);
}

#if 0
void
module_aerodyn_update(LoadableElem* pEl, 
	const VectorHandler& X, const VectorHandler& XP)
{
	DEBUGCOUTFNAME("module_aerodyn_update");
}
#endif

void 
AeroDynModule::AfterConvergence(
	const VectorHandler& X,
	const VectorHandler& XP)
{
	dOldTime = Time.dGet();
}

unsigned int
AeroDynModule::iGetInitialNumDof(void) const
{
	return 0;
}

void
AeroDynModule::InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 6*nblades*nelems;
	*piNumCols = 1;
}

VariableSubMatrixHandler&
AeroDynModule::InitialAssJac(
	VariableSubMatrixHandler& WorkMat, 
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler& 
AeroDynModule::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkVec.ResizeReset(0);

	return WorkVec;
}

void
AeroDynModule::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	bFirst = true;
}

unsigned int
AeroDynModule::iGetNumPrivData(void) const
{
	return 0;
}

int
AeroDynModule::GetNumConnectedNodes() const
{
	return nodes.size() + 2;
}

void
AeroDynModule::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	/*
	 * set args according to element connections
	 */
	connectedNodes.resize(nodes.size() + 2);

	unsigned n;
	for (n = 0; n < nodes.size(); n++) {
		connectedNodes[n] = nodes[n].pNode;
	}

	connectedNodes[n++] = pNacelle;
	connectedNodes[n++] = pHub;
}

// module_aerodyn->pNacelle
const StructNode *
AeroDynModule::pGetNacelleNode(void) const
{
	return pNacelle;
}

// module_aerodyn->pHub
const StructNode *
AeroDynModule::pGetHubNode(void) const
{
	return pHub;
}

// module_aerodyn->pRotorFurl
const StructNode *
AeroDynModule::pGetRotorFurlNode(void) const
{
	return pRotorFurl;
}
// module_aerodyn->pHub
const StructNode *
AeroDynModule::pGetTowerNode(void) const
{
	return pTower;
}
// module_aerodyn->Hub_Tower_xy_distance;
doublereal
AeroDynModule::dGetHubTowerXYDistance(void) const
{
	return Hub_Tower_xy_distance;
}

// module_aerodyn->nblades
F_INTEGER
AeroDynModule::iGetNumBlades(void) const
{
	return nblades;
}

// module_aerodyn->c_blade
F_INTEGER
AeroDynModule::iGetCurrBlade(void) const
{
	ASSERT(c_blade >= 0);
	ASSERT(c_blade < nblades);

	return c_blade;
}

// module_aerodyn->bladeR
const Mat3x3&
AeroDynModule::GetCurrBladeR(void) const
{
	ASSERT(c_blade >= 0);
	ASSERT(c_blade < nblades);

	return bladeR[c_blade - 1];
}

// module_aerodyn->nelems
F_INTEGER
AeroDynModule::iGetNumBladeElems(void) const
{
	return nelems;
}

// module_aerodyn->elem
F_INTEGER
AeroDynModule::iGetCurrBladeElem(void) const
{
	ASSERT(elem >= 0);
	ASSERT(elem < nblades*nelems);

	return elem;
}

// module_aerodyn->nodes
const StructNode *
AeroDynModule::pGetCurrBladeNode(void) const
{
	ASSERT(elem >= 0);
	ASSERT(elem < nblades*nelems);

	return nodes[elem].pNode;
}

const StructNode *
AeroDynModule::pGetCurrBladerootNode(void) const
{
	ASSERT(blade >= 0);
	ASSERT(blade < nblades);

	return pBladeroot[blade].pNode;
}

// module_aerodyn->nodes[::module_aerodyn->elem].dBuiltInTwist
doublereal
AeroDynModule::dGetCurrBladeNodeBuiltinTwist(void) const
{
	ASSERT(elem >= 0);
	ASSERT(elem < nblades*nelems);

	return nodes[elem].dBuiltInTwist;
}


// module_aerodyn->nodes[::module_aerodyn->elem].Ra
const Mat3x3&
AeroDynModule::GetCurrBladeNodeRa(void) const
{
	ASSERT(elem >= 0);
	ASSERT(elem < nblades*nelems);

	return nodes[elem].Ra;
}

doublereal 
AeroDynModule::dGetWindSpeed(void) const
{
	F_REAL HorWindV;
	F_REAL time = Time.dGet();
	__FC_DECL__(get_horizonal_wind_speed)(&time,&HorWindV);
	return HorWindV;
}
// Functions called by AeroDyn

int
__FC_DECL__(usrmes)(
	F_LOGICAL *Logical,
	F_CHAR msg[],
	F_INTEGER *code,
	F_CHAR level[])
{
#if 0
	// can't work, unless we know the size of the arrays!
	silent_cerr("module-aerodyn: msg=\"" << msg << "\" "
		"code=" << *code
		<< " level=" << level
		<< std::endl);
#endif
	silent_cerr("module-aerodyn: diagnostics from AeroDyn, "
		"code=" << *code << std::endl);

	return 0;
}

extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
	UserDefinedElemRead *rf = new UDERead<AeroDynModule>;

	if (!SetUDE("aerodyn", rf)) {
		delete rf;

		silent_cerr("module-aerodyn: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	UserDefinedElemRead *rfd = new UDERead<DisconModule>;

	if (!SetUDE("discon", rfd)) {
		delete rfd;

		silent_cerr("module-discon: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}

