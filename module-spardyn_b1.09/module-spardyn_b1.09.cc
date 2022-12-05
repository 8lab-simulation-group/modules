/* $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/modules/module-ModuleSpardyn/module-ModuleSpardyn.cc,v 1.11 2017/01/12 14:47:25 masarati Exp $ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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
/*
 AUTHOR: Reinhard Resch <r.resch@secop.com>
        Copyright (C) 2011(-2017) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
*/



#include "mbconfig.h"           // This goes first in every *.c,*.cc file

#include <cassert>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <iostream>
#include <iomanip>
#include <limits>

#include "dataman.h"
#include "userelem.h"
#include "module-spardyn_b1.09.h"
#include "irragularwave.h"

//static IrragularWave Wave;

class ModuleSpardyn
: virtual public Elem, public UserDefinedElem
{
public:
	ModuleSpardyn(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~ModuleSpardyn(void);
	virtual void Output(OutputHandler& OH) const;
	virtual unsigned int iGetNumDof(void) const;
	virtual DofOrder::Order GetDofType(unsigned int i) const;
	virtual std::ostream& DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const;
	virtual std::ostream& DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const;
	virtual unsigned int iGetNumPrivData(void) const;
	virtual unsigned int iGetPrivDataIdx(const char *s) const;
	virtual doublereal dGetPrivData(unsigned int i) const;
	virtual void SetInitialValue(VectorHandler& X);
	virtual void Update(const VectorHandler& XCurr,const VectorHandler& XPrimeCurr);
	virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);
	virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
	VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);
	SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);
	virtual int iGetNumConnectedNodes(void) const;
	virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
	virtual void SetValue(DataManager *pDM, VectorHandler& X, VectorHandler& XP,
			SimulationEntity::Hints *ph);
	virtual std::ostream& Restart(std::ostream& out) const;
	virtual unsigned int iGetInitialNumDof(void) const;
	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
   	virtual VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		      const VectorHandler& XCurr);
   	virtual SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);

private:

    StructNode              *pNode;			// pointer of nodes
	Mat3x3		            Ra;				// relative orientation

    doublereal          	length;			// length of spar
    doublereal          	CmCoef;			// coefficient of inertia Force of morison equation
    doublereal          	CaCoef;			// coefficient of addmass Force of morison equation
    doublereal          	CdCoef;			// coefficient of drag Force of morison equation
	doublereal          	Buoency;		// buoency
	doublereal				SubCoef;		// percentage submerged of spar

    Vec3                    Force;
    Vec3                    Moment;       

	bool					bHaveSurface = false;
	doublereal				offset_GS;
    doublereal          	area;
    doublereal          	volume;
    doublereal          	norm;
	
    doublereal          	CavCoef;
    doublereal          	CdvCoef;

	doublereal              	waterDensity;
    doublereal              	gravity; 
	const doublereal        	pi = M_PI;

    doublereal                  CurrTime;
    doublereal                  OldTime ;
    doublereal                  timestep;
    const DataManager*          m_pDM;
	DriveOwner                  FSF;
	bool 						bFirst = true;;
	mutable std::ofstream       out;

	IrragularWave				Wave;

	Mat3x3 dF_dx,  dF_dg,  dF_dv,  dF_dw,  dM_dx,  dM_dg,  dM_dv,  dM_dw;
	Mat3x3 dF_dxP, dF_dgP, dF_dvP, dF_dwP, dM_dxP, dM_dgP, dM_dvP, dM_dwP;

	Vec3	XPP;
	Vec3	WP;

	Vec3 	r1;
	Vec3	r2;
	Mat3x3  omega_ref_tild;

private:
	void CalcForce( const Vec3& XPP, const Vec3& WP);
	void CalcForceJac(const Vec3& vP, const Vec3& wP);
	Mat3x3 MultV1V2T(const Vec3 & v,const Vec3& w);

};


ModuleSpardyn::ModuleSpardyn(
	unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: 	Elem(uLabel, flag(0)),
	UserDefinedElem(uLabel, pDO),
    m_pDM(pDM)
{
	if (HP.IsKeyWord("help")) {
		silent_cout(
			"\n"
			"Module: 	ModuleSpardyn\n"
			"\n"
			<< std::endl);

		if (!HP.IsArg()) {
			// Exit quietly if nothing else is provided
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	// get node
	pNode = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));
	if (!pNode) {
		silent_cerr("spardyn(" << GetLabel() << "): structural node expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// getting rilative orientation if needed.
    if (HP.IsKeyWord("orientation")) {
		Ra = HP.GetRotRel(ReferenceFrame(pNode));
	} else {
		Ra = Eye3;
	}

	// get length
	if (!HP.IsKeyWord("length")) {
		silent_cerr("spardyn(" << GetLabel() << "): keyword \"length\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	length = HP.GetReal();

	// get diameter
	if (!HP.IsKeyWord("diameter")) {
		silent_cerr("spardyn(" << GetLabel() << "): keyword \"length\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	doublereal diameter = HP.GetReal();

	// get ca
	if (!HP.IsKeyWord("ca")) {
		silent_cerr("spardyn(" << GetLabel() << "): keyword \"ca\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	doublereal ca = HP.GetReal();

	// get cd
	if (!HP.IsKeyWord("cm")) {
		silent_cerr("spardyn(" << GetLabel() << "): keyword \"cm\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	doublereal cm = HP.GetReal();

	// get cm
	if (!HP.IsKeyWord("cd")) {
		silent_cerr("spardyn(" << GetLabel() << "): keyword \"cd\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	doublereal cd = HP.GetReal();

	// get rho
	if (!HP.IsKeyWord("rho")) {
		silent_cerr("spardyn(" << GetLabel() << "): keyword \"rho\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	waterDensity = HP.GetReal();

	// get g
	if (!HP.IsKeyWord("gravity")) {
		silent_cerr("spardyn(" << GetLabel() << "): keyword \"gravity\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	gravity = HP.GetReal();
	
	// calculate coefficient
	CmCoef  = 0.25 * cm * waterDensity * pi * diameter * diameter * length;
	CaCoef  = 0.25 * ca * waterDensity * pi * diameter * diameter * length;
	CdCoef  = 0.50 * cd * waterDensity * pi * diameter            * length;
	Buoency = 0.25 * diameter * diameter * pi * length * waterDensity * gravity;	

	// get surface information if needed.
	if (HP.IsKeyWord("surface")) { 
		bHaveSurface = true;
		// get offset
		if (!HP.IsKeyWord("offset")) {
			silent_cerr("spardyn(" << GetLabel() << "): keyword \"offset\" expected at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		offset_GS = HP.GetReal();

		// get cav
		if (!HP.IsKeyWord("cav")) {
			silent_cerr("spardyn(" << GetLabel() << "): keyword \"cav\" expected at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		doublereal cav = HP.GetReal();

		// get cdv
		if (!HP.IsKeyWord("cdv")) {
			silent_cerr("spardyn(" << GetLabel() << "): keyword \"cdv\" expected at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		doublereal cdv = HP.GetReal();

		// get area
		if (!HP.IsKeyWord("area")) {
			silent_cerr("spardyn(" << GetLabel() << "): keyword \"area\" expected at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		area = HP.GetReal();

		// get volume
		if (!HP.IsKeyWord("volume")) {
			silent_cerr("spardyn(" << GetLabel() << "): keyword \"volume\" expected at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		volume = HP.GetReal();

		// get norm
		if (!HP.IsKeyWord("norm")) {
			silent_cerr("spardyn(" << GetLabel() << "): keyword \"norm\" expected at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		norm = HP.GetReal();

		//calculate coefficients 
        CavCoef =     cav *waterDensity * volume;
        CdvCoef = 0.5*cdv *waterDensity * area;
	}

	// getting Force scale factor if needed.
	if (HP.IsKeyWord("Force" "scale" "factor")) {
		FSF.Set(HP.GetDriveCaller());

	} else {
		FSF.Set(new OneDriveCaller);
	}

	// getting  wave input file name 
	std::string ipFileName;
    if (HP.IsKeyWord("wave" "input" "file")) {
        const char *ifname = HP.GetFileName();
        ipFileName = std::string(ifname);
    }
	Wave.input(ipFileName);

	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));

	pDM->GetLogFile() << "spardyn: "
		<< uLabel << " "
		<< std::endl;
}


ModuleSpardyn::~ModuleSpardyn(void)
{
	// destroy private data
}


// the number of private data of this element, translational and rotational velocity for getting accelerations
// The number of private degrees of freedom is determined by this member function.
unsigned int
ModuleSpardyn::iGetNumDof(void) const
{
	return 6;
}

// the type of equations of the private degrees of freedom
DofOrder::Order
ModuleSpardyn::GetDofType(unsigned int i) const
{
	return DofOrder::DIFFERENTIAL;
}

// this fuction is called befor the integration starts in order to set the initial values for the private dgrees of freedom
void
ModuleSpardyn::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	const integer iElementIndex = iGetFirstIndex();

	const Vec3 v  = pNode->GetVCurr();
	const Vec3 w  = pNode->GetWCurr();

	const Vec3 a  = Zero3;
	const Vec3 wp = Zero3;

	X.Put( iElementIndex + 1, v );
	X.Put( iElementIndex + 4, w );
	XP.Put(iElementIndex + 1, a );
	XP.Put(iElementIndex + 4, wp);
}

// save the current state of private degree of freedom
// this is needed for implementation of the Output() and dGetPrivData functions
void
ModuleSpardyn::Update(const VectorHandler& X,const VectorHandler& XP)
{
	const integer iElementIndex = iGetFirstIndex();

	XPP = Vec3(XP(iElementIndex + 1), XP(iElementIndex + 2), XP(iElementIndex + 3));
	WP  = Vec3(XP(iElementIndex + 4), XP(iElementIndex + 5), XP(iElementIndex + 6));
}

void
ModuleSpardyn::AfterPredict(VectorHandler& X,VectorHandler& XP)
{
	bFirst = true;
	Update(X, XP);

}

void
ModuleSpardyn::Output(OutputHandler& OH) const
{
	if (bToBeOutput()) {

		if (OH.UseText(OutputHandler::LOADABLE)) {

			OH.Loadable() << GetLabel()
				<< " " << SubCoef
				<< " " << XPP
				<< " " << WP
				<< " " << Force
				<< " " << Moment
				<< std::endl;
		}
	}
}

// the dimention of subvector and submatrix which contribute to residual and jacobian matrix
void
ModuleSpardyn::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 12;
	*piNumCols = 12;
}

// calculate morison equation
// Note : Rotation matrix  v_global = R * V_local
void 
ModuleSpardyn::CalcForce( const Vec3& vP, const Vec3& wP)
{
	const Vec3 	 x = pNode->GetXCurr();
	const Mat3x3 R = pNode->GetRCurr()*Ra;
	const Vec3 	 v = pNode->GetVCurr();
	const Vec3	 w = pNode->GetWCurr();
	const Vec3  e3 = R.GetCol(3);

	CurrTime = m_pDM->dGetTime();
	
	// Force scale factor to "ramp up" loads
	doublereal dFSF = FSF.dGet();

	// sanity check
	ASSERT(dFSF >= 0.);
	ASSERT(dFSF <= 1.);

	Force  = Zero3;
	Moment = Zero3;
	SubCoef = 0.0;

	Wave.calc_uvw(CurrTime, x);
	const doublereal zeta = Wave.get_water_elevation();
	if(zeta > (x.dGet(3) - length/2.0) ) {

		if(zeta > (x.dGet(3) + length/2.0)) {
			SubCoef = 1.0;
		} else {
			SubCoef = (zeta - (x.dGet(3)-length/2.0) )/length;
		}

		const Vec3 r_GB = R * Vec3(0., 0., 0.5*(SubCoef-1.0)*length );

		const Vec3 ur    = Wave.get_velocity() - v;
		const Vec3 ur_N  = e3.Cross(ur.Cross(e3));
		const Vec3 uP    = Wave.get_acceleration();
		const Vec3 uP_N  = e3.Cross(uP.Cross(e3));

		const Vec3 vP_N = e3.Cross(vP.Cross(e3));

		const Vec3 F_n  =( (-vP_N * CaCoef) + (uP_N * CmCoef) + (ur_N * (ur_N.Norm() * CdCoef) ) )*SubCoef;
		const Vec3 M_n  = r_GB.Cross(F_n);

		Vec3	B(0.,0.,Buoency*SubCoef);
		Vec3	F_bn = e3.Cross(B.Cross(e3));
		const Vec3 M_bn = r_GB.Cross(F_bn);

		Vec3 F_v = Zero3;
		Vec3 F_bv = Zero3;

		if (bHaveSurface) {
			const Vec3 r_GS = R * Vec3( 0., 0., offset_GS );
			const Vec3 x_S = x + r_GS;
			const Vec3 v_S = v + w.Cross(r_GS);
			const Vec3 vP_S = vP + wP.Cross(r_GS) + w.Cross(w.Cross(r_GS));

			if(zeta > x_S.dGet(3)) {
				Wave.calc_uvw(CurrTime,x_S);
				const doublereal p_w = Wave.get_dynamic_pressure();

				const Vec3 urS = Wave.get_velocity() - v_S;
				const Vec3 urS_V = e3 * ( ur.Dot(e3) );
				const Vec3 uPS = Wave.get_acceleration();
				const Vec3 uPS_V = e3 * ( uP.Dot(e3) );
				const Vec3 vPS_V = e3 * ( vP_S.Dot(e3) );

				F_v = -e3 * (norm * p_w * area) - vPS_V * CavCoef + urS_V * ( urS_V.Norm() * CdvCoef); 
				F_bv = e3 * (norm *(area * waterDensity * gravity * x_S.dGet(3)));
			}
		}
	
		Force  = (F_n + F_v + F_bn + F_bv)*dFSF;
		Moment = (M_n + M_bn ) * dFSF;
	}
}

void
ModuleSpardyn::CalcForceJac(const Vec3& vP, const Vec3& wP)
{
	const Vec3&   x = pNode->GetXCurr();
	const Mat3x3& R = pNode->GetRCurr()*Ra;
	const Mat3x3& R_ref = pNode->GetRRef()*Ra;
	const Vec3&   v = pNode->GetVCurr();
	const Vec3&	  w = pNode->GetWCurr();
	const Vec3&  e3 = R.GetCol(3);
	const Vec3&  e3_ref = R_ref.GetCol(3);

	CurrTime = m_pDM->dGetTime();
	doublereal dFSF = FSF.dGet();

	dF_dx = Zero3x3; 
	dF_dg = Zero3x3; 
	dF_dv = Zero3x3; 
	dF_dw = Zero3x3; 

	dM_dx = Zero3x3; 
	dM_dg = Zero3x3; 
	dM_dv = Zero3x3; 
	dM_dw = Zero3x3;

	dF_dxP = Zero3x3; 
	dF_dgP = Zero3x3; 
	dF_dvP = Zero3x3; 
	dF_dwP = Zero3x3;

	dM_dxP = Zero3x3; 
	dM_dgP = Zero3x3; 
	dM_dvP = Zero3x3; 
	dM_dwP = Zero3x3;

	Wave.calc_uvw(CurrTime, x);
	const doublereal zeta = Wave.get_water_elevation();

	// check if bottom surface is under the water
	if(zeta > (x.dGet(3) - length/2.0))
	{
		// calulate parameters 
		const Vec3 ur    = Wave.get_velocity() - v;
		const Vec3 ur_N  = e3.Cross(ur.Cross(e3));
		const Vec3 uP    = Wave.get_acceleration();
		const Vec3 uP_N  = e3.Cross(uP.Cross(e3));
		const Vec3 vP_N = e3.Cross(vP.Cross(e3));

		// jacobian of SubCoefficient		
		Vec3 	dCsub_dx;
		if(zeta > (x.dGet(3) + length/2.0)) {
			SubCoef = 1.0;
			dCsub_dx = Zero3;
		} else {
			const Vec3 dzeta_dx = Wave.get_dzeta_dx();
			SubCoef = (zeta - (x.dGet(3)-length/2.0) )/length;
			dCsub_dx = (dzeta_dx - Vec3(0.,0., 1.))/length;
		}

		// the vector r_GB and Jacobians 
		const Vec3 r_GB = R * Vec3(0., 0., 0.5*(SubCoef-1.0)*length );
		const Vec3 drGB_dCsub = R * Vec3(0., 0., 0.5*length);
		const Mat3x3 drGB_dx = MultV1V2T(drGB_dCsub, dCsub_dx);
		const Mat3x3 drGB_dg = -Mat3x3(MatCross, R_ref * Vec3( 0. , 0. , 0.5*(SubCoef-1.)*length ));
		const Mat3x3 rGB_tild = Mat3x3(MatCross, r_GB);

		// jabobian and tild matrix of e3
		const Mat3x3 e3_tild = Mat3x3(MatCross, e3);
		// d(e3)/dg
		const Mat3x3 de3_dg = (Mat3x3(MatCross, e3_ref))*(-1.0);

		// jacobian of normal component of buoency force
		const Vec3 F_bn = e3.Cross(Vec3(0.,0.,Buoency).Cross(e3));
		const Mat3x3 Fbn_tild = Mat3x3(MatCross, F_bn);
		const Mat3x3 B_tild = Mat3x3(MatCross, Vec3(0.,0.,Buoency) );

		const Mat3x3 dFbn_dx = MultV1V2T( F_bn , dCsub_dx );
		const Mat3x3 dMbn_dx = rGB_tild * dFbn_dx - Fbn_tild * drGB_dx;

		const Mat3x3 dFbn_dg = ( e3_tild * B_tild - Mat3x3(MatCross, Vec3(0.,0.,Buoency).Cross(e3)) ) * de3_dg * (SubCoef);
		const Mat3x3 dMbn_dg = rGB_tild * dFbn_dg - Fbn_tild * drGB_dx;

		// jacobian of morison addmass force
		const Mat3x3 Fna_tild = Mat3x3(MatCross, (vP_N*(-1.0*CaCoef*SubCoef)));

		const Mat3x3 dFna_dx = MultV1V2T( (-vP_N*CaCoef), dCsub_dx);
		const Mat3x3 dFna_dg = ( e3_tild * Mat3x3(MatCross, vP) - Mat3x3(MatCross, vP.Cross(e3)) )* de3_dg * (-1.0*CaCoef*SubCoef);
		const Mat3x3 dFna_dvP = e3_tild * e3_tild * (CaCoef * SubCoef);
		const Mat3x3 dMna_dx = rGB_tild * dFna_dx - Fna_tild * drGB_dx;
		const Mat3x3 dMna_dg = rGB_tild * dFna_dg - Fna_tild * drGB_dg;
		const Mat3x3 dMna_dvP = rGB_tild * dFna_dvP;

		// jacobian of morison inertia force
		const Mat3x3 Fni_tild = Mat3x3(MatCross, (uP_N*(CmCoef*SubCoef)));

		const Mat3x3 dFni_dx =  MultV1V2T( ( uP_N*CmCoef), dCsub_dx);
		const Mat3x3 dFni_dg = ( e3_tild * Mat3x3(MatCross, uP) - Mat3x3(MatCross, uP.Cross(e3)))* de3_dg * (CmCoef*SubCoef);
		const Mat3x3 dMni_dx = rGB_tild * dFni_dx - Fni_tild * drGB_dx;
		const Mat3x3 dMni_dg = rGB_tild * dFni_dg - Fni_tild * drGB_dg;

		// jacobian of morison drag force
		const Mat3x3 Fnd_tild = Mat3x3(MatCross, (ur_N*(ur_N.Norm()*CdCoef*SubCoef)));
		const Mat3x3 urN_urNT = MultV1V2T(ur_N, ur_N);
		const doublereal urN_abs = ur_N.Norm();
		// d(ur_N)/dg = ( tild(e3)*tild(ur) - tild(ur x e3) )d(e3)/dg
		const Mat3x3 dur_N_dg = ( e3_tild * Mat3x3(MatCross, ur) - Mat3x3(MatCross, ur.Cross(e3)) ) * de3_dg;
		// d(ur_N)/dv = - tild(e3)*tild(e3) * d(ur)/dv ( d(ur)/dv = -I  : ur = u - v )
		const Mat3x3 dur_N_dv = e3_tild*e3_tild;

		const Mat3x3 dFnd_dx = MultV1V2T( (ur_N*(ur_N.Norm()*CdCoef)), dCsub_dx);
		const Mat3x3 dFnd_dg = ( (urN_urNT * dur_N_dg)/urN_abs + dur_N_dg*urN_abs) * (CdCoef*SubCoef);
		const Mat3x3 dFnd_dv = ( (urN_urNT * dur_N_dv)/urN_abs + dur_N_dv*urN_abs) * (CdCoef*SubCoef);
		const Mat3x3 dMnd_dx = rGB_tild * dFnd_dx - Fnd_tild * drGB_dx;
		const Mat3x3 dMnd_dg = rGB_tild * dFnd_dg - Fnd_tild * drGB_dg;
		const Mat3x3 dMnd_dv = rGB_tild * dFnd_dv;

		Mat3x3 dFva_dg = Zero3x3; 
		Mat3x3 dFva_dvP = Zero3x3;
		Mat3x3 dFvi_dg = Zero3x3;
		Mat3x3 dFvd_dg = Zero3x3;
		Mat3x3 dFvd_dv = Zero3x3;

		Mat3x3 dFbv_dx = Zero3x3;
		Mat3x3 dFbv_dg = Zero3x3;

		if(bHaveSurface) {

			// calculate surface parameters
			const Vec3 r_GS = R * Vec3( 0., 0., offset_GS );
			const Vec3 x_S = x + r_GS;
			const Vec3 v_S = v + w.Cross(r_GS);
			const Vec3 vP_S = vP + wP.Cross(r_GS) + w.Cross(w.Cross(r_GS));

			// check if the surface is under the water
			if(zeta > x_S.dGet(3)) {

				Wave.calc_uvw(CurrTime,x_S);
				const doublereal p_w = Wave.get_dynamic_pressure();

				const Vec3 urS = Wave.get_velocity() - v_S;
				const Vec3 urS_V = e3 * ( ur.Dot(e3) );
				const Vec3 uPS = Wave.get_acceleration();
				const Vec3 uPS_V = e3 * ( uP.Dot(e3) );
				const Vec3 vPS_V = e3 * ( vP_S.Dot(e3) );

				// jacobian of virtical morison addmass force
				dFva_dg = ( MultV1V2T(e3, vP) * de3_dg  +  de3_dg *( e3.Dot(vP) ) ) * (-1.0*CaCoef);
				dFva_dvP = MultV1V2T(e3, e3)*(-1.0*CavCoef);

				// jacobian of virtical morison inertia force
				dFvi_dg = de3_dg * (norm*area*p_w);

				// jacobian of virtical morison drag force
				const Mat3x3 urV_urVT = MultV1V2T(urS_V, urS_V);
				const doublereal urV_abs = urS_V.Norm();
				// d(ur_v)/dg
				const Mat3x3 dur_V_dg = MultV1V2T(e3, urS)*de3_dg + de3_dg * (e3.Dot(urS));
				// d(ur_v)/dv = e3 e3T d(ur)/dv ( d(ur)/dv = -I )
				const Mat3x3 dur_V_dv = -MultV1V2T(e3,e3);
				dFvd_dg = ( (urV_urVT*dur_V_dg)/urV_abs + dur_N_dg * urV_abs ) * (CdvCoef*SubCoef);
				dFvd_dv = ( (urV_urVT*dur_V_dv)/urV_abs + dur_N_dv * urV_abs ) * (CdvCoef*SubCoef);

				// jacobian of virtical component of buoency force
				dFbv_dx = MultV1V2T( e3 , Vec3(0.,0.,1.) ) * (area * norm * waterDensity * gravity);
				const Vec3 dzs_dg = drGB_dg.MulTV( Vec3(0.,0.,1.) );
				dFbv_dg = ( MultV1V2T( e3, dzs_dg ) + de3_dg * x_S.dGet(3) ) * (area * norm * waterDensity * gravity);
			}
		}
		//           horizona morison                  virtical morison             buoency
		dF_dx  = ( (dFna_dx  + dFni_dx + dFnd_dx)                                  + dFbn_dx + dFbv_dx)*dFSF;
		dF_dg  = ( (dFna_dg  + dFni_dg + dFnd_dg) + (dFva_dg + dFvi_dg + dFvd_dg)  + dFbn_dg + dFbv_dg)*dFSF;
		dF_dv  = (                       dFnd_dv                       + dFvd_dv                      )*dFSF;
		dF_dvP = (  dFna_dvP                      +  dFva_dvP                                         )*dFSF;
	
		dM_dx  = ( (dMna_dx  + dMni_dx + dMnd_dx)                                  + dMbn_dx          )*dFSF;
		dM_dg  = ( (dMna_dg  + dMni_dg + dMnd_dg)                                  + dMbn_dg          )*dFSF;
		dM_dv  = (                       dMnd_dv                                                      )*dFSF;
		dM_dvP = (  dMna_dvP                                                                          )*dFSF;
	}	
}
SubVectorHandler& 
ModuleSpardyn::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);

	WorkVec.ResizeReset(iNumRows);

	const integer iElementIndex = iGetFirstIndex();
	const integer iMomentumIndex = pNode->iGetFirstMomentumIndex();

	for(int iCnt = 1; iCnt <=6; iCnt++){
		WorkVec.PutRowIndex(iCnt, 	iMomentumIndex + iCnt);
		WorkVec.PutRowIndex(iCnt+6, iElementIndex  + iCnt);
	}

	// Geting pravate data of this element from XCurr and XPrimeCurr.
	const Vec3 y1(XCurr(iElementIndex + 1) , XCurr(iElementIndex + 2), XCurr(iElementIndex + 3));
	const Vec3 y2(XCurr(iElementIndex + 4) , XCurr(iElementIndex + 5), XCurr(iElementIndex + 6));

	// Geting derivative private data of this element, i,e, acceleration
	const Vec3 y1_dot(XPrimeCurr(iElementIndex + 1) , XPrimeCurr(iElementIndex + 2), XPrimeCurr(iElementIndex + 3));
	const Vec3 y2_dot(XPrimeCurr(iElementIndex + 4) , XPrimeCurr(iElementIndex + 5), XPrimeCurr(iElementIndex + 6));

	const Vec3& v = pNode->GetVCurr();
	const Vec3& w = pNode->GetWCurr();

	CalcForce(y1_dot, y2_dot);
	r1 =  y1 - v;
	r2 =  y2 - w ;


	WorkVec.Put( 1,  Force   );
	WorkVec.Put( 4,  Moment  );
	WorkVec.Put( 7,  r1  );
	WorkVec.Put(10,  r2  );
	

	return WorkVec;
}

VariableSubMatrixHandler& 
ModuleSpardyn::AssJac(VariableSubMatrixHandler& WorkMat_,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	/*
	 *                   [ x  ]                     [   F(y', y)               ]
	 * private data, y = [ g  ],   formulation, r = [   r_GB x F(y', y)        ],  
	 *                   [ v  ]                     [ v - x'                   ]    
	 *                   [ w  ]                     [ w - (g' + R_1d(g)*w_ref) ]   
	 *
	 *  x,g,are private data of node
	 *  v,w are private data of this element, and stored translational and rotatioal velocity
	 *  derivation private data XPrimeCurr is store translational and rotatinal acceleration
	 *            1        4       7     10
	 *          [ dF/dx  dF/dg  dF/dv dF/dw ]   1 
	 *  dr/dy = [ dM/dx  dM/dg  dM/dv dM/dw ]   4
	 *          [ 0     0       I     0     ]   7
	 *          [ 0    <w_ref>  0     I     ]  10
	 * 
	 *  d(R_1d(g)*w_ref)/dg = -<w_ref> , <> is tilda operator
	 * 			  1        4       7     10
	 *          [ dF/dx' dF/dg' dF/dv' dF/dw'] 1
	 *  dr/dy'= [ dM/dx' dM/dg' dM/dv' dM/dw'] 4
	 *          [  -I      0      0      0   ] 7
	 *          [   0     -I      0      0   ] 10
	 * 
	 *  contribution to Jacobian Matrix is  J = -( dr/dy' + dCoef * dr/dy )
	 */
	/*  
	 * set indices to submatrix WrokMat
	 *    [ x : position index of node + 1~3]  <= private data of node
	 * y =[ g : position index of node + 4~6] 
	 *    [ y1: element index + 1~3]           <= private data of element  
	 *    [ y2: element index + 4~6] 
	 * 
	 *    [ r1 : momentum index of node + 1~3] Force imposed to the node
	 * r =[ r2 : momentum index of node + 4~6] Moment imposed to the node
	 *    [ r3 : element index + 1~3]          <= adding eqation by this element 
	 *    [ r4 : element indes + 4~6]
	 * 
	 *            x   g  y1  y2 (index)     
	 *           [              ] r1
	 * WorkMat = [              ] r2
	 *           [              ] r3
	 *           [              ] r4 (index)
	 */ 
	integer iNumRows = 0;
	integer iNumCols = 0;

	WorkSpaceDim(&iNumRows, &iNumCols);

	FullSubMatrixHandler& WorkMat = WorkMat_.SetFull();

	WorkMat.ResizeReset(iNumRows, iNumCols);

	const integer iElementIndex = iGetFirstIndex();
	const integer iPositonIndex = pNode->iGetFirstPositionIndex();
	const integer iMomentumIndex = pNode->iGetFirstMomentumIndex();

	for(int iCnt = 1; iCnt<=6; iCnt++){
		WorkMat.PutRowIndex(iCnt, iMomentumIndex + iCnt);
		WorkMat.PutColIndex(iCnt, iPositonIndex  + iCnt);

		WorkMat.PutRowIndex(iCnt + 6, iElementIndex + iCnt);
		WorkMat.PutColIndex(iCnt + 6, iElementIndex + iCnt);
	}

	const Vec3& omega_ref = pNode->GetWRef(); // predicted angular velocity of node

	const Mat3x3 I = Eye3;
	//const Mat3x3 omega_ref_tild = Mat3x3(MatCross, omega_ref);

	// Geting pravate data of this element from XCurr and XPrimeCurr.
	const Vec3 y1(XCurr(iElementIndex + 1) , XCurr(iElementIndex + 2), XCurr(iElementIndex + 3));
	const Vec3 y2(XCurr(iElementIndex + 4) , XCurr(iElementIndex + 5), XCurr(iElementIndex + 6));

	// Geting derivative private data of this element, i,e, acceleration
	const Vec3 y1_dot(XPrimeCurr(iElementIndex + 1) , XPrimeCurr(iElementIndex + 2), XPrimeCurr(iElementIndex + 3));
	const Vec3 y2_dot(XPrimeCurr(iElementIndex + 4) , XPrimeCurr(iElementIndex + 5), XPrimeCurr(iElementIndex + 6));

	
	CalcForceJac(y1_dot, y2_dot);
	omega_ref_tild = Mat3x3(MatCross, omega_ref);

	WorkMat.Put( 1,  1, -dF_dxP   - (dF_dx          * dCoef) );
	WorkMat.Put( 1,  4, -dF_dgP   - (dF_dg          * dCoef) );
	WorkMat.Put( 1,  7, -dF_dvP   - (dF_dv          * dCoef) );
	WorkMat.Put( 1, 10, -dF_dwP   - (dF_dw          * dCoef) );

	WorkMat.Put( 4,  1, -dM_dxP   - (dM_dx          * dCoef) );
	WorkMat.Put( 4,  4, -dM_dgP   - (dM_dg          * dCoef) );
	WorkMat.Put( 4,  7, -dM_dvP   - (dM_dv          * dCoef) );
	WorkMat.Put( 4, 10, -dM_dwP   - (dM_dw          * dCoef) );

	WorkMat.Put( 7,  1,    I                                 );
	WorkMat.Put( 7,  7,           - (I              * dCoef) );

	WorkMat.Put(10,  4,    I      - (omega_ref_tild * dCoef) );
	WorkMat.Put(10, 10,           - (I              * dCoef) );
	
	return WorkMat_;
}




// This member function is called if the statement "print: dof description;" is specified in the input file.
std::ostream&
ModuleSpardyn::DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const
{
	/*
	integer iIndex = iGetFirstIndex();

	out
		<< prefix << iIndex + 1 << ": motor torque derivative [MP]" << std::endl
		<< prefix << iIndex + 2 << ": motor torque [M]" << std::endl
		<< prefix << iIndex + 3 << ": relative angular velocity [omega]" << std::endl;
	*/
	return out;
}

// This member function is called if the statement "print: equation description;" is specified in the input file.
std::ostream&
ModuleSpardyn::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{
	/*
	integer iIndex = iGetFirstIndex();

	out
		<< prefix << iIndex + 1 << ": motor DAE [f1]" << std::endl
		<< prefix << iIndex + 2 << ": motor torque derivative [f2]" << std::endl
		<< prefix << iIndex + 3 << ": angular velocity derivative [f9]" << std::endl;
	*/
	return out;
}

//  This member function is called when a bind statement or a element drive is used to access private data of this element.
unsigned int
ModuleSpardyn::iGetNumPrivData(void) const
{
	return 6;
}

// The following string values can be specified in a bind statement or in an element drive caller.
unsigned int
ModuleSpardyn::iGetPrivDataIdx(const char *s) const
{
	static const struct {
		int index;
		char name[6];
	}

	data[] = {
			{ 1, "XPPx"},		
			{ 2, "XPPy"},	
			{ 3, "XPPz"},	
			{ 4, "WPx"},
			{ 5, "WPy"},
			{ 6, "WPz"}	
	};

	for (unsigned i = 0; i < sizeof(data) / sizeof(data[0]); ++i ) {
		if (0 == strcmp(data[i].name,s)) {
			return data[i].index;
		}
	}

	silent_cerr("acceleration(" << GetLabel() << "): no private data \"" << s << "\"" << std::endl);

	return 0;
}

// This member function is called whenever a bind statement or a element drive is used to access private data.
doublereal
ModuleSpardyn::dGetPrivData(unsigned int i) const
{
	/*
	switch (i) {
	case 1:
		return XPP(1);
	case 2:
		return XPP(2);
	case 3:
		return XPP(3);
	case 4:
		return WP(1);
	case 5:
		return WP(2);
	case 6:
		return WP(3);
	}

	throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	*/
}

// This member function is called if the statement "print: node connection;" is specified in the input file.
int
ModuleSpardyn::iGetNumConnectedNodes(void) const
{
	return 0;
}

// This member function is called if the statement "print: node connection;" is specified in the input file.
void
ModuleSpardyn::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(0);
}

// This member function is called if "make restart file;" is specified in the input file.
std::ostream&
ModuleSpardyn::Restart(std::ostream& out) const
{
	return out << "# ModuleSpardyn: not implemented" << std::endl;
}

/*
 * The initial assembly phase is needed only if initial values are provided which are not consistent with the algebraic constraints.
 * Since this element does not use algebraic constraints we do not have to implement
 * iGetInitialNumDof(), InitialWorkSpaceDim(), InitialAssJac(), InitialAssRes() and SetInitialValue().
 */
void
ModuleSpardyn::SetInitialValue(VectorHandler& X)
{
	return;
}

unsigned int
ModuleSpardyn::iGetInitialNumDof(void) const
{
	return 0;
}

void 
ModuleSpardyn::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
ModuleSpardyn::InitialAssJac(
	VariableSubMatrixHandler& WorkMat, 
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler& 
ModuleSpardyn::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkVec.ResizeReset(0);

	return WorkVec;
}

/* 
 *                     [v1]              [v1*g1   v1*g2   v1*g3]
 * v(3x1) g^T(1x3) =   [v2] [g1 g2 g3] = [v2*g1   v2*g2   v2*g3] (3x3)
 *                     [v3]              [v3*g1   v3*g2   v3*g3]
 * 
 */
Mat3x3
ModuleSpardyn::MultV1V2T(const Vec3& v, const Vec3& w )
{
/* the one of the constructor of Mat3x3 class
   Mat3x3(const Vec3& v1, const Vec3& v2, const Vec3& v3) {
      pdMat[M11] = v1.pdVec[V1];
      pdMat[M21] = v1.pdVec[V2];
      pdMat[M31] = v1.pdVec[V3];

      pdMat[M12] = v2.pdVec[V1];
      pdMat[M22] = v2.pdVec[V2];
      pdMat[M32] = v2.pdVec[V3];

      pdMat[M13] = v3.pdVec[V1];
      pdMat[M23] = v3.pdVec[V2];
      pdMat[M33] = v3.pdVec[V3];
   };
*/
	return Mat3x3((v*w.dGet(1)), (v*w.dGet(2)), (v*w.dGet(3)) );
}


bool
spardyn_set(void)
{
#ifdef DEBUG
	std::cerr << __FILE__ <<":"<< __LINE__ << ":"<< __PRETTY_FUNCTION__ << std::endl;
#endif

	UserDefinedElemRead *rf = new UDERead<ModuleSpardyn>;

	if (!SetUDE("spardyn", rf)) {
		delete rf;
		return false;
	}

	return true;
}

#ifndef STATIC_MODULES

extern "C" 
{

int
module_init(const char *module_name, void *pdm, void *php)
{
	if (!spardyn_set()) {
		silent_cerr("spardyn: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);
		return -1;
	}

	return 0;
}

} // extern "C"

#endif // ! STATIC_MODULES

