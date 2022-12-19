/* $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/modules/module-asynchronous_machine/module-asynchronous_machine.h,v 1.7 2017/01/12 14:47:25 masarati Exp $ */
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
#ifndef ___MODULE_PITCH_CONTROL_H__INCLUDED___
#define ___MODULE_PITCH_CONTROL_H__INCLUDED___


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

extern "C" {

extern int
__FC_DECL__(call_controller)( float *HSS_Spd, float *BlPitch_in, int *NumBl, float *ZTime, float *GenEff,  float *GenTrq, float *ElecPwr, float *BlPitchCom);

}

class PitchControlModule
: virtual public Elem, public UserDefinedElem {
public:
	PitchControlModule(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~PitchControlModule(void);
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
	void AfterConvergence(const VectorHandler& X, 
		const VectorHandler& XP);
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

    StructNode          *pHSSNode;
	StructNode			*pLSSNode;
	StructNode			*pRFrlNode;
	Mat3x3		         Ra;
    Vec3                 Torque;                   

	struct PitchBearing {
		StructNode		*pBLootNode;
		StructNode		*pBottomNode;
	};

	integer				numBl;
	std::vector<PitchBearing>	PitchBr;
     	 
	float				HorWindV;		 
	doublereal 			TimeVSOn; 				// Time to enable active variable speed control
	doublereal			TimePCOn;				// Time to enable active pitch control
	float				PitchAngle[3];			// Pitch angle each blade
	doublereal			PitchVel[3];			// Pitch angular velocity
	float				InitPitch[3];
	float				PitchCom[3];			// Commanded blade pitch angel
	doublereal			PitchTorqe[3];			// Pitch Moment [N-m]
	doublereal			demandAngle;			// (option demand angle)
	float				GenEff;					// generator efficiency
	doublereal			GBRatio;				// Gear Box Ratio
	float 				GenSpeed;				// Generator Speed [rad/s]
	float				GenTorque;				// Generator Torque [N-m]
	float				GenPower;				// Generator Power [W]
	float				InitGenSpd;				// Initial Generator Speed [rad/s]

	doublereal 			PitchSpringConst;
	doublereal			PitchDumperConst;

    doublereal                      CurrTime;
    const DataManager*              m_pDM;
	DriveOwner                      FSF;      
    mutable std::ofstream           out;

	bool 			        bFirstAssRes;
	bool					bFirstAssJac;
	bool					bDemandAngle;

	const Vec3 			e1;
	const Vec3 			e2;
	const Vec3 			e3;

	// jacobian
	Mat3x3 dfb_dgb; 
	Mat3x3 dfb_dgh;
	Mat3x3 dfh_dgb; 
	Mat3x3 dfh_dgh;

	Mat3x3 dfb_dgb_dot; 
	Mat3x3 dfb_dgh_dot;
	Mat3x3 dfh_dgb_dot; 
	Mat3x3 dfh_dgh_dot;

	doublereal printTime;

	void call_discon();
	Mat3x3 MultV1V2T(const Vec3 & v,const Vec3& w);
};

PitchControlModule::PitchControlModule(
	unsigned uLabel, 
    const DofOwner *pDO,
	DataManager* pDM, 
    MBDynParser& HP)
:   Elem(uLabel, flag(0)),
    UserDefinedElem(uLabel, pDO),
    m_pDM(pDM),
    bFirstAssJac(true),
	bFirstAssRes(true),
	e1(1.,0.,0.),
	e2(0.,1.,0.),
	e3(0.,0.,1.),
	printTime(0.0)
{
	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
        "									\n"
        "Module: 	discon				    \n"
        "Author: 	            		    \n"
        "Organization:	kyushu university	\n"
        "		                    	    \n"
        "					                \n"
        "									\n"
        "	All rights reserved				\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	// get HSS node
	pHSSNode = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));
	if (!pHSSNode) {
		silent_cerr("discon(" << GetLabel() << "): structural node expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// get LSS node
	pLSSNode = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));
	if (!pLSSNode) {
		silent_cerr("discon(" << GetLabel() << "): structural node expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// get Rotor Furl node
	pRFrlNode = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));
	if (!pRFrlNode) {
		silent_cerr("discon(" << GetLabel() << "): structural node expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// getting rilative orientation if needed.
    if (HP.IsKeyWord("orientation")) {
		Ra = HP.GetRotRel(ReferenceFrame(pRFrlNode));
	} else {
		Ra = Eye3;
	}

	numBl = HP.GetInt();
	PitchBr.resize(numBl);
	for(int i=0; i<numBl; i++) {
		PitchBr[i].pBLootNode = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));
		PitchBr[i].pBottomNode = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));
	}

	// getting the generator speed and blade pitch, and initialize controller 
	const Vec3& HSS_rotvel = pHSSNode->GetWCurr();						// reference to global coordinate system
	const Vec3& RFrl_rotvel = pRFrlNode->GetWCurr();					// reference to global coordinate system
	const Mat3x3& RFrl_R	= pRFrlNode->GetRCurr()*Ra;	    			// Rotation matrix from RotorFurl system to global system

	GenSpeed = RFrl_R.MulTV(HSS_rotvel - RFrl_rotvel).dGet(1)*GBRatio;

	for(int i=0; i<numBl; i++){
		const Mat3x3& R_blade = PitchBr[i].pBLootNode->GetRCurr();    // Rotation Matrix from BladeRoot fo Global
		const Mat3x3& R_hub	= PitchBr[i].pBottomNode->GetRCurr();	// Rotation Matrix from PitchBrBottom to Global
		PitchAngle[i] = std::atan2( (R_hub * e1).Dot(R_blade * e2) , (R_hub * e1).Dot(R_blade * e1) ) ; 
		InitPitch[i] = PitchAngle[i];
	}

	TimePCOn		   = HP.GetReal();
	GenEff             = HP.GetReal()/100;
	GBRatio            = HP.GetReal();
	PitchSpringConst   = HP.GetReal();
	PitchDumperConst   = HP.GetReal();

	if(HP.IsKeyWord("demand" "angle")) {
		demandAngle = HP.GetReal();
		bDemandAngle = true;
	} else {
		bDemandAngle = false;
	}

	float ZTime = 0.0;
	int NBlades = numBl;
	HorWindV = 0.0;
	InitGenSpd = 0.0;
	// initialize the subroutine call_controller 
	// DISCON stores the generator speed on the first call. If you call it here, the initial speed of 0 will be stored.

	__FC_DECL__(call_controller)( &InitGenSpd, InitPitch, &NBlades, &ZTime,  &GenEff, &GenTorque, &GenPower, PitchCom);

	// getting Force scale factor if needed.
	if (HP.IsKeyWord("Force" "scale" "factor")) {
		FSF.Set(HP.GetDriveCaller());

	} else {
		FSF.Set(new OneDriveCaller);
	}


	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));

	pDM->GetLogFile() << "discon: "
		<< uLabel << " "
		<< std::endl;
	
}


PitchControlModule::~PitchControlModule(void)
{
	// destroy private data
}

// The number of private degrees of freedom is determined by this member function.
unsigned int
PitchControlModule::iGetNumDof(void) const
{
	return 0;
}

// the type of equations of the private degrees of freedom
DofOrder::Order
PitchControlModule::GetDofType(unsigned int i) const
{
	return DofOrder::DIFFERENTIAL;
}

// this fuction is called befor the integration starts in order to set the initial values for the private dgrees of freedom
void
PitchControlModule::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

// save the current state of private degree of freedom
// this is needed for implementation of the Output() and dGetPrivData functions
void
PitchControlModule::Update(const VectorHandler& X,const VectorHandler& XP)
{
	NO_OP;
}

void
PitchControlModule::AfterPredict(VectorHandler& X,VectorHandler& XP)
{
	bFirstAssJac = true;
	bFirstAssRes = true;
	Update(X, XP);
}

void 
PitchControlModule::AfterConvergence(const VectorHandler& X, const VectorHandler& XP)
{
	CurrTime    = m_pDM->dGetTime();
	if(CurrTime >=printTime) {
		std::cout<<" "<<CurrTime<<std::endl;
		if(printTime<=10) {
			printTime += 2.0;
		}else {
			printTime += 10.0;
		}
	}
}
void
PitchControlModule::Output(OutputHandler& OH) const
{
	if (bToBeOutput()) {

		if (OH.UseText(OutputHandler::LOADABLE)) {

			OH.Loadable() << GetLabel()
				<< " " << PitchAngle[0]
				<< " " << PitchAngle[1]
				<< " " << PitchAngle[2]
				<< " " << PitchCom[0]
				<< " " << PitchCom[1]
				<< " " << PitchCom[2]
				<< " " << PitchTorqe[0]
				<< " " << PitchTorqe[1]
				<< " " << PitchTorqe[2]
				<< std::endl;
		}
	}
}

// the dimention of subvector and submatrix which contribute to residual and jacobian matrix
void
PitchControlModule::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = numBl * 6;
	*piNumCols = numBl * 6;
}

void 
PitchControlModule::call_discon() {

	CurrTime    = m_pDM->dGetTime();
	
	//caluclate generator torque
	// Varible Speed Control
	// calculate the generator torque
	const Vec3& HSS_rotvel = pHSSNode->GetWCurr();						// reference to global coordinate system
	const Vec3& RFrl_rotvel = pRFrlNode->GetWCurr();		// reference to global coordinate system
	const Mat3x3& RFrl_R	= pRFrlNode->GetRCurr()*Ra;	    // Rotation matrix from RotorFurl system to global system
	GenSpeed = RFrl_R.MulTV(HSS_rotvel - RFrl_rotvel).dGet(1)*GBRatio;

	for(int i=0; i<numBl; i++){
		const Mat3x3& R_blade = PitchBr[i].pBLootNode->GetRCurr();    // Rotation Matrix from BladeRoot fo Global
		const Mat3x3& R_hub	= PitchBr[i].pBottomNode->GetRCurr();	// Rotation Matrix from PitchBrBottom to Global
		PitchAngle[i] = std::atan2( (R_hub * e1).Dot(R_blade * e2) , (R_hub * e1).Dot(R_blade * e1) ) ; 
	}

	int nbladesf = numBl;
	float c_time = CurrTime;
	__FC_DECL__(call_controller)( &GenSpeed, PitchAngle, &nbladesf, &c_time, &GenEff, &GenTorque, &GenPower, PitchCom);

	if(c_time<TimePCOn){
		for(int i=0; i<numBl; i++){
			PitchCom[i] = InitPitch[i];
		}
	}
}
VariableSubMatrixHandler& 
PitchControlModule::AssJac(VariableSubMatrixHandler& WorkMat_,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	WorkMat_.SetNullMatrix();

#if 0
	if(bFirstAssJac) {

		call_discon();

		for(int iBld = 0; iBld<numBl; iBld++) {
			// getting variables
			const Mat3x3& R_b     = PitchBr[iBld].pBLootNode->GetRCurr();  // Rotation Matrix from BladeRoot fo Global
			const Mat3x3& R_h	  = PitchBr[iBld].pBottomNode->GetRCurr();	// Rotation Matrix from PitchBrBottom to Global
			const Mat3x3& R_b_ref = PitchBr[iBld].pBLootNode->GetRRef();  // Rotation Matrix from BladeRoot fo Global
			const Mat3x3& R_h_ref = PitchBr[iBld].pBottomNode->GetRRef();	// Rotation Matrix from PitchBrBottom to Global

			const Vec3& W_b = PitchBr[iBld].pBLootNode->GetWCurr();    // angular velocity of blade root in global frame
			const Vec3& W_h = PitchBr[iBld].pBottomNode->GetWCurr();     // angular velocity of bearing bottom in global frame
			const Vec3& W_b_ref = PitchBr[iBld].pBLootNode->GetWRef(); 
			const Vec3& W_h_ref = PitchBr[iBld].pBottomNode->GetWRef(); 

			//temporary variable
			doublereal ey = (R_h * e1).Dot(R_b * e2);			// 	blade 2 axis component of hub 1 axis	
			doublereal ex = (R_h * e1).Dot(R_b * e1);           //  blade 1 axis component of hub 1 axis
			const Vec3 dey_dgb = (-Mat3x3(MatCross, R_b_ref*e2)).MulTV(R_h*e1);
			const Vec3 dey_dgh = (-Mat3x3(MatCross, R_h_ref*e1)).MulTV(R_b*e2);
			const Vec3 dex_dgb = (-Mat3x3(MatCross, R_b_ref*e1)).MulTV(R_h*e1);
			const Vec3 dex_dgh = (-Mat3x3(MatCross, R_h_ref*e1)).MulTV(R_b*e1);

			// temporary variable, the tangent value of pitch angle
			doublereal t = ey/ex;
			const Vec3 dt_dgb = dey_dgb*(1./ex) - dex_dgb *(ey/(ex*ex));
			const Vec3 dt_dgh = dey_dgh*(1./ex) - dex_dgh *(ey/(ex*ex));

			// pitch angle
			doublereal theta = std::atan(t);
			const Vec3 dtheta_dgb = dt_dgb * (1./(1.+t*t));
			const Vec3 dtheta_dgh = dt_dgh * (1./(1.+t*t));

			// pitch angular velocity
			const doublereal omega = (R_h.MulTV(W_b - W_h)).Dot(e3); // relative angular velocity
			const Vec3 domega_dgb = ( R_h.MulTM( -Mat3x3(MatCross, W_b_ref) ) ) * (e3);
			const Vec3 domega_dgh = ( R_h.MulTM( Mat3x3(MatCross,(W_b - W_h)) ) ).MulTV(e3)
			                      + ( R_h.MulTM( Mat3x3(MatCross, W_h_ref) ) ) * (e3);

			const Vec3 domega_dgb_dot = R_h.MulTV(e3);
			const Vec3 domega_dgh_dot =-R_h.MulTV(e3);

			// pitch torqe
			PitchTorqe[iBld] = -PitchSpringConst * (PitchCom[iBld] - PitchAngle[iBld]) - PitchDumperConst * omega;
			if(bDemandAngle) {
				PitchTorqe[iBld] = -PitchSpringConst * (demandAngle - PitchAngle[iBld]) - PitchDumperConst * omega;
			} 		

			const Vec3 dTorqe_dgb = dtheta_dgb * PitchSpringConst - domega_dgb * PitchDumperConst;
			const Vec3 dTorqe_dgh = dtheta_dgh * PitchSpringConst - domega_dgh * PitchDumperConst;
			const Vec3 dTorqe_dgb_dot =                           - domega_dgb_dot * PitchDumperConst;
			const Vec3 dTorqe_dgh_dot =                           - domega_dgh_dot * PitchDumperConst;

			/* formulation           private data
			 *  fb = R_h*e3*torqe     gb
			 *  fh =-R_h*e3*torqe     gh
			 */

			// calulate jacobian matirix
			dfb_dgb = MultV1V2T( R_h*e3, dTorqe_dgb);
			dfb_dgh = MultV1V2T( R_h*e3, dTorqe_dgh) - Mat3x3(MatCross, R_h_ref * e3) * PitchTorqe[iBld];
			dfh_dgb = -dfb_dgb;
			dfh_dgh = -dfb_dgh;

			dfb_dgb_dot = MultV1V2T( R_h*e3, dTorqe_dgb_dot);
			dfb_dgh_dot = MultV1V2T( R_h*e3, dTorqe_dgh_dot);
			dfh_dgb_dot = -dfb_dgb_dot;
			dfh_dgh_dot = -dfb_dgh_dot;

			bFirstAssRes = false;
		}
	}

	integer iNumRows = 0;
	integer iNumCols = 0;

	WorkSpaceDim(&iNumRows, &iNumCols);

	FullSubMatrixHandler& WorkMat = WorkMat_.SetFull();

	WorkMat.ResizeReset(iNumRows, iNumCols);

	for(int iBld = 0; iBld<numBl; iBld++) {
		// set indices in WorkMat
		integer iMomentumIndexBLoot = PitchBr[iBld].pBLootNode->iGetFirstMomentumIndex();
		integer iPositionIndexBLoot = PitchBr[iBld].pBLootNode->iGetFirstPositionIndex();
		integer iMomentumIndexPBottom = PitchBr[iBld].pBottomNode->iGetFirstMomentumIndex();
		integer iPositionIndexPBottom = PitchBr[iBld].pBottomNode->iGetFirstPositionIndex();

		for(int iCnt = 1; iCnt <= 3; iCnt++) {
			WorkMat.PutRowIndex( iBld*6 +     iCnt, iMomentumIndexBLoot   + 3 + iCnt);
			WorkMat.PutColIndex( iBld*6 +     iCnt, iPositionIndexBLoot   + 3 + iCnt);
			WorkMat.PutRowIndex( iBld*6 + 3 + iCnt, iMomentumIndexPBottom + 3 + iCnt);
			WorkMat.PutColIndex( iBld*6 + 3 + iCnt, iPositionIndexPBottom + 3 + iCnt);
		}


		WorkMat.Put( iBld*6 + 1, iBld*6 +  1, -dfb_dgb_dot   - (dfb_dgb * dCoef) );
		WorkMat.Put( iBld*6 + 1, iBld*6 +  4, -dfb_dgh_dot   - (dfb_dgh * dCoef) );
		WorkMat.Put( iBld*6 + 4, iBld*6 +  1, -dfh_dgb_dot   - (dfh_dgb * dCoef) );
		WorkMat.Put( iBld*6 + 4, iBld*6 +  4, -dfh_dgh_dot   - (dfh_dgh * dCoef) );

	}

#endif
	return WorkMat_;
}

SubVectorHandler& 
PitchControlModule::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	if(bFirstAssRes) {
		// calulate demand pitch angle for calling DISCON
		call_discon();

		// calculate pitch moment
		for(int i=0; i<numBl; i++){

			PitchTorqe[i] = -PitchSpringConst * PitchCom[i];
			if(bDemandAngle) {
				PitchTorqe[i] = -PitchSpringConst *demandAngle;
			} 
		}
		bFirstAssRes = false;
	}

	doublereal dFSF = FSF.dGet();

	integer iNumRows = 0;
	integer iNumCols = 0;

	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	// setting indeces in WorkVec
	for(int iBld = 0; iBld < numBl; iBld++) {

		integer iMomentumIndexBLoot = PitchBr[iBld].pBLootNode->iGetFirstMomentumIndex();
		integer iMomentumIndexPBottom = PitchBr[iBld].pBottomNode->iGetFirstMomentumIndex();

		for(int iCnt = 1; iCnt <= 3; iCnt++) {

			WorkVec.PutRowIndex( iBld*6 +     iCnt, iMomentumIndexBLoot   + 3 + iCnt);
			WorkVec.PutRowIndex( iBld*6 + 3 + iCnt, iMomentumIndexPBottom + 3 + iCnt);
		}
	}

	// setting force in WorkVec
	for(int iBld = 0; iBld < numBl; iBld++) {
		const Vec3 Torque = PitchBr[iBld].pBottomNode->GetRCurr()*Vec3(0.,0.,PitchTorqe[iBld])*dFSF;

		WorkVec.Put( iBld*6 + 1 ,  Torque);
		WorkVec.Put( iBld*6 + 4 , -Torque);
	}
	
	return WorkVec;
}


// This member function is called if the statement "print: dof description;" is specified in the input file.
std::ostream&
PitchControlModule::DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const
{
	return out;
}

// This member function is called if the statement "print: equation description;" is specified in the input file.
std::ostream&
PitchControlModule::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{
	return out;
}

//  This member function is called when a bind statement or a element drive is used to access private data of this element.
unsigned int
PitchControlModule::iGetNumPrivData(void) const
{
	return 0;
}

// The following string values can be specified in a bind statement or in an element drive caller.
unsigned int
PitchControlModule::iGetPrivDataIdx(const char *s) const
{
	return 0;
}

// This member function is called whenever a bind statement or a element drive is used to access private data.
doublereal
PitchControlModule::dGetPrivData(unsigned int i) const
{
	return 0.0;
}

// This member function is called if the statement "print: node connection;" is specified in the input file.
int
PitchControlModule::iGetNumConnectedNodes(void) const
{
	return 0;
}

// This member function is called if the statement "print: node connection;" is specified in the input file.
void
PitchControlModule::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(0);
}

// This member function is called if "make restart file;" is specified in the input file.
std::ostream&
PitchControlModule::Restart(std::ostream& out) const
{
	return out << "# PitchControlModule: not implemented" << std::endl;
}

/*
 * The initial assembly phase is needed only if initial values are provided which are not consistent with the algebraic constraints.
 * Since this element does not use algebraic constraints we do not have to implement
 * iGetInitialNumDof(), InitialWorkSpaceDim(), InitialAssJac(), InitialAssRes() and SetInitialValue().
 */
void
PitchControlModule::SetInitialValue(VectorHandler& X)
{
	return;
}

unsigned int
PitchControlModule::iGetInitialNumDof(void) const
{
	return 0;
}

void 
PitchControlModule::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
PitchControlModule::InitialAssJac(
	VariableSubMatrixHandler& WorkMat, 
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkMat.SetNullMatrix();

	return WorkMat;
}

/* 
 *                     [v1]              [v1*g1   v1*g2   v1*g3]
 * v(3x1) g^T(1x3) =   [v2] [g1 g2 g3] = [v2*g1   v2*g2   v2*g3] (3x3)
 *                     [v3]              [v3*g1   v3*g2   v3*g3]
 * 
 */
Mat3x3
PitchControlModule::MultV1V2T(const Vec3& v, const Vec3& w )
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

SubVectorHandler& 
PitchControlModule::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkVec.ResizeReset(0);

	return WorkVec;
}

extern bool pitchControl_set(void);

#endif // ___MODULE_DISCON_H__INCLUDED___
