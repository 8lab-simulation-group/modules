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
#ifndef ___MODULE_DISCON_H__INCLUDED___
#define ___MODULE_DISCON_H__INCLUDED___


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

class DisconModule
: virtual public Elem, public UserDefinedElem {
public:
	DisconModule(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~DisconModule(void);
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
	float				InitPitch[3];
	float				PitchCom[3];			// Commanded blade pitch angel
	doublereal			PitchMoment[3];			// Pitch Moment [N-m]
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

	bool 			        bFirst;

};

DisconModule::DisconModule(
	unsigned uLabel, 
    const DofOwner *pDO,
	DataManager* pDM, 
    MBDynParser& HP)
:   Elem(uLabel, flag(0)),
    UserDefinedElem(uLabel, pDO),
    m_pDM(pDM),
    bFirst(true)
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
		const Mat3x3& BlRoot_R	= PitchBr[i].pBLootNode->GetRCurr();		// Rotation Matrix from BladeRoot fo Global
		const Mat3x3& PitchBr_R	= PitchBr[i].pBottomNode->GetRCurr();	// Rotation Matrix from PitchBrBottom to Global
		Mat3x3 relativeR		= BlRoot_R.MulTM(PitchBr_R);			// Relative rotation Matrix from Bearing bottom to blade root
		// Rotation angle of blade root around z-axis with respect to bearing bottom
		// However, the pitch angle rotates clockwise about z-axis, so multiply -1
		PitchAngle[i] = -1*atan(relativeR.dGet(1,2)/relativeR.dGet(1,1)); 
		InitPitch[i] = PitchAngle[i];
	}

	TimeVSOn           = HP.GetReal();
	TimePCOn		   = HP.GetReal();
	GenEff             = HP.GetReal();
	GBRatio            = HP.GetReal();
	InitGenSpd		   = HP.GetReal()*GBRatio;
	PitchSpringConst   = HP.GetReal();
	PitchDumperConst   = HP.GetReal();

	if(HP.IsKeyWord("demand angle")) {
		demandAngle = HP.GetReal();
	} else {
		demandAngle = -1;
	}

	float ZTime = 0.0;
	int NBlades = numBl;
	HorWindV = 0.0;
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


DisconModule::~DisconModule(void)
{
	// destroy private data
}

// The number of private degrees of freedom is determined by this member function.
unsigned int
DisconModule::iGetNumDof(void) const
{
	return 0;
}

// the type of equations of the private degrees of freedom
DofOrder::Order
DisconModule::GetDofType(unsigned int i) const
{
	return DofOrder::DIFFERENTIAL;
}

// this fuction is called befor the integration starts in order to set the initial values for the private dgrees of freedom
void
DisconModule::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

// save the current state of private degree of freedom
// this is needed for implementation of the Output() and dGetPrivData functions
void
DisconModule::Update(const VectorHandler& X,const VectorHandler& XP)
{
	NO_OP;
}

void
DisconModule::AfterPredict(VectorHandler& X,VectorHandler& XP)
{
	bFirst = true;
	Update(X, XP);
}

void
DisconModule::Output(OutputHandler& OH) const
{
	if (bToBeOutput()) {

		if (OH.UseText(OutputHandler::LOADABLE)) {

			OH.Loadable() << GetLabel()
				<< " " << GenSpeed
				<< " " << GenTorque
				<< " " << GenPower
				<< " " << PitchAngle[0]
				<< " " << PitchAngle[1]
				<< " " << PitchAngle[2]
				<< " " << PitchCom[0]
				<< " " << PitchCom[1]
				<< " " << PitchCom[2]
				<< " " << PitchMoment[0]
				<< " " << PitchMoment[1]
				<< " " << PitchMoment[2]
				<< std::endl;
		}
	}
}

// the dimention of subvector and submatrix which contribute to residual and jacobian matrix
void
DisconModule::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 6 + numBl * 6;
	*piNumCols = 6 + numBl * 6;
}

VariableSubMatrixHandler& 
DisconModule::AssJac(VariableSubMatrixHandler& WorkMat_,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	/*
	 *                   [ g1  ]                     [   M(y', y) ]
	 * private data, y = [ g2  ],   formulation, r = [  -M(y', y) ],  
	 *
	 *  x,g,are private data of node
	 *  v,w are private data of this element, and stored translational and rotatioal velocity
	 *  derivation private data XPrimeCurr is store translational and rotatinal acceleration
	 *            1        4   
	 *          [ dF/dx  dF/dg]   1 
	 *  dr/dy = [ dM/dx  dM/dg]   4
	 * 
	 * 			  1        4      
	 *          [ dF/dx' dF/dg'] 1
	 *  dr/dy'= [ dM/dx' dM/dg'] 4
	 * 
	 *  contribution to Jacobian Matrix is  J = -( dr/dy' + dCoef * dr/dy )
	 */

	WorkMat_.SetNullMatrix();

#if 0
	integer iNumRows = 0;
	integer iNumCols = 0;

	WorkSpaceDim(&iNumRows, &iNumCols);

	FullSubMatrixHandler& WorkMat = WorkMat_.SetFull();

	WorkMat.ResizeReset(iNumRows, iNumCols);

	const integer iPositonIndex = pHSSNode->iGetFirstPositionIndex();
	const integer iMomentumIndex = pHSSNode->iGetFirstMomentumIndex();

	/*  
	 * set indices to submatrix WrokMat
	 *    [ x : position index of node + 1~3]  <= private data of node
	 * y =[ g : position index of node + 4~6] 
	 * 
	 *    [ r1 : momentum index of node + 1~3] Force imposed to the node
	 * r =[ r2 : momentum index of node + 4~6] Torque imposed to the node
	 * 
	 *            x   g   (index)     
	 *           [      ] r1
	 * WorkMat = [      ] r2
	 */ 

	for(int iCnt = 1; iCnt<=6; iCnt++){
		WorkMat.PutRowIndex(iCnt, iMomentumIndex + iCnt);
		WorkMat.PutColIndex(iCnt, iPositonIndex  + iCnt);
	}

	WorkMat.Put( 1,  1, -dF_dxP   - (dF_dx          * dCoef) );
	WorkMat.Put( 1,  4, -dF_dgP   - (dF_dg          * dCoef) );

	WorkMat.Put( 4,  1, -dM_dxP   - (dM_dx          * dCoef) );
	WorkMat.Put( 4,  4, -dM_dgP   - (dM_dg          * dCoef) );
#endif

	return WorkMat_;
}

SubVectorHandler& 
DisconModule::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	CurrTime    = m_pDM->dGetTime();

    // force scale factor to "ramp up" loads
	doublereal dFSF = FSF.dGet();

	if(bFirst) {

		//caluclate generator torque
		// Varible Speed Control
		// calculate the generator torque
		const Vec3& HSS_rotvel = pHSSNode->GetWCurr();						// reference to global coordinate system
		const Vec3& RFrl_rotvel = pRFrlNode->GetWCurr();		// reference to global coordinate system
		const Mat3x3& RFrl_R	= pRFrlNode->GetRCurr()*Ra;	    // Rotation matrix from RotorFurl system to global system
		GenSpeed = RFrl_R.MulTV(HSS_rotvel - RFrl_rotvel).dGet(1)*GBRatio;

		for(int i=0; i<numBl; i++){
			const Mat3x3& BlRoot_R	= PitchBr[i].pBLootNode->GetRCurr();		// Rotation Matrix from BladeRoot fo Global
			const Mat3x3& PitchBr_R	= PitchBr[i].pBottomNode->GetRCurr();	// Rotation Matrix from PitchBrBottom to Global
			Mat3x3 relativeR		= BlRoot_R.MulTM(PitchBr_R);			// Relative rotation Matrix from Bearing bottom to blade root
			// Rotation angle of blade root around z-axis with respect to bearing bottom
			// However, the pitch angle rotates clockwise about z-axis, so multiply -1
			PitchAngle[i] = -1*atan(relativeR.dGet(1,2)/relativeR.dGet(1,1)); 
		}

		int nbladesf = numBl;
		float c_time = CurrTime;
		__FC_DECL__(call_controller)( &GenSpeed, PitchAngle, &nbladesf, &c_time, &GenEff, &GenTorque, &GenPower, PitchCom);

		if(c_time<TimePCOn){
			for(int i=0; i<numBl; i++){
				PitchCom[i] = InitPitch[i];
			}
		}

		doublereal dFSF = FSF.dGet();
		for(int i= 0; i<numBl; i++){
			PitchMoment[i] = -PitchSpringConst*(PitchCom[i] - PitchAngle[i]);
			if(demandAngle!=-1) {
				PitchMoment[i] = -PitchSpringConst*(demandAngle - PitchAngle[i]);
			}
		}

	}

	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);

	WorkVec.ResizeReset(iNumRows);

	integer iMomentumIndexLSS = pHSSNode->iGetFirstMomentumIndex();
	integer iMomentumIndexRFril = pRFrlNode->iGetFirstMomentumIndex();

	for(int iCnt = 1; iCnt <=3; iCnt++){
		WorkVec.PutRowIndex(iCnt, 	iMomentumIndexLSS   + 3 + iCnt);
		WorkVec.PutRowIndex(iCnt+3, iMomentumIndexRFril + 3 + iCnt);
	}

	for(int iBld = 0; iBld < numBl; iBld++) {

		integer iMomentumIndexBLoot = PitchBr[iBld].pBLootNode->iGetFirstMomentumIndex();
		integer iMomentumIndexPBottom = PitchBr[iBld].pBottomNode->iGetFirstMomentumIndex();

		for(int iCnt = 1; iCnt <= 3; iCnt++) {

			WorkVec.PutRowIndex(6 + iBld*6 + iCnt, iMomentumIndexBLoot + 3 + iCnt);
			WorkVec.PutRowIndex(6 + iBld*6 + 3 + iCnt, iMomentumIndexPBottom + 3 + iCnt);
		}
	}

	const Vec3 Torque_HSS = Vec3(-GenTorque*GBRatio, 0., 0.)*dFSF;
	const Vec3 Torque_RFrl = Vec3(GenTorque*GBRatio, 0., 0.)*dFSF;

	WorkVec.Put( 1,  pHSSNode->GetRCurr() * Torque_HSS  );
	WorkVec.Put( 4,  pRFrlNode->GetRCurr()*Ra*Torque_RFrl );

	for(int iBld = 0; iBld < numBl; iBld++) {
		const Vec3 Torque_BLoot = Vec3(0.,0.,PitchMoment[iBld])*dFSF;

		WorkVec.Put( 6 + iBld*6 + 1 , PitchBr[iBld].pBLootNode->GetRCurr()*Torque_BLoot);
		WorkVec.Put( 6 + iBld*6 + 4 , PitchBr[iBld].pBottomNode->GetRCurr()*(-Torque_BLoot));
	}
	
	return WorkVec;
}


// This member function is called if the statement "print: dof description;" is specified in the input file.
std::ostream&
DisconModule::DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const
{
	return out;
}

// This member function is called if the statement "print: equation description;" is specified in the input file.
std::ostream&
DisconModule::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{
	return out;
}

//  This member function is called when a bind statement or a element drive is used to access private data of this element.
unsigned int
DisconModule::iGetNumPrivData(void) const
{
	return 0;
}

// The following string values can be specified in a bind statement or in an element drive caller.
unsigned int
DisconModule::iGetPrivDataIdx(const char *s) const
{
	return 0;
}

// This member function is called whenever a bind statement or a element drive is used to access private data.
doublereal
DisconModule::dGetPrivData(unsigned int i) const
{
	return 0.0;
}

// This member function is called if the statement "print: node connection;" is specified in the input file.
int
DisconModule::iGetNumConnectedNodes(void) const
{
	return 0;
}

// This member function is called if the statement "print: node connection;" is specified in the input file.
void
DisconModule::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(0);
}

// This member function is called if "make restart file;" is specified in the input file.
std::ostream&
DisconModule::Restart(std::ostream& out) const
{
	return out << "# DisconModule: not implemented" << std::endl;
}

/*
 * The initial assembly phase is needed only if initial values are provided which are not consistent with the algebraic constraints.
 * Since this element does not use algebraic constraints we do not have to implement
 * iGetInitialNumDof(), InitialWorkSpaceDim(), InitialAssJac(), InitialAssRes() and SetInitialValue().
 */
void
DisconModule::SetInitialValue(VectorHandler& X)
{
	return;
}

unsigned int
DisconModule::iGetInitialNumDof(void) const
{
	return 0;
}

void 
DisconModule::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
DisconModule::InitialAssJac(
	VariableSubMatrixHandler& WorkMat, 
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler& 
DisconModule::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkVec.ResizeReset(0);

	return WorkVec;
}

extern bool discon_set(void);

#endif // ___MODULE_DISCON_H__INCLUDED___
