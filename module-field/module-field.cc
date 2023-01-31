/* $Header$ */
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
#include "module-field.h"


class ModuleField
: virtual public Elem, public UserDefinedElem {
public:
	ModuleField(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~ModuleField(void);
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

    StructNode              	*pNode1;
	StructNode					*pNode2;

	Vec3						F1;
	Vec3						M1;
	Vec3						F2;
	Vec3						M2;

	Mat3x3					    KMatrix[2][2];
	Mat3x3 						CMatrix[2][2];
	Vec3						PreLength[2];

	Vec3						e1,e2,e3;

	bool 						bFirst;
    const DataManager*          m_pDM;
};

ModuleField::ModuleField(
	unsigned uLabel, 
    const DofOwner *pDO,
	DataManager* pDM, 
    MBDynParser& HP)
:   Elem(uLabel, flag(0)),
    UserDefinedElem(uLabel, pDO),
    m_pDM(pDM),
	e1(1.,0.,0.),
	e2(0.,1.,0.),
	e3(0.,0.,1.),
	bFirst(true)
{
	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
        "									\n"
        "Module: 	template				\n"
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

	// get node
	pNode2 = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));
	if (!pNode2) {
		silent_cerr("template(" << GetLabel() << "): structural node expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	pNode1 = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));
	if (!pNode1) {
		silent_cerr("template(" << GetLabel() << "): structural node expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// get Kmatrix

	KMatrix[0][0] = Zero3x3;
	KMatrix[0][1] = Zero3x3;
	KMatrix[1][0] = Zero3x3;
	KMatrix[1][1] = Zero3x3;

	if (!HP.IsKeyWord("KMatrix")) {
		silent_cerr("spardyn(" << GetLabel() << "): keyword \"KMatrix\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}	
	for (int i=0; i<6; i++) {
		for ( int j = 0; j < 6; j ++) {
			doublereal temp = HP.GetReal();

			KMatrix[int(i/3)][int(j/3)].Put(int(i%3)+1,int(j%3)+1, temp);

		}
	}
	std::cout<<""<<KMatrix[0][0]<<std::endl;
	std::cout<<""<<KMatrix[0][1]<<std::endl;
	std::cout<<""<<KMatrix[1][0]<<std::endl;
	std::cout<<""<<KMatrix[1][1]<<std::endl;
	// get Cmatrix

	CMatrix[0][0] = Zero3x3;
	CMatrix[0][1] = Zero3x3;
	CMatrix[1][0] = Zero3x3;
	CMatrix[1][1] = Zero3x3;

	if (!HP.IsKeyWord("CMatrix")) {
		silent_cerr("spardyn(" << GetLabel() << "): keyword \"CMatrix\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}	
	for (int i=0; i<6; i++) {
		for ( int j = 0; j < 6; j ++) {
			doublereal temp = HP.GetReal();
			
			CMatrix[int(i/3)][int(j/3)].Put(int(i%3)+1,int(j%3)+1, temp);

		}
	}

	// calculate prelength

	Vec3 x1 = pNode1->GetXCurr();
	Vec3 x2 = pNode2->GetXCurr();
	Mat3x3 R1 = pNode1->GetRCurr();
	Mat3x3 R2 = pNode2->GetRCurr();

	Vec3 r = R1.MulTV(x2-x1);

	doublereal theta1 = std::atan2(-(R1*e3).Dot(R2*e2), (R1*e3).Dot(R2*e3) );
	doublereal theta2 = std::atan2( (R1*e3).Dot(R2*e1), (R1*e3).Dot(R2*e3) );
	doublereal theta3 = std::atan2( (R1*e1).Dot(R2*e2), (R1*e1).Dot(R2*e1) );

	PreLength[0] = r;
	PreLength[1] = Vec3(theta1, theta2, theta3);

	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));

	pDM->GetLogFile() << "template: "
		<< uLabel << " "
		<< std::endl;

}


ModuleField::~ModuleField(void)
{
	// destroy private data
}

// The number of private degrees of freedom is determined by this member function.
// this module use private data of node, so this element has no private data
unsigned int
ModuleField::iGetNumDof(void) const
{
	return 0;
}

// the type of equations of the private degrees of freedom
DofOrder::Order
ModuleField::GetDofType(unsigned int i) const
{
	return DofOrder::DIFFERENTIAL;
}

// this fuction is called befor the integration starts in order to set the initial values for the private dgrees of freedom
void
ModuleField::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

// save the current state of private degree of freedom
// this is needed for implementation of the Output() and dGetPrivData functions
void
ModuleField::Update(const VectorHandler& X,const VectorHandler& XP)
{
	NO_OP;
}

void
ModuleField::AfterPredict(VectorHandler& X,VectorHandler& XP)
{
	bFirst = true;
	Update(X, XP);
}

void
ModuleField::Output(OutputHandler& OH) const
{
	if (bToBeOutput()) {

		if (OH.UseText(OutputHandler::LOADABLE)) {

			OH.Loadable() << GetLabel()
				<< " " << F1
				<< " " << M1
				<< " " << F2
				<< " " << M2
				<< std::endl;
		}
	}
}

// the dimention of subvector and submatrix which contribute to residual and jacobian matrix
void
ModuleField::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 12;
	*piNumCols = 12;
}

SubVectorHandler& 
ModuleField::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	if(bFirst) {
		Vec3 x1 = pNode1->GetXCurr();
		Mat3x3 R1 = pNode1->GetRCurr();
		Vec3 v1 = pNode1->GetVCurr();
		Vec3 w1 = pNode1->GetWCurr();

		Vec3 x2 = pNode2->GetXCurr();
		Mat3x3 R2 = pNode2->GetRCurr();
		Vec3 v2 = pNode2->GetVCurr();
		Vec3 w2 = pNode2->GetWCurr();

		Vec3 r = R1.MulTV(x2-x1);

		doublereal theta1 = std::atan2(-(R1*e3).Dot(R2*e2), (R1*e3).Dot(R2*e3) );
		doublereal theta2 = std::atan2( (R1*e3).Dot(R2*e1), (R1*e3).Dot(R2*e3) );
		doublereal theta3 = std::atan2( (R1*e1).Dot(R2*e2), (R1*e1).Dot(R2*e1) );		

		Vec3 theta(theta1, theta2, theta3);

		std::cout<<""<<r<<std::endl;
		std::cout<<""<<theta<<std::endl;
		Vec3 v = R1.MulTV(v2-v1);
		Vec3 w = R1.MulTV(w2-w1);

		F2 = - (KMatrix[0][0] * (r - PreLength[0])  + KMatrix[0][1] * (theta - PreLength[1]) )  - (CMatrix[0][0] * v + CMatrix[0][1] * w);
		M2 = - (KMatrix[1][0] * (r - PreLength[0])  + KMatrix[1][1] * (theta - PreLength[1]) )  - (CMatrix[1][0] * v + CMatrix[1][1] * w);

		F1 = -F2;
		M1 = -M2 - r.Cross(F2);

		bFirst = false;
	}
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);

	WorkVec.ResizeReset(iNumRows);

	// set indices to WorkVec
	const integer iMomentumIndex1 = pNode1->iGetFirstMomentumIndex();
	const integer iMomentumIndex2 = pNode2->iGetFirstMomentumIndex();
	for(int iCnt = 1; iCnt <=6; iCnt++){
		WorkVec.PutRowIndex(iCnt, 	iMomentumIndex1 + iCnt);
		WorkVec.PutRowIndex(6+iCnt, iMomentumIndex2 + iCnt);
	}

	WorkVec.Put( 1,  pNode1->GetRCurr() * F1  );
	WorkVec.Put( 4,  pNode1->GetRCurr() * M1  );
	WorkVec.Put( 7,  pNode1->GetRCurr() * F2  );
	WorkVec.Put(10,  pNode1->GetRCurr() * M2  );
	
	return WorkVec;
}

VariableSubMatrixHandler& 
ModuleField::AssJac(VariableSubMatrixHandler& WorkMat_,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	/*
	 *                   [ x ]                     [   F(y', y, t) ]
	 * private data, y = [   ],   formulation, r = [               ],  
	 *                   [ g ]                     [   M(y', y, t) ]
	 *  size of above vectors is 6, (translational and rotational displacement)
	 * 
	 *  x,g  are private data of node
	 *            1~3     4~6   
	 *          [ dF/dx  dF/dg]   1~3 
	 *  dr/dy = [             ]
	 *          [dM/dx  dM/dg ]   4~6
	 * 
	 * 			  1~3     4~6      
	 *          [ dF/dx' dF/dg'] 1~3
	 *  dr/dy'= [              ]
	 *          [ dM/dx' dM/dg'] 4~6
	 * 
	 *  size of above matirices is 6x6
	 *  contribution to Jacobian Matrix is  J = -( dr/dy' + dCoef * dr/dy )
	 *
	 * set indices to submatrix WrokMat
	 *    [ x : position index of node + 1~3]  <= private data of node
	 * y =[ g : position index of node + 4~6] 
	 * 
	 *    [ r1 : momentum index of node + 1~3] Force imposed to the node
	 * r =[ r2 : momentum index of node + 4~6] Moment imposed to the node
	 * 
	 * so, r shoud be set momentun index of node
	 *            x   g   (index)     
	 *           [      ] r1
	 * WorkMat = [      ] r2
	 */ 
	WorkMat_.SetNullMatrix();

#if 0
	integer iNumRows = 0;
	integer iNumCols = 0;

	WorkSpaceDim(&iNumRows, &iNumCols);

	FullSubMatrixHandler& WorkMat = WorkMat_.SetFull();

	WorkMat.ResizeReset(iNumRows, iNumCols);

	const integer iPositonIndex = pNode1->iGetFirstPositionIndex();
	const integer iMomentumIndex = pNode1->iGetFirstMomentumIndex();

	for(int iCnt = 1; iCnt<=6; iCnt++){
		WorkMat.PutRowIndex(iCnt, iMomentumIndex + iCnt);
		WorkMat.PutColIndex(iCnt, iPositonIndex  + iCnt);
	}

	Mat3x3 dF_dx, dF_dg, dM_dx, dM_dg;
	Mat3x3 dF_dx_dot, dF_dg_dot, dM_dx_dot, dM_dg_dot;

	/*
	 *calclate jacobian matrices her
	*/

	// set contribution jacobian matrix to WorkMat 
	WorkMat.Put( 1,  1, -dF_dx_dot   - (dF_dx          * dCoef) );
	WorkMat.Put( 1,  4, -dF_dg_dot   - (dF_dg          * dCoef) );

	WorkMat.Put( 4,  1, -dM_dx_dot   - (dM_dx          * dCoef) );
	WorkMat.Put( 4,  4, -dM_dg_dot   - (dM_dg          * dCoef) );
#endif
	return WorkMat_;
}


// This member function is called if the statement "print: dof description;" is specified in the input file.
std::ostream&
ModuleField::DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const
{
	return out;
}

// This member function is called if the statement "print: equation description;" is specified in the input file.
std::ostream&
ModuleField::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{
	return out;
}

//  This member function is called when a bind statement or a element drive is used to access private data of this element.
unsigned int
ModuleField::iGetNumPrivData(void) const
{
	return 0;
}

// The following string values can be specified in a bind statement or in an element drive caller.
unsigned int
ModuleField::iGetPrivDataIdx(const char *s) const
{
	return 0;
}

// This member function is called whenever a bind statement or a element drive is used to access private data.
doublereal
ModuleField::dGetPrivData(unsigned int i) const
{
	return 0.0;
}

// This member function is called if the statement "print: node connection;" is specified in the input file.
int
ModuleField::iGetNumConnectedNodes(void) const
{
	return 0;
}

// This member function is called if the statement "print: node connection;" is specified in the input file.
void
ModuleField::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(0);
}

// This member function is called if "make restart file;" is specified in the input file.
std::ostream&
ModuleField::Restart(std::ostream& out) const
{
	return out << "# ModuleField: not implemented" << std::endl;
}

/*
 * The initial assembly phase is needed only if initial values are provided which are not consistent with the algebraic constraints.
 * Since this element does not use algebraic constraints we do not have to implement
 * iGetInitialNumDof(), InitialWorkSpaceDim(), InitialAssJac(), InitialAssRes() and SetInitialValue().
 */
void
ModuleField::SetInitialValue(VectorHandler& X)
{
	return;
}

unsigned int
ModuleField::iGetInitialNumDof(void) const
{
	return 0;
}

void 
ModuleField::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
ModuleField::InitialAssJac(
	VariableSubMatrixHandler& WorkMat, 
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler& 
ModuleField::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkVec.ResizeReset(0);

	return WorkVec;
}


bool
field_set(void)
{
	UserDefinedElemRead *rf = new UDERead<ModuleField>;

	if (!SetUDE("field", rf)) {
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
	if (!field_set()) {
		silent_cerr("field: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);
		return -1;
	}

	return 0;
}

} // extern "C"

#endif // ! STATIC_MODULES

