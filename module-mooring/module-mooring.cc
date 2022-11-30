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
#include "module-mooring.h"

#include "moor3d.h"
#include "moorsys.h"


class ModuleMooring
: virtual public Elem, public UserDefinedElem {
public:
	ModuleMooring(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~ModuleMooring(void);
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

    StructNode              	*pNode;
    Vec3                    	Force;
    Vec3                    	Moment;       
	Mat3x3		            	Ra;            

    enum class MoorType {
        moor3d,
        moorsys,
    };

    MoorType 					moortype;
    Moor3d                     	moor3d;
	Moorsysf					moorsys;

    doublereal                  CurrTime;
    const DataManager*          m_pDM;
	DriveOwner                  FSF;      
    mutable std::ofstream       out;

	bool 					    bFirst;

};

ModuleMooring::ModuleMooring(
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
        "Module: 	mooring				    \n"
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
	pNode = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));
	if (!pNode) {
		silent_cerr("mooring(" << GetLabel() << "): structural node expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// getting rilative orientation if needed.
    if (HP.IsKeyWord("orientation")) {
		Ra = HP.GetRotRel(ReferenceFrame(pNode));
	} else {
		Ra = Eye3;
	}

	// getting Force scale factor if needed.
	if (HP.IsKeyWord("Force" "scale" "factor")) {
		FSF.Set(HP.GetDriveCaller());

	} else {
		FSF.Set(new OneDriveCaller);
	}

    if (HP.IsKeyWord("mooring" "program" "name")) {
        const char *moorname = HP.GetString();
        if(moorname == "moorsys") {
			moortype = MoorType::moorsys;
		}else if(moorname == "moor3d") {
			moortype = MoorType::moor3d;
		}
    } else {
		moortype = MoorType::moorsys; //default
	}

	std::string  ipFileNameMoorsys;
	std::string  otFileNameMoorsys; 
	std::string  dout = ".out";

    if (HP.IsKeyWord("mooring" "input" "file" "name")) {
        const char *ifname = HP.GetFileName();
        ipFileNameMoorsys = std::string(ifname);
		std::cout<<ipFileNameMoorsys<<std::endl;
		otFileNameMoorsys = std::string(ifname).erase(ipFileNameMoorsys.find(".")) + dout;
    }

    switch (moortype) {

    	case MoorType::moorsys:
    	    moorsys.input(ipFileNameMoorsys, otFileNameMoorsys);
    	    break;   

    	case MoorType::moor3d:
    	    moor3d.input(ipFileNameMoorsys);
    	    break; 
	
    	default:
    	    break;
    }

	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));

	pDM->GetLogFile() << "mooring: "
		<< uLabel << " "
		<< std::endl;

}


ModuleMooring::~ModuleMooring(void)
{
	// destroy private data
}

// The number of private degrees of freedom is determined by this member function.
unsigned int
ModuleMooring::iGetNumDof(void) const
{
	return 0;
}

// the type of equations of the private degrees of freedom
DofOrder::Order
ModuleMooring::GetDofType(unsigned int i) const
{
	return DofOrder::DIFFERENTIAL;
}

// this fuction is called befor the integration starts in order to set the initial values for the private dgrees of freedom
void
ModuleMooring::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

// save the current state of private degree of freedom
// this is needed for implementation of the Output() and dGetPrivData functions
void
ModuleMooring::Update(const VectorHandler& X,const VectorHandler& XP)
{
	NO_OP;
}

void
ModuleMooring::AfterPredict(VectorHandler& X,VectorHandler& XP)
{
	bFirst = true;
	Update(X, XP);
}

void
ModuleMooring::Output(OutputHandler& OH) const
{
	if (bToBeOutput()) {

		if (OH.UseText(OutputHandler::LOADABLE)) {

			OH.Loadable() << GetLabel()
				<< " " << Force
				<< " " << Moment
				<< std::endl;
		}
	}
}

// the dimention of subvector and submatrix which contribute to residual and jacobian matrix
void
ModuleMooring::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 6;
	*piNumCols = 6;
}

VariableSubMatrixHandler& 
ModuleMooring::AssJac(VariableSubMatrixHandler& WorkMat_,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	/*
	 *                   [ x  ]                     [   F(y', y) ]
	 * private data, y = [ g  ],   formulation, r = [   M(y', y) ],  
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

	const integer iPositonIndex = pNode->iGetFirstPositionIndex();
	const integer iMomentumIndex = pNode->iGetFirstMomentumIndex();

	/*  
	 * set indices to submatrix WrokMat
	 *    [ x : position index of node + 1~3]  <= private data of node
	 * y =[ g : position index of node + 4~6] 
	 * 
	 *    [ r1 : momentum index of node + 1~3] Force imposed to the node
	 * r =[ r2 : momentum index of node + 4~6] Moment imposed to the node
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
ModuleMooring::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	CurrTime    = m_pDM->dGetTime();

    // force scale factor to "ramp up" loads
	doublereal dFSF = FSF.dGet();

	if(bFirst) {
		Force = Zero3;
		Moment = Zero3;

	    const Vec3   &NodeX   = pNode->GetXCurr();		                	
	    const Mat3x3 &NodeR   = pNode->GetRCurr()*Ra; 

		switch (moortype) {
		case MoorType::moorsys:
			moorsys.set_position(CurrTime, NodeX, NodeR);
			Force = moorsys.get_force()*dFSF;
			Moment = moorsys.get_moment()*dFSF;
			break;
		
		case MoorType::moor3d:
        	moor3d.set_position(NodeX, NodeR);
			Force = moor3d.get_force()*dFSF;
			Moment = moor3d.get_moment()*dFSF;
			break;
		default:
			break;
		}
		bFirst = false;
	}
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);

	WorkVec.ResizeReset(iNumRows);

	const integer iMomentumIndex = pNode->iGetFirstMomentumIndex();

	for(int iCnt = 1; iCnt <=6; iCnt++){
		WorkVec.PutRowIndex(iCnt, 	iMomentumIndex + iCnt);
	}

	WorkVec.Put( 1,  Force   );
	WorkVec.Put( 4,  Moment  );
	
	return WorkVec;
}


// This member function is called if the statement "print: dof description;" is specified in the input file.
std::ostream&
ModuleMooring::DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const
{
	return out;
}

// This member function is called if the statement "print: equation description;" is specified in the input file.
std::ostream&
ModuleMooring::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{
	return out;
}

//  This member function is called when a bind statement or a element drive is used to access private data of this element.
unsigned int
ModuleMooring::iGetNumPrivData(void) const
{
	return 0;
}

// The following string values can be specified in a bind statement or in an element drive caller.
unsigned int
ModuleMooring::iGetPrivDataIdx(const char *s) const
{
	return 0;
}

// This member function is called whenever a bind statement or a element drive is used to access private data.
doublereal
ModuleMooring::dGetPrivData(unsigned int i) const
{
	return 0.0;
}

// This member function is called if the statement "print: node connection;" is specified in the input file.
int
ModuleMooring::iGetNumConnectedNodes(void) const
{
	return 0;
}

// This member function is called if the statement "print: node connection;" is specified in the input file.
void
ModuleMooring::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(0);
}

// This member function is called if "make restart file;" is specified in the input file.
std::ostream&
ModuleMooring::Restart(std::ostream& out) const
{
	return out << "# ModuleMooring: not implemented" << std::endl;
}

/*
 * The initial assembly phase is needed only if initial values are provided which are not consistent with the algebraic constraints.
 * Since this element does not use algebraic constraints we do not have to implement
 * iGetInitialNumDof(), InitialWorkSpaceDim(), InitialAssJac(), InitialAssRes() and SetInitialValue().
 */
void
ModuleMooring::SetInitialValue(VectorHandler& X)
{
	return;
}

unsigned int
ModuleMooring::iGetInitialNumDof(void) const
{
	return 0;
}

void 
ModuleMooring::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
ModuleMooring::InitialAssJac(
	VariableSubMatrixHandler& WorkMat, 
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler& 
ModuleMooring::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkVec.ResizeReset(0);

	return WorkVec;
}


bool
mooring_set(void)
{
#ifdef DEBUG
	std::cerr << __FILE__ <<":"<< __LINE__ << ":"<< __PRETTY_FUNCTION__ << std::endl;
#endif

	UserDefinedElemRead *rf = new UDERead<ModuleMooring>;

	if (!SetUDE("mooring", rf)) {
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
	if (!mooring_set()) {
		silent_cerr("mooring: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);
		return -1;
	}

	return 0;
}

} // extern "C"

#endif // ! STATIC_MODULES

