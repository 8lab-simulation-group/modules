/* 
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

#include "mbconfig.h" 		/* This goes first in every *.c,*.cc file */

#include "Rot.hh"

#include "dataman.h"
#include "userelem.h"
#include "drive_.h"
#include "iostream"

// uncomment to enable debug output
// #define MODULE_AERODYN_DEBUG

// keep this consistent with AeroDyn build
// #define USE_DOUBLE_PRECISION
#define USE_SINGLE_PRECISION

#include "NREL_AeroDyn.h"
#include "module-discon.h"

class AeroDynModule
: virtual public Elem,
public UserDefinedElem
{
private:

	struct AeroNode {
		StructNode				*pNode;				// Aero Node
		Vec3					f;					// offset of the aero point wrt./ the node,
													// constant in the node's reference frame
		Mat3x3					Ra;					// aerodynamic orientation of the aero point
													// wrt./ the node

		Vec3					Force;				// acting force
		Vec3					Moment;				// acting moment
	};
	std::vector<AeroNode>		nodes;				// nodes into Aerodynamic element

	std::vector<StructNode*>	pBladeRoot;		    // Blade root nodes
	std::vector<StructNode*>	pConed;			    // Coned nodes
	StructNode					*pNacelle;          // nacelle node
	StructNode					*pHub;              // hub node
	StructNode					*pLSSTeet;			// TeeterPin
	StructNode					*pShaftCS;			// ShaftCS
	StructNode  				*pRotorFurl;		// RotorFurl node
	StructNode  				*pTower;			// Tower node
	StructNode					*pPlatfromRef;		// Platform Reference

	integer						NumBlades;			// the number of blades.
	integer						NumElements;		// the number of elements per blade.

	std::string ofname;
	mutable std::ofstream out;

	/*
	 * Total aerodynamic data
	 */
	Vec3			TF;			// Total Force on the rotor in the absolute frame
	Vec3			TM;			// Total Moment on the rotor in the absolute frame.
	Vec3			TF_h;		// Total Force on the rotor in the hub frame
	Vec3			TM_h;		// Total Moment on the rotor in the hub frame.
	doublereal      Thrust;   // Rotor Thrust.
	doublereal      Torque;   // Rotor Torque.
	doublereal      Rotor_speed;   // Rotor angular velocity.
	doublereal 		Wind_speed;		// horizonal wind speed at hub height.



	F_INTEGER	c_elem;  // use to identify the current element in AeroDyn!
	F_INTEGER	c_blade; // use to identify the current blade in AeroDyn!
	F_REAL		rlocal;  // use to identify the current element position 
	F_REAL		bladeLength;// use to identify the blade length. 

	F_REAL      c_time;		// current time
	
	bool        bFirst;
	DriveOwner	Time;		// time drive
	doublereal	dOldTime;	// old time
	doublereal  dCurTime;   // current time
	F_REAL      dDT;		// time step

	DriveOwner  FSF;

	F_REAL Position[3];
	F_REAL TranslationVel[3];
	F_REAL RotationVel[3];
	F_REAL Orientation[3][3];

public:
	AeroDynModule(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	~AeroDynModule(void);

	unsigned int iGetNumDof(void) const;
	DofOrder::Order GetDofType(unsigned int i) const;
	void Output(OutputHandler& OH) const;
	std::ostream& Restart(std::ostream& out) const;
	void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
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
	void AfterPredict(VectorHandler& X, VectorHandler& XP);
	void AfterConvergence(const VectorHandler& X, 
		const VectorHandler& XP);
	unsigned int iGetInitialNumDof(void) const;
	void InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
   	VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat, 
		const VectorHandler& XCurr);
   	SubVectorHandler& 
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);
	void SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints* h = 0);
	unsigned int iGetNumPrivData(void) const;
	int GetNumConnectedNodes(void) const;
	void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
	doublereal dGetWindSpeed(void) const;

	private:
	void SetAeroMarkers();
	void SetADInitialInterfaceComponent();
	void SetADInterfaceComponent();
	void ConvertNodeToMarker(const StructNode *x, const StructNode *R, const StructNode *v, const StructNode *w);

};

// static handler for AeroDyn callbacks
// only one element per simulation can be active (AeroDyn limitation)
static AeroDynModule *module_aerodyn;


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

	// read blade length
	if (!HP.IsKeyWord("blade" "length")) {
		silent_cerr("aerodyn module(" << GetLabel() << "): keyword \"blade length\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	bladeLength = HP.GetReal();

	// read number of blades
	if (!HP.IsKeyWord("blades")) {
		silent_cerr("aerodyn module(" << GetLabel() << "): keyword \"blades\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	NumBlades = HP.GetInt();

	// number of elements per blade
	if (!HP.IsKeyWord("elements")) {
		silent_cerr("aerodyn module(" << GetLabel() << "): keyword \"elements\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	NumElements = HP.GetInt();

	// read nacelle node 
	if (!HP.IsKeyWord("NacelleCS")) {
		silent_cerr("aerodyn module(" << GetLabel() << "): keyword \"NacelleCS\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	pNacelle = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));
	if (!pNacelle) {
		silent_cerr("aerodyn module(" << GetLabel() << "): structural node expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// read hub node 
	if (!HP.IsKeyWord("HabCS")) {
		silent_cerr("aerodyn module(" << GetLabel() << "): keyword \"HubCS\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	pHub = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));
	if (!pHub) {
		silent_cerr("aerodyn module(" << GetLabel() << "): structural node expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// read teeterpin node 
	if (!HP.IsKeyWord("LSSTeetPin")) {
		silent_cerr("aerodyn module(" << GetLabel() << "): keyword \"LSSTeetPin\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	pLSSTeet = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));
	if (!pLSSTeet) {
		silent_cerr("aerodyn module(" << GetLabel() << "): structural node expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// read shaft node 
	if (!HP.IsKeyWord("ShaftCS")) {
		silent_cerr("aerodyn module(" << GetLabel() << "): keyword \"ShaftSC\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	pShaftCS = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));
	if (!pShaftCS) {
		silent_cerr("aerodyn module(" << GetLabel() << "): structural node expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// read rotorfurl node 
	if (!HP.IsKeyWord("RotorFurlRef")) {
		silent_cerr("aerodyn module(" << GetLabel() << "): keyword \"RotorFurlRef\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	pRotorFurl = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));
	if (!pRotorFurl) {
		silent_cerr("aerodyn module(" << GetLabel() << "): structural node expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	
	// read Tower node 
	if (!HP.IsKeyWord("TowerBaseCS")) {
		silent_cerr("aerodyn module(" << GetLabel() << "): keyword \"TowerBaseCS\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	pTower = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));
	if (!pTower) {
		silent_cerr("aerodyn module(" << GetLabel() << "): structural node expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// read platfrom reference node 
	if (!HP.IsKeyWord("PlatformRef")) {
		silent_cerr("aerodyn module(" << GetLabel() << "): keyword \"PlatformRef\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	pPlatfromRef = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));
	if (!pPlatfromRef) {
		silent_cerr("aerodyn module(" << GetLabel() << "): structural node expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// read blade root node 
	if (!HP.IsKeyWord("BladeKCS")) {
		silent_cerr("aerodyn module(" << GetLabel() << "): keyword \"BladeKCS\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	pBladeRoot.resize(NumBlades);
	for (int iBld = 0; iBld < NumBlades; iBld++){
		pBladeRoot[iBld] = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));
		if (!pBladeRoot[iBld]) {
			silent_cerr("aerodyn module(" << GetLabel() << "): structural node expected at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	// read coned node 
	if (!HP.IsKeyWord("ConedKCS")) {
		silent_cerr("aerodyn module(" << GetLabel() << "): keyword \"ConedKCS\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	pConed.resize(NumBlades);
	for (int iBld = 0; iBld < NumBlades; iBld++){
		pConed[iBld] = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));
		if (!pConed[iBld]) {
			silent_cerr("aerodyn module(" << GetLabel() << "): structural node expected at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}
	
	// read aero nodes
	if (!HP.IsKeyWord("Aero" "Point")) {
		silent_cerr("aerodyn module(" << GetLabel() << "): keyword \"Aero Point\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
    nodes.resize(NumBlades*NumElements);
	for (int iElem = 0; iElem < NumBlades*NumElements; iElem++) {

		nodes[iElem].pNode = dynamic_cast<StructNode *>(pDM->ReadNode(HP, Node::STRUCTURAL));
		if (!nodes[iElem].pNode) {
			silent_cerr("aerodyn module(" << GetLabel() << "): structural node expected at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		// (optional) aerodynamics offset with respect to the node,
		// constant in the reference frame of the node
		if (HP.IsKeyWord("position")) {
			nodes[iElem].f = HP.GetPosRel(ReferenceFrame(nodes[iElem].pNode));

		} else {
			nodes[iElem].f = Zero3;
		}

		// (optional) relative orientation between the aerodynamics
		// and the node
		if (HP.IsKeyWord("orientation")) {
			nodes[iElem].Ra = HP.GetRotRel(ReferenceFrame(nodes[iElem].pNode));

		} else {
			nodes[iElem].Ra = Eye3;
		}
	}

	// read force scalse factor
	if (HP.IsKeyWord("force" "scale" "factor")) {
		FSF.Set(HP.GetDriveCaller());
	} else {
		FSF.Set(new OneDriveCaller);
	}

	// set output 
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
	
	// get input and output filename
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

	// allocate memories in fortran variables
	char Version[26 + 1];
	snprintf(Version, sizeof(Version), "(%s)", VERSION);
	for (unsigned i = strlen(Version); i < sizeof(Version); i++) {
		Version[i] = ' ';
	}
	Version[sizeof(Version) - 1] = '\0';
	F_INTEGER num_blade = NumBlades;
	F_INTEGER num_elements = NumElements;
	__FC_DECL__(mbdyn_init)(Version, &num_blade, &num_elements);

	// pass initial interface component to AeroDyn
	SetADInitialInterfaceComponent();
	
	// initialize AeroDyn
	F_INTEGER input_file_name_len = input_file_name.size();
	F_INTEGER elem_file_name_len = elem_file_name.size();
	__FC_DECL__(mbdyn_ad_inputgate)(&bladeLength, (F_CHAR *)input_file_name.c_str(), &input_file_name_len, (F_CHAR *)elem_file_name.c_str(), &elem_file_name_len);

	// set time driver
	Time.Set(new TimeDriveCaller(pDM->pGetDrvHdl()));
	dCurTime = Time.dGet();
	::module_aerodyn = this;

}

void
AeroDynModule::ConvertNodeToMarker(const StructNode *Xnode, const StructNode *Rnode, const StructNode *Vnode, const StructNode *Wnode){
	const Vec3& x = Xnode->GetXCurr();
	const Vec3& v = Vnode->GetVCurr();
	const Vec3& w = Wnode->GetWCurr();
	const Mat3x3& R = Rnode->GetRCurr();	

	// convert variables 

	for (int i = 0; i < 3; i++) {
		Position[i]       = x.dGet(i + 1);
		TranslationVel[i] = v.dGet(i + 1);
		RotationVel[i]    = w.dGet(i + 1);

		for (int j = 0; j < 3; j++){
			/* Node! 
			 * Rows and columns are reversed when passing multiple arrays from C++ to Fortran
			 * However, the rotation matrices input to Aerodyn and retrieved from MBDyn have a transpose relationship, so they are automatically consistent.      
			 */
			Orientation[i][j] = R.dGet(i + 1, j + 1);
		}
	}
}

void 
AeroDynModule::SetADInitialInterfaceComponent() {
	// set blade root component   number of component is 1
	F_INTEGER num_component = 1;
	c_blade = 1;
	for(int ibld = 0; ibld < NumBlades; ibld++) {
		ConvertNodeToMarker(pBladeRoot[ibld], pBladeRoot[ibld], pBladeRoot[ibld], pBladeRoot[ibld]);
		__FC_DECL__(set_ad_interfacecomponent)(&num_component, &c_blade, Position, Orientation, TranslationVel, RotationVel);
		c_blade++;
	}
	// set hub component          number of component is 2
	num_component = 2;
	ConvertNodeToMarker(pHub, pHub, pHub, pHub);
	__FC_DECL__(set_ad_interfacecomponent)(&num_component, &c_blade, Position, Orientation, TranslationVel, RotationVel);
}

void 
AeroDynModule::SetADInterfaceComponent() {
	// set blade root component   number of component is 1
	// should be marker 10000*K (BladeKCS marker), but the current version of AeroDyn calculates forces normal and tangential to the cone of rotation
	F_INTEGER num_component = 1;
	c_blade = 1;
	for(int ibld = 0; ibld < NumBlades; ibld++) {
		ConvertNodeToMarker(pConed[ibld], pConed[ibld], pConed[ibld], pConed[ibld]);
		__FC_DECL__(set_ad_interfacecomponent)(&num_component, &c_blade, Position, Orientation, TranslationVel, RotationVel);
		c_blade++;
	}

	// set hub component          number of component is 2
	// should be marker 4000 (HubCS), but the current version of AeroDyn treats teeter deflections like blade deflections
	num_component = 2;
	ConvertNodeToMarker(pLSSTeet, pLSSTeet, pLSSTeet, pLSSTeet);
	__FC_DECL__(set_ad_interfacecomponent)(&num_component, &c_blade, Position, Orientation, TranslationVel, RotationVel);

	// set rotor furl component   number of component is 3
	num_component = 3;
	//RotorFurlRef: should be marker 2150 (ShaftCS_M), but 2100 is needed for HubVDue2Yaw
	ConvertNodeToMarker(pRotorFurl, pShaftCS, pShaftCS, pShaftCS);
	__FC_DECL__(set_ad_interfacecomponent)(&num_component, &c_blade, Position, Orientation, TranslationVel, RotationVel);

	// set nacelle component      number of component is 4
	num_component = 4;
	ConvertNodeToMarker(pNacelle, pNacelle, pNacelle, pNacelle);
	__FC_DECL__(set_ad_interfacecomponent)(&num_component, &c_blade, Position, Orientation, TranslationVel, RotationVel);

	// set tower component        number of component is 5
	num_component = 5;
	//PlatformRef: should be marker 1100 (TowerBaseCS_M) but 1000 is needed for HubVDue2Yaw
	ConvertNodeToMarker(pPlatfromRef,pTower, pTower, pTower);
	__FC_DECL__(set_ad_interfacecomponent)(&num_component, &c_blade, Position, Orientation, TranslationVel, RotationVel);
}
void
AeroDynModule::SetAeroMarkers()
{
	c_blade = 0;
	c_elem = 0;
	for(int ielem=0; ielem<NumBlades*NumElements; ielem++) {

		// get node information
		const Mat3x3 R = nodes[ielem].pNode->GetRCurr() * nodes[ielem].Ra;
		const Vec3 x = nodes[ielem].pNode->GetXCurr() + R * nodes[ielem].f;
		const Vec3 w = nodes[ielem].pNode->GetWCurr();
		const Vec3 v = nodes[ielem].pNode->GetVCurr() + w.Cross(R*nodes[ielem].f);

		// find current blade and current element
		if(ielem%NumElements==0){
			c_blade++;
		}
		c_elem = ielem%NumElements + 1;
			
		for (int i = 0; i < 3; i++) {
			Position[i]       = x.dGet(i + 1);
			TranslationVel[i] = v.dGet(i + 1);
			RotationVel[i]    = w.dGet(i + 1);

			for (int j = 0; j < 3; j++){
				/* Node! 
				 * Rows and columns are reversed when passing multiple arrays from C++ to Fortran
				 * However, the rotation matrices input to Aerodyn and retrieved from MBDyn have a transpose relationship, so they are automatically consistent.      
				 */
				Orientation[i][j] = R.dGet(i + 1, j + 1);
			}
		}

		// call fortran subroutine to pass variables
		__FC_DECL__(set_aeromarkers)(&c_blade, &c_elem, Position, Orientation, TranslationVel, RotationVel);
	}
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
	*piNumRows = 6*NumBlades*NumElements;
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

	if (bFirst) {
		// force scale factor to "ramp up" loads
		doublereal dFSF = FSF.dGet();

		// sanity check
		ASSERT(dFSF >= 0.);
		ASSERT(dFSF <= 1.);

		SetADInterfaceComponent();
		SetAeroMarkers();

		c_time = Time.dGet();
		c_blade = 0;
		for(int iElem = 0; iElem <NumBlades*NumElements; iElem++) {

			// find current blade and current element
			if(iElem%NumElements==0){
				c_blade++;
			}
			c_elem = iElem%NumElements + 1;

			// calculate aerodynamic loads
			F_REAL	Force_Aero[3];
			F_REAL  Moment_Aero[3];
			__FC_DECL__(call_ad_calculateloads)(&c_time, &c_blade, &c_elem, Force_Aero, Moment_Aero);

			const Mat3x3& R = nodes[iElem].Ra;   // rotation matrix aero point to struct node
			nodes[iElem].Force = R * Vec3(Force_Aero[0],Force_Aero[1],Force_Aero[2]);
			nodes[iElem].Moment = R * Vec3(Moment_Aero[0],Moment_Aero[1],Moment_Aero[2]) + nodes[iElem].f.Cross(nodes[iElem].Force);
		}

		bFirst = false;
	}

	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);

	WorkVec.ResizeReset(iNumRows);

	TF = Vec3(Zero3);
	TM = Vec3(Zero3);

	for (int iElem = 0; iElem < NumBlades*NumElements; iElem++) {
		/*
		 * set indices where force/moment need to be put
		 */
		integer iFirstIndex = nodes[iElem].pNode->iGetFirstMomentumIndex();
		for (int i = 1; i <= 6; i++) {
			WorkVec.PutRowIndex(6*iElem + i, iFirstIndex + i);
		}
		/*
		 * add force/moment to residual, after rotating them
		 * into the global frame
		 */
		/*
		 * Transfer the force/moment to global frame. 
	     */	 
  		WorkVec.Add(6*iElem + 1, nodes[iElem].pNode->GetRCurr()*(nodes[iElem].Force));
  		WorkVec.Add(6*iElem + 4, nodes[iElem].pNode->GetRCurr()*(nodes[iElem].Moment));
		/*
		 * calculate the force and moment contributions on the rotor in absolute frame.
	     */
		const Vec3& Xp = nodes[iElem].pNode->GetXCurr();
		const Vec3& Xh = pHub->GetXCurr();

		TF = TF + nodes[iElem].pNode->GetRCurr()*nodes[iElem].Force;
		TM = TM + nodes[iElem].pNode->GetRCurr()*nodes[iElem].Moment
			+ (Xp-Xh).Cross(nodes[iElem].pNode->GetRCurr()*nodes[iElem].Ra*nodes[iElem].Force);
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
	Rotor_speed =pNacelle->GetRCurr().MulTV(Wh - Wn).dGet(1);
	Wind_speed = dGetWindSpeed();

	return WorkVec;
}

void
AeroDynModule::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
	bFirst = true;
}


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
	*piNumRows = 6*NumBlades*NumElements;
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

doublereal 
AeroDynModule::dGetWindSpeed(void) const
{
	F_REAL HorWindV;
	F_REAL time = Time.dGet();
	__FC_DECL__(get_horizonal_wind_speed)(&time,&HorWindV);
	return HorWindV;
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

	UserDefinedElemRead *rfd1 = new UDERead<VSControlModule>;

	if (!SetUDE("discon", rfd1)) {
		delete rfd1;

		silent_cerr("module-discon: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}
/*
	UserDefinedElemRead *rfd2 = new UDERead<PitchControlModule>;

	if (!SetUDE("pitch" "control", rfd2)) {
		delete rfd2;

		silent_cerr("module-discon: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}
*/
	return 0;
}

