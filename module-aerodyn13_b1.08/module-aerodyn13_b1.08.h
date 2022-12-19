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

class AeroDynModule
: virtual public Elem,
public UserDefinedElem
{
private:

	struct AeroNode {
		StructNode	*pNode;
		Vec3		f;	// offset of the aero point wrt./ the node,
						// constant in the node's reference frame
		Mat3x3		Ra;	// aerodynamic orientation of the aero point
						// wrt./ the node

		doublereal	dBuiltInTwist;

		Vec3		F;	// Force acting on Node
		Vec3		M;	// Moment acting on Node, with respect
						// to node's position

 		doublereal      FN;    // Normal Force on each blade element (Not being used now!).
 		doublereal      FT;    // Tangental force on each blade element.(Not being used now!)
 		doublereal      AM;    // Aerodynamic moment on each blade element.(Not being used now!)
 
		doublereal	PITNOW; // Node pitch angle
	};

	/**
	 * Nacelle node; requirements:
	 * - axis 3 is the shaft axis
	 * - axis 3 in wind direction
	 */
	StructNode	*pNacelle;
	StructNode	*pHub;
	StructNode  *pRotorFurl;		// RotorFurl node追加
	StructNode  *pTower;			// Tower node追加
	std::vector<AeroNode>	pBladeroot;		// Bladeroot node

	integer		nblades;	// the number of blades.
	integer		nelems;		// the number of elements per blade.

	doublereal	Hub_Tower_xy_distance;		//hub ~ tower yaw axis distance

	/*
	 * node data
	 */
	std::vector<AeroNode>	nodes;		// nodes into Aerodynamic element
	std::vector<Mat3x3>		bladeR;		// orientation matrix of each blade root in the hub reference frame

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
	/* 
	 * internal states to access the Variables which is defined in the 
	 * common MODULEs of AeroDyn
	 */
	F_LOGICAL	FirstLoop;
	F_INTEGER	elem;    // use to identify the current element in the interface module!
	
	F_INTEGER	blade;    // use to identify the current element in the interface module!

	F_INTEGER	c_elem;  // use to identify the current element in AeroDyn!
	F_INTEGER	c_blade; // use to identify the current blade in AeroDyn!
	F_REAL		rlocal;  // use to identify the current element position 
	F_REAL		b_length;// use to identify the blade length. 
	F_REAL		nacelle_position[3];
	F_REAL		nacelle_orientation[9];
    F_REAL		nacelle_transvel[3];
	F_REAL		nacelle_rotvel[3];
	F_REAL		hub_position[3];
	F_REAL		hub_orientation[9];
    F_REAL		hub_transvel[3];
	F_REAL		hub_rotvel[3];
	F_REAL		rotorfurl_position[3];
	F_REAL		rotorfurl_orientation[9];
    F_REAL		rotorfurl_transvel[3];
	F_REAL		rotorfurl_rotvel[3];
	F_REAL		tower_position[3];
	F_REAL		tower_orientation[9];
    F_REAL		tower_transvel[3];
	F_REAL		tower_rotvel[3];
	F_REAL		bladeroot_position[3] ;
	F_REAL		bladeroot_transvel[3];
	F_REAL		bladeroot_rotvel[3];
	F_REAL		bladeroot_orientation[9];
	F_REAL		bladeelem_position[3];
	F_REAL		bladeelem_transvel[3];
	F_REAL		bladeelem_rotvel[3];
	F_REAL		bladeelem_orientation[9];

	F_REAL      c_time;		// current time
	
	bool        bFirst;
	DriveOwner	Time;		// time drive
	doublereal	dOldTime;	// old time
	doublereal  dCurTime;   // current time
	F_REAL      dDT;		// time step

	DriveOwner  FSF;

public:
	// interface for AeroDyn callbacks

	// module_aerodyn->pNacelle
	const StructNode *pGetNacelleNode(void) const;
	// module_aerodyn->pHub
	const StructNode *pGetHubNode(void) const;
	// module_aerodyn->pRotorFurl
	const StructNode *pGetRotorFurlNode(void) const;
	// module_aerodyn->pTower
	const StructNode *pGetTowerNode(void) const;
	// module_aerodyn->nodes
	const StructNode *pGetCurrBladerootNode(void) const;
	// module_aerodyn->Hub_Tower_xy_distance;
	doublereal dGetHubTowerXYDistance(void) const;
	// module_aerodyn->nblades
	F_INTEGER iGetNumBlades(void) const;
	// module_aerodyn->c_blade
	F_INTEGER iGetCurrBlade(void) const;
	// module_aerodyn->bladeR
	const Mat3x3& GetCurrBladeR(void) const;
	// module_aerodyn->nelems
	F_INTEGER iGetNumBladeElems(void) const;
	// module_aerodyn->elem
	F_INTEGER iGetCurrBladeElem(void) const;
	// module_aerodyn->nodes
	const StructNode *pGetCurrBladeNode(void) const;
	// module_aerodyn->nodes[::module_aerodyn->elem].dBuiltInTwist
	doublereal dGetCurrBladeNodeBuiltinTwist(void) const;
	// module_aerodyn->nodes[::module_aerodyn->elem].PITNOW
	void SetCurrBladeNodePITNOW(doublereal PITNOW);
	doublereal dGetCurrBladeNodePITNOW(void) const;
	// module_aerodyn->nodes[::module_aerodyn->elem].Ra
    const Mat3x3& GetCurrBladeNodeRa(void) const;
	// module_aerodyn->WindSpeed
	doublereal dGetWindSpeed(void) const;

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
#if 0
	unsigned int iGetPrivDataIdx(const char *s) const;
	doublereal dGetPrivData(unsigned int i) const;
#endif
	int GetNumConnectedNodes(void) const;
	void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
};

// static handler for AeroDyn callbacks
// only one element per simulation can be active (AeroDyn limitation)
static AeroDynModule *module_aerodyn;

