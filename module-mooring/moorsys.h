#include <fstream>
#include <iostream>
#include <math.h>
#include <string>
#include <vector>

#include "mbconfig.h"  


#ifndef __MOORSYS_H_INCLUDE__
#define __MOORSYS_H_INCLUDE__

/*
 *  fortran77で作製された TU11-00075_moorsys.fをC++から呼ぶインターフェース
 *　gfortran -c TU11-00075_moorsys.f cpp2moorsys.f90 -lstdc++ でオブジェクトファイル　同名.oのオブジェクトファイルを作製
 *  g++ -c moorsysf.cc　
 *  g++ cpp2moorsys.o moorsysf.o -o out -lgfortran
 */

extern "C" 
{
    extern int
    __FC_DECL__(input_moorsys)(char *ifname, long int *ifnamelen, char *ofname, long int *ofnamelen);
    extern int
    __FC_DECL__(call_moorsys)(double *t, double *xjp, double *rjp);
}

class Moorsysf {

    private :
    Vec3       force;
    Vec3       moment;

    double      initZCoord;
    bool        bFirst;

    inline Vec3 GetEular321formR(Mat3x3& R);
    
    public :
    Moorsysf() {bFirst = true;};
    ~Moorsysf() {};
    void input(std::string inFileName, std::string outFileName);
    void set_position(double t, const Vec3 &x, const Mat3x3& R);
    inline Vec3 get_force();
    inline Vec3 get_moment();
};

void
Moorsysf::input(std::string inFileName, std::string outFileName) {
    long int inFileNameLen = inFileName.size();
    long int outFileNameLen = outFileName.size();

    __FC_DECL__(input_moorsys)((char *)inFileName.c_str(), &inFileNameLen, (char *)outFileName.c_str(), &outFileNameLen);
}
void
Moorsysf::set_position(double t, const Vec3 &x, const Mat3x3& R) {
    if(bFirst) {
        initZCoord = x.dGet(3);
        bFirst = false;
    }

    Mat3x3 Aob(R);
    Vec3 eular321 = GetEular321formR(Aob);

    double Xpj[6], Rjp[6];

    for(int i=0; i<3; i++) {
        Xpj[i] = x.dGet(i+1);
        Xpj[i+3] = eular321.dGet(i+1);
    }

    Xpj[2] -= initZCoord;

    __FC_DECL__(call_moorsys)(&t, Xpj, Rjp);

    force  = Vec3(Rjp[0], Rjp[1], Rjp[2])*(-1000.);
    moment = Vec3(Rjp[3], Rjp[4], Rjp[5])*(-1000.);

}

// fuction caluclating the Eular3-2-1 angle form Rotation matrix
Vec3 
Moorsysf::GetEular321formR(Mat3x3& R)
{
	//the rotation matrix from mbdyn is transforming to global form body coordinate
	double theta;  		// roll angle [rad]
	double fai;			// pitch angle [rad]
	double psi;			// yaw anble [rad]

	theta = std::atan(R.dGet(3,2)/R.dGet(3,3));
	fai = std::asin(-R.dGet(3,1));
	psi = std::atan(R.dGet(2,1)/R.dGet(1,1));

	return Vec3(theta,fai,psi);
}

inline Vec3
Moorsysf::get_force() {
    return force;
}

inline Vec3
Moorsysf::get_moment() {
    return moment;
}

#endif
