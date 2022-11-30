#include <fstream>
#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include <sstream>
#include <random>

#include "dataman.h"

#ifndef __IRRAGULARWAVE_H_INCLUDE__
#define __IRRAGULARWAVE_H_INCLUDE__


class IrragularWave {
private:
    struct WaveParam
    {
        double          zeta;
        double          omega;
        double          kappa;
        double          epsilon;
        double          theta;
    };

    int                 numWave;
    std::vector<WaveParam>  cWaves;

    double              waterDen;
    double              waterDepth;
    double              gravity;
    double              waveAngle;
    double              rampTime;
    double              rampFactor;

    double              zStretch;
    double              surfaceElev;
    double              waveVelH;
    double              waveVelV;
    double              waveAccelH;
    double              waveAccelV;
    double              wavePress;

    bool                bset = false;

    // jcobian dzeta_dx(3x1) dzs_dx (3x1)
    Vec3                dzeta_dx;
    Vec3                dzs_dx;

    Mat3x3              du_dx;
    Mat3x3              duP_dx;
    Vec3                dp_dx;

    const double pi = M_PI;

private:
    double calc_dispersion_f(double omega);
    std::string 
    trim(const std::string& string, const char* trimCharacterList = " \t\v\r\n"); 

public:
    IrragularWave(){};
    ~IrragularWave(){};

    void input(std::string inFileName);
    void calc_uvw(double t, const Vec3 &x);

    inline double get_water_elevation();
    inline Vec3 get_velocity();
    inline Vec3 get_acceleration();
    inline double get_dynamic_pressure();
    inline const Vec3& get_dzeta_dx();
    inline const Mat3x3& get_du_dx();
    inline const Mat3x3& get_duP_dx();
    inline const Vec3& get_dp_dx();
    inline bool bWave_set(); 
    inline double get_parameter(int i, std::string name);
};
void 
IrragularWave::input(std::string inFileName)
{
    inFileName = trim(inFileName);
    std::ifstream   reading_file;    //object of ifstrem
    reading_file.open(inFileName, std::ios::in);

    if(!reading_file.is_open())
    {
        std::cerr << "Could not open the wave input file - '"
                  << inFileName << "'" <<std::endl;
    }
   
    std::string line;
    std::vector<std::string>    lines;
    while(std::getline(reading_file,line))
    {

        std::string::size_type posComment = line.find_first_of("!");
        if(posComment!=std::string::npos)
        { 
            line.erase(posComment, line.length());
        }
        
        if(!line.empty())
        {
            lines.push_back(line);
        }
    }     
    double deg2rad = pi/180.0;

    waterDen    = std::stod(lines[0]);
    waterDepth  = std::stod(lines[1]);
    gravity     = std::stod(lines[2]);
    double H3   = std::stod(lines[3]);
    double T3   = std::stod(lines[4]);
    waveAngle   = std::stod(lines[5])*deg2rad;
    numWave     = std::stoi(lines[6]);
    double alfaBM = std::stod(lines[7]);
    double betaBM = std::stod(lines[8]);
    rampTime    = std::stod(lines[9]);
    int seed1   = std::stoi(lines[10]);
    int seed2   = std::stoi(lines[11]);

    cWaves.resize(numWave);

    std::vector<double> eps;
    if(lines.size() > 12 )
    {
        for(int i=0; i<lines.size()-12; i++)
        {
            std::istringstream cline(lines[12+i]);
            std::string separate_line;
            while (std::getline(cline, separate_line, ' '))
            {
                eps.push_back(std::stod(separate_line)*2.0*pi);
            }            
        }
        for(int i=0; i< numWave; i++) {
            cWaves[i].epsilon = eps[i];
        }
    }
    else
    {
        std::mt19937 mt;     
        mt.seed((seed1,seed2) ); 
        std::uniform_real_distribution<> rand1(0, 1);  
        for(auto &i : cWaves){
            i.epsilon = 2.0*pi*rand1(mt);
        }      
    }

    // calculate the wave parameters
    for(int i=0; i<numWave; i++){
        cWaves[i].zeta = std::sqrt(0.5*alfaBM/betaBM)*H3/std::sqrt(double(numWave));
        cWaves[i].omega = 2.*pi*std::pow(betaBM,0.25)/T3*std::pow( std::log(2.*numWave/(2.*double(i)+ 1.)) , (-0.25) );

        cWaves[i].kappa = calc_dispersion_f(cWaves[i].omega);
    }
    if (numWave == 1)
    {
        cWaves[0].zeta = H3*0.5;
        cWaves[0].omega = 2.*pi/T3;
        cWaves[0].epsilon = 0.;

        cWaves[0].kappa = calc_dispersion_f(cWaves[0].omega);
    }   

    bset = true;

}

// mat3x3 result = v(3x1) v^T(1x3)
Mat3x3
MultfgT(const Vec3& v, const Vec3& w )
{
	Mat3x3 vd(Zero3,v,Zero3);
	Mat3x3 wd(Zero3,w,Zero3);
	return vd.MulMT(wd);
}

void 
IrragularWave::calc_uvw(double t, const Vec3 &x)
{
    surfaceElev=0.0;
    waveVelH = 0.0;
    waveVelV = 0.0;
    waveAccelH = 0.0;
    waveAccelV = 0.0;
    wavePress = 0.0;

    for(int i=0; i < numWave; i++)
    {
        cWaves[i].theta = cWaves[i].kappa*(x.dGet(1)*std::cos(waveAngle)+x.dGet(2)*std::sin(waveAngle)) - cWaves[i].omega*t + cWaves[i].epsilon;
        surfaceElev += cWaves[i].zeta * std::cos(cWaves[i].theta); 
    }


    zStretch = (x.dGet(3) - surfaceElev)/(waterDepth + surfaceElev) * waterDepth;


    for (int i = 0; i < numWave; i++)
    {
        // temporary variable
        double COSH    = std::cosh( cWaves[i].kappa * waterDepth );
        double COSH_zs = std::cosh( cWaves[i].kappa * (waterDepth + zStretch  ) ) / COSH;
        double SINH_zs = std::sinh( cWaves[i].kappa * (waterDepth + zStretch  ) ) / COSH;
        double COSH_z  = std::cosh( cWaves[i].kappa * (waterDepth + x.dGet(3) ) ) / COSH;
        double SINH_z  = std::sinh( cWaves[i].kappa * (waterDepth + x.dGet(3) ) ) / COSH;

        double cosTH = std::cos(cWaves[i].theta);  
        double sinTH = std::sin(cWaves[i].theta);

        // temporary variable, coefficient
        double ci = cWaves[i].zeta * cWaves[i].kappa * gravity;
        
        waveVelH   +=  ci / cWaves[i].omega * COSH_zs *cosTH;                
        waveVelV   +=  ci / cWaves[i].omega * SINH_zs *sinTH;                 
        waveAccelH +=  ci                   * COSH_zs *sinTH;                
        waveAccelV -=  ci                   * SINH_zs *cosTH;

        wavePress  +=  cWaves[i].zeta       * COSH_z * cosTH; 

#if 0
        // calcuate jacobian
        dVH_dzs    += ci / cWaves[i].omega * cWaves[i].kappa * SINH_zs *cosTH;
        dVV_dzs    += ci / cWaves[i].omega * cWaves[i].kappa * COSH_zs *sinTH;
        dAH_dzs    += ci                   * cWaves[i].kappa * SINH_zs *sinTH;
        dAV_dzs    += ci                   * cWaves[i].kappa * COSH_zs *cosTH;

        dVH_dth    -= ci / cWaves[i].omega * cWaves[i].kappa * COSH_zs *sinTH;
        dVV_dth    += ci / cWaves[i].omega * cWaves[i].kappa * SINH_zs *cosTH;
        dAH_dth    += ci                   * cWaves[i].kappa * COSH_zs *cosTH;
        dAV_dth    -= ci                   * cWaves[i].kappa * SINH_zs *sinTH;  

        dp_dz      +=  cWaves[i].zeta * cWaves[i].kappa      * COSH_z * cosTH;
        dp_dth     -=  cWaves[i].zeta * cWaves[i].kappa      * SINH_z * sinTH;
#endif
    }  

    double rampFactor;
    if(t < rampTime) {
        rampFactor = std::sin(pi/(rampTime*2) * t);
    }
    else {
        rampFactor = 1.0;
    }

    surfaceElev *= rampFactor;
    waveVelH *= rampFactor;
    waveVelV *= rampFactor;
    waveAccelH *= rampFactor;
    waveAccelV *= rampFactor;
    wavePress *= rampFactor;

    /*
    du_dx =   MultfgT(Vec3(dVH_dzs*std::cos(waveAngle), dVH_dzs*std::sin(waveAngle), dVV_dzs), dzs_dx)
            + MultfgT(Vec3(dVH_dth*std::cos(waveAngle), dVH_dth*std::sin(waveAngle), dVV_dth), Vec3(std::cos(waveAngle), std::sin(waveAngle), 0.));

    duP_dx=   MultfgT(Vec3(dAH_dzs*std::cos(waveAngle), dAH_dzs*std::sin(waveAngle), dAV_dzs), dzs_dx)
            + MultfgT(Vec3(dAH_dth*std::cos(waveAngle), dAH_dth*std::sin(waveAngle), dAV_dth), Vec3(std::cos(waveAngle), std::sin(waveAngle), 0.));
    
    dp_dx = Vec3(dp_dth*std::cos(waveAngle),
                 dp_dth*std::sin(waveAngle),
                 dp_dz                      );   
    dzeta_dx *= rampFactor;
    du_dx  *= rampFactor;
    duP_dx *= rampFactor;
    dp_dx  *= rampFactor;
    */
}

double 
IrragularWave::get_water_elevation(){
    return surfaceElev;
}

Vec3 
IrragularWave::get_velocity()
{
    return Vec3( waveVelH * std::cos(waveAngle),
                 waveVelH * std::sin(waveAngle),  
                 waveVelV);
}

Vec3 
IrragularWave::get_acceleration() 
{
    return Vec3( waveAccelH * std::cos(waveAngle),
                 waveAccelH * std::sin(waveAngle),
                 waveAccelV );
}

double 
IrragularWave::get_dynamic_pressure() 
{
    return wavePress*waterDen*gravity;
}

const Vec3& 
IrragularWave::get_dzeta_dx() 
{
    return dzeta_dx;
}

const Mat3x3& 
IrragularWave::get_du_dx() 
{
    return du_dx;
} 

const Mat3x3& 
IrragularWave::get_duP_dx() 
{
    return duP_dx;
}

const Vec3& 
IrragularWave::get_dp_dx() 
{
    return dp_dx;
}

bool 
IrragularWave::bWave_set() 
{
    return bset;
}

double 
IrragularWave::calc_dispersion_f(double omega) {
    /* dispersion formulation of water wave
    *  oemga *2 = gk tanh(kh)
    *  This function takes omega as the argument and computes and returns k
    */
    double OME    = waterDepth / gravity *omega * omega;
    double F      = 1. + OME * ( 0.666 + OME * ( 0.445 + OME * ( -0.105 + OME * 0.272 ) ) );
    double T      = 2. * pi / omega;
    double Lambda = T * std::sqrt( gravity * waterDepth ) * std::sqrt( F / ( 1.0 + OME * F ) );

    return  2.*pi/Lambda;
}

std::string 
IrragularWave::trim(const std::string& string, const char* trimCharacterList) 
{

    std::string result;
    
    std::string::size_type left = string.find_first_not_of(trimCharacterList);

    if (left != std::string::npos)
    {
        std::string::size_type right = string.find_last_not_of(trimCharacterList);
        result = string.substr(left, right - left + 1);
    }
    return result;

}

double 
IrragularWave::get_parameter(int i, std::string name) {
    double result;
    if (name == "zeta") {
        result = cWaves[i].zeta;
    }
    else if (name == "omega") {
        result =  cWaves[i].omega;
    }
    else if (name == "kappa") {
        result = cWaves[i].kappa;
    }
    else if (name == "epsilon") {
        result = cWaves[i].epsilon;
    }
    return result;
}




#endif