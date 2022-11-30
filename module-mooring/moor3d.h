#include <fstream>
#include <iostream>
#include <math.h>
#include <string>
#include <vector>


#include "mbconfig.h" 

#ifndef __MOOR3D_H_INCLUDE__
#define __MOOR3D_H_INCLUDE__

class Moor3d{
private:
    int     numLines;
    Vec3   force;
    Vec3   moment;

    double  initZCoord;

    struct catenaryLine{
        double  initHightFP;
        double  weight;
        double  initHorzTen;
        Vec3   localPosFP;
        double  horzAngle;

        /* 係留張力が０となる状態(係留鎖がL字上に海底に寝そべっている)のFP位置を原点とし、そこから初期のFP位置までの水平移動量 */
        double  initDelta; 
        double  lenCatenary;
        double  hTension;
        double  vTension;
        double  tension;
    };

    std::vector<catenaryLine> mLines;

    
    inline double catenary_func(double x, double d);
    inline double catenary_derivfunc(double x, double d);
    
public:
    Moor3d(){};
    void input(std::string inpFileName);
    void set_position(const Vec3 &x, const Mat3x3 &R);
    inline Vec3 get_force();
    inline Vec3 get_moment();
    inline double get_htension(int n);
    inline double get_vtension(int n);
    inline double get_catenary_length(int n);

    double catenary(double d);
    
};

void
Moor3d::input(std::string inpFileName){
    std::ifstream   reading_file;    //object of ifstrem
    reading_file.open(inpFileName, std::ios::in);

    if(!reading_file.is_open()){
        std::cerr << "Could not open the moorsys input file - '"
                  << inpFileName << "'" <<std::endl;
    }

    std::string line;
    std::vector<std::string>  lines;
    while(std::getline(reading_file, line)) {
        std::string::size_type posComment = line.find_first_of("!");
        if(posComment!=std::string::npos){
            line.erase(posComment, line.length());
        }
        if(!line.empty()) {
            lines.push_back(line);
        }
    }

    numLines = std::stod(lines[0]);
    initZCoord = std::stod(lines[1]);
    mLines.resize(numLines);

    for(int i=0; i<numLines; i++){

        double temp1, temp2, temp3;
        std::istringstream ssr(lines[2+i]);
        ssr          >> mLines[i].initHightFP
                     >> mLines[i].weight
                     >> mLines[i].initHorzTen
                     >> temp1
                     >> temp2
                     >> temp3
                     >> mLines[i].horzAngle;
        
        mLines[i].localPosFP = Vec3(temp1, temp2, temp3 - initZCoord);
        mLines[i].horzAngle *= M_PI/180.0;
        /*                                    ________                  
         *             H         wh         /      H
         *   delta/h = —— Acosh( ——　+1) - √  1 + 2——   + 1
         *             wh        H                 wh
         */
        double a = mLines[i].initHorzTen / (mLines[i].weight * mLines[i].initHightFP); 
        mLines[i].initDelta = mLines[i].initHightFP*(a*std::acosh(1.0/a + 1.0) - std::sqrt(1.0 + 2.0*a) + 1.0);
    }
}

void
Moor3d::set_position(const Vec3 &nodeX, const Mat3x3 & nodeR){

    force = Vec3(0.,0.,0.);
    moment = Vec3(0.,0.,0.);

    for(int i=0; i<numLines; i++){
        /* calculate current positon each fairleaders in global coordinate system*/
        Vec3 fairLPos = nodeX + nodeR * mLines[i].localPosFP;

        /* calculate temprary variables*/
        double currDelta = mLines[i].initDelta - ( (fairLPos[0] - mLines[i].localPosFP[0])*std::cos(mLines[i].horzAngle)
                                                  +(fairLPos[1] - mLines[i].localPosFP[1])*std::sin(mLines[i].horzAngle));

        double currHight = mLines[i].initHightFP + (fairLPos[2] - (initZCoord + mLines[i].localPosFP[2]));

        double d = currDelta/currHight; // d = delta/h
        double f = catenary(d);         // f = H/wh
    
        /* calculate parameters */
        mLines[i].hTension = f * mLines[i].weight * currHight;

        mLines[i].lenCatenary = currHight * std::sqrt(1.0 + 2.0*f);
        mLines[i].vTension = mLines[i].weight * mLines[i].lenCatenary;
        mLines[i].tension = mLines[i].hTension + mLines[i].weight*currHight;

        Vec3 df(mLines[i].hTension*std::cos(mLines[i].horzAngle), mLines[i].hTension*std::sin(mLines[i].horzAngle), -mLines[i].vTension);
        Vec3 dm = (fairLPos - nodeX).Cross(df);

        /* calculate force actiong on node*/
        force += df;
        moment += dm; 
    }
}

double 
Moor3d::catenary(double d){
    double x1 = 0.0;
    double x2 = 100.0;
    double exp = std::pow(10, -6);
    double maxIter = 100;

    double xcurr;

    double xl, xh;

    if(d<=0.0){
        xcurr = 0.0;
    }
    else{
        /* calulation by newton-raphson method*/
        double fl = catenary_func(x1,d);
        double fh = catenary_func(x2,d);

        if( (fl<0.0 && fh<0.0) || (fl>0.0 && fh>0.0) ){
            std::cerr << "No solution exists in the given range"<<std::endl;
            xcurr = 0.0;        
        } 
        else if( fl == 0.0){
            xcurr = x1;
        }
        else if(fh == 0.0){
            xcurr = x2;
        }
        else{

            if(fl < 0.0){
                xl = x1;
                xh = x2;
            }
            else{
                xl = x2;
                xh = x1;
            }

            xcurr = (x1 + x2)/2.0;
            double dxold = std::abs(x2-x1);
            double dx = dxold;

            double fcurr = catenary_func(xcurr,d);
            double dfcurr = catenary_derivfunc(xcurr,d);

            for(int i=0; i<= maxIter; i++){
                double temp = ((xcurr - xh)*dfcurr - fcurr) * ((xcurr - xl)*dfcurr - fcurr);

                if( (temp >= 0.0) || ( std::abs(2.0*fcurr) >= std::abs(dxold * dfcurr) ) ){
                    dxold = dx;
                    dx = (xh-xl)/2.0;
                    xcurr = xl + dx;
                    if(xcurr == xl){
                        break;
                    }

                }
                else{
                    dxold = dx;
                    dx = fcurr/dfcurr;
                    double t = xcurr;
                    xcurr = xl + dx;
                    if (xcurr == t){
                        break;
                    }
                }

                if(std::abs(dx)<=exp){
                    break;
                }

                fcurr = catenary_func(xcurr,d);
                dfcurr = catenary_derivfunc(xcurr,d);

                if(fcurr<0.0){
                    xl = xcurr;
                }
                else{
                    xh = xcurr;
                }

                if(i==maxIter){
                    std::cerr << "exceeded maximum iterations"<<std::endl;
                }
            }
        }
    }
    return xcurr;
}

double 
Moor3d::catenary_func(double x, double d){
    double f;
    if(x<=0.0){
        f = -d;
    }
    else{
        f =x * std::acosh(1.0/x + 1.0) - std::sqrt(1.0 + 2.0*x) + 1.0 - d;
    }
    return f;
}

double 
Moor3d::catenary_derivfunc(double x, double d){
    /*
    *                  1            1                  1
    *   df/dx = acosh( — + 1) - ————————   -  —————————————————————
    *                  x         ______         　_____________
    *                           √1 + 2x        x√　(1+1/x)^2 -1
    */
    double df;
    if(x<=0.0){
        df=0.0;
    }
    else{
        df = std::acosh(1.0/x + 1.0) -1.0/std::sqrt(1.0+2.0*x) -1.0/( x * std::sqrt( std::pow( (1.0+1.0/x), 2 ) - 1.0 ));
    }
    return df;
}

inline Vec3
Moor3d::get_force(){
    return force;
}

inline Vec3
Moor3d::get_moment(){
    return moment;
}

double
Moor3d::get_htension(int n){
    if(n>numLines){
        std::cerr << "exceeded mLines number in function 'get_htension'"<<std::endl;
    }
    return mLines[n-1].hTension; 
}

double
Moor3d::get_vtension(int n){
    if(n>numLines){
        std::cerr << "exceeded mLines number in function 'get_vtension'"<<std::endl;
    }
    return mLines[n-1].vTension;   
}

double 
Moor3d::get_catenary_length(int n) {
    if(n>numLines){
        std::cerr << "exceeded mLines number in function 'get_vtension'"<<std::endl;
    }
    return mLines[n-1].lenCatenary;   
}

#endif

