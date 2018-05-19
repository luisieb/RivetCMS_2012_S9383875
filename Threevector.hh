/*
 * Threevector.hh
 *
 *  Created on: 19/5/2018
 *      Author: ignacio
 */

#ifndef Threevector_HH_
#define Threevector_HH_

#include <iostream>
#include <math.h>

using namespace std;

class Threevector
{
    public:
	//	data member
	double _M_x ,_M_y ,_M_z ;


	//	constructor
	Threevector( double x, double y , double z)
	{
		_M_x = x;
		_M_y = y;
		_M_z = z;
	}



	//	elements of access

	double mag2() const{ return _M_x*_M_x +_M_y*_M_y +_M_z*_M_z;}
	double perp2() const{ return _M_x*_M_x +_M_y*_M_y;}

	//	magnitude and transverse component
    double mag() const { return sqrt(this -> mag2());}
    double perp() const { return sqrt(this -> perp2());}

    //	azimuth and polar angles
    double phi() const {return _M_x == 0.0 && _M_y == 0.0 ? 0.0 : atan2(_M_x,_M_y);}

    double theta() const{
    	double p = this -> perp();
    	return p == 0.0 && _M_z == 0.0 ? 0.0 : atan2(p,_M_z);
    }





// Operadores
    // suma +

    Threevector operator + (const Threevector& v1)
    {
    	Threevector result(0.0,0.0,0.0);

    	result._M_x = this-> _M_x + v1._M_x;
    	result._M_y = this-> _M_y + v1._M_y;
    	result._M_z = this-> _M_z + v1._M_z;

    	return result;
    }

    Threevector operator - (const Threevector& v1)
    {
    	Threevector result(0.0,0.0,0.0);

    	result._M_x = this-> _M_x - v1._M_x;
    	result._M_y = this-> _M_y - v1._M_y;
    	result._M_z = this-> _M_z - v1._M_z;

    	return result;
    }

    double operator * (const Threevector& v1)
    {
    	double result;

    	result = this-> _M_x * v1._M_x + this-> _M_y * v1._M_y + this-> _M_z * v1._M_z;
    	return result;
    }

    Threevector operator * (double v1)
    {
    	Threevector result(0.0,0.0,0.0);

    	result._M_x = this-> _M_x *v1;
    	result._M_y = this-> _M_y *v1;
    	result._M_z = this-> _M_z *v1;

    	return result;
    }
    Threevector operator / (double v1)
    {
    	Threevector result(0.0,0.0,0.0);

    	result._M_x = this-> _M_x /v1;
    	result._M_y = this-> _M_y /v1;
    	result._M_z = this-> _M_z /v1;

    	return result;
    }
};




#endif /* Threevector_HH_ */
