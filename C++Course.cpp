//============================================================================
// Name        : C++Course.cpp
// Author      : Luis Ignacio Estevez Banos
// Version     :
// Copyright   : LIEB
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include "Trheevector.hh"
using namespace std;

int main() {
	cout << "!!!Hello World!!!" << endl; // prints !!!Hello World!!!

	Trheevector a(0.1,0.2,0.3);
    Trheevector b(0.5,0.6,0.2);


  Trheevector  c = a+b;
  Trheevector  d = a-b;
    cout << "c =" << "("<<c._M_x << ","<< c._M_y << ","<< c._M_z << ")"<<endl;
    cout << "d =" << "("<<d._M_x << ","<< d._M_y << ","<< d._M_z << ")"<<endl;
    cout << "a*b = "<< a*b<<endl;
    cout << "a*a = "<< a*a<<"?="<<a.mag2()<<endl;
    Trheevector r = a*2.0 + b*2.0;
    cout << "a*2.0 + b*2.0 =" << "("<<r._M_x << ","<< r._M_y << ","<< r._M_z << ")"<<endl;
    Trheevector s = a/2.0 + b/2.0;
        cout << "a/2.0 + b/2.0 =" << "("<<s._M_x << ","<< s._M_y << ","<< s._M_z << ")"<<endl;

	return 0;
}

