#pragma once

//#include <complex>

typedef std::complex<double> cdouble;
#define TOL 1.0e-15
#define PI 3.1415926535897932384
#define I cdouble(0.0, 1.0)

class WeierstrassElliptic
{
public:
	static cdouble JacobiTheta(
		int type,
		cdouble z,
		cdouble q);
};