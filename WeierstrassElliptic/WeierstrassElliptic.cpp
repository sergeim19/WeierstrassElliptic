#include <stdio.h>
#include <iostream>
#include <fstream>

#include <math.h>
#include <complex>

#include "WeierstrassElliptic.h"

void main()
{
	// Work in progress.
	return;
}

cdouble WeierstrassElliptic::JacobiTheta(
	int type,
	cdouble z,
	cdouble q)
{
	// Initialize.
	cdouble J = 0.0, Jnew = 0.0;
	double diff = 1.0e10;
	int n = 0;

	switch (type)
	{
	case(1):
		// Evaluate series.
		do
		{
			Jnew = J + 2.0 * pow(-1, n) * pow(q, pow(n + 0.5, 2)) * sin((2.0 * n + 1.0) * z);
			diff = abs(Jnew - J);
			J = Jnew;
			n++;
		} while(diff > TOL);

		break;

	case(2):
		J = 0.0 + 0.0 * I;
		n = 0;

		do
		{
			Jnew = J + 2.0 * pow(q, pow(n + 0.5, 2)) * cos((2.0 * n + 1.0) * z);
			diff = abs(Jnew - J);
			J = Jnew;
			n++;
		} while(diff > TOL);

		break;

	case(3):
		J = 1.0 + 0.0 * I;
		n = 1;

		do
		{
			Jnew = J + 2.0 * pow(q, pow(n, 2)) * cos(2.0 * n * z);
			diff = abs(Jnew - J);
			J = Jnew;
			n++;
		} while(diff > TOL);

		break;

	case(4):
		J = 1.0 * 0.0 * I;
		n = 1;

		do
		{
			Jnew = 2.0 * pow(-1, n) * pow(q, pow(n, 2)) * cos(2.0 * n * z);
			diff = abs(Jnew - J);
			J = Jnew;
			n++;
		} while(diff > TOL);

		break;
	default:// TODO: Wrong input.
		break;
	}

	return J;
}

cdouble WeierstrassElliptic::LatticeRoot(
	int type,
	cdouble q,
	cdouble omega1)
{
	// TODO: Validate input.
	cdouble root = 0.0;

	switch(type)
	{
	case(1):
		root = pow(PI, 2) / (12.0 * pow(omega1, 2)) * (pow(JacobiTheta(2, 0, q), 4) + 2.0 * pow(JacobiTheta(4, 0, q), 4));
		break;

	case(2):
		root = pow(PI, 2) / (12.0 * pow(omega1, 2)) * (pow(JacobiTheta(2, 0, q), 4) - pow(JacobiTheta(4, 0, q), 4));
		break;

	case(3):
		root = pow(PI, 2) / (12.0 * pow(omega1, 2)) * (2.0 * pow(JacobiTheta(2, 0, q), 4) + pow(JacobiTheta(4, 0, q), 4));
		break;

	default: // TODO: Wrong input.
		break;
	}

	return root;
}

cdouble WeierstrassElliptic::WeierstrassP(
	cdouble z,
	cdouble q,
	cdouble omega1)
{
	// TODO: Validate input.

	// Evaluate lattice root.
	cdouble LatticeRoot1 = LatticeRoot(1, q, omega1);
	cdouble Numerator = PI * JacobiTheta(3, 0.0, q) * JacobiTheta(4, 0.0, q) * JacobiTheta(2, (PI * z) / (2.0*omega1), q);
	cdouble Denominator = 2.0 * omega1 * JacobiTheta(1, (PI * z) / (2.0 * omega1), q);

	// Evaluate WeierstrassP.
	cdouble WeierstrassP = pow(Numerator / Denominator, 2) + LatticeRoot1;

	return WeierstrassP;
}