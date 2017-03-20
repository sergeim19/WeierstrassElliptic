#include <stdio.h>
#include <math.h>
#include <complex>
#include "WeierstrassElliptic.h"

void main()
{
	int type = 4;
	cdouble z(2.0, -1.0);
	cdouble q(0.0, exp(-PI));

	cdouble J = WeierstrassElliptic::JacobiTheta(type, z, q);
	printf("J = %.16f + %.16f i\n", J.real(), J.imag());
	std::getchar();
}

cdouble WeierstrassElliptic::JacobiTheta(
	int type,
	cdouble z,
	cdouble q)
{
	cdouble J, Jnew;
	double diff;
	int n;

	switch (type)
	{
	case(1):
		J = 0.0 + 0.0 * I;
		n = 0;

		// Evaluate series.
		do
		{
			Jnew = J + 2.0 * pow(-1, n) * pow(q, pow(n + 0.5, 2)) * sin((2.0 * n + 1.0) * z);
			diff = abs(Jnew - J);
			J = Jnew;
			n++;
		} while (diff > TOL);

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
		} while (diff > TOL);

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
		} while (diff > TOL);

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
		} while (diff > TOL);

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
	cdouble root;

	switch (type)
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

cdouble WeierstrassEllipticP(
	cdouble z,
	cdouble q,
	cdouble omega1)
{
	return 0;
}