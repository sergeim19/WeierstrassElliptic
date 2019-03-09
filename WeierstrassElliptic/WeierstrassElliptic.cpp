#include <stdio.h>
#include <math.h>
#include <complex>
#include "WeierstrassElliptic.h"

int main()
{
	int type = 1;
	cdouble z(2.0, -1.0);
	cdouble q(0.0, exp(-PI));

	cdouble J = WeierstrassElliptic::JacobiTheta(type, z, q);
	printf("J = %.16f + %.16f i\n", J.real(), J.imag());
  std::getchar();
	return 0;
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

cdouble WeierstrassElliptic::WeierstrassEllipticP(
	cdouble z,
	cdouble q,
	cdouble omega1)
{
	// TODO: Validate input

	// Evaluate lattice root 1
	cdouble e1 = LatticeRoot(1, q, omega1);

	// Evaluate Weierstrass Elliptic P function
	cdouble z0 = PI * z / (2.0 * omega1);
	cdouble P = pow(
		PI * JacobiTheta(3, 0.0, q) * JacobiTheta(4, 0.0, q) * JacobiTheta(2, z0, q)/
		(2.0 * omega1 * JacobiTheta(1, z0, q)), 2);
	P += e1;

	return P;
}
