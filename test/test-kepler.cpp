#include <cmath>
#include <stdio.h>
#include "orbitalelements.h"
#include "kepler.h"

int main()
{
	double grav_param = 1.0;

	OrbitalElements oe;
	oe.sma = 1.0;
	oe.ecc = 0.0;
	oe.inc = 0.0;
	oe.lan = 0.0;
	oe.aop = 0.0;
	oe.tra = 0.0;

	double period = 2*M_PI*sqrt(pow(oe.sma,3.0)/grav_param);

	FILE* log = stdout;
	double correct;
	double calculated;
	double error;
	double tolerance = 1E-10;

	// test timeUntilAnomaly

	fprintf(log,"Test timeUntilAnomaly, circular, f = 0 to f = pi/2:\n");
	correct = period/4.0;
	calculated = timeUntilAnomaly(grav_param,oe,M_PI/2.0);
	error = std::abs(calculated-correct)/correct;
	fprintf(log, "\tCorrect: %10f  Calculated: %10f  Error Fraction: %10f  (%s)\n",correct,calculated,error,(error < tolerance)?"PASS":"FAIL");

	fprintf(log,"Test timeUntilAnomaly, circular, f = 0 to f = pi:\n");
	correct = period/2.0;
	calculated = timeUntilAnomaly(grav_param,oe,M_PI);
	error = std::abs(calculated-correct)/correct;
	fprintf(log, "\tCorrect: %10f  Calculated: %10f  Error Fraction: %10f  (%s)\n",correct,calculated,error,(error < tolerance)?"PASS":"FAIL");

	fprintf(log,"Test timeUntilAnomaly, circular, f = 3*pi/2 to f = 0:\n");
	oe.tra = 3*M_PI/2;
	correct = period/4.0;
	calculated = timeUntilAnomaly(grav_param,oe,0.0);
	error = std::abs(calculated-correct)/correct;
	fprintf(log, "\tCorrect: %10f  Calculated: %10f  Error Fraction: %10f  (%s)\n",correct,calculated,error,(error < tolerance)?"PASS":"FAIL");
}
