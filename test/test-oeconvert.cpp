#include <cmath>
#include <cstdio>
#include <vector>
#include "orbitalelements.h"
#include "oeconvert.h"

int main()
{
	double grav_param = 2;

	std::vector<double> rv0(6);
	rv0[0] =   10;
	rv0[1] =    2;
	rv0[2] =    3;
	rv0[3] = 0.01;
	rv0[4] =  0.3;
	rv0[5] = 0.05;

	OrbitalElements oe0;
	oe0.sma = 7.04909391142183;
	oe0.ecc = 0.5561545963072888;
	oe0.inc = 17.294397104279813	/ 180.0 * M_PI;
	oe0.lan = 300.4342357534433	/ 180.0 * M_PI;
	oe0.aop = 264.641995659073	/ 180.0 * M_PI;
	oe0.tra = 167.0393678503284	/ 180.0 * M_PI;

	OrbitalElements oe = rv2oe(grav_param,rv0);
	fprintf(stdout,"sma   error: %g\n",(oe.sma-oe0.sma)/oe0.sma);
	fprintf(stdout,"ecc   error: %g\n",(oe.ecc-oe0.ecc)/oe0.ecc);
	fprintf(stdout,"inc   error: %g\n",(oe.inc-oe0.inc)/oe0.inc);
	fprintf(stdout,"lan   error: %g\n",(oe.lan-oe0.lan)/oe0.lan);
	fprintf(stdout,"aop   error: %g\n",(oe.aop-oe0.aop)/oe0.aop);
	fprintf(stdout,"tra   error: %g\n",(oe.tra-oe0.tra)/oe0.tra);

	std::vector<double> rv = oe2rv(grav_param,oe);
	fprintf(stdout,"rv[0] error: %g\n",(rv[0]-rv0[0])/rv0[0]);
	fprintf(stdout,"rv[1] error: %g\n",(rv[1]-rv0[1])/rv0[1]);
	fprintf(stdout,"rv[2] error: %g\n",(rv[2]-rv0[2])/rv0[2]);
	fprintf(stdout,"rv[3] error: %g\n",(rv[3]-rv0[3])/rv0[3]);
	fprintf(stdout,"rv[4] error: %g\n",(rv[4]-rv0[4])/rv0[4]);
	fprintf(stdout,"rv[5] error: %g\n",(rv[5]-rv0[5])/rv0[5]);
}
