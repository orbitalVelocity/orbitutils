#include <cmath>
#include <cstdio>
#include <vector>
#include "orbitalelements.h"
#include "oeconvert.h"
#include "kepler.h"
#include "lambert.h"

int main()
{
	OrbitalElements oe0;
	oe0.sma = 1.5;
	oe0.ecc = 0.5;
	oe0.inc = 0.0 / 180.0 * M_PI;
	oe0.lan = 0.0 / 180.0 * M_PI;
	oe0.aop = 0.0 / 180.0 * M_PI;
	oe0.tra = 0.0 / 180.0 * M_PI;

	double grav_param = 1.0;
	double dt = 0.2;

	OrbitalElements oef = oe0;
	oef.tra = anomalyAfterTime(grav_param,oe0,dt);

	std::vector<double> rv0 = oe2rv(grav_param,oe0);
	std::vector<double> rvf = oe2rv(grav_param,oef);

	std::vector<double> r0 = oe2r(grav_param,oe0);
	std::vector<double> rf = oe2r(grav_param,oef);

	fprintf(stdout,"Initial position: [%.10g,%.10g,%.10g]\n",r0[0],r0[1],r0[2]);
	fprintf(stdout,"Final position: [%.10g,%.10g,%.10g]\n",rf[0],rf[1],rf[2]);

	std::vector<double> velocities = boundingVelocities(grav_param,r0,rf,dt,false);

	std::vector<double> err(6);
	err[0] = velocities[0]-rv0[3];
	err[1] = velocities[1]-rv0[4];
	err[2] = velocities[2]-rv0[5];
	err[3] = velocities[3]-rvf[3];
	err[4] = velocities[4]-rvf[4];
	err[5] = velocities[5]-rvf[5];

	fprintf(stdout,"Error vector: [%g,%g,%g,%g,%g,%g]\n",err[0],err[1],err[2],err[3],err[4],err[5]);
}
