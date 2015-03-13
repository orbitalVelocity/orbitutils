#ifndef LAMBERT
#define LAMBERT

#include <vector>

/**
 * Computes the initial velocity and final velocity for the orbit boundary value problem.
 * @param	grav_param	The gravitational parameter of the two-body system.
 * @param	r0	The initial position of the secondary body.
 * @param	rf	The final position of the secondary body.
 * @param	dt	The transfer time.
 * @param	long_way	Set to true if the long-way transfer is desired.
 * @returns	The initial velocity and final velocities packed into a single vector.
 */
std::vector<double> boundingVelocities(double grav_param, std::vector<double> r0,
	std::vector<double> rf, double dt, bool long_way);

#endif
