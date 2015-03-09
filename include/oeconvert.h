#ifndef OECONVERT
#define OECONVERT

#include <vector>

/**
 * Converts orbital elements to position.
 * @param	grav_param	The gravitational parameter of the two-body system.
 * @param	oe	The orbital elements to be converted from.
 * @returns	A vector corresponding to the position.
 */
std::vector<double> oe2r(double grav_param, OrbitalElements oe);

/**
 * Converts orbital elements to position and velocity.
 * @param	grav_param	The gravitational parameter of the two-body system.
 * @param	oe	The orbital elements to be converted from.
 * @returns	A vector corresponding to the concatenated position and velocity.
 */
std::vector<double> oe2rv(double grav_param, OrbitalElements oe);

/**
 * Converts position and velocity to orbital elements
 * @param	grav_param	The gravitational parameter of the two-body system.
 * @param	rv	The concatenated position and velocity to convert from.
 * @returns	The orbital elements.
 */
OrbitalElements rv2oe(double grav_param, std::vector<double> rv);

#endif
