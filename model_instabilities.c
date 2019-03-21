#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

/* @brief: A recipe to treat local Toomre instabilities. Each galaxy is included inside a do-while loop which finishes either when all rings of the
* galaxy are stable (determined by sum which equals 0 (1) if galaxy stable (unstable)) or when a limit (imax - free parameter) is reached. Inside the
* ring stable (unstable)). For each ring we calculate the Toomre Q parameter for gas (Qgas) and stars (Qstars) (if both components are present,
* else solve single component instabilities)  and combine them (Qtot) based on Romeo & Wiegert 2011. if Qtot < 1 the rings is considered to be
* unstable them mass is tranferred to adjucent rings in proportion (1/3 to j+1 ring and 2/3 to j-1 in order to conserve angular momentum). Then we
* identify three different instability scenarios:
* 1) Stable stellar ring - unstable mass is in the form of gas: invoke fix_gas_instabilities where a fraction (fmig - free parameter) of the unstable
* gas migrates to adjacent gas rings (along with metals) and the remaining is converted into stars. If the unstable ring is the innermost then all
* unstable gas feeds the SMBH. TODO add SF feedback
* 2) Stable gaseous ring - unstable mass is in the form of stars: invoke fix_stellar_instabilities where unstable stars migrate to adjacent
* stellar rings (along with metals). For rings with j < threshold (BulgeThresh - free parameter) all unstable stars are transferred to the bulge .
* 3) Unstable stellar and gaseous rings: invoke fix_gas_instabilities, recalculate SigmaStarsRing, invoke fix_stellar_instabilities. */

void check_disk_instability_DI (int p, double dt) {
	double RdStars, Vmax, kappa, sigmaGas, sigmaStars, SigmaCritGas, SigmaCritStars, SigmaStarRings, SigmaGasRings, q, s, W, Qstars, Qgas, Qtot,
	       Qmin, Qcrit, UnstableMass;
	double G      = 43.02; // Gravitational constant in Mpc 1e-10 Msun^-1 km^2 s^-2
	double fmig   = 0.3; // Fraction of unstable gas migrated to adjacent annuli.
	int    k, j;
	char   *phase = "one";

	// Calculate the scale length and Vmax.
	RdStars = get_stellar_disk_radius (p) / 3.0;
	if (Gal[p].Type != 0) {
		Vmax = Gal[p].InfallVmax;
	} else {
		Vmax = Gal[p].Vmax;
	}

	for (k  = 1; k < 2; k ++) {
		printf ("%d\n", k);

		// Star a loop from ring 10 to ring 0 and check the stability of each ring.
		for (j = 0; j < RNUM -1; j ++) {
			if (j != 0) { // Calculate surface densities and convert them from (1e10 M_sun/h / (Mpc/h)^2) to (1e10 M_sun / Mpc^2)
				SigmaGasRings  = (Gal[p].ColdGasRings[j] * Hubble_h) /
				                 (M_PI * (RingRadius[j] * RingRadius[j] - RingRadius[j - 1] * RingRadius[j - 1]));
				SigmaStarRings = (Gal[p].DiskMassRings[j] * Hubble_h) /
				                 (M_PI * (RingRadius[j] * RingRadius[j] - RingRadius[j - 1] * RingRadius[j - 1]));
			} else { // Follow the same process for the central ring.
				SigmaGasRings  = (Gal[p].ColdGasRings[j] * Hubble_h) / (M_PI * RingRadius[j] * RingRadius[j]);
				SigmaStarRings = (Gal[p].DiskMassRings[j] * Hubble_h) / (M_PI * RingRadius[j] * RingRadius[j]);
			}

			// Calculate the required quantities for the local Toomre Q parameter.
			sigmaGas   = 11.0; // Gas radial velocity dispersion in km s^-1
			kappa      = (sqrt (2) * Vmax * Hubble_h) / RingRadius[j]; // Epicyclic frequency in km s^-1 Mpc^-1
			sigmaStars = 0.87 * sqrt (M_PI * G * SigmaStarRings * RdStars / Hubble_h); // Leroy+08 velocity dispersion in km s^-1
//			sigmaStars = 0.5 * Vmax * exp ((- 1 * RingRadius[j]) / (2 * RdStars)); // Bottema+93 velocity dispersion in km s^-1

			// Calculate the Toomre Q parameters
			Qgas   = (kappa * sigmaGas) / (M_PI * G * SigmaGasRings);
			Qstars = (kappa * sigmaStars) / (3.36 * G * SigmaStarRings);

			if (Gal[p].ColdGasRings[j] > 1.0e-6 && Gal[p].DiskMassRings[j] > 1.0e-6) {

				// Calculate Qtot if both components are present. Else solve single component instabilities.
				q    = Qgas / Qstars;
				s    = sigmaGas / sigmaStars;
				W    = (2 * s) / (1 + (s * s)); // Weight factor.
				Qmin = 1.0 + W; // The minimum Q value for stellar and gas rings. For Q > Qmin the component is stable.

				// Calculate the conditional function for total Q
				if (q <= 1.0) {
					Qcrit = Qstars / (Qstars - W);
					Qtot  = (Qstars * Qgas) / (W * Qgas + Qstars);
				} else if (q > 1.0) {
					Qcrit = Qgas / (Qgas - W);
					Qtot  = (Qstars * Qgas) / (W * Qstars + Qgas);
				}

				// Check if the ring is unstable
				if (Qtot < 1.0 && Qtot > 0.0) {
					if (Qstars > Qmin && Qgas < Qmin) { // Then q must be <= 1.0 so the unstable mass is in the form of gas.
						// Since this ring is unstable we will transfer mass to the previous one, so we need to re-check its stability
						// Calculate the critical surface mass density of the ring.
						SigmaCritGas = (kappa * sigmaGas) / (M_PI * G * Qcrit);
						// Invoke a function to deal with instabilities in the gas ring.
						fix_gas_instabilities (p, j, SigmaGasRings, SigmaCritGas, fmig);

					} else if (Qgas > Qmin && Qstars < Qmin) { // Then q must be > 1.0 so the unstable mass is in the form of stars.
						// Calculate the critical surface mass density of the ring.
						SigmaCritStars = (kappa * sigmaStars) / (3.36 * G * Qcrit);
						// Invoke a function to deal with instabilities in the stellar ring.
						fix_stellar_instabilities (p, j, SigmaStarRings, SigmaCritStars, "one");
						fix_stellar_instabilities (p, j, SigmaStarRings, SigmaCritStars, "two");

					} else if (Qstars < Qmin && Qgas < Qmin) { // First raise Qgas to Qmin and then Qstars.
						// Calculate the critical surface mass density of the ring.
						SigmaCritGas = (kappa * sigmaGas) / (M_PI * G * Qmin);
						// Invoke a function to deal with instabilities in the gas ring.
						fix_gas_instabilities (p, j, SigmaGasRings, SigmaCritGas, fmig);

						// Re-calculate SigmaStarRings and raise as necessary.
						SigmaStarRings = (Gal[p].DiskMassRings[j] * Hubble_h) /
						                 (M_PI * (RingRadius[j] * RingRadius[j] - RingRadius[j - 1] * RingRadius[j - 1]));
						// Calculate the critical surface mass density of the ring.
						SigmaCritStars = (kappa * sigmaStars) / (3.36 * G * Qmin);
						// Invoke a function to deal with instabilities in the stellar ring.
						fix_stellar_instabilities (p, j, SigmaStarRings, SigmaCritStars, "one");
						fix_stellar_instabilities (p, j, SigmaStarRings, SigmaCritStars, "two");
					}
				} // if (Qtot < 1.0 && Qtot > 0.0)
			} else { // if (Gal[p].ColdGasRings[j] > 1.0e-6 && Gal[p].DiskMassRings[j] > 1.0e-6)
				Qmin = 1.0;
				if (Gal[p].DiskMassRings[j] == 0.0 && Gal[p].ColdGasRings[j] > 1e-6 && Qgas < Qmin && Qgas > 0.0) { // The unstable mass is in the
					// form of gas.
					// Calculate the critical surface mass density of the ring.
					SigmaCritGas = (kappa * sigmaGas) / (M_PI * G);
					// Invoke a function to deal with instabilities in the gas ring.
					fix_gas_instabilities (p, j, SigmaGasRings, SigmaCritGas, fmig);
				} else if (Gal[p].ColdGasRings[j] == 0.0 && Gal[p].DiskMassRings[j] > 1e-6 && Qstars < Qmin && Qstars > 0.0) { // The unstable
					// mass is in the form of stars.
					// Calculate the critical mass of the ring
					SigmaCritStars = (kappa * sigmaStars) / (3.36 * G);
					// Invoke a function to deal with instabilities in the stellar ring.
					fix_stellar_instabilities (p, j, SigmaStarRings, SigmaCritStars, "one");
					fix_stellar_instabilities (p, j, SigmaStarRings, SigmaCritStars, "two");
				}
			}
		} // for (j  = RNUM - 2; j >= 0; j --)
		strcpy (phase, "two");
	}
	printf ("%wefs\n");
} // end check_disk_instability_DI

void fix_gas_instabilities (int p, int j, double SigmaGasRings, double SigmaCritGas, double fmig) {
	double UnstableMass;

	if (j != 0) {
		// Convert unstable mass from (1e10 M_sun/h^2) to (1e10 M_sun/h)
		UnstableMass = (SigmaGasRings - SigmaCritGas) * M_PI * (RingRadius[j] * RingRadius[j] - RingRadius[j - 1] * RingRadius[j - 1]) / Hubble_h;
	} else {
		UnstableMass = (SigmaGasRings - SigmaCritGas) * M_PI * RingRadius[j] * RingRadius[j] / Hubble_h;
	}

	if (UnstableMass > 0.0) {
		if (UnstableMass > Gal[p].ColdGasRings[j]) {
			UnstableMass = Gal[p].ColdGasRings[j];
		}

		if (j != 0) { // fmig*UnstableMass is moved to adjacent annuli and the remaining is transformed into stars.
			Gal[p].ColdGasRings[j - 1] += (2 * fmig * UnstableMass) / 3.;
			Gal[p].ColdGasRings[j + 1] += (fmig * UnstableMass) / 3.;
			Gal[p].DiskMassRings[j] += (1.0 - fmig) * UnstableMass;
			Gal[p].DiskMass += (1.0 - fmig) * UnstableMass;
			Gal[p].ColdGasRings[j] -= UnstableMass;
			Gal[p].ColdGas -= (1.0 - fmig) * UnstableMass;

		} else { // If you remove mass from the last ring (i.e., 0), then transfer it to the black hole.
			Gal[p].BlackHoleMass += UnstableMass;
			Gal[p].ColdGasRings[j] -= UnstableMass;
			Gal[p].ColdGas -= UnstableMass;
		}
	}
	mass_checks (p, "model_starformation_and_feedback.c", __LINE__);
} // fix_gas_instabilities

void fix_stellar_instabilities (int p, int j, double SigmaStarRings, double SigmaCritStars, char phase[]) {
	double UnstableMass, fbulge, fcarry = 0.3, fmig = 0.3;
	fbulge = 1 / (2 - fcarry);

	if (j != 0) {
		// Convert unstable mass from (1e10 M_sun/h^2) to (1e10 M_sun/h)
		UnstableMass = (SigmaStarRings - SigmaCritStars) * M_PI * (RingRadius[j] * RingRadius[j] - RingRadius[j - 1] * RingRadius[j - 1]) / Hubble_h;
	} else {
		UnstableMass = (SigmaStarRings - SigmaCritStars) * M_PI * RingRadius[j] * RingRadius[j] / Hubble_h;
	}

	if (UnstableMass > 0.0) {
		if (UnstableMass > Gal[p].DiskMassRings[j]) {
			UnstableMass = Gal[p].DiskMassRings[j];
		}
		if (strcmp (phase, "one") == 0) {
			// Move 2/3 of the unstable mass to the inner and 1/3 to the outer ring along with the associated metals.
			Gal[p].DiskMassRings[j - 1] += fmig * ((2 * UnstableMass) / 3.);
			Gal[p].DiskMassRings[j + 1] += fmig * (UnstableMass / 3.);
			Gal[p].DiskMassRings[j] -= fmig * UnstableMass;

			Gal[p].BulgeMass += UnstableMass;
			Gal[p].DiskMass -= UnstableMass;
			Gal[p].DiskMassRings[j] -= UnstableMass;
		} else {

		}
	}
	mass_checks (p, "model_starformation_and_feedback.c", __LINE__);
} // end fix_stellar_instabilities

void update_sizes_from_instabilities (int p, int j, double UnstableMass, double fcarry, char process[]) {
	double frac, bulgesize, fint = 4.0;
	int    ii;

	if (strcmp (process, "Add") == 0) {
		if (Gal[p].DiskMass + UnstableMass > 1.e-8) {
			for (ii = 0; ii < 3; ii ++) {
				Gal[p].DiskSpin[ii] = ((Gal[p].DiskSpin[ii]) * (Gal[p].DiskMass) + UnstableMass * Gal[p].ColdGasSpin[ii]) /
				                      (Gal[p].DiskMass + UnstableMass);
			}
		}
	} // if (strcmp (process, "Add") == 0)

	// Instabilities are assumed to NOT remove angular momentum. Since the mass of the disk changes, the spin is changed by the same amount to
	// keep angular momentum constant.
	if (strcmp (process, "Subtract") == 0) {
		frac = UnstableMass / Gal[p].DiskMass;

		for (ii = 0; ii < 3; ii ++) {
			if (frac == 1) { //everything transferred to the bulge
				Gal[p].DiskSpin[ii] = 0;
			} else {
				Gal[p].DiskSpin[ii] = Gal[p].DiskSpin[ii] / (1 - frac);
			}
		}
		// Size of newly formed bulge, which consists of the stellar mass transfered from the disk. This is calculated from eq 35 in Guo+11.
		bulgesize = bulge_from_disk (frac) * Gal[p].DiskRadius / 3.;
		if (Gal[p].BulgeMass < TINY_MASS) {
			// if previous Bulge Mass = 0 -> bulge size is given directly from newly formed bulge */
			Gal[p].BulgeSize = bulgesize;
		} else {
			Gal[p].BulgeSize = (Gal[p].BulgeMass + UnstableMass) * (Gal[p].BulgeMass + UnstableMass) /
			                   (Gal[p].BulgeMass * Gal[p].BulgeMass / Gal[p].BulgeSize + UnstableMass * UnstableMass / bulgesize +
			                    fint * Gal[p].BulgeMass * UnstableMass / (Gal[p].BulgeSize + bulgesize));
		}
	} // if (strcmp (process, "Subtract") == 0)

	// Update the disc scale length
	if (DiskRadiusModel == 0) {
		Gal[p].DiskRadius = get_stellar_disk_radius (p);
	}
} // end update_sizes_from_instabilities