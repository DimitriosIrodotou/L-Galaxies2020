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
* do-while loop there is a  for loop which runs from ring 10 to ring 0 and checks the local stability of each ring (determined by n[j] = 0 (1) if
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
	       Qmin, Qcrit, UnstableMass, fRing[RNUM - 1];
	double G    = 43.02; // Gravitational constant in Mpc 1e-10 Msun^-1 km^2 s^-2
	double fmig = 0.3; // Fraction of unstable gas migrated to adjacent annuli.
	int    i    = 0, imax = 100, k, j, n[RNUM - 1], sum;

	do {
		// Set the instability flag for this galaxy to zero.
		sum = 0;

		// Calculate the scale length and Vmax.
		RdStars = get_stellar_disk_radius (p) / 3.0;
		if (Gal[p].Type != 0) {
			Vmax = Gal[p].InfallVmax;
		} else {
			Vmax = Gal[p].Vmax;
		}

		// Star a loop from ring 10 to ring 0 and check the stability of each ring.
		for (j = RNUM - 2; j >= 0; j --) {
			// Set the instability flag for this ring to zero.
			n[j]     = 0;
			fRing[j] = 0.0;

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
						n[j] = 1;
						// Calculate the critical surface mass density of the ring.
						SigmaCritGas = (kappa * sigmaGas) / (M_PI * G * Qcrit);
						// Invoke a function to deal with instabilities in the gas ring.
						fix_gas_instabilities (p, j, SigmaGasRings, SigmaCritGas, fmig);

					} else if (Qgas > Qmin && Qstars < Qmin) { // Then q must be > 1.0 so the unstable mass is in the form of stars.
						n[j] = 1;
						// Calculate the critical surface mass density of the ring.
						SigmaCritStars = (kappa * sigmaStars) / (3.36 * G * Qcrit);
						// Invoke a function to deal with instabilities in the stellar ring.
						fix_stellar_instabilities (p, j, SigmaStarRings, SigmaCritStars);

					} else if (Qstars < Qmin && Qgas < Qmin) { // First raise Qgas to Qmin and then Qstars.
						n[j] = 1;
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
						fix_stellar_instabilities (p, j, SigmaStarRings, SigmaCritStars);
					}
				} // if (Qtot < 1.0 && Qtot > 0.0)
			} else { // if (Gal[p].ColdGasRings[j] > 1.0e-6 && Gal[p].DiskMassRings[j] > 1.0e-6)
				Qmin = 1.0;
				if (Gal[p].DiskMassRings[j] == 0.0 && Gal[p].ColdGasRings[j] > 1e-6 && Qgas < Qmin && Qgas > 0.0) { // The unstable mass is in the
					// form of gas.
					n[j] = 1;
					// Calculate the critical surface mass density of the ring.
					SigmaCritGas = (kappa * sigmaGas) / (M_PI * G);
					// Invoke a function to deal with instabilities in the gas ring.
					fix_gas_instabilities (p, j, SigmaGasRings, SigmaCritGas, fmig);
				} else if (Gal[p].ColdGasRings[j] == 0.0 && Gal[p].DiskMassRings[j] > 1e-6 && Qstars < Qmin && Qstars > 0.0) { // The unstable
					// mass is in the form of stars.
					n[j] = 1;
					// Calculate the critical mass of the ring
					SigmaCritStars = (kappa * sigmaStars) / (3.36 * G);
					// Invoke a function to deal with instabilities in the stellar ring.
					fix_stellar_instabilities (p, j, SigmaStarRings, SigmaCritStars);
				}
			}
		} // for (j  = RNUM - 2; j >= 0; j --)
		// Once you finished the entire disk (11 rings) check if each ring is stable i.e., n[j] = 0 for all j. If not repeat the whole process
		// until n[j] = 0 for all j but no more than 30 times. Note: in the above for loop j goes from 10 to 0, so after the last iteration its
		// value is -1.
		if (j < 0) {
			for (k = 0; k < RNUM - 1; k ++) {
				sum += n[k];
			}
		}
		i += 1;
	} while (i < imax && sum > 0);
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
			update_sizes_from_instabilities (p, j, ((1.0 - fmig) * UnstableMass), "Add");
//			update_metals_from_instabilities (p, j, ((1.0 - fmig) * UnstableMass), "ColdGas", "DiskMass");

//			// Deal with the associated feedback.
//			update_h2fraction (p);
//			sfe    = SfrEfficiency * UnitTime_in_years / Hubble_h;
//			if (SFRtdyn == 1) {
//				sfe = (sfe / tdyn) / UnitTime_in_years * Hubble_h * 1e7;
//			} // for star formation rate proportional to 1/t_dyn
//			strdot += strdotr[j];
//			for (j = 0; j < RNUM; j ++) {
//				if (strdotr[j] < 0.0) {
//					strdotr[j] = 0.;
//				}
//				starsRings[j] = strdotr[j] * dt;
//				if (starsRings[j] < 0.0) {
//					starsRings[j] = 0.0;
//				}
//			}
//#ifdef COMPUTE_SPECPHOT_PROPERTIES
//#ifndef POST_PROCESS_MAGS
//			if (Gal[p].ColdGas > 0.) {
//				metallicitySF = metals_total (Gal[p].MetalsColdGas) / Gal[p].ColdGas;
//			} else {
//				metallicitySF = 0.;
//			}
//#endif // COMPUTE_SPECPHOT_PROPERTIES
//#endif // POST_PROCESS_MAGS
//#ifndef FEEDBACK_COUPLED_WITH_MASS_RETURN
//			update_stars_due_to_reheat (p, centralgal, &UnstableMass, starsRings);
//#endif //FEEDBACK_COUPLED_WITH_MASS_RETURN
//			Gal[p].SfrRings[j] += starsRings[j] / (dt * STEPS);
//			update_from_star_formation (p, UnstableMass, starsRings, "insitu", nstep);
//
//#ifdef COMPUTE_SPECPHOT_PROPERTIES
//#ifndef POST_PROCESS_MAGS
//			//  Update the luminosities due to the stars formed
//			add_to_luminosities (p, UnstableMass, time, dt, metallicitySF);
//#endif // POST_PROCESS_MAGS
//#endif // COMPUTE_SPECPHOT_PROPERTIES
//			update_massweightage (p, UnstableMass, time);
//
//#ifndef FEEDBACK_COUPLED_WITH_MASS_RETURN
//			SN_feedback (p, centralgal, UnstableMass, starsRings, "ColdGas");
//#endif // FEEDBACK_COUPLED_WITH_MASS_RETURN

		} else { // If you remove mass from the last ring (i.e., 0), then transfer it to the black hole.
			Gal[p].BlackHoleMass += UnstableMass;
			Gal[p].ColdGasRings[j] -= UnstableMass;
			Gal[p].ColdGas -= UnstableMass;
//			update_metals_from_instabilities (p, j, UnstableMass, "ColdGas", "BlackHoleMass"); // removes only
		}
	}
	mass_checks (p, "model_starformation_and_feedback.c", __LINE__);
} // fix_gas_instabilities

void fix_stellar_instabilities (int p, int j, double SigmaStarRings, double SigmaCritStars) {
	double UnstableMass, fractionRings[RNUM];
	int    BulgeThresh = 2;

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
		if (j > BulgeThresh) { // Move 2/3 of the unstable mass to the inner and 1/3 to the outer ring along with the associated metals.
			Gal[p].DiskMassRings[j - 1] += (2 * UnstableMass) / 3.;
			Gal[p].DiskMassRings[j + 1] += UnstableMass / 3.;
			Gal[p].DiskMassRings[j] -= UnstableMass;
//			update_metals_from_instabilities (p, j, UnstableMass, "DiskMass", "DiskMass");
		} else { // If you remove mass from the last ring (i.e., 0), then transfer it to the bulge.
			Gal[p].BulgeMass += UnstableMass;
			Gal[p].DiskMass -= UnstableMass;
			Gal[p].DiskMassRings[j] -= UnstableMass;
			update_sizes_from_instabilities (p, j, UnstableMass, "Subtract");
//			update_metals_from_instabilities (p, j, UnstableMass, "DiskMass", "BulgeMass");
		}
	}
	mass_checks (p, "model_starformation_and_feedback.c", __LINE__);
} // end fix_stellar_instabilities

void update_sizes_from_instabilities (int p, int j, double UnstableMass, char process[]) {
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

//void update_metals_from_instabilities (int p, int j, double UnstableMass, char from[], char to[]) {
//
//	// Variables for TOTAL quantities
//	double Mass, fraction, Metals[NUM_METAL_CHANNELS];
//	int    mm, jj;
//#ifdef DETAILED_METALS_AND_MASS_RETURN
//#ifdef INDIVIDUAL_ELEMENTS
//	int ee;
//	double Yield[NUM_ELEMENTS];
//#endif // INDIVIDUAL_ELEMENTS
//#endif // DETAILED_METALS_AND_MASS_RETURN
//
//	// Variables for SFH
//#ifdef STAR_FORMATION_HISTORY
//	int ii;
//  double sfh_Mass[SFH_NBIN], sfh_MassRings[RNUM][SFH_NBIN], sfh_Metals[SFH_NBIN][NUM_METAL_CHANNELS],
//  sfh_MetalsRings[RNUM][SFH_NBIN][NUM_METAL_CHANNELS];
//#ifdef DETAILED_METALS_AND_MASS_RETURN
//#ifdef INDIVIDUAL_ELEMENTS
//  double sfh_Elements[SFH_NBIN][NUM_ELEMENTS],sfh_ElementsRings[RNUM][SFH_NBIN][NUM_ELEMENTS];
//#endif // INDIVIDUAL_ELEMENTS
//#endif // DETAILED_METALS_AND_MASS_RETURN
//#endif // STAR_FORMATION_HISTORY
//
//	// Variables for LUMINOSITIES
//#ifdef COMPUTE_SPECPHOT_PROPERTIES
//#ifndef POST_PROCESS_MAGS
//	int ll, outputbin;
//#ifdef OUTPUT_REST_MAGS
//	double Lum[NMAG][NOUT],YLum[NMAG][NOUT];
//#endif // OUTPUT_REST_MAGS
//#ifdef COMPUTE_OBS_MAGS
//	double ObsLum[NMAG][NOUT], ObsYLum[NMAG][NOUT];
//#ifdef OUTPUT_MOMAF_INPUTS
//	double dObsLum[NMAG][NOUT],dObsYLum[NMAG][NOUT];
//#endif // OUTPUT_MOMAF_INPUTS
//#endif // COMPUTE_OBS_MAGS
//#endif // POST_PROCESS_MAGS
//#endif // COMPUTE_SPECPHOT_PROPERTIES
//
//	// Variables for RINGS
//	double MassRings[RNUM], MetalsRings[RNUM][NUM_METAL_CHANNELS];
//#ifdef DETAILED_METALS_AND_MASS_RETURN
//#ifdef INDIVIDUAL_ELEMENTS
//	double YieldRings[RNUM][NUM_ELEMENTS];
//#endif // INDIVIDUAL_ELEMENTS
//#endif //DETAILED_METALS_AND_MASS_RETURN
//
//
//	// Subtract from component
//	if (strcmp (from, "ColdGas") == 0) {
//		fraction = UnstableMass / Gal[p].ColdGasRings[j];
//		Metals   = fraction * Gal[p].MetalsColdGasRings[jj][mm];
//
//#ifdef INDIVIDUAL_ELEMENTS
//		for(ee=0;ee<NUM_ELEMENTS;ee++)
//	 {
//	   Yield[ee] += fractionRings[jj] * Gal[p].ColdGasRings_elements[jj][ee];
//	   YieldRings[jj][ee] = fractionRings[jj] * Gal[p].ColdGasRings_elements[jj][ee];
//	 }
//#endif
//
//		// if there is SF, gas goes to stars into the last sfh bin
//#ifdef STAR_FORMATION_HISTORY
//		for (ii=0; ii<Gal[p].sfh_ibin; ii++)
//	   {
//		   sfh_Mass[ii]=0.;
//		   for (jj=0;jj<RNUM;jj++)
//		   {
//			   sfh_MassRings[jj][ii]=0.;
//			   for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
//				   sfh_MetalsRings[jj][ii][mm] = 0.;
//#ifdef INDIVIDUAL_ELEMENTS
//			   for(ee=0;ee<NUM_ELEMENTS;ee++)
//				   sfh_ElementsRings[jj][ii][ee]=0.;
//#endif
//		   }
//		   for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
//			   sfh_Metals[ii][mm] = 0.;
//#ifdef INDIVIDUAL_ELEMENTS
//		   for(ee=0;ee<NUM_ELEMENTS;ee++)
//			   sfh_Elements[ii][ee]=0.;
//#endif
//	   }
//
//	   for (jj=0;jj<RNUM;jj++)
//	   {
//		   sfh_Mass[Gal[p].sfh_ibin]+=fractionRings[jj]*Gal[p].ColdGasRings[jj];
//		   sfh_MassRings[jj][Gal[p].sfh_ibin]=fractionRings[jj]*Gal[p].ColdGasRings[jj];
//		   for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
//		   {
//			   sfh_Metals[Gal[p].sfh_ibin][mm] += fractionRings[jj]*Gal[p].MetalsColdGasRings[jj][mm];
//			   sfh_MetalsRings[jj][Gal[p].sfh_ibin][mm] = fractionRings[jj]*Gal[p].MetalsColdGasRings[jj][mm];
//		   }
//#ifdef INDIVIDUAL_ELEMENTS
//		   for(ee=0;ee<NUM_ELEMENTS;ee++)
//		   {
//			   sfh_Elements[Gal[p].sfh_ibin][ee]+=fractionRings[jj]*Gal[p].ColdGasRings_elements[jj][ee];
//			   sfh_ElementsRings[jj][Gal[p].sfh_ibin][ee]=fractionRings[jj]*Gal[p].ColdGasRings_elements[jj][ee];
//		   }
//#endif
//	   }
//#endif //STAR_FORMATION_HISTORY
//	} else if (strcmp (cq, "DiskMass") == 0) {
//		Mass += fractionRings[jj] * Gal[p].DiskMassRings[jj];
//		MassRings[jj] = fractionRings[jj] * Gal[p].DiskMassRings[jj];
//		for (mm = 0; mm < NUM_METAL_CHANNELS; mm ++) {
//			Metals[mm] += fractionRings[jj] * Gal[p].MetalsDiskMassRings[jj][mm];
//			MetalsRings[jj][mm] = fractionRings[jj] * Gal[p].MetalsDiskMassRings[jj][mm];
//		}
//#ifdef INDIVIDUAL_ELEMENTS
//		for(ee=0;ee<NUM_ELEMENTS;ee++)
//		{
//		Yield[ee] += fractionRings[jj]*Gal[p].DiskMassRings_elements[jj][ee];
//		YieldRings[jj][ee] = fractionRings[jj]*Gal[p].DiskMassRings_elements[jj][ee];
//		}
//#endif
//
//#ifdef STAR_FORMATION_HISTORY
//		for (ii=0; ii<=Gal[p].sfh_ibin; ii++)
//		{
//		if(Gal[p].sfh_DiskMass[ii]>0.)
//		 {
//		   for (jj=0;jj<RNUM;jj++)
//			 {
//			   sfh_Mass[ii]+=fractionRings[jj]*Gal[p].sfh_DiskMassRings[jj][ii];
//			   sfh_MassRings[jj][ii]=fractionRings[jj]*Gal[p].sfh_DiskMassRings[jj][ii];
//
//			   for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
//			   {
//				   sfh_Metals[ii][mm] += fractionRings[jj]*Gal[p].sfh_MetalsDiskMassRings[jj][ii][mm];
//				   sfh_MetalsRings[jj][ii][mm] = fractionRings[jj]*Gal[p].sfh_MetalsDiskMassRings[jj][ii][mm];
//			   }
//#ifdef INDIVIDUAL_ELEMENTS
//			   for(ee=0;ee<NUM_ELEMENTS;ee++)
//			   {
//				   sfh_Elements[ii][ee]+=fractionRings[jj]*Gal[p].sfh_DiskMass_elementsRings[jj][ii][ee];
//				   sfh_ElementsRings[jj][ii][ee]=fractionRings[jj]*Gal[p].sfh_DiskMass_elementsRings[jj][ii][ee];
//			   }
//#endif
//			 }
//		 }
//		}
//#endif //STAR_FORMATION_HISTORY
//
//#ifdef COMPUTE_SPECPHOT_PROPERTIES
//#ifndef POST_PROCESS_MAGS
//		for(outputbin = 0; outputbin < NOUT; outputbin++)
//		{
//		for(ll = 0; ll < NMAG; ll++)
//		  {
//#ifdef OUTPUT_REST_MAGS
//			Lum[ll][outputbin]=fraction*(Gal[p].Lum[ll][outputbin]-Gal[p].LumBulge[ll][outputbin]);
//			YLum[ll][outputbin]=fraction*(Gal[p].YLum[ll][outputbin]-Gal[p].YLumBulge[ll][outputbin]);
//#endif
//#ifdef COMPUTE_OBS_MAGS
//			ObsLum[ll][outputbin]=fraction*(Gal[p].ObsLum[ll][outputbin]-Gal[p].ObsLumBulge[ll][outputbin]);
//			ObsYLum[ll][outputbin]=fraction*(Gal[p].ObsYLum[ll][outputbin]-Gal[p].ObsYLumBulge[ll][outputbin]);
//#ifdef OUTPUT_MOMAF_INPUTS
//			dObsLum[ll][outputbin]=fraction*(Gal[p].dObsLum[ll][outputbin]-Gal[p].dObsLumBulge[ll][outputbin]);
//			dObsYLum[ll][outputbin]=fraction*(Gal[p].dObsYLum[ll][outputbin]-Gal[p].dObsYLumBulge[ll][outputbin]);
//#endif
//#endif
//		  }
//		}
//#endif //POST_PROCESS_MAGS
//#endif //COMPUTE_SPECPHOT_PROPERTIES
//
//	}
//
//// Add to component
//	if (
//			strcmp (cp,
//			        "ColdGas") == 0) {
//		Gal[p].ColdGas +=
//				Mass;
//		for (
//				mm = 0;
//				mm < NUM_METAL_CHANNELS; mm ++) {
//			Gal[p].MetalsColdGas[mm] += Metals[mm];
//		}
//#ifdef INDIVIDUAL_ELEMENTS
//		for(ee=0;ee<NUM_ELEMENTS;ee++)
//		Gal[p].ColdGas_elements[ee] += Yield[ee];
//#endif
////RINGS
//		for (
//				jj = 0;
//				jj < RNUM;
//				jj ++) {
//			Gal[p].ColdGasRings[jj] += MassRings[jj];
//			for (
//					mm = 0;
//					mm < NUM_METAL_CHANNELS; mm ++) {
//				Gal[p].MetalsColdGasRings[jj][mm] += MetalsRings[jj][mm];
//			}
//#ifdef INDIVIDUAL_ELEMENTS
//			for(ee=0;ee<NUM_ELEMENTS;ee++)
//			Gal[p].ColdGasRings_elements[jj][ee] += YieldRings[jj][ee];
//#endif
//		}
//	} else if (
//			strcmp (cp,
//			        "DiskMass") == 0) {
//		Gal[p].DiskMass +=
//				Mass;
//		for (
//				mm = 0;
//				mm < NUM_METAL_CHANNELS; mm ++) {
//			Gal[p].MetalsDiskMass[mm] += Metals[mm];
//		}
//#ifdef INDIVIDUAL_ELEMENTS
//		for(ee=0;ee<NUM_ELEMENTS;ee++)
//		Gal[p].DiskMass_elements[ee] += Yield[ee];
//#endif
////RINGS
//		for (
//				jj = 0;
//				jj < RNUM;
//				jj ++) {
//			Gal[p].DiskMassRings[jj] += MassRings[jj];
//			for (
//					mm = 0;
//					mm < NUM_METAL_CHANNELS; mm ++) {
//				Gal[p].MetalsDiskMassRings[jj][mm] += MetalsRings[jj][mm];
//			}
//#ifdef INDIVIDUAL_ELEMENTS
//			for(ee=0;ee<NUM_ELEMENTS;ee++)
//			Gal[p].DiskMassRings_elements[jj][ee] += YieldRings[jj][ee];
//#endif
//		}
//
//#ifdef STAR_FORMATION_HISTORY
//		for (ii=0; ii<=Gal[p].sfh_ibin; ii++)
//		 if(sfh_Mass[ii]>0.)
//		   {
//			 Gal[p].sfh_DiskMass[ii] += sfh_Mass[ii];
//			 for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
//				 Gal[p].sfh_MetalsDiskMass[ii][mm] += sfh_Metals[ii][mm];
//#ifdef INDIVIDUAL_ELEMENTS
//			 for(ee=0;ee<NUM_ELEMENTS;ee++)
//				 Gal[p].sfh_DiskMass_elements[ii][ee] += sfh_Elements[ii][ee];
//#endif
//
//			 for(jj=0;jj<RNUM;jj++)
//			 {
//				 Gal[p].sfh_DiskMassRings[jj][ii] += sfh_MassRings[jj][ii];
//				 for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
//					 Gal[p].sfh_MetalsDiskMassRings[jj][ii][mm] += sfh_MetalsRings[jj][ii][mm];
//#ifdef INDIVIDUAL_ELEMENTS
//				 for(ee=0;ee<NUM_ELEMENTS;ee++)
//					 Gal[p].sfh_DiskMass_elementsRings[jj][ii][ee] += sfh_ElementsRings[jj][ii][ee];
//#endif
//			 }
//		   }
//#endif
//
//#ifdef COMPUTE_SPECPHOT_PROPERTIES
//#ifndef POST_PROCESS_MAGS
//		for(outputbin = 0; outputbin < NOUT; outputbin++)
//		{
//		for(ll = 0; ll < NMAG; ll++)
//		  {
//#ifdef OUTPUT_REST_MAGS
//			Gal[p].Lum[ll][outputbin]+=Lum[ll][outputbin];
//			Gal[p].YLum[ll][outputbin]+=YLum[ll][outputbin];
//#endif
//#ifdef COMPUTE_OBS_MAGS
//			Gal[p].ObsLum[ll][outputbin]+=ObsLum[ll][outputbin]; // Mass += fractionRings[jj]/RNUM*Gal[p].HotGas;
//			Gal[p].ObsYLum[ll][outputbin]+=ObsYLum[ll][outputbin];
//#ifdef OUTPUT_MOMAF_INPUTS
//			Gal[p].dObsLum[ll][outputbin]+=dObsLum[ll][outputbin];
//			Gal[p].dObsYLum[ll][outputbin]+=dObsYLum[ll][outputbin];
//#endif
//#endif
//		  }
//		}
//#endif //POST_PROCESS_MAGS
//#endif //COMPUTE_SPECPHOT_PROPERTIES
//	} else if (
//			strcmp (cp,
//			        "BulgeMass") == 0) {
//		Gal[p].BulgeMass +=
//				Mass;
//		for (
//				mm = 0;
//				mm < NUM_METAL_CHANNELS; mm ++) {
//			Gal[p].MetalsBulgeMass[mm] += Metals[mm];
//		}
//#ifdef INDIVIDUAL_ELEMENTS
//		for(ee=0;ee<NUM_ELEMENTS;ee++)
//		Gal[p].BulgeMass_elements[ee] += Yield[ee];
//#endif
//
////RINGS
//		for (
//				jj = 0;
//				jj < RNUM;
//				jj ++) {
//			Gal[p].BulgeMassRings[jj] += MassRings[jj];
//			for (
//					mm = 0;
//					mm < NUM_METAL_CHANNELS; mm ++) {
//				Gal[p].MetalsBulgeMassRings[jj][mm] += MetalsRings[jj][mm];
//			}
//#ifdef INDIVIDUAL_ELEMENTS
//			for(ee=0;ee<NUM_ELEMENTS;ee++)
//			   Gal[p].BulgeMassRings_elements[jj][ee] += YieldRings[jj][ee];
//#endif
//		}
//
//#ifdef STAR_FORMATION_HISTORY
//		for (ii=0; ii<=Gal[p].sfh_ibin; ii++)
//		 if(sfh_Mass[ii]>0.)
//		   {
//			 Gal[p].sfh_BulgeMass[ii] += sfh_Mass[ii];
//			 for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
//				 Gal[p].sfh_MetalsBulgeMass[ii][mm] += sfh_Metals[ii][mm];
//#ifdef INDIVIDUAL_ELEMENTS
//			 for(ee=0;ee<NUM_ELEMENTS;ee++)
//				 Gal[p].sfh_BulgeMass_elements[ii][ee] += sfh_Elements[ii][ee];
//#endif
//
//			 for(jj=0;jj<RNUM;jj++)
//			 {
//				 Gal[p].sfh_BulgeMassRings[jj][ii] += sfh_MassRings[jj][ii];
//				 for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
//					 Gal[p].sfh_MetalsBulgeMassRings[jj][ii][mm] += sfh_MetalsRings[jj][ii][mm];
//#ifdef INDIVIDUAL_ELEMENTS
//				 for(ee=0;ee<NUM_ELEMENTS;ee++)
//					 Gal[p].sfh_BulgeMass_elementsRings[jj][ii][ee] += sfh_ElementsRings[jj][ii][ee];
//#endif
//			 }
//		   }
//#endif
//
//#ifdef COMPUTE_SPECPHOT_PROPERTIES
//#ifndef POST_PROCESS_MAGS
//		for(outputbin = 0; outputbin < NOUT; outputbin++)
//		{
//		for(ll = 0; ll < NMAG; ll++)
//		  {
//#ifdef OUTPUT_REST_MAGS
//			Gal[p].Lum[ll][outputbin]+=Lum[ll][outputbin];
//			Gal[p].YLum[ll][outputbin]+=YLum[ll][outputbin];
//			Gal[p].LumBulge[ll][outputbin]+=Lum[ll][outputbin];
//			Gal[p].YLumBulge[ll][outputbin]+=YLum[ll][outputbin];
//#endif
//#ifdef COMPUTE_OBS_MAGS
//			Gal[p].ObsLum[ll][outputbin]+=ObsLum[ll][outputbin];
//			Gal[p].ObsYLum[ll][outputbin]+=ObsYLum[ll][outputbin];
//			Gal[p].ObsLumBulge[ll][outputbin]+=ObsLum[ll][outputbin];
//			Gal[p].ObsYLumBulge[ll][outputbin]+=ObsYLum[ll][outputbin];
//#ifdef OUTPUT_MOMAF_INPUTS
//			Gal[p].dObsLum[ll][outputbin]+=dObsLum[ll][outputbin];
//			Gal[p].dObsYLum[ll][outputbin]+=dObsYLum[ll][outputbin];
//			Gal[p].dObsLumBulge[ll][outputbin]+=dObsLum[ll][outputbin];
//			Gal[p].dObsYLumBulge[ll][outputbin]+=dObsYLum[ll][outputbin];
//#endif
//#endif
//		  }
//		}
//#endif //POST_PROCESS_MAGS
//#endif //COMPUTE_SPECPHOT_PROPERTIES
//	}
////Subtract from galaxy q;
//	if (
//			strcmp (cq,
//			        "ColdGas") == 0) {
//		Gal[p].ColdGas -=
//				Mass;
//		for (
//				mm = 0;
//				mm < NUM_METAL_CHANNELS; mm ++) {
//			Gal[p].MetalsColdGas[mm] -= Metals[mm];
//		}
//#ifdef INDIVIDUAL_ELEMENTS
//		for(ee=0;ee<NUM_ELEMENTS;ee++)
//		Gal[p].ColdGas_elements[ee] -= Yield[ee];
//#endif
//		for (
//				jj = 0;
//				jj < RNUM;
//				jj ++) {
//			Gal[p].ColdGasRings[jj] -= MassRings[jj];
//			for (
//					mm = 0;
//					mm < NUM_METAL_CHANNELS; mm ++) {
//				Gal[p].MetalsColdGasRings[jj][mm] -= MetalsRings[jj][mm];
//			}
//#ifdef INDIVIDUAL_ELEMENTS
//			for(ee=0;ee<NUM_ELEMENTS;ee++)
//			Gal[p].ColdGasRings_elements[jj][ee] -= YieldRings[jj][ee];
//#endif
//		}
//	} else if (
//			strcmp (cq,
//			        "DiskMass") == 0) {
//		Gal[p].DiskMass -=
//				Mass;
//		for (
//				mm = 0;
//				mm < NUM_METAL_CHANNELS; mm ++) {
//			Gal[p].MetalsDiskMass[mm] -= Metals[mm];
//		}
//#ifdef INDIVIDUAL_ELEMENTS
//		for(ee=0;ee<NUM_ELEMENTS;ee++)
//		Gal[p].DiskMass_elements[ee] -= Yield[ee];
//#endif
////RINGS
//		for (
//				jj = 0;
//				jj < RNUM;
//				jj ++) {
//			Gal[p].DiskMassRings[jj] -= MassRings[jj];
//			for (
//					mm = 0;
//					mm < NUM_METAL_CHANNELS; mm ++) {
//				Gal[p].MetalsDiskMassRings[jj][mm] -= MetalsRings[jj][mm];
//			}
//#ifdef INDIVIDUAL_ELEMENTS
//			for(ee=0;ee<NUM_ELEMENTS;ee++)
//			Gal[p].DiskMassRings_elements[jj][ee] -= YieldRings[jj][ee];
//#endif
//		}
//
//#ifdef STAR_FORMATION_HISTORY
//		for (ii=0; ii<=Gal[p].sfh_ibin; ii++)
//		 if(sfh_Mass[ii]>0.)
//		   {
//			 Gal[p].sfh_DiskMass[ii] -= sfh_Mass[ii];
//			 for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
//				 Gal[p].sfh_MetalsDiskMass[ii][mm] -= sfh_Metals[ii][mm];
//#ifdef INDIVIDUAL_ELEMENTS
//			 for(ee=0;ee<NUM_ELEMENTS;ee++)
//				 Gal[p].sfh_DiskMass_elements[ii][ee] -= sfh_Elements[ii][ee];
//#endif
//
//			 for(jj=0;jj<RNUM;jj++)
//			   {
//				 Gal[p].sfh_DiskMassRings[jj][ii] -= sfh_MassRings[jj][ii];
//				 for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
//					 Gal[p].sfh_MetalsDiskMassRings[jj][ii][mm] -= sfh_MetalsRings[jj][ii][mm];
//
//#ifdef INDIVIDUAL_ELEMENTS
//				 for(ee=0;ee<NUM_ELEMENTS;ee++)
//					 Gal[p].sfh_DiskMass_elementsRings[jj][ii][ee] -= sfh_ElementsRings[jj][ii][ee];
//#endif
//			   }
//		   }
//#endif
//
//#ifdef COMPUTE_SPECPHOT_PROPERTIES
//#ifndef POST_PROCESS_MAGS
//		for(outputbin = 0; outputbin < NOUT; outputbin++)
//		{
//		for(ll = 0; ll < NMAG; ll++)
//		  {
//#ifdef OUTPUT_REST_MAGS
//			Gal[p].Lum[ll][outputbin]-=Lum[ll][outputbin];
//			Gal[p].YLum[ll][outputbin]-=YLum[ll][outputbin];
//#endif
//#ifdef COMPUTE_OBS_MAGS
//			Gal[p].ObsLum[ll][outputbin]-=ObsLum[ll][outputbin];
//			Gal[p].ObsYLum[ll][outputbin]-=ObsYLum[ll][outputbin];
//#ifdef OUTPUT_MOMAF_INPUTS
//			Gal[p].dObsLum[ll][outputbin]-=dObsLum[ll][outputbin];
//			Gal[p].dObsYLum[ll][outputbin]-=dObsYLum[ll][outputbin];
//#endif
//#endif
//		  }
//		}
//#endif //POST_PROCESS_MAGS
//#endif// COMPUTE_SPECPHOT_PROPERTIES
//	} else if (
//			strcmp (cq,
//			        "BulgeMass") == 0) {
//		Gal[p].BulgeMass -=
//				Mass;
//		for (
//				mm = 0;
//				mm < NUM_METAL_CHANNELS; mm ++) {
//			Gal[p].MetalsBulgeMass[mm] -= Metals[mm];
//		}
//#ifdef INDIVIDUAL_ELEMENTS
//		for(ee=0;ee<NUM_ELEMENTS;ee++) // Mass += fractionRings[jj]/RNUM*Gal[p].HotGas;
//		Gal[p].BulgeMass_elements[ee] -= Yield[ee];
//#endif
//
////RINGS
//		for (
//				jj = 0;
//				jj < RNUM;
//				jj ++) {
//			Gal[p].BulgeMassRings[jj] -= MassRings[jj];
//			for (
//					mm = 0;
//					mm < NUM_METAL_CHANNELS; mm ++) {
//				Gal[p].MetalsBulgeMassRings[jj][mm] -= MetalsRings[jj][mm];
//			}
//#ifdef INDIVIDUAL_ELEMENTS
//			for(ee=0;ee<NUM_ELEMENTS;ee++)
//			Gal[p].BulgeMassRings_elements[jj][ee] -= YieldRings[jj][ee];
//#endif
//		}
//
//#ifdef STAR_FORMATION_HISTORY
//		for (ii=0; ii<=Gal[p].sfh_ibin; ii++)
//		 if(sfh_Mass[ii]>0.)
//		   {
//			 Gal[p].sfh_BulgeMass[ii] -= sfh_Mass[ii];
//			 for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
//				 Gal[p].sfh_MetalsBulgeMass[ii][mm] -= sfh_Metals[ii][mm];
//#ifdef INDIVIDUAL_ELEMENTS
//			 for(ee=0;ee<NUM_ELEMENTS;ee++)
//				 Gal[p].sfh_BulgeMass_elements[ii][ee] -= sfh_Elements[ii][ee];
//#endif
//
//			 for(jj=0;jj<RNUM;jj++)
//			   {
//				 Gal[p].sfh_BulgeMassRings[jj][ii] -= sfh_MassRings[jj][ii];
//				 for(mm=0;mm<NUM_METAL_CHANNELS;mm++)
//					 Gal[p].sfh_MetalsBulgeMassRings[jj][ii][mm] -= sfh_MetalsRings[jj][ii][mm];
//
//#ifdef INDIVIDUAL_ELEMENTS
//				 for(ee=0;ee<NUM_ELEMENTS;ee++)
//					 Gal[p].sfh_BulgeMass_elementsRings[jj][ii][ee] -= sfh_ElementsRings[jj][ii][ee];
//#endif
//			   }
//		   }
//#endif
//
//#ifdef COMPUTE_SPECPHOT_PROPERTIES
//#ifndef POST_PROCESS_MAGS
//		for(outputbin = 0; outputbin < NOUT; outputbin++)
//		{
//		for(ll = 0; ll < NMAG; ll++)
//		  {
//#ifdef OUTPUT_REST_MAGS
//			   Gal[p].Lum[ll][outputbin]-=Lum[ll][outputbin];
//			   Gal[p].YLum[ll][outputbin]-=YLum[ll][outputbin];
//			   Gal[p].LumBulge[ll][outputbin]-=Lum[ll][outputbin];
//			   Gal[p].YLumBulge[ll][outputbin]-=YLum[ll][outputbin];
//#endif
//#ifdef COMPUTE_OBS_MAGS
//			   Gal[p].ObsLum[ll][outputbin]-=ObsLum[ll][outputbin];
//			   Gal[p].ObsYLum[ll][outputbin]-=ObsYLum[ll][outputbin];
//			   Gal[p].ObsLumBulge[ll][outputbin]-=ObsLum[ll][outputbin];
//			   Gal[p].ObsYLumBulge[ll][outputbin]-=ObsYLum[ll][outputbin];
//#ifdef OUTPUT_MOMAF_INPUTS
//			   Gal[p].dObsLum[ll][outputbin]-=dObsLum[ll][outputbin];
//			   Gal[p].dObsYLum[ll][outputbin]-=dObsYLum[ll][outputbin];
//			   Gal[p].dObsLumBulge[ll][outputbin]-=dObsLum[ll][outputbin];
//			   Gal[p].dObsYLumBulge[ll][outputbin]-=dObsYLum[ll][outputbin];
//#endif
//#endif
//		  }
//		}
//#endif //POST_PROCESS_MAGS
//#endif//
//
//	}
//}