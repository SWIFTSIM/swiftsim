/*
! FORTRAN MODULE FOR CALCULATING THE CLUSTER FORMATION EFFICIENCY (LOCAL
VERSION) ! Copyright (C) 2012  Diederik Kruijssen
!
! NAME:
!       CFElocalmod (module)
!       F_CFElocal (main function)
!
! PURPOSE:
!       Calculate the fraction of star formation occurring in bound clusters,
!       i.e. the cluster formation efficiency or CFE
!
! TERMS OF USE:
!       If you use this routine while preparing a paper, please cite the paper
!       in which this model was presented:
!       Kruijssen, J. M. D., 2012, MNRAS 426, 3008
!
! CALLING SEQUENCE:
!       f_cfelocal(rholoc,sigmaloc,csloc)
!
! INPUT PARAMETERS:
!       NOTE: ALL INPUT SHOULD BE IN SI UNITS!
!       rholoc - local gas volume density
!       sigmaloc - local 1D gas velocity dispersion
!       csloc - local gas sound speed
!
! OPTIONAL INPUT IS SET IN THE FILE "PARAMETERS"
!       sflaw - star formation law: 0=Elmegreen (2002, default)
!                                   1=Krumholz & McKee (2005)
!       qvir - giant molecular cloud virial parameter (default=1.3)
!       tsn - time of the first supernova (default=3.*Myr)
!       tview - time at which CFE is determined (default=10.*Myr)
!       surfGMC - giant molecular cloud surface density
(default=100.*Msun/pc**2.) !       ecore - maximum (protostellar core) star
formation efficiency (default=.5) !       beta0 - turbulent-to-magnetic pressure
ratio (default=1.e10) !       radfb - feedback: 0=supernova, 1=radiative, 2=both
(default=0)
!
! OUTPUT:
!       cfe - fraction of star formation occurring in bound clusters
!
! EXAMPLE IS GIVEN IN THE FILE testCFE.F90:
!       Calculate the CFE in the solar neighbourhood
!
!       PROGRAM testCFElocal
!           use CFElocalmod
!           REAL pc,msun,myr,rhoMW,sigmaMW,csMW,cfe(1:4)
!
!           pc=3.1e16 !parsec in meters
!           Msun=2.e30 !solar mass in kg
!           Myr=3.16e13 !million years in seconds
!           rhoMW=0.03*msun/pc^3. !volume density
!           sigmaMW=7.e3 !velocity dispersion
!           csMW=0.2e3 !sound speed
!           cfe=f_cfelocal(rhoMW,sigmaMW,csMW) !CFE
!           PRINT*,cfe(1)*100.
!       END PROGRAM
!
!       > g95 f_cfelocal.f90 testCFElocal.f90 -o testCFElocal
!       > ./testCFElocal
!          8.827569
!
! MODULE STRUCTURE:
!       Auxiliary functions come first, the actual CFE function is located at
!       the bottom of this file
!
! REVISION HISTORY:
!       Written by Diederik Kruijssen, August 2012
!
!       Translated to C:
!       JLP (2015) -- j.l.pfeffer@ljmu.ac.uk
!
*/

#ifndef SWIFT_MOSAICS_CFELOCAL_H
#define SWIFT_MOSAICS_CFELOCAL_H

#include <math.h>

/* free-fall time as a function of density */
INLINE static double f_tff(double rho) {
  const double G = 6.67e-11;               // gravitational constant
  return sqrt(3. * M_PI / 32. / G / rho);  // free-fall time
}

/* Mach number as a function of surface density, angular velocity and Toomre Q,
 * from Krumholz & McKee (2005) */
INLINE static double f_mach(double sigmaloc, double csloc) {
  return sigmaloc / csloc;  // Mach number
}

/* dispersion of overdensity PDF as a function of  Mach number and magnetic
 * pressure ratio, based on e.g. Padoan & Nordlund (2011) */
INLINE static double f_sigrho(double mach, double beta0) {
  const double b = 0.5;  // constant
  return sqrt(
      log(1. + 3. * b * b * mach * mach * beta0 / (beta0 + 1.)));  // dispersion
}

/* critical overdensity for SF in the KM05 sSFR_ff as a function of Mach number
 * and GMC virial ratio, from Krumholz & McKee (2005) */
INLINE static double f_xcrit(double qvir, double mach) {
  const double phix = 1.12;  // constant
  return M_PI * M_PI * phix * phix / 15. * qvir * mach *
         mach;  // critical overdensity
}

/* specific star formation rate per free-fall time (sSFR_ff) for Elmegreen
 * (2002, sflaw=0) or Krumholz & McKee (2005, sflaw=1) */
INLINE static double f_sfrff(double qvir, double mach, double beta0,
                             double ecore, int sflaw) {
  double xcrit, sigrho, f;
  const double phit = 1.91;        // constant
  xcrit = f_xcrit(qvir, mach);     // critical overdensity for star formation in
                                   // the KM05 model, see above
  sigrho = f_sigrho(mach, beta0);  // overdensity PDF dispersion, see above
  if (sflaw == 1)
    f = .5 * ecore / phit *
        (1. +
         erf((-2. * log(xcrit) + sigrho * sigrho) /
             (2. * M_SQRT2 *
              sigrho)));  // specific star formation rate per free-fall time
  else
    f = 0.012;  // specific star formation rate per free-fall time
  return f;
}

/* overdensity PDF of the ISM as a function of overdensity x and its logarithmic
 * mean and dispersion */
INLINE static double f_dpdx(double x, double mulnx, double sig) {
  double dif = log(x) - mulnx;  // dummy
  return 1. / (sqrt(2. * M_PI * sig * sig) * x) *
         exp(-.5 * dif * dif / (sig * sig));  // overdensity PDF
}

/* rough relation between gas surface density and volume density and velocity
 * dispersion, assuming an equilibrium disc */
INLINE static double f_surfg(double rholoc, double sigmaloc) {
  const double G = 6.67e-11;  // gravitational constant
  const double phiP = 3.;     // constant
  return sqrt(2. * rholoc * sigmaloc * sigmaloc /
              (M_PI * G * phiP));  // gas surface density
}

/* naturally bound fraction of star formation */
INLINE static double f_fstar(double rholoc, double sigmaloc, double csloc,
                             double x, double ecore, double beta0, double qvir,
                             double tsn, double tview, double surfGMC,
                             int sflaw, int radfb) {
  double surfg, surffb, surffb3, mach, sfrff, rhog, tff, efb, einc, efbrad,
      fstar;

  const double G = 6.67e-11;     // gravitational constant
  const double sigSB = 5.67e-8;  // Stefan-Boltzmann constant
  const double c = 299792458.;   // speed of light
  const double phifb = 1.6e-5;   // feedback efficiency
  const double kappa0 = 2.4e-5;  // opacity constant
  const double psi = .3;         // light-to-mass ratio
  const double phitrap = .2;     // trapping ratio

  surfg = f_surfg(rholoc, sigmaloc);  // estimate of gas surface density
  surffb =
      fmax(surfGMC, surfg);  // surface density on which radiative feedback acts
  surffb3 = surffb * surffb * surffb;
  mach = f_mach(sigmaloc, csloc);  // Mach number, see above
  sfrff = f_sfrff(
      qvir, mach, beta0, ecore,
      sflaw);  // specific star formation rate per free-fall time, see above
  rhog = x * rholoc;  // gas volume density
  tff = f_tff(rhog);  // free-fall time, see above

  // Star Formation Efficiencies
  if (radfb == 0. || radfb == 2)
    efb =
        0.5 * sfrff * tsn / tff *
        (1. + sqrt(1. + 4. * tff * sigmaloc * sigmaloc /
                            (phifb * sfrff * tsn * tsn * x)));  //  SN feedback
  else
    efb = 1.;
  einc = sfrff * tview / tff;  // star formation is incomplete/still ongoing
  if (radfb > 0.)
    efbrad = 2. * sigSB / (phitrap * kappa0 * kappa0 * psi * surffb3) *
             (sqrt(1. + 2. * M_PI * c * G * phitrap * kappa0 * kappa0 * surffb *
                            surffb3 / (1. * sigSB)) -
              1.);  // radiative feedback
  else
    efbrad = 1.;
  // epsilons=[ecore,efb,efbrad,einc]; //SFEs for
  // [maximum,SNfeedback,radiativefeedback, incomplete]
  // fstar=MIN(epsilons(1),epsilons(2),epsilons(3),epsilons(4)); //local SFE is
  // the minimum of those
  fstar = fmin(fmin(ecore, efb), fmin(efbrad, einc));
  return fstar;
}

/* obtain fractions from integrating the overdensity PDFs */
INLINE static double f_integrate(double xsurv, double mulnx, double sig,
                                 double rholoc, double sigmaloc, double csloc,
                                 double ecore, double beta0, double qvir,
                                 double tsn, double tview, double surfGMC,
                                 int cce, int sflaw, int radfb) {
  double xmin1, xmin2, xmax, xarr[1000], f1, dx, xg, fstar, bound, integral,
      dpdx, f2, frac;
  int nx, ix;

  nx = 1000;  // number of integration steps (checked to be sufficient for
              // convergence)
  xmin1 = exp(mulnx - 5. * sig);  // minimum overdensity
  xmax = exp(mulnx + 10. * sig);  // maximum overdensity

  if (cce > 0 && xsurv < xmin1)
    xmin1 = xsurv;  // if calculating the cruel cradle effect and critical
                    // overdensity below minimum, then adjust

  if (cce > 0 && xsurv > xmax)
    xmax = xsurv;  // if calculating the cruel cradle effect and critical
                   // overdensity below maximum, then adjust

  for (ix = 0; ix < nx; ix++)
    xarr[ix] = xmin1 * pow(xmax / xmin1, (ix - 0.5) / nx);  // integration array

  f1 = 0.;                       // initialize integral
  for (ix = 0; ix < nx; ix++) {  // denominator integral
    dx = xarr[ix] * (pow(xmax / xmin1, 1. / (2. * nx)) -
                     pow(xmax / xmin1, -1. / (2. * nx)));  // step size
    xg = xarr[ix];                                         // overdensity
    fstar = f_fstar(rholoc, sigmaloc, csloc, xg, ecore, beta0, qvir, tsn, tview,
                    surfGMC, sflaw, radfb);  // local SFE
    bound = fstar / ecore;                   // local bound fraction
    if (cce == 0)
      bound = 1.;  // if not calculating the cruel cradle effect but the
                   // naturally bound fraction of SF, the denominator should
                   // contain all SF
    if (cce == 2)
      bound = 1.;  // if calculating the cruel cradle effect with respect to all
                   // SF, the denominator should contain all SF
    integral = bound * fstar * xarr[ix];  // integral part 1
    dpdx = f_dpdx(xg, mulnx, sig);   // overdensity PDF, i.e. integral part 2
    f1 = f1 + integral * dpdx * dx;  // integral
  }

  if (cce > 0)  // if calculating the cruel cradle effect set minimum
                // overdensity to critical overdensity
    xmin2 = xsurv;
  else
    xmin2 = xmin1;

  for (ix = 0; ix < nx; ix++)
    xarr[ix] = xmin2 * pow(xmax / xmin2, (ix - 0.5) / nx);  // integration array

  f2 = 0.;                       // initialize integral
  for (ix = 0; ix < nx; ix++) {  // numerator integral
    dx = xarr[ix] * (pow(xmax / xmin2, 1. / (2. * nx)) -
                     pow(xmax / xmin2, -1. / (2. * nx)));  // step size
    xg = xarr[ix];                                         // overdensity
    fstar = f_fstar(rholoc, sigmaloc, csloc, xg, ecore, beta0, qvir, tsn, tview,
                    surfGMC, sflaw, radfb);  // local SFE
    bound = fstar / ecore;                   // local bound fraction
    if (cce == 2)
      bound = 1.;  // if calculating the cruel cradle effect with respect to all
                   // SF, the numerator should contain all SF
    integral = bound * fstar * xarr[ix];  // integral part 1
    dpdx = f_dpdx(xg, mulnx, sig);   // overdensity PDF, i.e. integral part 2
    f2 = f2 + integral * dpdx * dx;  // integral
  }

  if (f1 == 0)
    frac = 0.;
  else
    frac = f2 / f1;  // numerator divided by denominator
  return frac;
}

/* ratio of encounter timescale to energy dissipation timescale as a function of
 * cloud virial ratio and overdensity */
INLINE static double f_phit(double qvir, double x) {
  return 3.1 * sqrt((qvir / 1.3) * (x / 1.e4));  // ratio of encounter timescale
                                                 // to energy dissipation
                                                 // timescale
}

/* adiabatic correction as a function of cloud virial ratio and overdensity */
INLINE static double f_phiad(double qvir, double x) {
  double phit;
  phit = f_phit(
      qvir, x);  // ratio of encounter timescale to energy dissipation timescale
  return exp(-2. * phit);  // adiabatic correction
}

/* critical overdensity to remain bound despite the cruel cradle effect (also
 * see Kruijssen et al. 2011) */
INLINE static double f_xcce(double sigmaloc, double surfGMC, double qvir,
                            double tview) {
  double xmin, xmax, xfit, accuracy, xfit0, xarr[101], xarr2[101], phiad, diff,
      diffx;
  int nx, niter, itermax, ix, ixfit;

  const double G = 6.67e-11;   // gravitational constant
  const double g_close = 1.5;  // close encounter correction
  const double phish = 2.8;    // higher-order energy loss correction
  const double f =
      0.7;  // fraction of injected energy that is used for unbinding the region
  // const double eta = 2.*1.305*3.*M_PI/64.; //for Plummer
  // const double rh2r2av = .25; //for Plummer

  // SOLVE implicit relation for x_cce
  xmin = 1.e-4;               // minimum x
  xmax = 1.e8;                // maximum x
  nx = 101;                   // length of x array
  niter = 0;                  // number of elapsed iterations
  itermax = 10;               // maximum number of iterations
  xfit = xmin * xmin / xmax;  // initialisation of fitted x
  accuracy = 1.e-6;           // desired logarithmic accuracy
  xfit0 = xfit * accuracy;    // initialisation of previously fitted x

  while (fabs(log10(xfit / xfit0)) > accuracy &&
         niter < itermax) {        // while iteration     does not give desired
                                   // convergence, do
    xfit0 = xfit;                  // previously fitted x
    for (ix = 0; ix < nx; ix++) {  // for all x
      xarr[ix] = pow(10., (ix - 1) / (nx - 1.) * log10(xmax / xmin) +
                              log10(xmin));  // x array
      phiad = f_phiad(qvir, xarr[ix]);       // adiabatic correction, see above
      xarr2[ix] = 87.5 * sqrt(M_PI) * f * g_close * G * phish * surfGMC *
                  phiad * tview / sigmaloc;  // right-hand side of equation
    }
    diff = 1.e30;
    for (ix = 0; ix < nx; ix++) {  // for all x
      diffx = fabs(xarr[ix] - xarr2[ix]);
      if (diffx < diff && ix > 1) {
        diff = diffx;
        ixfit = ix;  // index where x equals right-hand side of equation
      }
    }
    xfit = xarr[ixfit];      // solution for x_cce
    xmin = xarr[ixfit - 1];  // new minimum x
    xmax = xarr[ixfit + 1];  // new maximum x
    niter++;                 // increase number of elapsed iterations by 1
  }

  /** TODO hide this within a DEBUG def
      if (niter == itermax){ //if we stopped due to reaching maximum number of
     iterations, then stop printf("f_xcce: no convergence, increase nx or
     itermax in the f_xcce subroutine"); endrun(932);
      }
  */
  return xfit;
}

/* FUNCTION TO CALCULATE THE CLUSTER FORMATION EFFICIENCY */
INLINE static double f_cfelocal(double rholoc, double sigmaloc, double csloc,
                                const struct stars_props* stars_properties) {
  double surfg, mach, sigrho, mulnx, fbound;

  // SET CONSTANTS
  const double pc = 3.086e16;                 // parsec in meters
  const double Msun = 1.989e30;               // solar mass in kg
  const double Myr = 1.e6 * 86400. * 365.25;  // million years in seconds

  const int sflaw = stars_properties->sflaw;
  const double qvir = stars_properties->qvir;
  const double tsn = stars_properties->tsn * Myr;
  const double tview = stars_properties->tview * Myr;
  const double ecore = stars_properties->ecore;
  const double beta0 = stars_properties->beta0;
  const int radfb = stars_properties->radfb;

  double surfGMC = stars_properties->surfGMC * Msun / (pc * pc);
  surfg = f_surfg(rholoc, sigmaloc);  // estimate of gas surface density
  if (surfg > surfGMC) surfGMC = surfg;

  // CALCULATE DERIVED PARAMETERS
  mach = f_mach(sigmaloc, csloc);  // Mach number
  sigrho = f_sigrho(mach, beta0);  // dispersion of overdensity PDF
  mulnx = -.5 * sigrho * sigrho;   // logarithmic mean of overdensity PDF

  // CALCULATE F_BOUND
  fbound = f_integrate(0., mulnx, sigrho, rholoc, sigmaloc, csloc, ecore, beta0,
                       qvir, tsn, tview, surfGMC, 0, sflaw,
                       radfb);  // naturally bound part of star formation

  /** NOTE We never need CCE in this version
      //CALCULATE F_CCE
  #ifdef MOSAICS_CCE
      double xsurv, fcce;
      xsurv=f_xcce(sigmaloc,surfGMC,qvir,tview); //critical overdensity to
  remain bound despite the cruel cradle effect
      fcce=f_integrate(xsurv,mulnx,sigrho,rholoc,sigmaloc,csloc,ecore,beta0,qvir,tsn,tview,surfGMC,1,sflaw,radfb);
  //part of bound SF surviving the cruel cradle effect #else const double fcce
  = 1.; #endif

      //CALCULATE CFE
      return fbound*fcce;
  */

  // The CFE
  return fbound;
}

/* local supernovae feedback timescale */
INLINE static double feedback_timescale(
    double rholoc, double sigmaloc, double csloc,
    const struct stars_props* stars_properties) {
  double mach, sfrff, rhog, tff;

  // SET CONSTANTS
  const double Myr = 1.e6 * 86400. * 365.25;  // million years in seconds

  const int sflaw = stars_properties->sflaw;
  const double qvir = stars_properties->qvir;
  const double tsn = stars_properties->tsn * Myr;
  const double ecore = stars_properties->ecore;
  const double beta0 = stars_properties->beta0;
  const double phifb = 1.6e-5;  // feedback efficiency

  const double x = 1.;

  mach = f_mach(sigmaloc, csloc);  // Mach number, see above
  sfrff = f_sfrff(
      qvir, mach, beta0, ecore,
      sflaw);  // specific star formation rate per free-fall time, see above
  rhog = x * rholoc;  // gas volume density
  tff = f_tff(rhog);  // free-fall time, see above

  // SN feedback
  return 0.5 * tsn *
         (1. + sqrt(1. + 4. * tff * sigmaloc * sigmaloc /
                             (phifb * sfrff * tsn * tsn * x)));
}

#endif /* SWIFT_MOSAICS_CFELOCAL_H */
