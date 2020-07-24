/* sample.c
 * written by Sanghyuk Moon
 * 
 * do rejection sampling using the halton sequence
 * density is assumed to be normalized; rho ~ [0,1]
 */

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_qrng.h>
#include <gsl/gsl_rng.h>

const unsigned int dim=2;

double density(double R, double z) {
  double rho, scaleh, Rmin;
  Rmin = 0.1;
  scaleh = 0.1; // Gaussian scale height of the disk
  rho = Rmin/R*exp(-z*z/scaleh/scaleh/2);
  return rho;
}

void sample_particles(long int ngas, double x1min, double x1max, double x3min, double x3max) {
  long id=0;
  double pos[dim], u4;
  double x1,x2,x3; // particle position in cylindrical coordinates
  double den; // desired gas density
  gsl_qrng *qrng;
  gsl_rng *rng;

  // GSL random number generator setup
  qrng = gsl_qrng_alloc(gsl_qrng_halton, dim);
  rng = gsl_rng_alloc(gsl_rng_ranlxd2);

  while (id < ngas) {
    gsl_qrng_get(qrng, pos);
    x1 = x1min + (x1max-x1min)*pos[0];
    x3 = x3min + (x3max-x3min)*pos[1];
    u4 = x1max*gsl_rng_uniform(rng); // uniform random number in [0, x1max]
    den = density(x1, x3);

    // accept if u4 < R*density.
    // additional weighting by R compensates the R*d\phi volume factor.
    if (u4 <= x1*den) {
      x2 = 2.0*M_PI*(double)id/(double)ngas;
      printf("%le %le %le\n", x1*cos(x2), x1*sin(x2), x3);
      id++;
    }
  }
  gsl_qrng_free(qrng);
  gsl_rng_free(rng);
}

int main() {
  long int ngas; // number of particles
  double x1min, x1max, x3min, x3max;

  // set domain info
  ngas = pow(10,6);
  x1min = 0.1;
  x1max = 1.0;
  x3min = -0.5;
  x3max = 0.5;

  // run rejection samplig
  sample_particles(ngas, x1min, x1max, x3min, x3max);
  return 0;
}
