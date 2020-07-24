/* sample.c
 * written by Sanghyuk Moon
 * read density profile from scf.py and sample points from it
 */

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_qrng.h>
#include <gsl/gsl_rng.h>

void* malloc(size_t);
void read_data();
void sample_particles();

static double *rho_ring;
static double *pos_r;
static double *pos_mu;

gsl_rng *RNG; /* global generator */

static long Npart;

static double Rmax;
static double zmax;
static double rho_max;
// To achieve smooth resolution gradient, we need to define transition region.
// eta_0 is a thickness of the ring, which is defined in icgen.py.
// Transition region is from eta = eta_0 to eta = eta_plus.
// mh is a desired mass of a halo gas particle, mr is of a ring mass particle
// alpha is a power law index for the mass transition
static double eta_0, eta_plus, alpha;
static int N;

double max(double a[], long N) {
  double res = 0.0;
  for (int i=0; i<N; ++i) {
    if (a[i] > res) {
      res = a[i];
    }   
  }
  return res;
}

void read_data() {
  int i,j;
  FILE* fp_ring;
  FILE* fp_prms;

  pos_r = malloc(sizeof *pos_r * N);
  pos_mu = malloc(sizeof *pos_mu * N);
  rho_ring = malloc(sizeof *rho_ring * N * N);

  fp_ring = fopen("rho_ring.out","r");
  fp_prms = fopen("icgen_to_sample","r");
  for (i=0; i<N; i++) {
    for (j=0; j<N; j++) {
      fscanf(fp_ring, "%le %le %le", &pos_r[i], &pos_mu[j], &rho_ring[j*N + i]);
    }
  }
  fscanf(fp_prms, "%le %le %le", &eta_0, &eta_plus, &alpha);
  fclose(fp_ring);
  fclose(fp_prms);
}

void sample_particles() {
  int dim=2;
  long id=0;
  double pos[dim], u4;
  double theta;
  double r, mu;
  double eta, probability;
  FILE *fp;
  /* GSL quasi random number & interpolation */
  gsl_qrng *q = gsl_qrng_alloc(gsl_qrng_halton, 2);
  gsl_interp_accel *xacc = gsl_interp_accel_alloc();
  gsl_interp_accel *yacc = gsl_interp_accel_alloc();
  const gsl_interp2d_type *T = gsl_interp2d_bilinear;
  gsl_interp2d *interp = gsl_interp2d_alloc(T, N, N);
  gsl_interp2d_init(interp, pos_r, pos_mu, rho_ring, N, N);
  printf("Npart = %ld\n", Npart);

  fp = fopen("ring.sample","w");
  while (id < Npart) {
    gsl_qrng_get(q, pos);  
    pos[0] *= Rmax;
    pos[1] = zmax*(2.0*pos[1] - 1.0);
    eta = sqrt((pos[0] - 1.0)*(pos[0] - 1.0) + pos[1]*pos[1]);
    r = sqrt(pos[0]*pos[0] + pos[1]*pos[1]);
    mu = pos[1] / r;
    u4 = rho_max * gsl_rng_uniform(RNG);
    if (fabs(mu) > max(pos_mu, N)) continue; // sample inside data box

    if ((0 < eta)&&(eta < eta_0)) {
      probability = pos[0]*gsl_interp2d_eval(interp, pos_r, pos_mu, rho_ring, r, fabs(mu), xacc, yacc);
    }
    else if ((eta_0 < eta)&&(eta < eta_plus)) {
      probability = pos[0]*gsl_interp2d_eval(interp, pos_r, pos_mu, rho_ring, r, fabs(mu), xacc, yacc) / pow(eta/eta_0, alpha);
    }
    else if (eta > eta_plus) {
      probability = pos[0]*gsl_interp2d_eval(interp, pos_r, pos_mu, rho_ring, r, fabs(mu), xacc, yacc) / pow(eta_plus/eta_0, alpha);
    }
    else {
      printf("strange eta!\n");
    }

    if (u4 <= probability)
    {
      theta = 2.0*M_PI*(double)id/(double)Npart;
      fprintf(fp, "%le %le %le\n",
              pos[0]*cos(theta), pos[0]*sin(theta), pos[1]);
      id++;
    }

//    if (u4 <= pos[0]*gsl_interp2d_eval(interp, pos_r, pos_mu, rho_ring, r,
//      fabs(mu), xacc, yacc))
//    {
//      theta = 2.0*M_PI*(double)id/(double)Npart;
//      fprintf(fp, "%le %le %le\n",
//              pos[0]*cos(theta), pos[0]*sin(theta), pos[1]);
//      id++;
//    }
  }
  fclose(fp);
  gsl_qrng_free(q);
  gsl_interp2d_free(interp);
}

int main() {
  Npart = pow(2,21);
  Rmax = 2.0;   // should be consistent with icgen.py (halo mass calc.)
  zmax = 0.5;   // should be consistent with icgen.py (halo mass calc.)
  rho_max = 1.2;
  N = 1025;

  /* GSL random number generator setup */
  const gsl_rng_type *T;
  T = gsl_rng_ranlxd2;
  RNG = gsl_rng_alloc(T);
  printf("generator type: %s\n", gsl_rng_name(RNG));
  printf("seed = %lu\n", gsl_rng_default_seed);
  printf("first value = %le\n", gsl_rng_uniform(RNG));
  read_data();
  sample_particles();
  gsl_rng_free(RNG);
  return 0;
}
