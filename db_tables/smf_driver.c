// smf_driver.c -- machine-readable output of glfunc (edgarcosta/lfunctions)
//
// Reads the same one-line spec as examples/rational.c,
//   label:degree:conductor:weight:mus:euler_factors
// (weight = motivic weight w, normalisation = w/2; mus = Gamma_R shifts in the
// algebraic normalisation; euler_factors = ascending integer coefficients of the
// local polynomial at the primes 2, 3, 5, ... in order, each of length degree+1),
// and prints one JSON object per input line with everything an LMFDB
// lfunc_lfunctions / lfunc_plot row needs.
//
// Usage: smf_driver.exe input.txt [output.json] [--self-dual] [--plot-points N]
//                                 [--target-prec BITS]
//   --self-dual    tell glfunc the L-function is self-dual (Lparams.self_dual=YES);
//                  otherwise DK, as in rational.c.
//   --plot-points  number of Z(t) samples requested on [0, max_t] (default 2048;
//                  glfunc caps max_t at 512/(8*degree) and the density at its grid).
//
// If fewer Euler factors are given than glfunc needs (Lfunc_nmax), the driver
// lowers nmax with Lfunc_reduce_nmax, as ~/lfun/rational_r.c does, and reports
// "certified": false.  Everything else in the output is then not rigorous.
//
// Balls are printed as {"mid": "<40 significant digits>", "rad": <double>}.
// Each zero also carries "x2e100" = round(mid * 2^100), the CMF Lhash convention.

#include <flint/acb_poly.h>
#include <flint/arb.h>
#include <flint/arf.h>
#include <flint/fmpz.h>
#include <ctype.h>
#include <inttypes.h>
#include <primesieve.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "glfunc.h"

#define MID_DIGITS 40

static const struct { uint64_t bit; const char *name; } ERRNAMES[] = {
  {ERR_NO_DATA, "ERR_NO_DATA"}, {ERR_ZERO_ERROR, "ERR_ZERO_ERROR"},
  {ERR_OOM, "ERR_OOM"}, {ERR_UPSAMPLE, "ERR_UPSAMPLE"}, {ERR_MU_HALF, "ERR_MU_HALF"},
  {ERR_M_ERROR, "ERR_M_ERROR"}, {ERR_STAT_POINT, "ERR_STAT_POINT"},
  {ERR_SPEC_VALUE, "ERR_SPEC_VALUE"}, {ERR_G_INFILE, "ERR_G_INFILE"},
  {ERR_BAD_DEGREE, "ERR_BAD_DEGREE"}, {ERR_SPEC_NZ, "ERR_SPEC_NZ"},
  {ERR_G_EXTENT, "ERR_G_EXTENT"},
  {ERR_SOME_DATA, "ERR_SOME_DATA"}, {ERR_ZERO_PREC, "ERR_ZERO_PREC"},
  {ERR_RH_ERROR, "ERR_RH_ERROR"}, {ERR_INSUFF_EULER, "ERR_INSUFF_EULER"},
  {ERR_NO_RANK, "ERR_NO_RANK"}, {ERR_CONFLICT_RANK, "ERR_CONFLICT_RANK"},
  {ERR_DBL_ZERO, "ERR_DBL_ZERO"}, {ERR_SPEC_PREC, "ERR_SPEC_PREC"},
  {ERR_G_OUTFILE, "ERR_G_OUTFILE"},
};

static void print_errors_json(FILE *f, Lerror_t e) {
  int first = 1;
  fprintf(f, "[");
  for (size_t i = 0; i < sizeof(ERRNAMES) / sizeof(ERRNAMES[0]); i++)
    if (e & ERRNAMES[i].bit) {
      fprintf(f, "%s\"%s\"", first ? "" : ", ", ERRNAMES[i].name);
      first = 0;
    }
  fprintf(f, "]");
}

static void print_arb_json(FILE *f, const arb_t x) {
  char *s = arf_get_str(arb_midref(x), MID_DIGITS);
  fprintf(f, "{\"mid\": \"%s\", \"rad\": %.6e}", s, mag_get_d(arb_radref(x)));
  flint_free(s);
}

static void print_zero_json(FILE *f, const arb_t z) {
  arf_t t;
  fmpz_t n;
  arf_init(t);
  fmpz_init(n);
  arf_mul_2exp_si(t, arb_midref(z), 100);
  arf_get_fmpz(n, t, ARF_RND_NEAR);
  char *s = arf_get_str(arb_midref(z), MID_DIGITS);
  char *h = fmpz_get_str(NULL, 10, n);
  fprintf(f, "{\"mid\": \"%s\", \"rad\": %.6e, \"x2e100\": \"%s\"}", s,
          mag_get_d(arb_radref(z)), h);
  flint_free(s);
  flint_free(h);
  fmpz_clear(n);
  arf_clear(t);
}

// parse "[a,b,c]" of doubles; returns count or -1
static int parse_doubles(double *v, int n, const char *s) {
  while (isspace(*s)) s++;
  if (*s != '[') return -1;
  s++;
  int i = 0;
  while (i < n) {
    char *e;
    v[i++] = strtod(s, &e);
    if (e == s) return -1;
    s = e;
    while (isspace(*s)) s++;
    if (*s == ']') break;
    if (*s != ',') return -1;
    s++;
  }
  return i;
}

// parse "[[..],[..],..]" of int64 with rows of length d; returns rows or -1
static int64_t parse_factors(int64_t **out, int d, const char *s) {
  size_t cap = 1024, n = 0;
  int64_t *v = malloc(cap * sizeof(int64_t));
  while (*s && *s != '[') s++;
  if (!*s) return -1;
  s++;
  for (;;) {
    while (*s && *s != '[' && *s != ']') s++;
    if (*s != '[') break; // end of outer list
    s++;
    for (int j = 0; j < d; j++) {
      char *e;
      if (n == cap) { cap *= 2; v = realloc(v, cap * sizeof(int64_t)); }
      v[n++] = strtoll(s, &e, 10);
      if (e == s) { free(v); return -1; }
      s = e;
      while (isspace(*s)) s++;
      if (j < d - 1) {
        if (*s != ',') { free(v); return -1; }
        s++;
      }
    }
    while (isspace(*s)) s++;
    if (*s != ']') { free(v); return -1; } // row of wrong length
    s++;
  }
  *out = v;
  return (int64_t)(n / d);
}

static int run_line(FILE *out, char *line, int self_dual, uint64_t plot_points,
                    int64_t target_prec) {
  char *tok[6];
  int nt = 0;
  char *p = line;
  tok[nt++] = p;
  while (nt < 6 && (p = strchr(p, ':')) != NULL) {
    *p++ = 0;
    tok[nt++] = p;
  }
  if (nt != 6) {
    fprintf(stderr, "bad line (%d fields)\n", nt);
    return -1;
  }
  const char *label = tok[0];
  int degree = atoi(tok[1]);
  uint64_t conductor = strtoull(tok[2], NULL, 10);
  int weight = atoi(tok[3]);
  double mus[MAX_DEGREE];
  if (degree < 1 || degree > MAX_DEGREE || parse_doubles(mus, degree, tok[4]) != degree) {
    fprintf(stderr, "%s: bad degree or mus\n", label);
    return -1;
  }
  int64_t *ef = NULL;
  int64_t nef = parse_factors(&ef, degree + 1, tok[5]);
  if (nef < 0) {
    fprintf(stderr, "%s: bad euler factors\n", label);
    return -1;
  }

  Lerror_t ecode = ERR_SUCCESS;
  Lparams_t Lp;
  Lp.degree = degree;
  Lp.conductor = conductor;
  Lp.normalisation = 0.5 * weight;
  Lp.mus = mus;
  Lp.target_prec = target_prec;
  Lp.wprec = 0;
  Lp.gprec = 0;
  Lp.rank = DK;
  Lp.self_dual = self_dual ? YES : DK;
  Lp.cache_dir = ".";
  Lfunc_t L = Lfunc_init_advanced(&Lp, &ecode);
  if (L == NULL || fatal_error(ecode)) {
    fprintf(stderr, "%s: Lfunc_init failed\n", label);
    fprint_errors(stderr, ecode);
    free(ef);
    return -1;
  }

  uint64_t nmax_required = Lfunc_nmax(L);
  size_t nprimes;
  uint64_t *primes = (uint64_t *)primesieve_generate_primes(0, nmax_required, &nprimes, UINT64_PRIMES);
  int certified = 1;
  uint64_t nmax_used = nmax_required;
  size_t nuse = nprimes;
  if ((size_t)nef < nprimes) {
    // not enough data: lower nmax to just below the first missing prime
    certified = 0;
    nmax_used = primes[nef] - 1;
    if (!Lfunc_reduce_nmax(L, nmax_used)) {
      fprintf(stderr, "%s: Lfunc_reduce_nmax(%" PRIu64 ") failed\n", label, nmax_used);
      free(ef);
      primesieve_free(primes);
      Lfunc_clear(L);
      return -1;
    }
    nuse = nef;
  }
  acb_poly_t lp;
  acb_poly_init(lp);
  for (size_t i = 0; i < nuse; i++) {
    acb_poly_zero(lp);
    for (int j = 0; j <= degree; j++) {
      acb_t c;
      acb_init(c);
      acb_set_si(c, ef[i * (degree + 1) + j]);
      acb_poly_set_coeff_acb(lp, j, c);
      acb_clear(c);
    }
    Lfunc_use_lpoly(L, primes[i], lp);
  }
  acb_poly_clear(lp);

  ecode |= Lfunc_compute(L);
  int fatal = fatal_error(ecode);
  if (ecode & (ERR_RH_ERROR | ERR_INSUFF_EULER | ERR_NO_RANK | ERR_CONFLICT_RANK | ERR_ZERO_PREC))
    certified = 0;
  if (fatal) certified = 0;

  fprintf(out, "{\"label\": \"%s\", \"degree\": %d, \"conductor\": %" PRIu64
          ", \"motivic_weight\": %d, \"normalisation\": %.1f, \"mus\": [",
          label, degree, conductor, weight, 0.5 * weight);
  for (int i = 0; i < degree; i++) fprintf(out, "%s%g", i ? ", " : "", mus[i]);
  fprintf(out, "], \"self_dual_hint\": %s, \"target_prec\": %" PRId64,
          self_dual ? "true" : "false", target_prec);
  fprintf(out, ", \"nmax_required\": %" PRIu64 ", \"nmax_used\": %" PRIu64
          ", \"primes_required\": %zu, \"euler_factors_given\": %" PRId64
          ", \"euler_factors_used\": %zu, \"last_prime_given\": %" PRIu64,
          nmax_required, nmax_used, nprimes, nef, nuse,
          nef > 0 && (size_t)nef <= nprimes ? primes[nef - 1] : (uint64_t)0);
  fprintf(out, ", \"ecode\": %" PRIu64 ", \"errors\": ", (uint64_t)ecode);
  print_errors_json(out, ecode);
  fprintf(out, ", \"fatal\": %s, \"certified\": %s", fatal ? "true" : "false",
          certified ? "true" : "false");
  if (!fatal) {
    fprintf(out, ", \"wprec\": %" PRId64 ", \"rank\": %" PRId64, Lfunc_wprec(L), Lfunc_rank(L));
    acb_srcptr sg = Lfunc_sign(L);
    fprintf(out, ", \"sign\": {\"re\": ");
    print_arb_json(out, acb_realref(sg));
    fprintf(out, ", \"im\": ");
    print_arb_json(out, acb_imagref(sg));
    acb_srcptr sq = Lfunc_sqrt_sign(L);
    fprintf(out, "}, \"sqrt_sign\": {\"re\": ");
    print_arb_json(out, acb_realref(sq));
    fprintf(out, ", \"im\": ");
    print_arb_json(out, acb_imagref(sq));
    fprintf(out, "}, \"taylor\": ");
    print_arb_json(out, Lfunc_Taylor(L));
    int nsides = self_dual ? 1 : 2;
    for (int side = 0; side < nsides; side++) {
      arb_srcptr z = Lfunc_zeros(L, side);
      fprintf(out, ", \"zeros%s\": [", side ? "_dual" : "");
      for (int i = 0; i < MAX_ZEROS && !arb_is_zero(z + i); i++) {
        if (i) fprintf(out, ", ");
        print_zero_json(out, z + i);
      }
      fprintf(out, "]");
      Lplot_t *pl = Lfunc_plot_data(L, side, 64.0, plot_points);
      if (pl) {
        fprintf(out, ", \"plot%s\": {\"spacing\": %.17g, \"values\": [", side ? "_dual" : "",
                pl->spacing);
        for (uint64_t i = 0; i < pl->n_points; i++)
          fprintf(out, "%s%.17g", i ? ", " : "", pl->points[i]);
        fprintf(out, "]}");
        Lfunc_clear_plot(pl);
      }
    }
  }
  fprintf(out, "}\n");
  fflush(out);
  if (ecode) fprint_errors(stderr, ecode);

  free(ef);
  primesieve_free(primes);
  Lfunc_clear(L);
  return fatal ? -1 : 0;
}

int main(int argc, char **argv) {
  const char *in = NULL, *outname = NULL;
  int self_dual = 0;
  uint64_t plot_points = 2048;
  int64_t target_prec = DEFAULT_TARGET_PREC;
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "--self-dual")) self_dual = 1;
    else if (!strcmp(argv[i], "--plot-points") && i + 1 < argc) plot_points = strtoull(argv[++i], NULL, 10);
    else if (!strcmp(argv[i], "--target-prec") && i + 1 < argc) target_prec = strtoll(argv[++i], NULL, 10);
    else if (!in) in = argv[i];
    else if (!outname) outname = argv[i];
  }
  if (!in) {
    fprintf(stderr, "usage: %s input [output.json] [--self-dual] [--plot-points N] [--target-prec B]\n", argv[0]);
    return 2;
  }
  FILE *input = fopen(in, "r");
  if (!input) { fprintf(stderr, "cannot open %s\n", in); return 1; }
  FILE *out = outname ? fopen(outname, "w") : stdout;
  if (!out) { fprintf(stderr, "cannot open %s\n", outname); return 1; }
  char *line = NULL;
  size_t len = 0;
  int status = 0;
  while (getline(&line, &len, input) != -1) {
    size_t l = strlen(line);
    while (l && (line[l - 1] == '\n' || line[l - 1] == '\r')) line[--l] = 0;
    if (!l || line[0] == '#') continue;
    if (run_line(out, line, self_dual, plot_points, target_prec)) status = 1;
  }
  free(line);
  fclose(input);
  if (out != stdout) fclose(out);
  return status;
}
