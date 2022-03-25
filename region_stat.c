#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>


// note that in order to use memcpy we will need to use
// int instead of size_t..
struct regions {
  size_t length;
  size_t capacity;
  int *start;
  int *end;
  double *score;
};

struct regions init_regions(size_t capacity){
  struct regions reg;
  reg.length = 0;
  reg.capacity = capacity;
  reg.start = malloc(sizeof(int) * capacity);
  reg.end = malloc(sizeof(int) * capacity);
  reg.score = malloc(sizeof(double) * capacity);
  return(reg);
}

void free_regions(struct regions *reg){
  free(reg->start);
  free(reg->end);
  free(reg->score);
}

void grow_regions(struct regions *reg){
  size_t capacity = reg->capacity * 2;
  int *start = malloc(sizeof(int) * capacity);
  int *end = malloc(sizeof(int) * capacity);
  double *score = malloc(sizeof(double) * capacity);
  memcpy( reg->start, start, sizeof(int) * reg->capacity );
  memcpy( reg->end, end, sizeof(int) * reg->capacity );
  memcpy( reg->score, score, sizeof(double) * reg->capacity );
  free_regions(reg);
  reg->capacity = capacity;
  reg->start = start;
  reg->end = end;
  reg->score = score;
}

void push_regions(struct regions *reg, int start, int end, double score){
  if(reg->length == reg->capacity)
    grow_regions(reg);
  size_t i = reg->length;
  reg->start[i] = start;
  reg->end[i] = end;
  reg->score[i] = score;
  reg->length = i + 1;
}



// Given a vector of ranks (ranks), the population size (N)
// and a windowed size (w) will calculate a positive
// enrichment score for windows of size w along the
// the r. Fills the scores into the vector scores, which
// must be preallocated and of the correct size (1 + r_n - w)
void window_scores(double *ranks, size_t N, size_t w, size_t r_n,
		   double *scores){
  double score = 0;
  for(size_t i=0; i < w && i < r_n; ++i){
    score += log( ranks[i] / (1 + N - ranks[i]) );
  }
  scores[0] = score;
  for(size_t i=w; i < r_n; ++i){
    score -= log( ranks[i-w] / (1 + N - ranks[i-w]) );
    score += log( ranks[i] / (1 + N - ranks[i]) );
    scores[1+i-w] = score;
  }
}

struct regions define_regions(double *ranks, size_t N, size_t r_n, size_t min_wl){
  struct regions reg = init_regions(100);
  size_t b = 0;
  size_t e = 0;
  size_t max_i = 0;
  double score = 0;
  double last_score = 0;
  double max_score = 0;
  while(b < r_n && e < r_n){
    score = last_score + log( ranks[e] / (1 + N - ranks[e]) );
    if((last_score > 0 && score <= 0) || e == r_n - 1){ // end of
      if(1 + max_i - b >= min_wl)
	push_regions(&reg, b, max_i, max_score);
      score = last_score = max_score = 0;
      b = e = max_i = max_i + 1;
      continue;
    }
    if(score < 0){
      score = max_score = last_score = 0;
      b = e = e + 1;
      continue;
    }
    if(max_score < score){
      max_score = score;
      max_i = e;
    }
    last_score = score;
    ++e;
  }
  if(1 + max_i - b >= min_wl)
    push_regions(&reg, b, max_i, max_score);
  return(reg);
}

SEXP r_window_scores(SEXP ranks_r, SEXP N_r, SEXP w_r){
  if(!isReal(ranks_r) || !isInteger(N_r) || !isInteger(w_r))
    error("Arguments should be ranks (double), N (integer), w (integer)");
  int r_n = length(ranks_r);
  int N_n = length(N_r);
  int w_n = length(w_r);
  if(N_n != 1 || w_n != 1)
    error("N and w should contain a single value each");
  int N = asInteger(N_r);
  int w = asInteger(w_r);
  if(w > r_n || w < 1)
    error("w must be a positive integer that is not longer than ranks");
  double *ranks = REAL(ranks_r);
  SEXP scores_r = PROTECT(allocVector(REALSXP, 1 + r_n - w));
  double *scores = REAL(scores_r);
  window_scores( ranks, N, w, r_n, scores );
  UNPROTECT(1);
  return(scores_r);
}


SEXP r_define_regions(SEXP ranks_r, SEXP N_r, SEXP min_wl_r){
  if(!isReal(ranks_r) || !isInteger(N_r) || !isInteger(min_wl_r))
    error("Arguments should be ranks (double), N (integer), min_wl (integer)");
  int r_n = length(ranks_r);
  int N_n = length(N_r);
  int min_wl_n = length(min_wl_r);
  if(N_n != 1 || min_wl_n != 1)
    error("N and min_wl should be of length 1");
  int N = asInteger(N_r);
  int min_wl = asInteger(min_wl_r);
  if(N < min_wl || N < 1)
    error("N must be positive and larger than min_wl");
  double *ranks = REAL(ranks_r);
  // Return a list of 3 components, starts, ends, and scores
  // This can be turned into a dataframe in R
  struct regions reg = define_regions( ranks, N, r_n, min_wl );
  SEXP regions_r = PROTECT( allocVector( VECSXP, 3 ));
  SET_VECTOR_ELT( regions_r, 0, allocVector(INTSXP, reg.length) );
  SET_VECTOR_ELT( regions_r, 1, allocVector(INTSXP, reg.length) );
  SET_VECTOR_ELT( regions_r, 2, allocVector(REALSXP, reg.length) );
  int *reg_start = INTEGER( VECTOR_ELT(regions_r, 0));
  int *reg_end = INTEGER( VECTOR_ELT(regions_r, 1));
  double *reg_score = REAL( VECTOR_ELT(regions_r, 2));
  memcpy( reg_start, reg.start, sizeof(int) * reg.length );
  memcpy( reg_end, reg.end, sizeof(int) * reg.length );
  memcpy( reg_score, reg.score, sizeof(double) * reg.length );
  free_regions(&reg);
  UNPROTECT(1);
  return(regions_r);
}
