/*
 * @file eigen.h
 * @author Karima SADYKOVA
 * @date 2022-12-29
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include "../include/eigen.h"
#include "../include/givens.h"
#include "../include/hessenberg.h"
#include <time.h>

int main(void) {

/* Test for Given's algorithm */
  double	**A;
  double	**Q;
  double	**R;
 
  A = su_rand_mat(N, N);

  Q = su_init_mat(M, M);
  R = su_init_mat(M, N);

  su_given(A, Q, R, M, N);
  printf("final Q : \n");
  su_print_matrice(Q, M, M);
  printf("final R : \n");
  su_print_matrice(R, M, N);
  su_free_mat(M, Q);
  su_free_mat(M, R);

  return (0);
}