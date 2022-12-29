/*
 * @file eigen.h
 * @author Karima SADYKOVA
 * @date 2022-12-29
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#ifndef GIVENS_H
#define GIVENS_H
#include "hessenberg.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
//#include <gmp.h>
//#include <mpfr.h>
#define N 5
#define M 4

/* Given */
void su_given(double **A, double **Q, double **R, int m, int n);
double **su_transpose(int m, int n, double **mat);
void su_make_given(int m, int i, int j, double **R, double **G);
void su_coeff_given(int i, int j, double *c, double *s, double **r);
double **su_transpose(int m, int n, double **mat);
/* QAQ algorithm */
void su_gen_threshold_mat(double **A, int size, int n, double threshold);
void su_gen_mat_check_eigen(double **A, int size, int n, double threshold);
/* Given's utils */
void su_print_matrice(double **mat, int m, int n);
void su_free_mat(int m, double **mat);
double **su_init_mat(int m, int n);
void su_make_id(int m, int n, double **mat);
void su_copy_matrix(double **dst, double **src, int m, int n);
double **su_mul_mat(int m, int n, double **mat1, double **mat2);
double **su_rand_mat(int m, int n);
#endif