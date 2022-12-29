/*
 * @file eigen.h
 * @author Karima SADYKOVA
 * @date 2022-12-29
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#ifndef EIGEN_H
#define EIGEN_H
#include "hessenberg.h"
#include "givens.h"
#define RADIX 2.0

void balanc(double **a, int n);
void elmhes(double **a, int n);
void hqr(double **a, int n, double *wr, double *wi);
void n_eigen(double *_a, int n, double *wr, double *wi);
#endif