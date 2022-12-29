/*
 * @file eigen.h
 * @author Karima SADYKOVA
 * @date 2022-12-29
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include "../include/givens.h"
#include "../include/eigen.h"

/*************************
 * balance a real matrix *
 *************************/
void balanc(double **a, int n) {
  int i, j, last = 0;
  double s, r, g, f, c, sqrdx;
  sqrdx = RADIX * RADIX;
  while (last == 0) {
    last = 1;
    for (i = 0; i < n; i++) {
      r = c = 0.0;
      for (j = 0; j < n; j++)
        if (j != i) {
          c += fabs(a[j][i]);
          r += fabs(a[i][j]);
        }
      if (c != 0.0 && r != 0.0) {
        g = r / RADIX;
        f = 1.0;
        s = c + r;
        while (c < g) {
          f *= RADIX;
          c *= sqrdx;
        }
        g = r * RADIX;
        while (c > g) {
          f /= RADIX;
          c /= sqrdx;
        }
        if ((c + r) / f < 0.95 * s) {
          last = 0;
          g = 1.0 / f;
          for (j = 0; j < n; j++)
            a[i][j] *= g;
          for (j = 0; j < n; j++)
            a[j][i] *= f;
        }
      }
    }
  }
}

/*****************************************************
 * convert a non-symmetric matrix to Hessenberg form *
 *****************************************************/
void elmhes(double **a, int n) {
  int i, j, m;
  double y, x;
  for (m = 1; m < n - 1; m++) {
    x = 0.0;
    i = m;
    for (j = m; j < n; j++) {
      if (fabs(a[j][m - 1]) > fabs(x)) {
        x = a[j][m - 1];
        i = j;
      }
    }
    if (i != m) {
      for (j = m - 1; j < n; j++)
        SWAP(a[i][j], a[m][j]);
      for (j = 0; j < n; j++)
        SWAP(a[j][i], a[j][m]);
    }
    if (x != 0.0) {
      for (i = m + 1; i < n; i++) {
        if ((y = a[i][m - 1]) != 0.0) {
          y /= x;
          a[i][m - 1] = y;
          for (j = m; j < n; j++)
            a[i][j] -= y * a[m][j];
          for (j = 0; j < n; j++)
            a[j][m] += y * a[j][i];
        }
      }
    }
  }
}
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
/**************************************
 * QR algorithm for Hessenberg matrix *
 **************************************/
void hqr(double **a, int n, double *wr, double *wi) {
  int nn, m, l, k, j, its, i, mmin;
  double z, y, x, w, v, u, t, s, r, q, p, anorm;
  p = q = r = 0.0;
  anorm = 0.0;
  for (i = 0; i < n; i++)
    for (j = i - 1 > 0 ? i - 1 : 0; j < n; j++)
      anorm += fabs(a[i][j]);
  nn = n - 1;
  t = 0.0;
  while (nn >= 0) {
    its = 0;
    do {
      for (l = nn; l > 0; l--) {
        s = fabs(a[l - 1][l - 1]) + fabs(a[l][l]);
        if (s == 0.0)
          s = anorm;
        if (fabs(a[l][l - 1]) + s == s) {
          a[l][l - 1] = 0.0;
          break;
        }
      }
      x = a[nn][nn];
      if (l == nn) {
        wr[nn] = x + t;
        wi[nn--] = 0.0;
      } else {
        y = a[nn - 1][nn - 1];
        w = a[nn][nn - 1] * a[nn - 1][nn];
        if (l == nn - 1) {
          p = 0.5 * (y - x);
          q = p * p + w;
          z = sqrt(fabs(q));
          x += t;
          if (q >= 0.0) {
            z = p + SIGN(z, p);
            wr[nn - 1] = wr[nn] = x + z;
            if (z != 0.0)
              wr[nn] = x - w / z;
            wi[nn - 1] = wi[nn] = 0.0;
          } else {
            wr[nn - 1] = wr[nn] = x + p;
            wi[nn - 1] = -(wi[nn] = z);
          }
          nn -= 2;
        } else {
          if (its == 30) {
            fprintf(stderr, "[hqr] too many iterations.\n");
            break;
          }
          if (its == 10 || its == 20) {
            t += x;
            for (i = 0; i < nn + 1; i++)
              a[i][i] -= x;
            s = fabs(a[nn][nn - 1]) + fabs(a[nn - 1][nn - 2]);
            y = x = 0.75 * s;
            w = -0.4375 * s * s;
          }
          ++its;
          for (m = nn - 2; m >= l; m--) {
            z = a[m][m];
            r = x - z;
            s = y - z;
            p = (r * s - w) / a[m + 1][m] + a[m][m + 1];
            q = a[m + 1][m + 1] - z - r - s;
            r = a[m + 2][m + 1];
            s = fabs(p) + fabs(q) + fabs(r);
            p /= s;
            q /= s;
            r /= s;
            if (m == l)
              break;
            u = fabs(a[m][m - 1]) * (fabs(q) + fabs(r));
            v = fabs(p) *
                (fabs(a[m - 1][m - 1]) + fabs(z) + fabs(a[m + 1][m + 1]));
            if (u + v == v)
              break;
          }
          for (i = m; i < nn - 1; i++) {
            a[i + 2][i] = 0.0;
            if (i != m)
              a[i + 2][i - 1] = 0.0;
          }
          for (k = m; k < nn; k++) {
            if (k != m) {
              p = a[k][k - 1];
              q = a[k + 1][k - 1];
              r = 0.0;
              if (k + 1 != nn)
                r = a[k + 2][k - 1];
              if ((x = fabs(p) + fabs(q) + fabs(r)) != 0.0) {
                p /= x;
                q /= x;
                r /= x;
              }
            }
            if ((s = SIGN(sqrt(p * p + q * q + r * r), p)) != 0.0) {
              if (k == m) {
                if (l != m)
                  a[k][k - 1] = -a[k][k - 1];
              } else
                a[k][k - 1] = -s * x;
              p += s;
              x = p / s;
              y = q / s;
              z = r / s;
              q /= p;
              r /= p;
              for (j = k; j < nn + 1; j++) {
                p = a[k][j] + q * a[k + 1][j];
                if (k + 1 != nn) {
                  p += r * a[k + 2][j];
                  a[k + 2][j] -= p * z;
                }
                a[k + 1][j] -= p * y;
                a[k][j] -= p * x;
              }
              mmin = nn < k + 3 ? nn : k + 3;
              for (i = l; i < mmin + 1; i++) {
                p = x * a[i][k] + y * a[i][k + 1];
                if (k != (nn)) {
                  p += z * a[i][k + 2];
                  a[i][k + 2] -= p * r;
                }
                a[i][k + 1] -= p * q;
                a[i][k] -= p;
              }
            }
          }
        }
      }
    } while (l + 1 < nn);
  }
}
/*********************************************************
 * calculate eigenvalues for a non-symmetric real matrix *
 *********************************************************/
void n_eigen(double *_a, int n, double *wr, double *wi) {
  int i;
  double **a = (double **)calloc(n, sizeof(void *));
  for (i = 0; i < n; ++i)
    a[i] = _a + i * n;
  balanc(a, n);
  elmhes(a, n);
  hqr(a, n, wr, wi);
  free(a);
}