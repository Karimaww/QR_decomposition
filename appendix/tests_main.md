## Test for eigenvalues
```c
	int             i, j;
    static double   u[5], v[5];
    static double   a[5][5] = {{1.0, 6.0, -3.0, -1.0, 7.0},
    {8.0, -15.0, 18.0, 5.0, 4.0}, {-2.0, 11.0, 9.0, 15.0, 20.0},
    {-13.0, 2.0, 21.0, 30.0, -6.0}, {17.0, 22.0, -5.0, 3.0, 6.0}};

    for (i = 0; i <= 4; i++) {
        for (j = 0; j <= 4; j++)
            printf("%13.7e ", a[i][j]);
        printf("\n");
    }
    printf("\n");
    n_eigen(a[0], 5, u, v);
    for (i = 0; i <= 4; i++)
        printf("%13.7e +J %13.7e\n", u[i], v[i]);
    printf("\n");
```
Correct output:
```sh
	4.2961008e+01 +J 0.0000000e+00
	-6.6171383e-01 +J 0.0000000e+00
	-1.5339638e+01 +J -6.7556929e+00
	-1.5339638e+01 +J 6.7556929e+00
	1.9379983e+01 +J 0.0000000e+00
```

## Test for recursive implementation of Hessenberg's algorithm using threshold
```c
  double	**A;
  A = su_rand_mat(N, N);
  printf("A init : \n");
  su_print_matrice(A, N, N);

  su_hessenberg_elimination(A, N);
  printf("A hessenberg : \n");
  su_print_matrice(A, N, N);

  su_hessenberg_rec(A, N, N - 1, 0);

  printf("A finale : \n");
  su_print_matrice(A, N, N);

  su_free_mat(N, A);
```

## Test for Given's algorithm
```c
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
```

## Test to generate matrices
```c
  double	**A;
  A = su_rand_mat(N, N);
  printf("A init : \n");
  su_print_matrice(A, N, N);

  su_hessenberg_elimination(A, N);
  printf("A hessenberg : \n");
  su_print_matrice(A, N, N);

  su_gen_threshold_mat(A, N, N - 1, 0);

  su_free_mat(N, A);
```

## Test to check eigenvalues of generated matrices
```c
	double	**A;
	A = su_rand_mat(N, N);
	printf("A init : \n");
	su_print_matrice(A, N, N);

	su_hessenberg_elimination(A, N);
	printf("A hessenberg : \n");
	su_print_matrice(A, N, N);

	static double	u[50], v[50];

	printf("Eigenvalues of matrix A after Hessenberg: \n");
	n_eigen(A[0], N, u, v);
	for (int i = 0; i < N; i++)
		printf("%13.7e +J %13.7e\n", u[i], v[i]);
	printf("\n");

	su_gen_mat_check_eigen(A, N, N - 1, 0);

	su_free_mat(N, A);
```