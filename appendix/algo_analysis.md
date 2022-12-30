## CPU time for su_given() and su_hessenberg_rec_threshold() functions which implement Given's and Hessenberg's algorithms
```c
	clock_t temps_initial;
	clock_t temps_final;
	double	temps_cpu1, temps_cpu2;
	FILE *fd = fopen("givens_hessenberg_times.txt", "w");
	double	**A;
	double	**Q;
	double	**R;
	for (int i = 5; i <= 50; i += 5){
			A = su_rand_mat(i, i);
			Q = su_init_mat(i, i);
			R = su_init_mat(i, i);
			/* time for Given's algorithm */
			temps_initial = clock () ;
			su_given(A, Q, R, i, i);
			temps_final = clock () ;
			temps_cpu1 = ((double)(temps_final - temps_initial)) / CLOCKS_PER_SEC;
			su_free_mat(i, Q);
			su_free_mat(i, R);
			su_free_mat(i, A);

			A = su_rand_mat(i, i);
			/* time for Hessenberg's algorithm */
			temps_initial = clock () ;
			su_hessenberg_rec(A, i, i - 1, 0);
			temps_final = clock () ;
			temps_cpu2 = ((double)(temps_final - temps_initial)) / CLOCKS_PER_SEC;
			su_free_mat(i, A);

			fprintf(fd, "%d %lf %lf\n", i, temps_cpu1, temps_cpu2);
	}
	fclose(fd);
```

## CPU time for su_gen_threshold_mat() function which generates Ai matrices using Q*A*Qt formula
```c
	clock_t temps_initial;
    clock_t temps_final;
    double temps_cpu1;
    FILE *fd = fopen("generate_mat_times.txt", "w");
    double **A;

    for (int i = 2; i <= 50; i += 2) {
      A = su_rand_mat(i, i);

      /* generating matices algorithm */
      temps_initial = clock();
      su_gen_threshold_mat(A, i, i, 0);
      temps_final = clock();
      temps_cpu1 = ((double)(temps_final - temps_initial)) / CLOCKS_PER_SEC;
      su_free_mat(i, A);
      fprintf(fd, "%d %lf\n", i, temps_cpu1);
    }
    fclose(fd);
```

## CPU time for su_gen_mat_check_eigen() function which generates Ai matrices using Q*A*Qt formula and calculates eigenvalues of a matrix
```c
	clock_t temps_initial;
    clock_t temps_final;
    double temps_cpu1;
    FILE *fd = fopen("check_eigen_times.txt", "w");
    double **A;

    for (int i = 2; i <= 50; i += 2) {
		A = su_rand_mat(i, i);
		su_hessenberg_elimination(A, i);

      /* check eigen values */
      temps_initial = clock();
      su_gen_mat_check_eigen(A, i, i - 1, 0);
      temps_final = clock();
      temps_cpu1 = ((double)(temps_final - temps_initial)) / CLOCKS_PER_SEC;
      su_free_mat(i, A);
      fprintf(fd, "%d %lf\n", i, temps_cpu1);
    }
    fclose(fd);
```