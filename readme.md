# README - Given's and Hessenberg's algorithm #
## _Supervised by Mr. Vincent Neiger and Mr. Jeremy Berthomieu_

This project was carried out by Karima SADYKOVA and Kimia SADREDDINI in the context of the EU MODEL: Numerical and Symbolic Algorithms Modeling. It was supervised by Mr. Jeremy Berthomieu in lectures and Vincent Neiger in tutorials.

## Global variables you can change at compilation

- M : the number of rows
- N : the number of colums

To change these variables change in either givens.h as such

```c
# define N <the number of columns of your matrice>
# define M <the number of rows of your matrice>
```

or directly in the makefile change the line

```sh
CFLAGS		= -Wall -Wextra -Werror -lm
```
to

```sh
CFLAGS		= -Wall -Wextra -Werror -lm -D N=<the number of columns of your matrice> M=<the number of rows of your matrice> 
```

## Files and implemented functions

- givens.c : implements a QR decomposition using Given's method
- givens_utils.c : useful functions used in givens.c
- hessenberg.c : implements upper Hessenberg
- eigen.c : calculates eigen values of a non-symmetric matrix
- main.c : tester, to test with known matrices uncomment the two versions of matrice A and change the macros M and N accordingly
- eigen.c : compute the eigenvalues of matrixes
- *.h : header files

## Run the program

Run:

```sh
make && gcc givens.a -lm && ./a.out 
```

after making the correct changes in the main according to what you want to test.

All tests are provided in **appendix/tests_main.md**

You can also run:

```sh
make re # delete and make the project again
make clean # delete all the compiled files
make fclean # delete everything
```
## Gnuplot for algorithm's analysis

Algorithms analyses used are in **appendix/algo_analyse.md**. For each analysis we used gnuplot to create graphes.

Use the command below for Given's and Hessenberg CPU time:
```sh
gnuplot -p < appendix/given_hessenberg.txt
```

Use the command below to show Generating Matrices Algorithm's CPU time:
```sh
gnuplot -p < appendix/gen_mat.txt
```

Use the command below to show Matrix Eigenvalues Calculating Algorithm's CPU time:
```sh
gnuplot -p < appendix/check_eigen.txt
```