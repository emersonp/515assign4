//------------------------------------------------------------------------- 
// This is supporting software for CS415/515 Parallel Programming.
// Copyright (c) Portland State University
//------------------------------------------------------------------------- 

// Jacobi method for solving a Laplace equation.  
//
// Usage: ./jacobi [N]
// 
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define EPSILON 0.001 	// convergence tolerance
#define VERBOSE 1 	// printing control

// Initialize the mesh with a fixed set of boundary conditions.
// 
void init_array(int n, double a[n][n])  {
  int i, j;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) 
      a[i][j] = 0;
  }
  for (i = 1; i < n-1; i++) {
    a[n-1][i] = 1.0;
    a[i][n-1] = 1.0;
  }
}

// Display the whole mesh.
// 
void print_array(int n, double a[n][n])  {
  int i, j;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++)
      printf("%8.4f ", a[i][j]);
    printf("\n");
  }
}

// Jacobi iteration -- return the iteration count.
// 
int jacobi(int n, double x[n][n], double epsilon) {
  double xnew[n][n];	// buffer for new values    
  double delta;		// measure of convergence   
  int cnt = 0;		// iteration counter              
  int i, j;

  do {	
    delta = 0.0;
    for (i = 1; i < n-1; i++) {
      for (j = 1; j < n-1; j++) {
	xnew[i][j] = (x[i-1][j] + x[i][j-1] + x[i+1][j] + x[i][j+1]) / 4.0;
	delta = fmax(delta, fabs(xnew[i][j] - x[i][j]));
      }
    }	
    for (i = 1; i < n-1; i++) {
      for (j = 1; j < n-1; j++) {
	x[i][j] = xnew[i][j];
      }
    }
    cnt++;
    if (VERBOSE) {
      printf("Iter %d: (delta=%6.4f)\n", cnt, delta);
      print_array(n, x);
    }
  } while (delta > epsilon);
  return cnt;
}

// Gauss-Seidel Natural Index
//
int natind(int n, double x[n][n], double epsilon) {
  double delta, old_x;
  int count = 0;
  int i, j;

  do {
    delta = 0.0;
    for (i = 1; i < n-1; i++) {
      for (j = 1; j < n-1; j++) {
        old_x = x[i][j];
        x[i][j] = (x[i-1][j] + x[i][j-1] + x[i+1][j] + x[i][j+1]) / 4.0;
        delta = fmax(delta, fabs(x[i][j] - old_x));
      }
    }
    count++;
    if (VERBOSE) {
      printf("Iter %d: (delta=%6.4f)\n", count, delta);
      print_array(n, x);
    }
  } while (delta > epsilon);
  return count;
}

// Gauss-Seidel Red Black
//
int redblack(int n, double x[n][n], double epsilon) {
  double delta, old_x;
  int count = 0;
  int i, j;

  do {
    delta = 0.0;
    for (i = 1; i < n-1; i++) {
      for (j = (1 + ((i + count) % 2)); j < n-1; j += 2) {
        old_x = x[i][j];
        x[i][j] = (x[i-1][j] + x[i][j-1] + x[i+1][j] + x[i][j+1]) / 4.0;
        delta = fmax(delta, fabs(x[i][j] - old_x));
      }
    }
    count++;
    if (VERBOSE) {
      printf("Iter %d: (delta=%6.4f)\n", count, delta);
      print_array(n, x);
    }
  } while (delta > epsilon);
  return count / 2 + count % 2; // Divided by two to acknowledge red vs black iterations.
}

// Main routine.
//
int main(int argc, char **argv) {

  int n = 8;  	   	// mesh size, default 8 x 8
  if (argc > 1) {  	// check command line for overwrite
    if ((n = atoi(argv[1])) < 2) {
      printf("Mesh size must must be greater than 2, use default\n");
      n = 8;
    }
  }

  double a[n][n];	// mesh array
  double b[n][n];
  double c[n][n];
  init_array(n, a);
  init_array(n, b);
  init_array(n, c);
  

  // Jacobi iteration, return value is the total iteration number
  int cnt = jacobi(n, a, EPSILON);
  int nat_ind_count = natind(n, b, EPSILON);
  int red_black_count = redblack(n, c, EPSILON);
  printf("Mesh size: %d x %d, epsilon=%6.4f, total Jacobi iterations: %d\n", 
	 n, n, EPSILON, cnt);
  if (VERBOSE) 
    print_array(n, a);
  
  printf("\n\nMesh size: %d x %d, epsilon=%6.4f, total Gauss-Seidel natural index iterations: %d\n", 
	 n, n, EPSILON, nat_ind_count);
  if (VERBOSE) 
    print_array(n, b);

  printf("\n\nMesh size: %d x %d, epsilon=%6.4f, total Gauss-Seidel redblack iterations: %d\n", 
	 n, n, EPSILON, red_black_count);
  if (VERBOSE) 
    print_array(n, c);
}
