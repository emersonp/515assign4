//------------------------------------------------------------------------- 
// This is supporting software for CS415/515 Parallel Programming.
// Copyright (c) Portland State University
//------------------------------------------------------------------------- 

// Jacobi method for solving a Laplace equation.  
//
// Usage: ./jacobi-shm -nl <#locales>
// 
//

config const epsilon = 0.001;	// convergence tolerance
config const verbose = false; 	// printing control
config const n = 8; 	        // mesh size (including boundary)

// Jacobi iteration -- return the iteration count.
// 
proc jacobi(D: domain(2), x: [D] real, epsilon: real) { 
  const ID = D.expand(-1,-1); 	// domain for interior points
  var xnew: [D] real;           // buffer for new values
  var delta: real; 		// measure of convergence 
  var cnt = 0;			// iteration counter

  do {
    forall ij in ID do
      xnew(ij) = (x(ij+(0,1)) + x(ij+(0,-1)) 
                   + x(ij+(1,0)) + x(ij+(-1,0))) / 4.0;
    delta = max reduce abs(xnew[ID] - x[ID]);
    x[ID] = xnew[ID];

    cnt += 1;
    if (verbose) {
      writeln("Iter: ", cnt, " (delta=", delta, ")\n");
      writeln(x);
    }
  } while (delta > epsilon);

  return cnt;
}

// Gauss-Seidel Natural Index iteration -- return the iteration count.
// 
proc natind(D: domain(2), x: [D] real, epsilon: real) { 
  const ID = D.expand(-1,-1); 	// domain for interior points
  var delta: real; 		// measure of convergence 
  var cnt = 0;			// iteration counter
  var proc_delta: [n] real;

  do {
    proc_delta = 0;
    forall ij in ID do {
      var old_x: real;
      old_x = x(ij);
      x(ij) = (x(ij+(0,1)) + x(ij+(0,-1)) 
                   + x(ij+(1,0)) + x(ij+(-1,0))) / 4.0;
      proc_delta[here.id] = max(proc_delta[here.id], abs(x(ij) - old_x));
    }
    delta = max reduce proc_delta;

    cnt += 1;
    if (verbose) {
      writeln("Iter: ", cnt, " (delta=", delta, ")\n");
      writeln(x);
    }
  } while (delta > epsilon);

  return cnt;
}
 
//Gauss-Seidel Red Black iteration -- return the iteration count.
// 
proc redblack(D: domain(2), x: [D] real, epsilon: real) { 
  const ID = D.expand(-1,-1); 	// domain for interior points
  var delta: real; 		// measure of convergence 
  var cnt = 0;			// iteration counter
  var proc_delta: [n] real;

  do {
    proc_delta = 0;
    var proc_delta: [n] real;
    forall ij in ID do {
      if (ij[1] % 2 == ij[2] % 2) {
        var old_x: real;
        old_x = x(ij);
        x(ij) = (x(ij+(0,1)) + x(ij+(0,-1)) 
                     + x(ij+(1,0)) + x(ij+(-1,0))) / 4.0;
        proc_delta[here.id] = max(proc_delta[here.id], abs(x(ij) - old_x));
      }
    }

    forall ij in ID do {
      if (ij[1] % 2 != ij[2] % 2) {
        var old_x: real;
        old_x = x(ij);
        x(ij) = (x(ij+(0,1)) + x(ij+(0,-1)) 
                     + x(ij+(1,0)) + x(ij+(-1,0))) / 4.0;
        proc_delta[here.id] = max(proc_delta[here.id], abs(x(ij) - old_x));
      }
    }

    delta = max reduce proc_delta;

    cnt += 1;
    if (verbose) {
      writeln("Iter: ", cnt, " (delta=", delta, ")\n");
      writeln(x);
    }
  } while (delta > epsilon);

  return cnt;
}

// Main routine.
//
proc main() {
  const D = {0..n-1, 0..n-1};   // domain including boundary points
  var a: [D] real = 0.0;	// mesh array
  a[n-1, 0..n-1] = 1.0;         // - setting boundary values
  a[0..n-1, n-1] = 1.0;

  var b: [D] real = 0.0;	// mesh array
  b[n-1, 0..n-1] = 1.0;         // - setting boundary values
  b[0..n-1, n-1] = 1.0;

  var c: [D] real = 0.0;	// mesh array
  c[n-1, 0..n-1] = 1.0;         // - setting boundary values
  c[0..n-1, n-1] = 1.0;
  
  var cnt = jacobi(D, a, epsilon);
  var natindcount = natind(D, b, epsilon);
  var redblackcount = redblack(D, c, epsilon);
  
  writeln("Mesh size: ", n, " x ", n, ", epsilon=", epsilon, 
            ", total Jacobi iterations: ", cnt);
  writeln("Mesh size: ", n, " x ", n, ", epsilon=", epsilon, 
            ", total nat index iterations: ", natindcount);
  writeln("Mesh size: ", n, " x ", n, ", epsilon=", epsilon, 
            ", total red black iterations: ", redblackcount);
}
