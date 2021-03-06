//------------------------------------------------------------------------- 
// This is supporting software for CS415/515 Parallel Programming.
// Copyright (c) Portland State University
//------------------------------------------------------------------------- 

// Domain Map Examples
//
// Usage: ./domMap -nl <#locales>
// 
//

use BlockDist, CyclicDist;

config const n = 8;
const D = {1..n, 1..n};

// A 2D Block-distributed domain BlockD and a Block-distributed 
// array BA declared over the domain.
//
const BlockD = D dmapped Block(boundingBox=D);
//var BA = 0; 
var BA: [BlockD] int;
var InnerBA: subdomain(BlockD) = BlockD[2..n-1, 2..n-1];

// A 2D Cyclic-distributed domain CyclicD and a Cyclic-distributed 
// array CA declared over the domain.
//
const CyclicD = D dmapped Cyclic(startIdx=D.low);
var CA: [CyclicD] int;
var InnerCA: subdomain(CyclicD) = CyclicD[2..n-1, 2..n-1];

// To illustrate how the index set is distributed across locales,
// we'll use forall loop to initialize each array element to the
// locale ID that stores that element.
//
forall e in InnerBA do
  BA[e] = here.id;

forall e in InnerCA do
  CA[e] = here.id;

// Output the arrays to visually see how the elements are
// partitioned across the locales.
//
writeln("Block Array:");
writeln(BA);

writeln("Cyclic Array:");
writeln(CA);
