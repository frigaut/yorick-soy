/*
  SOY: Sparse Operations with Yorick
  Copyright (C) 2004 Ralf Flicker
  Copyright (C) 2010 Marcos van Dam (marcos@flatwavefronts.com)

  This work free software; you can redistribute it and/or
  modify it under the terms of the Creative Commons License
  Attribution-NonCommercial-ShareAlike 2.0. This software
  is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY. Follow the CC links on the distribution
  site for more details: http://homepage.mac.com/rflicker/soy.htm
*/


//==================================================================

plug_in, "soy";
soy_version = "1.3.2";
write,format="SOY %s plugin loaded\n",soy_version;

struct rco
{
  long r;
  long c;
  long n;
  pointer ix;
  pointer jx;
  pointer xn;
  float t;
}

struct rco_d
{
  long r;
  long c;
  long n;
  pointer ix;
  pointer jx;
  pointer xn;
  double t;
}

struct ruo
{
  long r;
  long n;
  pointer ix;
  pointer jx;
  pointer xn;
  pointer xd;
  float t;
}

struct ruo_d
{
  long r;
  long n;
  pointer ix;
  pointer jx;
  pointer xn;
  pointer xd;
  double t;
}

//==================================================================

func sprco(x,ur=,un=)
  /* DOCUMENT func sprco(x,ur=,un=)
     Store a 2D matrix on sparse RCO form (see SOY documentation
     for more info on this and other formats used).
     SEE ALSO: sprco_float.c, sprco_double.c (soy.c)
  */
{
  argc = 2;
  if (ur == []) ur = MR;
  if (un == []) un = MN;
  if (typeof(x) == "float") {
    s = rco();
    s.c = long((dimsof(x))(2));
    s.r = long((dimsof(x))(3));
    s.ix = &(array(long,ur));
    s.jx = &(array(long,un));
    s.xn = &(array(float,un));
    a = [&s,&x];
    tmp = sprco_float(argc,&a);
    return s;}
  else if (typeof(x) == "double") {
    s = rco_d();
    s.c = long((dimsof(x))(2));
    s.r = long((dimsof(x))(3));
    s.ix = &(array(long,ur));
    s.jx = &(array(long,un));
    s.xn = &(array(double,un));
    a = [&s,&x];
    tmp = sprco_double(argc,&a);
    return s;}
  else error,"Unsupported Data Type";
}

extern sprco_float;
/* PROTOTYPE
   long sprco_float(int argc, pointer a)
*/

extern sprco_double;
/* PROTOTYPE
   long sprco_double(int argc, pointer a)
*/

//==================================================================
func spruo(x,ur=,un=)
  /* DOCUMENT func spruo(x)
     Store a square matrix on sparse RUO format.
     SEE ALSO: spruo_float.c, spruo_double.c (soy.c)
  */
{
  argc = 2;
  if (ur == []) ur = MR;
  if (un == []) un = MN;
  if (typeof(x) == "float") {
    s = ruo();
    s.r = long((dimsof(x))(3));
    s.ix = &(array(long,ur));
    s.jx = &(array(long,un));
    s.xn = &(array(float,un));
    s.xd = &(array(float,ur));
    a = [&s,&x];
    tmp = spruo_float(argc,&a);
    return s;}
  else if (typeof(x) == "double") {
    s = ruo_d();
    s.r = long((dimsof(x))(3));
    s.ix = &(array(long,ur));
    s.jx = &(array(long,un));
    s.xn = &(array(double,un));
    s.xd = &(array(double,ur));
    a = [&s,&x];
    tmp = spruo_double(argc,&a);
    return s;}
  else error,"Unsupported Data Type";
}

extern spruo_float;
/* PROTOTYPE
   long spruo_float(int argc, pointer a)
*/

extern spruo_double;
/* PROTOTYPE
   long spruo_double(int argc, pointer a)
*/

//==================================================================
func rcoxv(a,v)
  /* DOCUMENT func rcoxv(a,v)
     Matrix-vector multiplication with RCO matrix "a" and real "v".
     SEE ALSO: rcoxv.c (soy.c)
  */
{
  argc = 3;
  if (typeof(v) == "float" && typeof(*a.xn) == "float") {
    u = array(float,a.r);
    s = [&a,&v,&u];
    tmp = rcoxv_float(argc,&s);
    return u;}
  else if (typeof(v) == "double" && typeof(*a.xn) == "double") {
    u = array(double,a.r);
    s = [&a,&v,&u];
    tmp = rcoxv_double(argc,&s);
    return u;}
  else error,"Unsupported/mixed Data Types";
}

extern rcoxv_float;
/* PROTOTYPE
   int rcoxv_float(int argc, pointer s)
*/

extern rcoxv_double;
/* PROTOTYPE
   int rcoxv_double(int argc, pointer s)
*/

//==================================================================
func ruoxv(a,v)
  /* DOCUMENT func ruoxv(a,v)
     Matrix-vetor multiplication with RUO matrix "a" and real "v".
     SEE ALSO: ruoxv_float.c, ruoxv_double.c (soy.c)
  */
{
  argc = 4;
  if (typeof(v) == "float" && typeof(*a.xn) == "float") {
    u = array(float,a.r);
    w = array(float,a.r);
    s = [&a,&v,&u,&w];
    tmp = ruoxv_float(argc, &s);
    return u;}
  else if (typeof(v) == "double" && typeof(*a.xn) == "double") {
    u = array(double,a.r);
    w = array(double,a.r);
    s = [&a,&v,&u,&w];
    tmp = ruoxv_double(argc, &s);
    return u;}
  else error,"Unsupported/mixed Data Types";
}

extern ruoxv_float;
/* PROTOTYPE
   long ruoxv_float(int argc, pointer s)
*/

extern ruoxv_double;
/* PROTOTYPE
   long ruoxv_double(int argc, pointer s)
*/

//==================================================================
func rcoadd(a,b,ur=,un=)
  /* DOCUMENT func rcoadd(a,b)
     Addition of two RCO matrices "a" and "b".
     SEE ALSO: rcoadd_float.c, rcoadd_double.c (soy.c)
  */
{
  if (typeof(*a.xn) != typeof(*b.xn)) error,"Mixed Data Types";
  if (a.r != b.r || a.c != b.c) error,"Matrices have incompatible dimensions!";
  argc = 5;
  if (ur == []) ur = MR;
  if (un == []) un = MN;
  ss = array(long,ur);
  if (typeof(*a.xn) == "float") {
    c = rco();
    c.r = a.r;
    c.c = a.c;
    c.t = a.t;
    c.ix = &(array(long,ur));
    c.jx = &(array(long,un));
    c.xn = &(array(float,un));
    t = array(float,ur);
    s = [&a,&b,&c,&t,&ss];
    tmp = rcoadd_float(argc, &s);
    return c;}
  else if (typeof(*a.xn) == "double") {
    c = rco_d();
    c.r = a.r;
    c.c = a.c;
    c.t = a.t;
    c.ix = &(array(long,ur));
    c.jx = &(array(long,un));
    c.xn = &(array(double,un));
    t = array(double,ur);
    s = [&a,&b,&c,&t,&ss];
    tmp = rcoadd_double(argc, &s);
    return c;}
  else error,"Unsupported Data Type";
}

extern rcoadd_float;
/* PROTOTYPE
   long rcoadd_float(int argc, pointer s)
*/

extern rcoadd_double;
/* PROTOTYPE
   long rcoadd_double(int argc, pointer s)
*/

//==================================================================
func ruoadd(a,b,ur=,un=)
  /* DOCUMENT func ruoadd(a,b)
     Addition of two square symmetric RUO matrices "a" and "b".
     SEE ALSO: ruoadd_float.c, ruoadd_double.c (soy.c)
  */
{
  if (typeof(*a.xn) != typeof(*b.xn)) error,"Mixed Data Types";
  if (a.r != b.r) error,"Matrices have incompatible dimensions!";
  argc = 5;
  if (ur == []) ur = MR;
  if (un == []) un = MN;
  s = array(long,a.r);
  if (typeof(*a.xn) == "float") {
    c = ruo();
    c.r = a.r;
    c.t = a.t;
    c.ix = &(array(long,ur));
    c.jx = &(array(long,un));
    c.xn = &(array(float,un));
    c.xd = &(array(float,ur));
    tt = array(float,a.n+b.n);
    t = [&a,&b,&c,&tt,&s];
    tmp = ruoadd_float(argc, &t);
    return c;}
  else if (typeof(*a.xn) == "double") {
    c = ruo_d();
    c.r = a.r;
    c.t = a.t;
    c.ix = &(array(long,ur));
    c.jx = &(array(long,un));
    c.xn = &(array(double,un));
    c.xd = &(array(double,ur));
    tt = array(double,a.n+b.n);
    t = [&a,&b,&c,&tt,&s];
    tmp = ruoadd_double(argc, &t);
    return c;}
  else error,"Unsupported Data Type";
}

extern ruoadd_float;
/* PROTOTYPE
   long ruoadd_float(int argc, pointer t)
*/

extern ruoadd_double;
/* PROTOTYPE
   long ruoadd_double(int argc, pointer t)
*/

//==================================================================
func rcoata(a,ur=,un=)
  /* DOCUMENT func rcoata(a)
     Matrix mutiplication of RCO matrix "a" with its transpose from
     the left, e.g. returns transpose(a)#a.
     SEE ALSO: rcoata_float.c, rcoata_double.c (soy.c)
  */
{
  argc = 2;
  if (ur == []) ur = MR;
  if (un == []) un = MN;
  if (typeof(*a.xn) == "float") {
    b = ruo();
    b.r = a.r;
    b.t = a.t;
    b.ix = &(array(long,ur));
    b.jx = &(array(long,un));
    b.xn = &(array(float,un));
    b.xd = &(array(float,ur));
    s = [&a,&b];
    tmp = rcoata_float(argc, &s);
    return b;}
  else if (typeof(*a.xn) == "double") {
    b = ruo_d();
    b.r = a.r;
    b.t = a.t;
    b.ix = &(array(long,ur));
    b.jx = &(array(long,un));
    b.xn = &(array(double,un));
    b.xd = &(array(double,ur));
    s = [&a,&b];
    tmp = rcoata_double(argc, &s);
    return b;}
  else error,"Unsupported Data Type";
}

extern rcoata_float;
/* PROTOTYPE
   long rcoata_float(int argc, pointer s)
*/

extern rcoata_double;
/* PROTOTYPE
   long rcoata_double(int argc, pointer s)
*/

//==================================================================
func rcoatb(a,b,ur=,un=)
  /* DOCUMENT func rcoatb(a)
     Matrix mutiplication of two RCO matrices "a" and transpose "b".
     SEE ALSO: rcoatb_float.c (soy.c)
  */
{
  argc = 3;
  if (ur == []) ur = MR;
  if (un == []) un = MN;
  if (typeof(*a.xn) == "float" && typeof(*b.xn) == "float") {
    c = rco();
    c.r = a.r;
    c.c = b.r;
    c.t = min([a.t,b.t]);
    c.ix = &(array(long,ur));
    c.jx = &(array(long,un));
    c.xn = &(array(float,un));
    s = [&a,&b,&c];
    tmp = rcoatb_float(argc, &s);
    return c;}
  else if (typeof(*a.xn) == "double" && typeof(*b.xn) == "double") {
    c = rco_d();
    c.r = a.r;
    c.c = b.r;
    c.t = min([a.t,b.t]);
    c.ix = &(array(long,ur));
    c.jx = &(array(long,un));
    c.xn = &(array(double,un));
    s = [&a,&b,&c];
    tmp = rcoatb_double(argc, &s);
    return c;}
  else error,"Unsupported or mixed data type(s)";
}

extern rcoatb_float;
/* PROTOTYPE
   long rcoatb_float(int argc, pointer s)
*/

extern rcoatb_double;
/* PROTOTYPE
   long rcoatb_double(int argc, pointer s)
*/

//==================================================================
func rcotr(a)
  /* DOCUMENT func rcotr(a)
     Transposes an RCO matrix. This is actually quite complicated,
     so I use one builtin Yorick function first (sort) to make life
     just a little bit easier.
     Modified 2010/07/29 by Marcos van Dam
     Fixed problem with zero matrices not transposing
  */
{
  argc = 5;
  ur = dimsof(*a.ix);
  un = dimsof(*a.xn);
  if (typeof(*a.xn) == "float") {
    at = rco();
    at.xn = &(array(float,un));
  }
  else if (typeof(*a.xn) == "double") {
    at = rco_d();
    at.xn = &(array(double,un));
  }
  else error,"Unsupported data type";
  at.n = a.n;
  at.r = a.c;
  at.c = a.r;
  at.t = a.t;
  at.ix = &(array(long,ur));
  at.jx = &(array(long,un));
  if (a.n > 0){    
    sjx = long(sort((*a.jx)(1:a.n)));
    hjx = (*a.jx)(sjx);
    ax = array(long,a.c);
    acx = array(long,a.c+1);
    rind = array(long,at.n);
    s = [&ax,&acx,&hjx,&rind,&a];
    if (typeof(*a.xn) == "float") tmp = rcotr_float(argc, &s);
    if (typeof(*a.xn) == "double") tmp = rcotr_double(argc, &s);
    (*at.ix)(1:at.r+1) = acx;
    (*at.jx)(1:at.n) = rind(sjx);
    (*at.xn)(1:at.n) = (*a.xn)(sjx);
  }
  return at;
}

extern rcotr_float;
/* PROTOTYPE
   long rcotr_float(int argc, pointer s)
*/

extern rcotr_double;
/* PROTOTYPE
   long rcotr_double(int argc, pointer s)
*/

//==================================================================
func ruosgs(a,b,d,n,x,v)
  /* DOCUMENT func ruosgs(a,b,n,x,v)
     Symmetric Gauss-Seidel iterations.
  */
{
  argc = 6;
  u = array(float,a.r);
  s = [&a,&b,&d,&u,&v,&x];
  for (nn=1; nn<=n; nn++) {
    tmp = ruosgs_float(argc, &s);
    if (nn < n) u = array(float,a.r);
  }
  return x;
}

extern ruosgs_float;
/* PROTOTYPE
   long ruosgs_float(int argc, pointer s)
*/

//==================================================================
//func rcoxv_mt(a,v)
//  /* DOCUMENT func rcoxv_mt(a,v)
//     Multi-threaded sparse matrix-vector multiplication
//     with RCO matrix "a" and real vector "v".
//     SEE ALSO: rcoxv_threads_float.c (soy.c)
//  */
//{
//  argc = 3;
//  u = array(float,a.r);
//  s = [&a,&v,&u];
//  tmp = rcoxv_thread(argc,&s);
//  return u;
//}
//
//extern rcoxv_thread;
///* PROTOTYPE
//   int rcoxv_thread(int argc, pointer s)
//*/

//==================================================================
func rcoinf(a)
  /* DOCUMENT func rcoinf(a)
     Inflates RCO matrix "a" to its full form, the explicit 2D matrix.
  */
{
  if (typeof(*a.xn) == "float") x = array(float,a.c,a.r);
  else if (typeof(*a.xn) == "double") x = array(double,a.c,a.r);
  else error,"Unsupported Data Type";
  for (i=1; i<=a.r; i++) {
    if ((*a.ix)(i+1) > (*a.ix)(i)) {
      for (j=(*a.ix)(i)+1; j<=(*a.ix)(i+1); j++) {
        x((*a.jx)(j)+1,i)=(*a.xn)(j);}}}
  return x;
}

func ruoinf(a)
  /* DOCUMENT func ruoinf(a)
     Inflates RUO matrix "a" to full form.
  */
{
  if (typeof(*a.xn) == "float") x = array(float,a.r,a.r);
  else if (typeof(*a.xn) == "double") x = array(double,a.r,a.r);
  else error,"Unsupported Data Type";  
  for (i=1; i<=a.r; i++) x(i,i) = (*a.xd)(i);
  for (i=1; i<a.r; i++) {
    if ((*a.ix)(i+1) > (*a.ix)(i)) {
      for (j=(*a.ix)(i)+1; j<=(*a.ix)(i+1); j++) {
        x(i,(*a.jx)(j)+1)=x((*a.jx)(j)+1,i)=(*a.xn)(j);}}}
  return x;
}

func rcox(a,c)
  /* DOCUMENT func rcox(a,c)
     Multiplies RCO matrix "a" with a scalar "c".
  */
{
  if (typeof(*a.xn) != typeof(c)) error,"Mixed Data Types";
  (*a.xn)(1:a.n) = c*(*a.xn)(1:a.n);
}

func ruox(a,c)
  /* DOCUMENT func ruox(a,c)
     Multiplies RUO matrix "a" with scalar "c".
  */
{
  if (typeof(*a.xn) != typeof(c)) error,"Mixed Data Types";
  (*a.xn)(1:a.n) = c*(*a.xn)(1:a.n);
  (*a.xd)(1:a.r) = c*(*a.xd)(1:a.r);
}

func rcos(a)
  /* DOCUMENT func rcos(a)
     Computes the sparseness (or rather, the "fill") of a matrix on RCO format.
  */
{
  if (a.r != 0) {xfill = float(a.n)/float(a.c*a.r);}
  else {xfill = 0;}
  return xfill;
}

func ruos(a)
  /* DOCUMENT func ruos(a)
     Computes the sparseness (or rather, the "fill") of a matrix on RUO format.
  */
{
  if (a.r != 0) {xfill = float(a.n*2+a.r)/float(a.r^2.);}
  else {xfill = 0;}
  return xfill;
}

func spinfo(a)
  /* DOCUMENT func spinfo(a)
     Prints information about a sparse matrix in RCO or RUO format.
  */
{
  if (typeof(a) != "struct_instance") {error,"Argument not a RCO or RUO structure!";}
  members= strtok(strtok(print(structof(a))(2:-1))(2,)," (;")(1,);
  if (numberof(members) != 7) {error,"Argument not a RCO or RUO structure!";}
  if (members(2) == "c") {
    sptype = "RCO";
    xfill = rcos(a)*100.;
  } else {
    sptype = "RUO";
    xfill = ruos(a)*100.;
  }
  ur = numberof(*a.ix);
  un = numberof(*a.xn);
  sfilln = float(a.n)/float(un)*100.;
  sfillr = float(a.r)/float(ur)*100.;
  write,"["+sptype+"]       stored      max.   usage";
  write,format="no. rows : %8i %9i   %4.2f%%\n",a.r,ur,sfillr;
  write,format="elements : %8i %9i   %4.2f%%\n",a.n,un,sfilln;
  if (sptype == "RCO") write,format="   (cols : %8i)\n",a.c;
  write,format="matrix fill : %4.2f%%\n",xfill;
}

//==================================================================
func ruopcg(a,b,x0,&nit,tol=,itmax=,sgs=)
  /* DOCUMENT func ruopcg(a,b,x0,&nit,tol=,itmax=,sgs=)
     Preconditioned conjugate gradient solver for a symmetric positive
     definite sparse linear system, with Jacobi preconditioner. This
     algorithm is implemented straight out of Numerical Recipes, with
     the matrix-vector multiplications carried out sparsely by the
     ruoxv(a,v) function.
     Optionally one may invoke symmetric Gauss-Seidel iterations upon
     the Jacobi preconditioning, by setting the keyword sgs=#iters.
     (at least 1). SGS requires the U, D and L to be externally defined.
     SEE ALSO: ruoxv, aotest
  */
{
  prec = typeof(*a.xn);
  if (sum(prec == [typeof(b),typeof(x0)]) != 2) {
    error,"Inconsistent Data Types a/b/x0";}
  if (itmax == []) itmax = a.r;
  if (tol == []) tol = 1.e-4;
  x = array(double,a.r);
  bnrm = sum(b^2.);
  u = ruoxv(a,x0);
  r = b-u;
  if (sgs == []) {z = r/(*a.xd)(1:a.r);}  // Jacobi preconditioner
  //  else {z = ruosgs_Y(x0,r,sgs,U,D,L);}  // Gauss-Seidel iterations
  else {z = ruosgs(U,L,D,sgs,x0*0.f,r);}  // Gauss-Seidel iterations
  k = 0;
  err = 1.;
  while ((k <= itmax) && (err > tol)) {
    k += 1;
    bknum = sum(z*r);
    if (k > 1) {
      bk = bknum/bkden;
      p = bk*p+z;}
    else p = z;
    bkden = bknum;
    if (prec == "float") u = ruoxv(a,float(p));
    else u = ruoxv(a,double(p));
    z = u;
    akden = sum(z*p);
    ak = bknum/akden;
    x += ak*p;
    r -= ak*z;
    if (sgs == []) {z = r/(*a.xd)(1:a.r);}  // Jacobi preconditioner
    //    else {z = ruosgs(x0,r,sgs,U,D,L);}  // Gauss-Seidel iterations
    else {z = ruosgs(U,L,D,sgs,x0*0.f,r);}  // Gauss-Seidel iterations
    err = sum(r^2.)/bnrm;
  }
  nit = k;
  if (prec == "float") return float(x);
  else return double(x);
}

//==================================================================
func rcobuild(&a,v,t,ur=,un=)
  /* DOCUMENT func rcobuild(a,v,r,t,ur=,un=)
     Appends row-vectors v to an RCO matrix a, at threshold t.
  */
{
  if (ur == []) ur = MR;
  if (un == []) un = MN;
  if (a == [] && typeof(v) == "float") a = rco();
  if (a == [] && typeof(v) == "double") a = rco_d();
  if (*a.xn == [] && typeof(v) == "float") a.xn = &(array(float,un));
  if (*a.xn == [] && typeof(v) == "double") a.xn = &(array(float,un));
  if (*a.ix == []) a.ix = &(array(long,ur));
  if (*a.jx == []) a.jx = &(array(long,un));
  if (a.c == 0) {
    a.c = numberof(v);
    a.r = 0;
    a.t = t;
    a.n = 0;
  }
  a.r += 1;
  tmp0 = (abs(v) >= t);
  if (anyof(tmp0)) {
    tmp = where(tmp0);
    n = numberof(tmp);
    (*a.jx)(1+a.n:a.n+n) = tmp-1;
    (*a.xn)(1+a.n:a.n+n) = v(tmp);
    (*a.ix)(a.r+1) = (*a.ix)(a.r)+n;
    a.n += n;
  } else {
    (*a.ix)(a.r+1) = (*a.ix)(a.r);
  }
}

//==================================================================
func Laplace_FDA(nact,aind,ur=,un=)
  /* DOCUMENT Laplace_FDA(CS,nact,aind)
     Builds the discrete 2D Laplacian operator directly on RCO
     format. Requires a bit of pre-processing (computing "aind"
     - see aotest.i for how to do that).
  */
{
  if (ur == []) ur = nact+2;
  if (un == []) un = nact*5+1;
  w = array(long,5);
  v = [&float([-1.,0.25,0.25,0.25,0.25]),\
       &float([-0.75,0.25,0.25,0.25]),\
       &float([-2./3.,1./3.,1./3.])];
  c = rco();
  c.r = nact;
  c.c = nact;
  c.ix = &(array(long,ur));
  c.jx = &(array(long,un));
  c.xn = &(array(float,un));
  c.t = 0.;
  cp = 0;
  for (i=1; i<=nact; i++) {
    w *=0;
    w(1) = i; // center grid point
    actr = (aind(i,1)+1 == aind(,1)) & (aind(i,2) == aind(,2));
    actu = (aind(i,1) == aind(,1)) & (aind(i,2)+1 == aind(,2));
    actl = (aind(i,1)-1 == aind(,1)) & (aind(i,2) == aind(,2));
    actd = (aind(i,1) == aind(,1)) & (aind(i,2)-1 == aind(,2));
    if (anyof(actr)) w(2) = where(actr); // neighbor right
    if (anyof(actu)) w(3) = where(actu); // neighbor up
    if (anyof(actl)) w(4) = where(actl); // neighbor left
    if (anyof(actd)) w(5) = where(actd); // neighbor down
    z = where(w > 0);
    nz = numberof(z);
    (*c.xn)(cp+1:cp+nz) = *v(6-nz);
    (*c.jx)(cp+1:cp+nz) = w(z)-1;
    cp += nz;
    (*c.ix)(i+1) = cp;
  }
  c.n = cp;
  return c;
}

//==================================================================
func save_rco(a,fn)
  /* DOCUMENT save_rco(a,fn)
     Saves an RCO structure a to the binary file fn by converting
     all of its elements to float (double) and putting them into a
     single vector.
  */
{
  r = a.r;
  c = a.c;
  n = a.n;
  if (typeof(*a.xn) == "float") {
    v = array(float,n*2+r+5);
    v(1:3) = float([n,r,c]);
    v(4:n+3) = (*a.xn)(1:n);
    v(n+4:2*n+3) = float((*a.jx)(1:n));
    v(2*n+4:2*n+r+5) = float((*a.ix)(1:r+2));
  }
  else if (typeof(*a.xn) == "double") {
    v = array(double,n*2+r+5);
    v(1:3) = double([n,r,c]);
    v(4:n+3) = (*a.xn)(1:n);
    v(n+4:2*n+3) = double((*a.jx)(1:n));
    v(2*n+4:2*n+r+5) = double((*a.ix)(1:r+2));
  }
  else error,"Unsupported data type";
  save,createb(fn),v;
}

//==================================================================
func save_ruo(a,fn)
  /* DOCUMENT save_ruo(a,fn)
     Saves an RUO structure a to the binary file fn by converting
     all of its elements to float (double) and putting them into a
     single vector.
  */
{
  r = a.r;
  n = a.n;
  if (typeof(*a.xn) == "float") {
    v = array(float,n*2+r*2+4);
    v(1:2) = float([n,r]);
    v(3:n+2) = (*a.xn)(1:n);
    v(n+3:2*n+2) = float((*a.jx)(1:n));
    v(2*n+3:2*n+r+4) = float((*a.ix)(1:r+2));
    v(2*n+r+5:2*n+2*r+4) = (*a.xd)(1:r);
  }
  else if (typeof(*a.xn) == "double") {
    v = array(double,n*2+r*2+4);
    v(1:2) = double([n,r]);
    v(3:n+2) = (*a.xn)(1:n);
    v(n+3:2*n+2) = double((*a.jx)(1:n));
    v(2*n+3:2*n+r+4) = double((*a.ix)(1:r+2));
    v(2*n+r+5:2*n+2*r+4) = (*a.xd)(1:r);
  }
  else error,"Unsupported data type";
  save,createb(fn),v;
}

//==================================================================
func restore_rco(fn, ur=, un=)
  /* DOCUMENT restore_rco(fn, ur=, un=)
     Returns the RCO structure saved in the file fn by save_rco.
     Modified by Marcos van Dam, August 2010.
     Now pads out the pointers, so that the matrices can be manipulated
  */
{
  if (ur == []) ur = MR;
  if (un == []) un = MN;
  restore,openb(fn),v;
  if (typeof(v) == "float") {a = rco();}
  else if (typeof(v) == "double") {a = rco_d();}
  else error,"Unsupported data type";
  a.n = long(v(1));
  a.r = long(v(2));
  a.c = long(v(3));

  a.ix = &(array(long,ur));
  a.jx = &(array(long,un));
  a.xn = &(array(float,un));

  xn = v(4:a.n+3);
  jx = long(v(a.n+4:2*a.n+3));
  ix = long(v(2*a.n+4:2*a.n+a.r+5));
  
  (*a.xn)(1:a.n) = xn;
  (*a.jx)(1:numberof(jx)) = jx;
  (*a.ix)(1:numberof(ix)) = ix;
  
  return a;
}

//==================================================================
func restore_ruo(fn, ur=, un=)
  /* DOCUMENT restore_ruo(fn, ur=, un=)
     Returns the RUO structure saved in the file fn by save_rco.
     Modified by Marcos van Dam, August 2010.
     Now pads out the pointers, so that the matrices can be manipulated
  */
{
  if (ur == []) ur = MR;
  if (un == []) un = MN;
  restore,openb(fn),v;
  if (typeof(v) == "float") a = ruo();
  else if (typeof(v) == "double") a = ruo_d();
  else error,"Unsupported data type in";
  a.n = long(v(1));
  a.r = long(v(2));

  a.ix = &(array(long,ur));
  a.jx = &(array(long,un));
  a.xn = &(array(float,un));
  a.xd = &(array(float,ur));

  xn = v(3:a.n+2);
  jx = long(v(a.n+3:2*a.n+2));
  ix = long(v(2*a.n+3:2*a.n+a.r+4));
  xd = v(2*a.n+a.r+5:2*a.n+2*a.r+4);

  (*a.xn)(1:numberof(xn)) = xn;
  (*a.jx)(1:numberof(jx)) = jx;
  (*a.ix)(1:numberof(ix)) = ix;
  (*a.xd)(1:numberof(xd)) = xd;

  return a;
}

//==================================================================
func rcodr(&a,r)
  /* DOCUMENT rcodr(&a,r)
     Delete a specific row from an RCO structure.
  */
{
  nel = (*a.ix)(r+1)-(*a.ix)(r);
  if (nel == 0) {
    if (r == 1) {
      (*a.ix)(2:a.r) = (*a.ix)(3:a.r+1);
    } else if (r < a.r) {
      (*a.ix)(r+1:a.r) = (*a.ix)(r+2:a.r+1);
    }
    (*a.ix)(a.r+1) = 0;
  } else {
    if (r == a.r) {
      (*a.jx)(a.n-nel+1:a.n) *= 0;
      (*a.xn)(a.n-nel+1:a.n) *= 0.0f;
      (*a.ix)(a.r+1) = 0; 
    } else if (r == 1) {
      (*a.jx)(1:a.n-nel) = (*a.jx)(nel+1:a.n);
      (*a.xn)(1:a.n-nel) = (*a.xn)(nel+1:a.n);
      (*a.ix)(2:a.r) = (*a.ix)(3:a.r+1)-nel;
      (*a.ix)(a.r+1) = 0;
    } else {
      // modification made by Marcos van Dam, 28 May 2010
      //(*a.jx)((*a.ix)(r):a.n-nel-1) = (*a.jx)((*a.ix)(r+1):a.n-1); //orig
      //(*a.xn)((*a.ix)(r):a.n-nel-1) = (*a.xn)((*a.ix)(r+1):a.n-1); //orig
      (*a.jx)((*a.ix)(r)+1:a.n) = (*a.jx)((*a.ix)(r)+1+nel:a.n+nel); //rev
      (*a.xn)((*a.ix)(r)+1:a.n) = (*a.xn)((*a.ix)(r)+1+nel:a.n+nel); //rev
      
      (*a.ix)(r+1:a.r) = (*a.ix)(r+2:a.r+1)-nel;
      (*a.ix)(a.r+1) = 0;
    }
  }
  a.r -= 1;
  a.n -= nel;
}

//==================================================================
func spcon(&a,b,diag=,ruo=)
  /* DOCUMENT spcon(&a,b,diag=,ruo=)
     Concatenate two RCO matrices row-wise (default), diagonally
     (diag=1), or two RUO matrices diagonally (ruo=1).

     Modified 2010/7/29, Marcos van Dam
     Fixed bugs with concatenating ruo matrices

     Modified 2010/09/04, Marcos van Dam
     Fixed bug concatenating a matrix with a matrix of zeros
  */
{
  if (b.n > 0){
    (*a.jx)(a.n+1:a.n+b.n) = (*b.jx)(1:b.n);
    if (diag == 1 && ruo != 1) (*a.jx)(a.n+1:a.n+b.n) += a.c;
    if (ruo == 1) (*a.jx)(a.n+1:a.n+b.n) += a.r;
    (*a.xn)(a.n+1:a.n+b.n) = (*b.xn)(1:b.n);
  }
  if (ruo == 1) (*a.xd)(a.r+1:a.r+b.r) = (*b.xd)(1:b.r);
  if (ruo == 1) (*a.ix)(a.r+1:a.r+b.r) = (*b.ix)(1:b.r)+(*a.ix)(a.r);
  else (*a.ix)(a.r+2:a.r+b.r+1) = (*b.ix)(2:b.r+1)+(*a.ix)(a.r+1);
  a.r += b.r;
  a.n += b.n;
  if (diag == 1 && ruo != 1) a.c += b.c;
  else if (diag != 1 && ruo != 1) a.c = max([a.c,b.c]);
  // it is possible for two RCO matrices of different widths to
  // be concatennated row-wise - the narrower one is automatically
  // padded with zeroes to the right.
}

//==================================================================
func ruo_UDL(a)
  /* DOCUMENT func ruo_UDL(a)
     Splits an RUO matrix into its upper triangular (U),
     diagonal (D), and lower tringular (L) parts.
  */
{
  D = (*a.xd)(1:a.r);
  U = rco();
  U.jx = &((*a.jx)(1:a.n));
  U.xn = &((*a.xn)(1:a.n));
  U.ix = &((*a.ix)(1:a.r+1));
  (*U.ix)(a.r+1) = (*U.ix)(a.r);
  U.r = a.r;
  U.c = U.r;
  U.n = a.n;
  L = rcotr(U);
  return [&U,&D,&L];
}
//==================================================================
func ruo2rco(a)
  /* DOCUMENT func ruo2rco(a)
     Converts an RUO matrix into RCO.
     Calls rcotr and rcoadd.
  */
{
  d = (*a.xd)(1:a.r);
  u = rco();
  u.jx = &((*a.jx)(1:a.n+2));
  u.xn = &((*a.xn)(1:a.n+2));
  u.ix = &((*a.ix)(1:a.r+3));
  (*u.ix)(a.r+1) = (*u.ix)(a.r);
  u.r = a.r;
  u.c = u.r;
  u.n = a.n;
  l = rcotr(u);

  // now embed the diagonal into u
  dix = (*u.ix)(1:u.r+2);
  dix(1:u.r+1) += long(span(1,u.r+1,u.r+1)-1);
  dix(u.r+2) = dix(u.r+1);
  djx = array(long,u.n+u.r+2);
  dxn = array(float,u.n+u.r+2);
  for (i=1; i<=u.r; i++) {
    djx(dix(i)+1) = i-1;
    dxn(dix(i)+1) = d(i);
    if ( (*u.ix)(i+1) >  (*u.ix)(i)  ) {
      djx(dix(i)+2:dix(i+1)+1) = (*u.jx)((*u.ix)(i)+1:(*u.ix)(i+1)+1);
      dxn(dix(i)+2:dix(i+1)+1) = (*u.xn)((*u.ix)(i)+1:(*u.ix)(i+1)+1);
    }
  }
  u.ix = &(dix);
  u.jx = &(djx);
  u.xn = &(dxn);
  u.n += u.r;
  b = rcoadd(u,l);
  
  return b;
}
//==================================================================
func ruosgs_Y(u,v,n,U,D,L)
  /* DOCUMENT func ruosgs_Y(u,v,n,U,D,L)
     Symmetric Gauss-Seidel iterations.
     SEE ALSO: ruo_UDL
  */
{
  u = v/D;
  for (i=1;i<=n;i++) {
    tmp = v-rcoxv(U,u);
    u = tmp/D;
    tmp = v-rcoxv(L,u);
    u = tmp/D;
  }
  return u;
}
//==================================================================
func intop(dims)
  /* DOCUMENT int(dims)
     Interpolating operators implemented as sparse matrices.
     "dims" is a dimension list with 2^i entries., example:
     test = binop([64,32,16]) will return a 2-element pointer
     array, each pointing to a 3-element array of RCO structures.
     SEE ALSO:
  */
{
  nd = numberof(dims);
  a = array(rco,nd);
  b = array(rco,nd);
  for (i=1;i<=nd;i++) {
    a.ix = &(array(long,MR));
    a.jx = &(array(long,MN));
    a.xn = &(array(float,MN));
    b.ix = &(array(long,MR));
    b.jx = &(array(long,MN));
    b.xn = &(array(float,MN));
  }
  v1 = array(0.5f,2);
  v2 = array(0.25f,4);

  /* // debugging the truth table...
  for (i=1;i<=dims(k);i++) {
    for (j=1;j<=dims(k);j++) {
      //ind = (i-1)*dims(k)+j-1;
      ind = (long((i+1)/2)-1)*dims(k)/2+long((j+1)/2)-1;
      i2 = i%2; j2 = j%2;
      write,format="%2d %2d %2d %d %2d %d %d %d %4d\n",j,i,j2,i2,\
        ((1-j2) && i2),((1-i2) && j2),((1-i2) && (1-j2)),(i2 && j2),ind;
    }
  }
  */

  for (k=1;k<=nd;k++) {
    ii = 1;
    n = 1;
    for (i=1;i<=dims(k);i++) {
      for (j=1;j<=dims(k);j++) {
        //        ind = (i-1)*dims(k)+j-1;
        ind = (long((i+1)/2)-1)*dims(k)/2+long((j+1)/2)-1;
        i2 = i%2; j2 = j%2;
        if ((1-j2) && i2) {             // x-interpolation
          (*a(k).xn)(n:n+1) = v1;
          inds = [ind,ind+1];
          if (j == dims(k)) inds(2) -= dims(k)/2;
          (*a(k).jx)(n:n+1) = inds;
          (*a(k).ix)(ii+1) = (*a(k).ix)(ii)+2;
          n += 2;
        } else if ((1-i2) && j2) {      // y-interpolation
          (*a(k).xn)(n:n+1) = v1;
          inds = [ind,ind+dims(k)/2];
          if (i == dims(k)) inds(2) -= (dims(k)/2)^2;
          (*a(k).jx)(n:n+1) = inds;
          (*a(k).ix)(ii+1) = (*a(k).ix)(ii)+2;
          n += 2;
        } else if ((1-i2) && (1-j2)) {    // xy-interpolation
          (*a(k).xn)(n:n+3) = v2;
          inds = [ind,ind+1,ind+dims(k)/2,ind+dims(k)/2+1];
          if (j == dims(k)) {
            inds(2) -= dims(k)/2;
            inds(4) -= dims(k)/2;
              }
          if (i == dims(k)) {
            inds(3) -= (dims(k)/2)^2;
            inds(4) -= (dims(k)/2)^2;
          }
          (*a(k).jx)(n:n+3) = inds;
          (*a(k).ix)(ii+1) = (*a(k).ix)(ii)+4;
          n += 4;
        } else if (i2 && j2) {        //existing point
          (*a(k).xn)(n) = 1.f;
          (*a(k).jx)(n) = ind;
          (*a(k).ix)(ii+1) = (*a(k).ix)(ii)+1;
          n += 1;
        }
        ii++;
      }
    }
    a(k).r = long(dims(k)^2);
    a(k).c = long(dims(k)^2/4);
    a(k).n = n-1;
    b(k) = rcotr(a(k));
  }
  return [&a,&b];
}

//==================================================================
func spunit(n, precision=)
  /* DOCUMENT spunit(n, precision=)
     Create an identity matrix with n rows and columns
     Marcos van Dam, July 2010
  */
{
  n = long(n);
  if (precision=="double"){
    spidentity = spruo(double(unit(1)));
  }
  else {
    spidentity = spruo(float(unit(1)));    
  }

  spidentity.r = n;
  (*spidentity.xd)(1:n) = 1;
  
  return spidentity;
}              


//==================================================================
func spzeros(m, n, precision=)
/* DOCUMENT  spzeros(n, m, precision=)
   Create an matrix of zeros with m rows and n columns
   Marcos van Dam, July 2010
*/

{
  if (precision=="double"){
    zero_row = sprco(array(double,[2,n,1]));
  }
  else {
    zero_row = sprco(array(float,[2,n,1]));
  }
  
  for (row_counter=1;row_counter<=m;row_counter++){
    if (row_counter == 1){
      zero_matrix = zero_row;
    }
    else {
      spcon, zero_matrix, zero_row;
    }
  }
  return zero_matrix;
}  
