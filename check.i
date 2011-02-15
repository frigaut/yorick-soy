// check-plug.i script for SOY, 18 Nov 2004

/*
  SOY: Sparse Operations with Yorick
  Copyright (C) 2004 Ralf Flicker (rflicker@mac.com)
  Copyright (C) 2010 MArcos van Dam (marcos@flatwavefronts.com)

  This work free software; you can redistribute it and/or
  modify it under the terms of the Creative Commons License
  Attribution-NonCommercial-ShareAlike 2.0. This software
  is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY. Follow the CC links on the distribution
  site for more details: http://homepage.mac.com/rflicker/soy.htm
*/
//plug_dir,".";

require,"soy.i";
require,"random.i";

extern MR,MN,errors,errflg;
MR = 50000;
MN = 500000;

errors = [];
errflg = [];

ftol = 1.e-6;
dtol = 1.e-12;

write,"\nSetting error tolerances : ";
write,format="  >> single precision: %5.0e\n",ftol;
write,format="  >> double precision: %5.0e\n\n",dtol;


// Random rectangular matrices
A = float(random([2,200,100])-0.5);
A = A*(abs(A)<0.2);
B = double(random([2,200,100])-0.5);
B = B*(abs(B)<0.2);

// Random symmetrix matrices
S = float(random([2,100,100])-0.5);
S = S(,+)*S(,+);
S = S*(S>0.5);
T = double(random([2,100,100])-0.5);
T = T(,+)*T(,+);
T = T*(T>0.5);

// Random vectors
v = float(random(200));
w = double(random(200));
s = float(random(100));
t = double(random(100));

// Test1: sprco/spruo and rcoinf/ruoinf
write,format="testing sprco(a) and rcoinf (float)...%s",".";
a = sprco(A);
AA = rcoinf(a);
err = max(abs(AA-A));
grow,errors,err;
grow,errflg,(err != 0);
if ((err != 0)) {
  write,format="ERROR OVERFLOW%s","!\n";
} else {
  write,format="OK%s\n",".";
}
write,format="testing sprco(a) and rcoinf (double)...%s",".";
a = sprco(B);
AA = rcoinf(a);
err = max(abs(AA-B));
grow,errors,err;
grow,errflg,(err != 0);
if ( (err != 0)) {
  write,format="ERROR OVERFLOW%s","!\n";
} else {
  write,format="OK%s\n",".";
}
write,format="testing sprco(a,ur=200,un=10000) and rcoinf (float)...%s",".";
a = sprco(A,ur=200,un=10000);
AA = rcoinf(a);
err = max(abs(AA-A));
grow,errors,err;
grow,errflg,(err != 0);
if ( (err != 0)) {
  write,format="ERROR OVERFLOW%s","!\n";
} else {
  write,format="OK%s\n",".";
}
write,format="testing sprco(a,ur=200,un=10000) and rcoinf (double)...%s",".";
a = sprco(B,ur=200,un=10000);
AA = rcoinf(a);
err = max(abs(AA-B));
grow,errors,err;
grow,errflg,(err != 0);
if ( (err != 0)) {
  write,format="ERROR OVERFLOW%s","!\n";
} else {
  write,format="OK%s\n",".";
}
write,format="testing spruo(s) and ruoinf (float)...%s",".";
a = spruo(S);
AA = ruoinf(a);
err = max(abs(AA-S));
grow,errors,err;
grow,errflg,(err != 0);
if ( (err != 0)) {
  write,format="ERROR OVERFLOW%s","!\n";
} else {
  write,format="OK%s\n",".";
}
write,format="testing sprco(a) and rcoinf (double)...%s",".";
a = spruo(T);
AA = ruoinf(a);
err = max(abs(AA-T));
grow,errors,err;
grow,errflg,(err != 0);
if ( (err != 0)) {
  write,format="ERROR OVERFLOW%s","!\n";
} else {
  write,format="OK%s\n",".";
}
write,format="testing sprco(a,ur=200,un=10000) and rcoinf (float)...%s",".";
a = spruo(S,ur=200,un=10000);
AA = ruoinf(a);
err = max(abs(AA-S));
grow,errors,err;
grow,errflg,(err != 0);
if ( (err != 0)) {
  write,format="ERROR OVERFLOW%s","!\n";
} else {
  write,format="OK%s\n",".";
}
write,format="testing sprco(a,ur=200,un=10000) and rcoinf (double)...%s",".";
a = spruo(T,ur=200,un=10000);
AA = ruoinf(a);
err = max(abs(AA-T));
grow,errors,err;
grow,errflg,(err != 0);
if ( (err != 0)) {
  write,format="ERROR OVERFLOW%s","!\n";
} else {
  write,format="OK%s\n",".";
}



// Test 2: rcoxv and ruoxv
write,format="testing rcoxv(a,v) (float)...%s",".";
vv = A(+,)*v(+);
a = sprco(A,ur=200,un=10000);
uu = rcoxv(a,v);
err = max(abs(vv-uu));
grow,errors,err;
grow,errflg,(err > ftol);
if ( (err > ftol)) {
  write,format="ERROR OVERFLOW%s","!\n";
} else {
  write,format="OK%s\n",".";
}
write,format="testing rcoxv(a,v) (double)...%s",".";
a = sprco(B,ur=200,un=10000);
vv = B(+,)*w(+);
uu = rcoxv(a,w);
err = max(abs(vv-uu));
grow,errors,err;
grow,errflg,(err > dtol);
if ( (err > dtol)) {
  write,format="ERROR OVERFLOW%s","!\n";
} else {
  write,format="OK%s\n",".";
}
write,format="testing ruoxv(a,v) (float)...%s",".";
vv = S(+,)*s(+);
a = spruo(S,ur=200,un=10000);
uu = ruoxv(a,s);
err = max(abs(vv-uu));
grow,errors,err;
grow,errflg,(err > ftol);
if ( (err > ftol)) {
  write,format="ERROR OVERFLOW%s","!\n";
} else {
  write,format="OK%s\n",".";
}
write,format="testing ruoxv(a,v) (double)...%s",".";
a = spruo(T,ur=200,un=10000);
vv = T(+,)*t(+);
uu = ruoxv(a,t);
err = max(abs(vv-uu));
grow,errors,err;
grow,errflg,(err > dtol);
if ( (err > dtol)) {
  write,format="ERROR OVERFLOW%s","!\n";
} else {
  write,format="OK%s\n",".";
}


// Test 3: rcoadd and ruoadd
write,format="testing rcoadd(a,b) (float)...%s",".";
C = float(B);
D = A+C;
a = sprco(A,ur=200,un=10000);
c = sprco(C,ur=200,un=10000);
d = rcoadd(a,c);
dd = rcoinf(d);
err = max(abs(D-dd));
grow,errors,err;
grow,errflg,(err > ftol);
if ( (err > ftol)) {
  write,format="ERROR OVERFLOW%s","!\n";
} else {
  write,format="OK%s\n",".";
}
write,format="testing rcoadd(a,b) (double)...%s",".";
C = double(A);
D = B+C;
b = sprco(B,ur=200,un=10000);
c = sprco(C,ur=200,un=10000);
d = rcoadd(b,c);
dd = rcoinf(d);
err = max(abs(D-dd));
grow,errors,err;
grow,errflg,(err > dtol);
if ( (err > dtol)) {
  write,format="ERROR OVERFLOW%s","!\n";
} else {
  write,format="OK%s\n",".";
}
write,format="testing ruoadd(a,b) (float)...%s",".";
U = float(T);
V = S+U;
a = spruo(S,ur=200,un=10000);
c = spruo(U,ur=200,un=10000);
d = ruoadd(a,c);
dd = ruoinf(d);
err = max(abs(V-dd));
grow,errors,err;
grow,errflg,(err > ftol);
if ( (err > ftol)) {
  write,format="ERROR OVERFLOW%s","!\n";
} else {
  write,format="OK%s\n",".";
}
write,format="testing ruoadd(a,b) (double)...%s",".";
U = double(S);
V = U+T;
a = spruo(T,ur=200,un=10000);
c = spruo(U,ur=200,un=10000);
d = ruoadd(a,c);
dd = ruoinf(d);
err = max(abs(V-dd));
grow,errors,err;
grow,errflg,(err > dtol);
if ( (err > dtol)) {
  write,format="ERROR OVERFLOW%s","!\n";
} else {
  write,format="OK%s\n",".";
}
write,format="testing rcoadd(a,b,ur=200,un=20000) (float)...%s",".";
C = float(B);
D = A+C;
a = sprco(A,ur=200,un=10000);
c = sprco(C,ur=200,un=10000);
d = rcoadd(a,c,ur=200,un=20000);
dd = rcoinf(d);
err = max(abs(D-dd));
grow,errors,err;
grow,errflg,(err > ftol);
if ( (err > ftol)) {
  write,format="ERROR OVERFLOW%s","!\n";
} else {
  write,format="OK%s\n",".";
}
write,format="testing rcoadd(a,b,ur=200,un=20000) (double)...%s",".";
C = double(A);
D = B+C;
b = sprco(B,ur=200,un=10000);
c = sprco(C,ur=200,un=10000);
d = rcoadd(b,c,ur=200,un=20000);
dd = rcoinf(d);
err = max(abs(D-dd));
grow,errors,err;
grow,errflg,(err > dtol);
if ( (err > dtol)) {
  write,format="ERROR OVERFLOW%s","!\n";
} else {
  write,format="OK%s\n",".";
}
write,format="testing ruoadd(a,b,ur=200,un=20000) (float)...%s",".";
U = float(T);
V = S+U;
a = spruo(S,ur=200,un=10000);
c = spruo(U,ur=200,un=10000);
d = ruoadd(a,c,ur=200,un=20000);
dd = ruoinf(d);
err = max(abs(V-dd));
grow,errors,err;
grow,errflg,(err > ftol);
if ( (err > ftol)) {
  write,format="ERROR OVERFLOW%s","!\n";
} else {
  write,format="OK%s\n",".";
}
write,format="testing ruoadd(a,b,ur=200,un=20000) (double)...%s",".";
U = double(S);
V = U+T;
a = spruo(T,ur=200,un=10000);
c = spruo(U,ur=200,un=10000);
d = ruoadd(a,c,ur=200,un=20000);
dd = ruoinf(d);
err = max(abs(V-dd));
grow,errors,err;
grow,errflg,(err > dtol);
if ( (err > dtol)) {
  write,format="ERROR OVERFLOW%s","!\n";
} else {
  write,format="OK%s\n",".";
}


// Test 4: rcoata
write,format="testing rcoata(a) (float)...%s",".";
ATA = A(+,)*A(+,);
a = sprco(A,ur=200,un=10000);
d = rcoata(a);
dd = ruoinf(d);
err = max(abs(ATA-dd));
grow,errors,err;
grow,errflg,(err > ftol);
if ( (err > ftol)) {
  write,format="ERROR OVERFLOW%s","!\n";
} else {
  write,format="OK%s\n",".";
}
write,format="testing rcoata(a) (double)...%s",".";
BTB = B(+,)*B(+,);
a = sprco(B,ur=200,un=10000);
d = rcoata(a);
dd = ruoinf(d);
err = max(abs(BTB-dd));
grow,errors,err;
grow,errflg,(err > dtol);
if ( (err > dtol)) {
  write,format="ERROR OVERFLOW%s","!\n";
} else {
  write,format="OK%s\n",".";
}
write,format="testing rcoata(a,ur=200,un=10000) (float)...%s",".";
ATA = A(+,)*A(+,);
a = sprco(A,ur=200,un=10000);
d = rcoata(a,ur=200,un=10000);
dd = ruoinf(d);
err = max(abs(ATA-dd));
grow,errors,err;
grow,errflg,(err > ftol);
if ( (err > ftol)) {
  write,format="ERROR OVERFLOW%s","!\n";
} else {
  write,format="OK%s\n",".";
}
write,format="testing rcoata(a,ur=200,un=10000) (double)...%s",".";
BTB = B(+,)*B(+,);
a = sprco(B,ur=200,un=10000);
d = rcoata(a,ur=200,un=10000);
dd = ruoinf(d);
err = max(abs(BTB-dd));
grow,errors,err;
grow,errflg,(err > dtol);
if ( (err > dtol)) {
  write,format="ERROR OVERFLOW%s","!\n";
} else {
  write,format="OK%s\n",".";
}



// Test 5: rcoatb
write,format="testing rcoatb(a,b) (float)...%s",".";
C = float(B);
ATC = C(+,)*A(+,);
a = sprco(A,ur=200,un=10000);
c = sprco(C,ur=200,un=10000);
d = rcoatb(a,c);
dd = rcoinf(d);
err = max(abs(ATC-dd));
grow,errors,err;
grow,errflg,(err > ftol);
if ( (err > ftol)) {
  write,format="ERROR OVERFLOW%s","!\n";
} else {
  write,format="OK%s\n",".";
}
write,format="testing rcoatb(a,b) (double)...%s",".";
C = double(A);
ATC = C(+,)*B(+,);
a = sprco(B,ur=200,un=10000);
c = sprco(C,ur=200,un=10000);
d = rcoatb(a,c);
dd = rcoinf(d);
err = max(abs(ATC-dd));
grow,errors,err;
grow,errflg,(err > dtol);
if ( (err > dtol)) {
  write,format="ERROR OVERFLOW%s","!\n";
} else {
  write,format="OK%s\n",".";
}
write,format="testing rcoatb(a,b,ur=200,un=10000) (float)...%s",".";
C = float(B);
ATC = C(+,)*A(+,);
a = sprco(A,ur=200,un=10000);
c = sprco(C,ur=200,un=10000);
d = rcoatb(a,c,ur=200,un=10000);
dd = rcoinf(d);
err = max(abs(ATC-dd));
grow,errors,err;
grow,errflg,(err > ftol);
if ( (err > ftol)) {
  write,format="ERROR OVERFLOW%s","!\n";
} else {
  write,format="OK%s\n",".";
}
write,format="testing rcoatb(a,b,ur=200,un=10000) (double)...%s",".";
C = double(A);
ATC = C(+,)*B(+,);
a = sprco(B,ur=200,un=10000);
c = sprco(C,ur=200,un=10000);
d = rcoatb(a,c,ur=200,un=10000);
dd = rcoinf(d);
err = max(abs(ATC-dd));
grow,errors,err;
grow,errflg,(err > dtol);
if ( (err > dtol)) {
  write,format="ERROR OVERFLOW%s","!\n";
} else {
  write,format="OK%s\n",".";
}



// Test 6: rcoatb
write,format="testing rcotr(a) (float)...%s",".";
AT = transpose(A);
a = sprco(A,ur=210,un=10000);
at = rcotr(a);
dd = rcoinf(at);
err = max(abs(AT-dd));
grow,errors,err;
grow,errflg,(err > ftol);
if ( (err > ftol)) {
  write,format="ERROR OVERFLOW%s","!\n";
} else {
  write,format="OK%s\n",".";
}
write,format="testing rcotr(a) (double)...%s",".";
BT = transpose(B);
a = sprco(B,ur=210,un=10000);
at = rcotr(a);
dd = rcoinf(at);
err = max(abs(BT-dd));
grow,errors,err;
grow,errflg,(err > ftol);
if ( (err > ftol)) {
  write,format="ERROR OVERFLOW%s","!\n";
} else {
  write,format="OK%s\n",".";
}


write,"\nFinished!\n";

if (anyof(errflg)) {
  write,format="The script encountered %i (recoverable) error(s).\n",sum(errflg);
  write,format="The largest error by absolute value was %5.3e.\n",max(errors);
  write,"If this is much greater than the single precision error";
  write,format=" tolerance (%4.0e), then something went seriously wrong.\n",ftol;
  write,"(if it is close, then no worries)\n"
  write,"To submit a bug report, email: rflicker@mac.com.\n";
} else {
  write,"All functions returned without errors.\n";
}
