/*
  SOY: Sparse Operations with Yorick
  Copyright (C) 2004 Ralf Flicker (rflicker@mac.com)
  Copyright (C) 2010-2012 Marcos van Dam (marcos@flatwavefronts.com)

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

#include <math.h>
#include <stdlib.h>

typedef struct {
  long r,c,n,*ix,*jx;
  float *xn,t;} rco; 

typedef struct {
  long r,c,n,*ix,*jx;
  double *xn,t;} rco_d;

typedef struct {
  long r,n,*ix,*jx;
  float *xn,*xd,t;} ruo; 

typedef struct {
  long r,n,*ix,*jx;
  double *xn,*xd,t;} ruo_d; 

long sprco_float(argc,argv)
int argc;
void *argv[];
{
  rco* s;
  long i,j,k;
  float *x;
  s = (rco *)argv[0];
  x = (float *)argv[1];
  k = 0;
  for (i=0; i < s->r; i++) {
    for (j=0; j < s->c; j++) {
      if (fabs(x[s->c*i+j]) > s->t) {
	(s->xn)[k] = x[s->c*i+j];
	(s->jx)[k] = j;
	k++;}}
    (s->ix)[i+1] = k;}
  s->n = k;
  return k;
}

long sprco_double(argc,argv)
int argc;
void *argv[];
{
  rco_d* s;
  long i,j,k;
  double *x;
  s = (rco_d *)argv[0];
  x = (double *)argv[1];
  k = 0;
  for (i=0; i < s->r; i++) {
    for (j=0; j < s->c; j++) {
      if (fabs(x[s->c*i+j]) > s->t) {
	(s->xn)[k] = x[s->c*i+j];
	(s->jx)[k] = j;
	k++;}}
    (s->ix)[i+1] = k;}
  s->n = k;
  return k;
}

long spruo_float(argc,argv)
int argc;
void *argv[];
{
  ruo* s;
  long i,j,k;
  float *x;
  s = (ruo *)argv[0];
  x = (float *)argv[1];
  k = 0;
  for (i=0; i < s->r; i++) (s->xd)[i] = x[s->r*i+i];
  for (i=0; i < s->r-1; i++) {
    for (j=i+1; j < s->r; j++) {
      if (fabs(x[s->r*i+j]) > s->t) {
	(s->xn)[k] = x[s->r*i+j];
	(s->jx)[k] = j;
	k++;}}
    (s->ix)[i+1] = k;}
  s->n = k;
  return k;
}

long spruo_double(argc,argv)
int argc;
void *argv[];
{
  ruo_d* s;
  long i,j,k;
  double *x;
  s = (ruo_d *)argv[0];
  x = (double *)argv[1];
  k = 0;
  for (i=0; i < s->r; i++) (s->xd)[i] = x[s->r*i+i];
  for (i=0; i < s->r-1; i++) {
    for (j=i+1; j < s->r; j++) {
      if (fabs(x[s->r*i+j]) > s->t) {
	(s->xn)[k] = x[s->r*i+j];
	(s->jx)[k] = j;
	k++;}}
    (s->ix)[i+1] = k;}
  s->n = k;
  return k; 
}

long rcoxv_float(argc,argv)
int argc;
void *argv[];
{
  rco* s;
  long i,j,k;
  float *v, *u;
  s = (rco *)argv[0];
  v = (float *)argv[1];
  u = (float *)argv[2];
  /* struct timeval tod0,tod1; 
     gettimeofday(&tod0,NULL); */ 
  for (i=0; i < s->r; i++) {
    if ((s->ix)[i+1]-(s->ix)[i] > 0) {
      for (j=(s->ix)[i]; j < (s->ix)[i+1]; j++) {
	u[i] += (s->xn)[j]*v[(s->jx)[j]];}}}
  /* gettimeofday(&tod1,NULL); 
     printf("Time X: %u\n",(unsigned long)tod1.tv_usec-(unsigned long)tod0.tv_usec); */
  return j;
}

long rcoxv_double(argc,argv)
int argc;
void *argv[];
{
  rco_d* s;
  long i,j;
  double *v, *u;
  s = (rco_d *)argv[0];
  v = (double *)argv[1];
  u = (double *)argv[2];
  for (i=0; i < s->r; i++) {
    if ((s->ix)[i+1]-(s->ix)[i] > 0) {
      for (j=(s->ix)[i]; j < (s->ix)[i+1]; j++) {
	u[i] += (s->xn)[j]*v[(s->jx)[j]];}}}
      return j;
}

long ruoxv_float(argc,argv)
int argc;
void *argv[];
{
  ruo* a;
  long i,j;
  float *v, *u, *w;
  a = (ruo *)argv[0];
  v = (float *)argv[1];
  u = (float *)argv[2];
  w = (float *)argv[3];
  for (i=0; i < a->r; i++) u[i] = (a->xd)[i]*v[i];
  for (i=0; i < a->r-1; i++) {
    if ((a->ix)[i+1] > (a->ix)[i]) {
      for (j=(a->ix)[i]; j < (a->ix)[i+1]; j++) {
	u[i] += (a->xn)[j]*v[(a->jx)[j]];
	w[(a->jx)[j]] += (a->xn)[j]*v[i];}}}
  for (i=0; i < a->r; i++) u[i] += w[i];
  return j;
}

long ruoxv_double(argc,argv)
int argc;
void *argv[];
{
  ruo_d* a;
  long i,j;
  double *v, *u, *w;
  a = (ruo_d *)argv[0];
  v = (double *)argv[1];
  u = (double *)argv[2];
  w = (double *)argv[3];
  for (i=0; i < a->r; i++) u[i] = (a->xd)[i]*v[i];
  for (i=0; i < a->r-1; i++) {
    if ((a->ix)[i+1] > (a->ix)[i]) {
      for (j=(a->ix)[i]; j < (a->ix)[i+1]; j++) {
	u[i] += (a->xn)[j]*v[(a->jx)[j]];
	w[(a->jx)[j]] += (a->xn)[j]*v[i];}}}
  for (i=0; i < a->r; i++) u[i] += w[i];
  return j;
}

long rcoadd_float(argc,argv)
int argc;
void *argv[];
{
  rco* a;
  rco* b;
  rco* c;
  long i,j,p,*s;
  float *t;
  a = (rco *)argv[0];
  b = (rco *)argv[1];
  c = (rco *)argv[2];
  t = (float *)argv[3];
  s = (long *)argv[4];
  p = 0;
  for (i=0; i < a->c; i++) s[i] = -1;
  for (i=0; i < a->r; i++) {
    (c->ix)[i] = p;
    if ((a->ix)[i+1] > (a->ix)[i]) {
      for (j=(a->ix)[i]; j < (a->ix)[i+1]; j++) {
	(c->jx)[p++] = (a->jx)[j];
	s[(a->jx)[j]] = i;}}
    if ((b->ix)[i+1] > (b->ix)[i]) {
      for (j=(b->ix)[i]; j < (b->ix)[i+1]; j++) {
	if (s[(b->jx)[j]] != i) (c->jx)[p++] = (b->jx)[j];}}}
  (c->ix)[a->r] = p;  
  for (i=0; i < a->r; i++) {
    if ((c->ix)[i+1] > (c->ix)[i]) {
      for (j=(c->ix)[i]; j < (c->ix)[i+1]; j++) t[(c->jx)[j]] = 0.0;
      if ((a->ix)[i+1] > (a->ix)[i]) {
	for (j=(a->ix)[i]; j < (a->ix)[i+1]; j++) t[(a->jx)[j]] = (a->xn)[j];}
      if ((b->ix)[i+1] > (b->ix)[i]) {
	for (j=(b->ix)[i]; j < (b->ix)[i+1]; j++) t[(b->jx)[j]] = t[(b->jx)[j]]+(b->xn)[j];}
      for (j=(c->ix)[i]; j < (c->ix)[i+1]; j++) (c->xn)[j] = t[(c->jx)[j]];}}
  c->n = p;
  return p;
}

long rcoadd_double(argc,argv)
int argc;
void *argv[];
{
  rco_d* a;
  rco_d* b;
  rco_d* c;
  long i,j,p,*s;
  double *t;
  a = (rco_d *)argv[0];
  b = (rco_d *)argv[1];
  c = (rco_d *)argv[2];
  t = (double *)argv[3];
  s = (long *)argv[4];
  p = 0;
  for (i=0; i < a->c; i++) s[i] = -1;
  for (i=0; i < a->r; i++) {
    (c->ix)[i] = p;
    if ((a->ix)[i+1] > (a->ix)[i]) {
      for (j=(a->ix)[i]; j < (a->ix)[i+1]; j++) {
	(c->jx)[p++] = (a->jx)[j];
	s[(a->jx)[j]] = i;}}
    if ((b->ix)[i+1] > (b->ix)[i]) {
      for (j=(b->ix)[i]; j < (b->ix)[i+1]; j++) {
	if (s[(b->jx)[j]] != i) (c->jx)[p++] = (b->jx)[j];}}}
  (c->ix)[a->r] = p;  
  for (i=0; i < a->r; i++) {
    if ((c->ix)[i+1] > (c->ix)[i]) {
      for (j=(c->ix)[i]; j < (c->ix)[i+1]; j++) t[(c->jx)[j]] = 0.0;
      if ((a->ix)[i+1] > (a->ix)[i]) {
	for (j=(a->ix)[i]; j < (a->ix)[i+1]; j++) t[(a->jx)[j]] = (a->xn)[j];}
      if ((b->ix)[i+1] > (b->ix)[i]) {
	for (j=(b->ix)[i]; j < (b->ix)[i+1]; j++) t[(b->jx)[j]] = t[(b->jx)[j]]+(b->xn)[j];}
      for (j=(c->ix)[i]; j < (c->ix)[i+1]; j++) (c->xn)[j] = t[(c->jx)[j]];}}
  c->n = p;
  return p;
}

long ruoadd_float(argc,argv)
int argc;
void *argv[];
{
  ruo* a;
  ruo* b;
  ruo* c;
  long i,j,p,*s;
  float *tt;
  a = (ruo *)argv[0];
  b = (ruo *)argv[1];
  c = (ruo *)argv[2];
  tt = (float *)argv[3];
  s = (long *)argv[4];
  for (i=0; i < a->r; i++) (c->xd)[i] = (a->xd)[i]+(b->xd)[i];
  for (i=0; i < a->r; i++) s[i] = -1;
  p = 0;
  for (i=0; i < a->r-1; i++) {
    (c->ix)[i] = p;
    if ((a->ix)[i+1] > (a->ix)[i]) {
      for (j=(a->ix)[i]; j < (a->ix)[i+1]; j++) {
	(c->jx)[p++] = (a->jx)[j];
	s[(a->jx)[j]] = i;}}
    if ((b->ix)[i+1] > (b->ix)[i]) {
      for (j=(b->ix)[i]; j < (b->ix)[i+1]; j++) {
	if (s[(b->jx)[j]] != i) (c->jx)[p++] = (b->jx)[j];}}}
  (c->ix)[a->r-1] = p;
  for (i=0; i < a->r-1; i++) {
    if ((c->ix)[i+1] > (c->ix)[i]) {
      for (j=(c->ix)[i]; j < (c->ix)[i+1]; j++) tt[(c->jx)[j]] = 0.0;
      if ((a->ix)[i+1] > (a->ix)[i]) {
	for (j=(a->ix)[i]; j < (a->ix)[i+1]; j++) tt[(a->jx)[j]] = (a->xn)[j];}
      if ((b->ix)[i+1] > (b->ix)[i]) {
	for (j=(b->ix)[i]; j < (b->ix)[i+1]; j++) tt[(b->jx)[j]] = tt[(b->jx)[j]]+(b->xn)[j];}
      for (j=(c->ix)[i]; j < (c->ix)[i+1]; j++) (c->xn)[j] = tt[(c->jx)[j]];}}
      c->n = p; 
  return p;
}

long ruoadd_double(argc,argv)
int argc;
void *argv[];
{
  ruo_d* a;
  ruo_d* b;
  ruo_d* c;
  long i,j,p,*s;
  double *tt;
  a = (ruo_d *)argv[0];
  b = (ruo_d *)argv[1];
  c = (ruo_d *)argv[2];
  tt = (double *)argv[3];
  s = (long *)argv[4];
  for (i=0; i < a->r; i++) (c->xd)[i] = (a->xd)[i]+(b->xd)[i];
  for (i=0; i < a->r; i++) s[i] = -1;
  p = 0;
  for (i=0; i < a->r-1; i++) {
    (c->ix)[i] = p;
    if ((a->ix)[i+1] > (a->ix)[i]) {
      for (j=(a->ix)[i]; j < (a->ix)[i+1]; j++) {
	(c->jx)[p++] = (a->jx)[j];
	s[(a->jx)[j]] = i;}}
    if ((b->ix)[i+1] > (b->ix)[i]) {
      for (j=(b->ix)[i]; j < (b->ix)[i+1]; j++) {
	if (s[(b->jx)[j]] != i) (c->jx)[p++] = (b->jx)[j];}}}
  (c->ix)[a->r-1] = p;
  for (i=0; i < a->r-1; i++) {
    if ((c->ix)[i+1] > (c->ix)[i]) {
      for (j=(c->ix)[i]; j < (c->ix)[i+1]; j++) tt[(c->jx)[j]] = 0.0;
      if ((a->ix)[i+1] > (a->ix)[i]) {
	for (j=(a->ix)[i]; j < (a->ix)[i+1]; j++) tt[(a->jx)[j]] = (a->xn)[j];}
      if ((b->ix)[i+1] > (b->ix)[i]) {
	for (j=(b->ix)[i]; j < (b->ix)[i+1]; j++) tt[(b->jx)[j]] = tt[(b->jx)[j]]+(b->xn)[j];}
      for (j=(c->ix)[i]; j < (c->ix)[i+1]; j++) (c->xn)[j] = tt[(c->jx)[j]];}}
      c->n = p; 
  return p;
}

long rcoata_float(argc,argv)
int argc;
void *argv[];
{
  rco* a;
  ruo* b;
  long i,j,k,l,p,ni,nj;
  float mb=0.0;
  a = (rco *)argv[0];
  b = (ruo *)argv[1];
  for (i=0; i<a->r; i++) {
    if ((a->ix)[i+1]-(a->ix)[i] > 0) {  
      for (k=(a->ix)[i]; k<(a->ix)[i+1]; k++) {
	(b->xd)[i] += (a->xn)[k]*(a->xn)[k];}}}
  p = 0;
  for (i=0; i < a->r-1; i++) {
    for (j=i+1; j < a->r; j++) {
      ni = (a->ix)[i+1]-(a->ix)[i];
      nj = (a->ix)[j+1]-(a->ix)[j];
      if (ni > 0 && nj > 0) {
	mb = 0.0;
	for (k=0; k < ni; k++) {
	  for (l=0; l < nj; l++) {
	    if ((a->jx)[(a->ix)[i]+k] == (a->jx)[(a->ix)[j]+l])
	      mb += (a->xn)[(a->ix)[i]+k]*((a->xn)[(a->ix)[j]+l]);}}
	if (fabs(mb) > b->t) {
	  (b->xn)[p] = mb;
	  (b->jx)[p] = j;
	  p++;}}}
    (b->ix)[i+1] = p;}
  b->n = p;
  return p;
}

long rcoata_double(argc,argv)
int argc;
void *argv[];
{
  rco_d* a;
  ruo_d* b;
  long i,j,k,l,p,ni,nj;
  double mb=0.0;
  a = (rco_d *)argv[0];
  b = (ruo_d *)argv[1];
  for (i=0; i<a->r; i++) {
    if ((a->ix)[i+1]-(a->ix)[i] > 0) {  
      for (k=(a->ix)[i]; k<(a->ix)[i+1]; k++) {
	(b->xd)[i] += (a->xn)[k]*(a->xn)[k];}}}
  p = 0;
  for (i=0; i < a->r-1; i++) {
    for (j=i+1; j < a->r; j++) {
      ni = (a->ix)[i+1]-(a->ix)[i];
      nj = (a->ix)[j+1]-(a->ix)[j];
      if (ni > 0 && nj > 0) {
	mb = 0.0;
	for (k=0; k < ni; k++) {
	  for (l=0; l < nj; l++) {
	    if ((a->jx)[(a->ix)[i]+k] == (a->jx)[(a->ix)[j]+l])
	      mb += (a->xn)[(a->ix)[i]+k]*((a->xn)[(a->ix)[j]+l]);}}
	if (fabs(mb) > b->t) {
	  (b->xn)[p] = mb;
	  (b->jx)[p] = j;
	  p++;}}}
    (b->ix)[i+1] = p;}
  b->n = p;
  return p;
}

long rcoatb_float(argc,argv)
int argc;
void *argv[];
{
  rco* a;
  rco* b;
  rco* c;
  long i,j,k,l,p,ni,nj;
  float mb=0.0;
  a = (rco *)argv[0];
  b = (rco *)argv[1];
  c = (rco *)argv[2];
  p = 0;
  for (i=0; i < a->r; i++) {
    for (j=0; j < b->r; j++) {
      ni = (a->ix)[i+1]-(a->ix)[i];
      nj = (b->ix)[j+1]-(b->ix)[j];
      if (ni > 0 && nj > 0) {
	mb = 0.0;
	for (k=0; k < ni; k++) {
	  for (l=0; l < nj; l++) {
	    if ((a->jx)[(a->ix)[i]+k] == (b->jx)[(b->ix)[j]+l])
	      mb += (a->xn)[(a->ix)[i]+k]*((b->xn)[(b->ix)[j]+l]);}}
	if (fabs(mb) > c->t) {
	  (c->xn)[p] = mb;
	  (c->jx)[p] = j;
	  p++;}}}
    (c->ix)[i+1] = p;}
  c->n = p;
  return p;
}

long rcoatb_double(argc,argv)
int argc;
void *argv[];
{
  rco_d* a;
  rco_d* b;
  rco_d* c;
  long i,j,k,l,p,ni,nj;
  double mb=0.0;
  a = (rco_d *)argv[0];
  b = (rco_d *)argv[1];
  c = (rco_d *)argv[2];
  p = 0;
  for (i=0; i < a->r; i++) {
    for (j=0; j < b->r; j++) {
      ni = (a->ix)[i+1]-(a->ix)[i];
      nj = (b->ix)[j+1]-(b->ix)[j];
      if (ni > 0 && nj > 0) {
	mb = 0.0;
	for (k=0; k < ni; k++) {
	  for (l=0; l < nj; l++) {
	    if ((a->jx)[(a->ix)[i]+k] == (b->jx)[(b->ix)[j]+l])
	      mb += (a->xn)[(a->ix)[i]+k]*((b->xn)[(b->ix)[j]+l]);}}
	if (fabs(mb) > c->t) {
	  (c->xn)[p] = mb;
	  (c->jx)[p] = j;
	  p++;}}}
    (c->ix)[i+1] = p;}
  c->n = p;
  return p;
}

long rcotr_float(argc,argv)
int argc;
void *argv[];
{
  rco* a;
  long i,j,nz;
  long *ax, *acx, *hjx, *rind;
  ax = (long *)argv[0];
  acx = (long *)argv[1];
  hjx = (long *)argv[2];
  rind = (long *)argv[3];
  a = (rco *)argv[4];
  nz = 0;
  for (i=0; i < a->n; i++) {ax[hjx[i]] += 1;}
  for (i=1; i < a->c+1; i++) {acx[i] = ax[i-1]+acx[i-1];}
  for (i=0; i < a->r; i++) {
    nz = (a->ix)[i+1]-(a->ix)[i];
    if (nz > 0) {
    for (j=0; j < nz; j++) {rind[(a->ix)[i]+j] = i;}}}
  return i;
}

long rcotr_double(argc,argv)

int argc;
void *argv[];
{
  rco_d* a;
  long i,j,nz;
  long *ax, *acx, *hjx, *rind;
  ax = (long *)argv[0];
  acx = (long *)argv[1];
  hjx = (long *)argv[2];
  rind = (long *)argv[3];
  a = (rco_d *)argv[4];
  nz = 0;
  for (i=0; i < a->n; i++) {ax[hjx[i]] += 1;}
  for (i=1; i < a->c+1; i++) {acx[i] = ax[i-1]+acx[i-1];}
  for (i=0; i < a->r; i++) {
    nz = (a->ix)[i+1]-(a->ix)[i];
    if (nz > 0) {
    for (j=0; j < nz; j++) {rind[(a->ix)[i]+j] = i;}}}
  return i;
}

long ruosgs_float(argc,argv)
int argc;
void *argv[];
{
  rco* a;
  rco* b;
  long i,j;
  float *d, *u, *v, *x;
  a = (rco *)argv[0];
  b = (rco *)argv[1];
  d = (float *)argv[2];
  u = (float *)argv[3];
  v = (float *)argv[4];
  x = (float *)argv[5];
  for (i=0; i < a->r; i++) {
    if ((a->ix)[i+1] > (a->ix)[i]) {
      for (j=(a->ix)[i]; j < (a->ix)[i+1]; j++) {
        u[i] = (a->xn)[j]*x[(a->jx)[j]];
      }
    }
    if (i > 0) {
      if ((b->ix)[i+1] > (b->ix)[i]) {
        for (j=(b->ix)[i]; j < (b->ix)[i+1]; j++) {
          u[i] += (b->xn)[j]*x[(b->jx)[j]];
        } 
      }
    }
    x[i] = (v[i]-u[i])/d[i];
  }
  return j;
} 

