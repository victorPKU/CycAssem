/*                              SDOCK2.0
 *   Global Protein-Protein Docking Using Stepwise Force-Field Potentials   
 *                Written by ZHANG Changsheng 
 *  Copyright (c) 2024 Molecular Design Lab at Peking University, Beijing, China  
 
 *    contact: changshengzhang@pku.edu.cn, lhlai@pku.edu.cn
 *
 * surface.c : Calculate solvent accessible surface area for every protein atom.
 *  The algorithm is from NSC2.0 (Frank Eisenhaber etc.)
 ***************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "surface.h"
#include "mem.h"
	
float*	xpunsp=NULL;
float	del_cube;
int*    ico_wk=NULL, *ico_pt=NULL;
int     n_dot, ico_cube, last_n_dot=0, last_densit=0, last_unsp=0;
int     last_cubus=0;
float	rg, rh;

void divarc(float x1, float y1, float z1,float x2, float y2, float z2,\
            int div1, int div2, float *xr, float *yr, float *zr)
{
	float xd, yd, zd, dd, d1, d2, s, x, y, z;
	float phi, sphi, cphi;

	xd = y1*z2-y2*z1;
	yd = z1*x2-z2*x1;
	zd = x1*y2-x2*y1;
	dd = sqrt(xd*xd+yd*yd+zd*zd);
	d1 = x1*x1+y1*y1+z1*z1;
	d2 = x2*x2+y2*y2+z2*z2;
 
	phi = asin(dd/sqrt(d1*d2));
	phi = phi*((float)div1)/((float)div2);
	sphi = sin(phi); cphi = cos(phi);
	s  = (x1*xd+y1*yd+z1*zd)/dd;

	x = xd*s*(1.-cphi)/dd + x1 * cphi + (yd*z1-y1*zd)*sphi/dd;
	y = yd*s*(1.-cphi)/dd + y1 * cphi + (zd*x1-z1*xd)*sphi/dd;
	z = zd*s*(1.-cphi)/dd + z1 * cphi + (xd*y1-x1*yd)*sphi/dd;
	dd = sqrt(x*x+y*y+z*z);
	*xr = x/dd; *yr = y/dd; *zr = z/dd;
}

static void icosaeder_vertices(float *xus) {
  rh = sqrt(1.-2.*cos(TORAD(72.)))/(1.-cos(TORAD(72.)));
  rg = cos(TORAD(72.))/(1.-cos(TORAD(72.)));
  xus[ 0] = 0.;                  xus[ 1] = 0.;                  xus[ 2] = 1.;
  xus[ 3] = rh*cos(TORAD(72.));  xus[ 4] = rh*sin(TORAD(72.));  xus[ 5] = rg;
  xus[ 6] = rh*cos(TORAD(144.)); xus[ 7] = rh*sin(TORAD(144.)); xus[ 8] = rg;
  xus[ 9] = rh*cos(TORAD(216.)); xus[10] = rh*sin(TORAD(216.)); xus[11] = rg;
  xus[12] = rh*cos(TORAD(288.)); xus[13] = rh*sin(TORAD(288.)); xus[14] = rg;
  xus[15] = rh;                  xus[16] = 0;                   xus[17] = rg;
  xus[18] = rh*cos(TORAD(36.));  xus[19] = rh*sin(TORAD(36.));  xus[20] = -rg;
  xus[21] = rh*cos(TORAD(108.)); xus[22] = rh*sin(TORAD(108.)); xus[23] = -rg;
  xus[24] = -rh;                 xus[25] = 0;                   xus[26] = -rg;
  xus[27] = rh*cos(TORAD(252.)); xus[28] = rh*sin(TORAD(252.)); xus[29] = -rg;
  xus[30] = rh*cos(TORAD(324.)); xus[31] = rh*sin(TORAD(324.)); xus[32] = -rg;
  xus[33] = 0.;                  xus[34] = 0.;                  xus[35] = -1.;
}

static int ico_dot_arc(int densit) { 
  int i, j, k, tl, tl2, tn, tess;
  float a, d, x, y, z, x2, y2, z2, x3, y3, z3;
  float xij, yij, zij, xji, yji, zji, xik, yik, zik, xki, yki, zki,
    xjk, yjk, zjk, xkj, ykj, zkj;
  float* xus=NULL;

a= sqrt (16.0 );  
a = sqrt((((float) densit)-2.)/10.);
  tess = (int) ceil(a);
  n_dot = 10*tess*tess+2;
  snew(xus,3*n_dot);
  xpunsp = xus;
  icosaeder_vertices(xus);

  if (tess > 1) {
    tn = 12;
    a = rh*rh*2.*(1.-cos(TORAD(72.)));
    for (i=0; i<11; i++) {
      for (j=i+1; j<12; j++) {
        x = xus[3*i]-xus[3*j];
        y = xus[1+3*i]-xus[1+3*j]; z = xus[2+3*i]-xus[2+3*j];
        d = x*x+y*y+z*z;
        if (fabs(a-d) > DP_TOL) continue;
        for (tl=1; tl<tess; tl++) {
            divarc(xus[3*i], xus[1+3*i], xus[2+3*i],
                 xus[3*j], xus[1+3*j], xus[2+3*j],
            tl, tess, &xus[3*tn], &xus[1+3*tn], &xus[2+3*tn]);
          tn++;
          }
        }
      }
/* calculate tessalation of icosaeder faces */
    for (i=0; i<10; i++) {
      for (j=i+1; j<11; j++) {
        x = xus[3*i]-xus[3*j];
        y = xus[1+3*i]-xus[1+3*j]; z = xus[2+3*i]-xus[2+3*j];
        d = x*x+y*y+z*z;
        if (fabs(a-d) > DP_TOL) continue;

        for (k=j+1; k<12; k++) {
          x = xus[3*i]-xus[3*k];
          y = xus[1+3*i]-xus[1+3*k]; z = xus[2+3*i]-xus[2+3*k];
          d = x*x+y*y+z*z;
          if (fabs(a-d) > DP_TOL) continue;
          x = xus[3*j]-xus[3*k];
          y = xus[1+3*j]-xus[1+3*k]; z = xus[2+3*j]-xus[2+3*k];
          d = x*x+y*y+z*z;
          if (fabs(a-d) > DP_TOL) continue;
          for (tl=1; tl<tess-1; tl++) {
            divarc(xus[3*j], xus[1+3*j], xus[2+3*j],
                   xus[3*i], xus[1+3*i], xus[2+3*i],
              tl, tess, &xji, &yji, &zji);
            divarc(xus[3*k], xus[1+3*k], xus[2+3*k],
                   xus[3*i], xus[1+3*i], xus[2+3*i],
              tl, tess, &xki, &yki, &zki);

            for (tl2=1; tl2<tess-tl; tl2++) {
              divarc(xus[3*i], xus[1+3*i], xus[2+3*i],
                     xus[3*j], xus[1+3*j], xus[2+3*j],
                tl2, tess, &xij, &yij, &zij);
              divarc(xus[3*k], xus[1+3*k], xus[2+3*k],
                     xus[3*j], xus[1+3*j], xus[2+3*j],
                tl2, tess, &xkj, &ykj, &zkj);
              divarc(xus[3*i], xus[1+3*i], xus[2+3*i],
                     xus[3*k], xus[1+3*k], xus[2+3*k],
                tess-tl-tl2, tess, &xik, &yik, &zik);
              divarc(xus[3*j], xus[1+3*j], xus[2+3*j],
                     xus[3*k], xus[1+3*k], xus[2+3*k],
                tess-tl-tl2, tess, &xjk, &yjk, &zjk);
              divarc(xki, yki, zki, xji, yji, zji, tl2, tess-tl,
                &x, &y, &z);
              divarc(xkj, ykj, zkj, xij, yij, zij, tl, tess-tl2,
                &x2, &y2, &z2);
              divarc(xjk, yjk, zjk, xik, yik, zik, tl, tl+tl2,
                &x3, &y3, &z3);
              x = x+x2+x3; y = y+y2+y3; z = z+z2+z3;
              d = sqrt(x*x+y*y+z*z);
              xus[3*tn] = x/d;
              xus[1+3*tn] = y/d;
              xus[2+3*tn] = z/d;
              tn++;
              }		/* cycle tl2 */
            }		/* cycle tl */
          }		/* cycle k */
        }		/* cycle j */
      }			/* cycle i */
     }		/* end of if (tess > 1) */
  return n_dot;
}		/* end of routine ico_dot_arc */
static int ico_dot_dod(int densit) { 
  int i, j, k, tl, tl2, tn, tess, j1, j2;
  float a, d, x, y, z, x2, y2, z2, x3, y3, z3, ai_d, adod;
  float xij, yij, zij, xji, yji, zji, xik, yik, zik, xki, yki, zki,
    xjk, yjk, zjk, xkj, ykj, zkj;
  float* xus=NULL;
/* calculate tesselation level */
  a = sqrt((((float) densit)-2.)/30.);
  tess = MAXIMUM((int) ceil(a), 1);
  n_dot = 30*tess*tess+2;

  snew(xus,3*n_dot);
  xpunsp = xus;
  icosaeder_vertices(xus);

  tn=12;
/* square of the edge of an icosaeder */
    a = rh*rh*2.*(1.-cos(TORAD(72.)));
/* dodecaeder vertices */
  for (i=0; i<10; i++) {
    for (j=i+1; j<11; j++) {
      x = xus[3*i]-xus[3*j];
      y = xus[1+3*i]-xus[1+3*j]; z = xus[2+3*i]-xus[2+3*j];
      d = x*x+y*y+z*z;
      if (fabs(a-d) > DP_TOL) continue;
      for (k=j+1; k<12; k++) {
        x = xus[3*i]-xus[3*k];
        y = xus[1+3*i]-xus[1+3*k]; z = xus[2+3*i]-xus[2+3*k];
        d = x*x+y*y+z*z;
        if (fabs(a-d) > DP_TOL) continue;
        x = xus[3*j]-xus[3*k];
        y = xus[1+3*j]-xus[1+3*k]; z = xus[2+3*j]-xus[2+3*k];
        d = x*x+y*y+z*z;
        if (fabs(a-d) > DP_TOL) continue;
        x = xus[  3*i]+xus[  3*j]+xus[  3*k];
        y = xus[1+3*i]+xus[1+3*j]+xus[1+3*k];
        z = xus[2+3*i]+xus[2+3*j]+xus[2+3*k];
        d = sqrt(x*x+y*y+z*z);
        xus[3*tn]=x/d; xus[1+3*tn]=y/d; xus[2+3*tn]=z/d;
        tn++;
        }
      }
    }

  if (tess > 1) {
    tn = 32;
/* square of the edge of an dodecaeder */
    adod = 4.*(cos(TORAD(108.))-cos(TORAD(120.)))/(1.-cos(TORAD(120.)));
/* square of the distance of two adjacent vertices of ico- and dodecaeder */
    ai_d = 2.*(1.-sqrt(1.-a/3.));

/* calculate tessalation of mixed edges */
    for (i=0; i<31; i++) {
      j1 = 12; j2 = 32; a = ai_d;
      if (i>=12) { j1=i+1; a = adod; }
      for (j=j1; j<j2; j++) {
        x = xus[3*i]-xus[3*j];
        y = xus[1+3*i]-xus[1+3*j]; z = xus[2+3*i]-xus[2+3*j];
        d = x*x+y*y+z*z;
        if (fabs(a-d) > DP_TOL) continue;
        for (tl=1; tl<tess; tl++) {
          divarc(xus[3*i], xus[1+3*i], xus[2+3*i],
                 xus[3*j], xus[1+3*j], xus[2+3*j],
            tl, tess, &xus[3*tn], &xus[1+3*tn], &xus[2+3*tn]);
          tn++;
          }
        }
      }
/* calculate tessalation of pentakisdodecahedron faces */
    for (i=0; i<12; i++) {
      for (j=12; j<31; j++) {
        x = xus[3*i]-xus[3*j];
        y = xus[1+3*i]-xus[1+3*j]; z = xus[2+3*i]-xus[2+3*j];
        d = x*x+y*y+z*z;
        if (fabs(ai_d-d) > DP_TOL) continue;

        for (k=j+1; k<32; k++) {
          x = xus[3*i]-xus[3*k];
          y = xus[1+3*i]-xus[1+3*k]; z = xus[2+3*i]-xus[2+3*k];
          d = x*x+y*y+z*z;
          if (fabs(ai_d-d) > DP_TOL) continue;
          x = xus[3*j]-xus[3*k];
          y = xus[1+3*j]-xus[1+3*k]; z = xus[2+3*j]-xus[2+3*k];
          d = x*x+y*y+z*z;
          if (fabs(adod-d) > DP_TOL) continue;
          for (tl=1; tl<tess-1; tl++) {
            divarc(xus[3*j], xus[1+3*j], xus[2+3*j],
                   xus[3*i], xus[1+3*i], xus[2+3*i],
              tl, tess, &xji, &yji, &zji);
            divarc(xus[3*k], xus[1+3*k], xus[2+3*k],
                   xus[3*i], xus[1+3*i], xus[2+3*i],
              tl, tess, &xki, &yki, &zki);

            for (tl2=1; tl2<tess-tl; tl2++) {
              divarc(xus[3*i], xus[1+3*i], xus[2+3*i],
                     xus[3*j], xus[1+3*j], xus[2+3*j],
                tl2, tess, &xij, &yij, &zij);
              divarc(xus[3*k], xus[1+3*k], xus[2+3*k],
                     xus[3*j], xus[1+3*j], xus[2+3*j],
                tl2, tess, &xkj, &ykj, &zkj);
              divarc(xus[3*i], xus[1+3*i], xus[2+3*i],
                     xus[3*k], xus[1+3*k], xus[2+3*k],
                tess-tl-tl2, tess, &xik, &yik, &zik);
              divarc(xus[3*j], xus[1+3*j], xus[2+3*j],
                     xus[3*k], xus[1+3*k], xus[2+3*k],
                tess-tl-tl2, tess, &xjk, &yjk, &zjk);
               divarc(xki, yki, zki, xji, yji, zji, tl2, tess-tl,
                &x, &y, &z);
              divarc(xkj, ykj, zkj, xij, yij, zij, tl, tess-tl2,
                &x2, &y2, &z2);
              divarc(xjk, yjk, zjk, xik, yik, zik, tl, tl+tl2,
                &x3, &y3, &z3);
              x = x+x2+x3; y = y+y2+y3; z = z+z2+z3;
              d = sqrt(x*x+y*y+z*z);
              xus[3*tn] = x/d;
              xus[1+3*tn] = y/d;
              xus[2+3*tn] = z/d;
              tn++;
              }		/* cycle tl2 */
            }		/* cycle tl */
          }		/* cycle k */
        }		/* cycle j */
      }			/* cycle i */
     }		/* end of if (tess > 1) */
  return n_dot;
}		/* end of routine ico_dot_dod */
static int unsp_type(int densit) {
  int i1, i2;
  i1 = 1;
  while (10*i1*i1+2 < densit) i1++;
  i2 = 1;
  while (30*i2*i2+2 < densit) i2++;
  if (10*i1*i1-2 < 30*i2*i2-2) return UNSP_ICO_ARC;
  else return UNSP_ICO_DOD;
}

static int make_unsp(int densit, int mode, int * num_dot, int cubus) {
  int ndot, ico_cube_cb, i, j, k, l, ijk, tn, tl, tl2;
  float* xus;
  int*    work;
  float x, y, z;

  if (xpunsp) sfree(xpunsp); if (ico_wk) sfree(ico_wk);

  k=1; if (mode < 0) { k=0; mode = -mode; }
  if (mode == UNSP_ICO_ARC)      { ndot = ico_dot_arc(densit); }
  else if (mode == UNSP_ICO_DOD)      { ndot = ico_dot_dod(densit); }
  else  return 1;
  last_n_dot = ndot; last_densit = densit; last_unsp = mode;
  *num_dot=ndot; if (k) return 0;

/* in the following the dots of the unit sphere may be resorted */
  last_unsp = -last_unsp;

/* determine distribution of points in elementary cubes */
  if (cubus) {
    ico_cube = cubus;
    }
  else {
    last_cubus = 0;
    i=1;
    while (i*i*i*2 < ndot) i++;
    ico_cube = MAXIMUM(i-1, 0);
    }
  ico_cube_cb = ico_cube*ico_cube*ico_cube;
  del_cube=2./((float)ico_cube);
  snew(work,ndot);
  xus = xpunsp;
  for (l=0; l<ndot; l++) {
    i = MAXIMUM((int) floor((1.+xus[3*l])/del_cube), 0);
    if (i>=ico_cube) i = ico_cube-1;
    j = MAXIMUM((int) floor((1.+xus[1+3*l])/del_cube), 0);
    if (j>=ico_cube) j = ico_cube-1;
    k = MAXIMUM((int) floor((1.+xus[2+3*l])/del_cube), 0);
    if (k>=ico_cube) k = ico_cube-1;
    ijk = i+j*ico_cube+k*ico_cube*ico_cube;
    work[l] = ijk;
    }

  snew(ico_wk,2*ico_cube_cb+1);
  ico_pt = ico_wk+ico_cube_cb;
  for (l=0; l<ndot; l++) {
    ico_wk[work[l]]++;   /* dots per elementary cube */
    }

/* reordering of the coordinate array in accordance with box number */
  tn=0;
  for (i=0; i<ico_cube; i++) {
    for (j=0; j<ico_cube; j++) {
      for (k=0; k<ico_cube; k++) {
        tl=0;
        tl2 = tn;
        ijk = i+ico_cube*j+ico_cube*ico_cube*k;
        *(ico_pt+ijk) = tn;
        for (l=tl2; l<ndot; l++) {
          if (ijk == work[l]) {
            x = xus[3*l]; y = xus[1+3*l]; z = xus[2+3*l];
            xus[3*l] = xus[3*tn];
            xus[1+3*l] = xus[1+3*tn]; xus[2+3*l] = xus[2+3*tn];
            xus[3*tn] = x; xus[1+3*tn] = y; xus[2+3*tn] = z;
            ijk = work[l]; work[l]=work[tn]; work[tn]=ijk;
            tn++; tl++;
            }
          }
        *(ico_wk+ijk) = tl;
        }		/* cycle k */
      }			/* cycle j */
    }			/* cycle i */
  sfree(work); return 0;
}

typedef struct _stwknb {
  float x;
  float y;
  float z;
  float dot;
} Neighb;
		

float surfarea(int an, float ax[][3],float *rad,float *asurf,int densit, float solvR)
{
	int i, iii, i_ac,i_pa;
  int j;
  int distribution;
  int l;
  int nat, maxnei, nnei, last;
  int* wkdot=NULL;
  Neighb  *wknb, *ctnb;
  int j_at;
  float dx, dy, dz, dd, ai, aisq, ajsq, aj, as, a;
  float dotarea, area;
  float* xus;
  float ddx[3];
 
  float ra2max;

  distribution = unsp_type(densit);
  if (distribution != -last_unsp || last_cubus != 4 ||
      (densit != last_densit && densit != last_n_dot)) {
    if (make_unsp(densit, (-distribution), &n_dot, 4)) 
      return 1;
  }
  xus = xpunsp;

  dotarea = FOURPI/(float) n_dot;
  area = 0.;

  ra2max = rad[0]+solvR;
  for(i_pa=1;i_pa<an;i_pa++)
	  ra2max = rad[i_pa]+solvR;
  ra2max = 2.*ra2max;
  
  nat = an;
	maxnei = nat;
   snew(wknb, maxnei);
  snew(wkdot, n_dot*2);
  
   /* calculate surface for all atoms, step cube-wise */
  for(i_pa=0;i_pa<an;i_pa++)
{
    ai   = rad[i_pa]+solvR; 
    aisq = ai*ai;
    
    for (i=0; (i<n_dot); i++) 
      wkdot[i] = 0;
    
    ctnb = wknb; 
    nnei = 0;
    for (j=0; (j<nat); j++) {
      if (j == i_pa)
	continue;
      j_at = j;
      aj   = rad[j_at]+solvR; 
      ajsq = aj*aj;
      
      /* DvdS 11/02/02 To be modified for periodicity */
	  for(iii=0;iii<3;iii++)
		  ddx[iii]=ax[i_pa][iii]-ax[j_at][iii];
      dx = ddx[0], dy = ddx[1], dz = ddx[2];
      dd = dx*dx+dy*dy+dz*dz;
      
      as = ai+aj; 
      if (dd > as*as) 
	continue;
      nnei++;
      ctnb->x = dx; 
      ctnb->y = dy; 
      ctnb->z = dz;
      ctnb->dot = (dd+aisq-ajsq)/(2*ai); /* reference dot product */
      ctnb++;
    }

    /* check points on accessibility */
    i_ac = 0;
    if (nnei) {
      last = 0; 
      for (l=0; (l<n_dot); l++) {
	if (xus[3*l]*wknb[last].x+
	    xus[1+3*l]*wknb[last].y+
	    xus[2+3*l]*wknb[last].z <= wknb[last].dot) {
	  for (j=0; (j<nnei); j++) {
	    if ((xus[3*l]*wknb[j].x + xus[1+3*l]*wknb[j].y + xus[2+3*l]*wknb[j].z) > 
		wknb[j].dot) {
	      last = j; 
	      break;
	    }
	  }
	  if (j >= nnei) { 
	    i_ac++; 
	    wkdot[l] = 1; 
	  }
	}     /* end of cycle j */
      }       /* end of cycle l */
    }
    else {
      i_ac  = n_dot;
      for (l=0; l < n_dot; l++) 
	wkdot[l] = 1;
    }
    
    a = aisq*dotarea* (float) i_ac;
    area = area + a;
	if(asurf!=NULL)
    	asurf[i_pa] = a;
 }
  /*free(wkatm); */
  sfree(wkdot);
  sfree(wknb);
  
  return area;
}	
