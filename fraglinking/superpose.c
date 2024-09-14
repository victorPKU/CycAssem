#include "superpose.h"
#include "mem.h"
#include <math.h>
#include <stdio.h>
#include "geometry.h"

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
  a[k][l]=h+s*(g-h*tau);
	
void jacobi(double **a,int n,double d[],double **v,int *nrot)
{
  int j,i;
  int iq,ip;
  double tresh,theta,tau,t,sm,s,h,g,c,*b,*z;

  snew(b,n);
  snew(z,n);
  for (ip=0; ip<n; ip++) {
    for (iq=0; iq<n; iq++) v[ip][iq]=0.0;
    v[ip][ip]=1.0;
 }
  for (ip=0; ip<n;ip++) {
    b[ip]=d[ip]=a[ip][ip];
    z[ip]=0.0;
  }
  *nrot=0;
  for (i=1; i<=50; i++) {
    sm=0.0;
    for (ip=0; ip<n-1; ip++) {
      for (iq=ip+1; iq<n; iq++)
        sm += fabs(a[ip][iq]);
    }
    if (sm == 0.0) {
  		sfree(b);
  		sfree(z);
      return;
    }
    if (i < 4)
      tresh=0.2*sm/(n*n);
    else
      tresh=0.0;
    for (ip=0; ip<n-1; ip++) {
      for (iq=ip+1; iq<n; iq++) {
        g=100.0*fabs(a[ip][iq]);
        if (i > 4 && fabs(d[ip])+g == fabs(d[ip])
            && fabs(d[iq])+g == fabs(d[iq]))
          a[ip][iq]=0.0;
        else if (fabs(a[ip][iq]) > tresh) {
          h=d[iq]-d[ip];
          if (fabs(h)+g == fabs(h))
            t=(a[ip][iq])/h;
          else {
            theta=0.5*h/(a[ip][iq]);
            t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
            if (theta < 0.0) t = -t;
          }
          c=1.0/sqrt(1+t*t);
          s=t*c;
          tau=s/(1.0+c);
          h=t*a[ip][iq];
          z[ip] -= h;
          z[iq] += h;
          d[ip] -= h;
          d[iq] += h;
          a[ip][iq]=0.0;
          for (j=0; j<ip; j++) {
            ROTATE(a,j,ip,j,iq)
	  }
          for (j=ip+1; j<iq; j++) {
            ROTATE(a,ip,j,j,iq)
            }
          for (j=iq+1; j<n; j++) {
            ROTATE(a,ip,j,iq,j)
            }
          for (j=0; j<n; j++) {
            ROTATE(v,j,ip,j,iq)
            }
          ++(*nrot);
        }
      }
    }
    for (ip=0; ip<n; ip++) {
      b[ip] +=  z[ip];
      d[ip]  =  b[ip];
      z[ip]  =  0.0;
    }
  }
  printf("ERROR:  Too many iterations in routine JACOBI in superposition!\n");
}

void calc_fit_R(int an,float xp[][3],float x[][3],float R[3][3])
{
  int    c,r,n,j,m,i,irot;
  double **omega,**om;
  double d[6],xnr,xpc;
  float vh[3][3],vk[3][3],u[3][3];
  float   mn;
  int    index;
  float   max_d;

  snew(omega,6);
  snew(om,6);
  for(i=0; i<6; i++) {
    snew(omega[i],6);
    snew(om[i],6);
  }
  
  for(i=0; i<6; i++) {
    d[i]=0;
    for(j=0; j<6; j++) {
      omega[i][j]=0;
      om[i][j]=0;
    }
  }
  
  /*calculate the matrix U*/
  for(i=0;i<3;i++)
	  for(j=0;j<3;j++)
		  u[i][j]=0;
  for(n=0;(n<an);n++)
      for(c=0; (c<3); c++) {
	xpc=xp[n][c];
	for(r=0; (r<3); r++) {
	  xnr=x[n][r];
	  u[c][r]+=xnr*xpc;
	}
      }
  
  /*construct omega*/
  /*omega is symmetric -> omega==omega' */
  for(r=0; r<6; r++)
    for(c=0; c<=r; c++)
      if (r>=3 && c<3) {
        omega[r][c]=u[r-3][c];
        omega[c][r]=u[r-3][c];
      } else {
        omega[r][c]=0;
        omega[c][r]=0;
      }
  
  /*determine h and k*/
  jacobi(omega,6,d,om,&irot);
  /*real   **omega = input matrix a[0..n-1][0..n-1] must be symmetric
   *int     natoms = number of rows and columns
   *real      NULL = d[0]..d[n-1] are the eigenvalues of a[][]
   *real       **v = v[0..n-1][0..n-1] contains the vectors in columns
   *int      *irot = number of jacobi rotations
   */
  
  
  index=0; /* For the compiler only */

  /* Copy only the first two eigenvectors */  
  for(j=0; j<2; j++) {
    max_d=-1000;
    for(i=0; i<6; i++)
      if (d[i]>max_d) {
        max_d=d[i];
        index=i;
      }
    d[index]=-10000;
    for(i=0; i<3; i++) {
      vh[j][i]=M_SQRT2*om[i][index];
      vk[j][i]=M_SQRT2*om[i+3][index];
    }
  }
  /* Calculate the last eigenvector as the outer-product of the first two.
   * This insures that the conformation is not mirrored and
   * prevents problems with completely flat reference structures.
   */  
  oprod(vh[0],vh[1],vh[2]);
  oprod(vk[0],vk[1],vk[2]);

  /*determine R*/
  for(r=0; r<3; r++)
    for(c=0; c<3; c++)
      R[r][c] = vk[0][r]*vh[0][c] +
	        vk[1][r]*vh[1][c] +
	        vk[2][r]*vh[2][c];

  for(i=0; i<6; i++) {
    sfree(omega[i]);
    sfree(om[i]);
  }
  sfree(omega);
  sfree(om);
}
void do_rot(int an,float x[][3],float xt[][3],float R[3][3])
{
  int    i,j,m,r,c;

  for(j=0; j<an; j++) {
    for(r=0; r<3; r++) {
      xt[j][r]=0.0;
      for(c=0; c<3; c++)
        xt[j][r]+=R[r][c]*x[j][c];
    }
  }
}

float calcrmsd(int an,float xp[][3],float x[][3])
{
  int i,d;
  float xd, rd;
  
  rd=0;
  for(i=0; i<an; i++) {
    for(d=0 ; d<3; d++) {
      xd = x[i][d] - xp[i][d];
      rd += xd*xd;
      }
  }
   return sqrt(rd/an);
}
void printterminalxyz2(float terxyz[4][3])
{
	int i;

	for(i=0;i<4;i++)  printf("%8.3f%8.3f%8.3f\n",terxyz[i][0],terxyz[i][1],terxyz[i][2]);
	printf("======================\n");
}

float fit_terminal(float lptatx[4][3],float tptermxyz[4][3],float trans[3], float rotR[3][3])
{
	float tptatx[4][3];
	float rmsd;

	get_center(4, tptermxyz, trans);
	irtranslate(4,tptermxyz ,tptatx,trans);
	calc_fit_R(4,lptatx,tptatx,rotR);
	do_rot(4,tptatx,tptermxyz,rotR);
	/*printterminalxyz2(tptermxyz);*/
	rmsd=calcrmsd(4,lptatx,tptermxyz);
	/*printf("rmsd= %8.3f\n",rmsd);*/
	
	return rmsd;
}

