/*                              SDOCK2.0
 *   Global Protein-Protein Docking Using Stepwise Force-Field Potentials   
 *                Written by ZHANG Changsheng 
 *  Copyright (c) 2024 Molecular Design Lab at Peking University, Beijing, China  
 
 *    contact: changshengzhang@pku.edu.cn, lhlai@pku.edu.cn
 *
 * GNM.c : GNM model for atom fluctuation calculation.
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "GNM.h"
#include "geometry.h"
#include "mem.h"
#include <math.h>


void genblocks(PROTEIN *pro)
{
	int blk;
	int ch,r,ia;
	int i;
	RES *res,*res0,*res1;
	CHN *chn;
	ATM *atm0,*atme;
	ATM *atm1,*atm2;
	int start,start0;
	int atend;
	int ai,aj;
      int ssc;
	int *ssb;

	blk=0;
	for(ch=0;ch<pro->chN;ch++){
		chn=pro->c+ch;
		for(r=0;r<chn->resN;r++){
			res=chn->r+r;
			/*printf("%d %s\n",r+1,res->name);*/
			if(res->ishet==SMOL){
				res->blockN=1;
				res->block_i[0]=blk++;
				res->block_at[0][0]=res->atmN;
				res->block_at[0][1]=(res->a-pro->a);
				continue;
			}
			atm0=res->a;
			atme=res->a+res->atmN-1;
			start=(res->a-pro->a);
			if(r!=0){
				res0=chn->r+r-1;
				start0=(res0->a-pro->a);
			}
			res1=chn->r+r+1;
			if(!strncmp(atme->abbrname,"OXT",3)) atend=res->atmN-1;
			else atend=res->atmN;
			if(r==0||res0->ishet==SMOL||distance2(res0->CA->xyz,res->CA->xyz)>BLDIS2){
				if(!strncmp(res->abbrname,"PRO",3)){
					res->block_i[0]=blk++;
					for(ia=0,i=1;ia<atend;ia++){
						if(ia==2||ia==3) continue;
						res->block_at[0][i++]=start+ia;
					}
					res->block_at[0][0]=i-1;
				}
				else{
					if(!strncmp(atm0->abbrname,"N  ",3)||!strncmp(atm0->abbrname,"NT ",3)){
						res->block_i[0]=blk++;
						res->block_at[0][0]=1;
						res->block_at[0][1]=start;
					}
				}
			}
			else{
				if(!strncmp(res->abbrname,"PRO",3)){
					res->block_i[0]=blk++;
					i=1;
					res->block_at[0][i++]=start0+2;
					res->block_at[0][i++]=start0+3;
					for(ia=0;ia<res->atmN;ia++){
						if(ia==2||ia==3) continue;
						res->block_at[0][i++]=start+ia;
					}
					res->block_at[0][0]=i-1;
				}
				else{
					if(!strncmp(atm0->abbrname,"N  ",3)||!strncmp(atm0->abbrname,"NT ",3)){
						res->block_i[0]=blk++;
						res->block_at[0][0]=3;
						res->block_at[0][1]=start0+2;
						res->block_at[0][2]=start0+3;
						res->block_at[0][3]=start;
					}
				}
			}
			res->blockN=1;
			if(!strncmp(res->abbrname,"GLY",3)||!strncmp(res->abbrname,"ALA",3)||!strncmp(res->abbrname,"SER",3)||!strncmp(res->abbrname,"CYS",3)||!strncmp(res->abbrname,"VAL",3)||!strncmp(res->abbrname,"THR",3)){
				for(ia=1,i=1;ia<atend;ia++){
					if(ia==2||ia==3) continue;
					res->block_at[1][i++]=start+ia;
				}
				if(i>1){
					res->blockN++;
					res->block_i[1]=blk++;
					res->block_at[1][0]=i-1;
				}
			}
			else if((!strncmp(res->abbrname,"LYS",3))||(!strncmp(res->abbrname,"ARG",3))){
				for(ia=1,i=1;ia<=5&&ia<atend;ia++){
					if(ia==2||ia==3) continue;
					res->block_at[1][i++]=start+ia;
				}
				if(i>1){
					res->blockN++;
					res->block_i[1]=blk++;
					res->block_at[1][0]=i-1;
				}
				for(i=1;ia<atend;ia++){
					res->block_at[2][i++]=start+ia;
				}
				if(i>1){
					res->blockN++;
					res->block_i[2]=blk++;
					res->block_at[2][0]=i-1;
				}
			}
			else if(strncmp(res->abbrname,"PRO",3)){
				for(ia=1,i=1;ia<=4&&ia<atend;ia++){
					if(ia==2||ia==3) continue;
					res->block_at[1][i++]=start+ia;
				}
				if(i>1){
					res->blockN++;
					res->block_i[1]=blk++;
					res->block_at[1][0]=i-1;
				}
				for(i=1;ia<atend;ia++){
					res->block_at[2][i++]=start+ia;
				}
				if(i>1){
					res->blockN++;
					res->block_i[2]=blk++;
					res->block_at[2][0]=i-1;
				}
			}
			if(r==chn->resN-1||res1->ishet==SMOL||distance2(res->CA->xyz,res1->CA->xyz)>BLDIS2){
				atm1=res->a+2;
				atm2=res->a+3;
				if((!strncmp(atm1->abbrname,"C  ",3)||!strncmp(atm1->abbrname,"CT ",3))&&(!strncmp(atm2->abbrname,"O  ",3)||!strncmp(atm2->abbrname,"OT ",3))){
					res->block_i[res->blockN]=blk++;
					res->block_at[res->blockN][0]=2;
					res->block_at[res->blockN][1]=start+2;
					res->block_at[res->blockN][2]=start+3;
					if(!strncmp(atme->abbrname,"OXT",3)){
						res->block_at[res->blockN][0]++;
						res->block_at[res->blockN][3]=start+res->atmN-1;
					}
					res->blockN++;
				}
			}
		}
	}
	pro->blockN=blk;
	/*printf("%d\n",blk);*/
	ssc=0;
      snew(ssb, pro->resN/4);
 	for(ai=0;ai<pro->atmN;ai++){
		for(aj=ai+1;aj<pro->atmN;aj++){
			if(!strncmp(pro->a[ai].abbrname,"SG ",3)&&!strncmp(pro->a[aj].abbrname,"SG ",3)){
				if(distance2(pro->a[ai].xyz,pro->a[aj].xyz)<2.2*2.2){
					ssb[ssc*2]=ai;
					ssb[ssc*2+1]=aj;
					ssc++;
				}
			}
		}
	}
	pro->ssbN=ssc;
	snew(pro->ssb,ssc*2);
	for(i=0;i<ssc*2;i++) pro->ssb[i]=ssb[i];
	sfree(ssb);
}

#define CONTACTDIS2 4.0*4.0
int iscontact(int *bat1, int *bat2,int het1, int het2, PROTEIN *pro)
{
	int ai,aj;
	int a1,a2;
	for(ai=0;ai<bat1[0];ai++){
		if(het1==SMOL) a1=bat1[1]+ai;
		else a1=bat1[ai+1];
		for(aj=0;aj<bat2[0];aj++){
			if(het2==SMOL) a2=bat2[1]+aj;
			else a2=bat2[aj+1];
			if(distance2(pro->a[a1].xyz,pro->a[a2].xyz)<CONTACTDIS2)
				return 1;
		}
	}
	return 0;
}
int getlink5(PROTEIN *pro, int FRAGLINK[][2])
{
	int linki=0;
	int ir,jr;
	int ifr,jfr;

	for(ir=0;ir<pro->resN;ir++){
		for(ifr=0;ifr<pro->r[ir].blockN;ifr++){
			for(jr=ir;jr<pro->resN;jr++){
				for(jfr=0;jfr<pro->r[jr].blockN;jfr++){
					if(ir==jr&&ifr>=jfr)
						continue;
					if(iscontact(pro->r[ir].block_at[ifr],pro->r[jr].block_at[jfr],pro->r[ir].ishet, pro->r[jr].ishet,pro)==1){
						FRAGLINK[linki][0]=pro->r[ir].block_i[ifr];
						FRAGLINK[linki][1]=pro->r[jr].block_i[jfr];
						/*printf("  link5 %d %d %d\n", linki+1,FRAGLINK[linki][0],FRAGLINK[linki][1]);*/
						linki++;
					}
				}
			}
		}
	}
	return linki;
}

int getlink6(PROTEIN *pro, int BONDLINK[][2])
{
	int linki=0;
	int ch,r;
	RES *res,*res1;
	CHN *chn;
	for(ch=0;ch<pro->chN;ch++){
		chn=pro->c+ch;
		for(r=0;r<chn->resN;r++){
			res=chn->r+r;
			if(res->blockN>1){
				BONDLINK[linki][0]=res->block_i[0];
				BONDLINK[linki][1]=res->block_i[1];
				linki++;
			}
			if(r<chn->resN-1){
				res1=chn->r+r+1;
				if(res->blockN>1){
					BONDLINK[linki][0]=res->block_i[1];
				}
				else BONDLINK[linki][0]=res->block_i[0];
				BONDLINK[linki][1]=res1->block_i[0];
				linki++;
			}
			if(r==chn->resN-1&&res->block_at[res->blockN-1][1]==(res->a-pro->a)+2){
				if(res->blockN>2)
					BONDLINK[linki][0]=res->block_i[1];
				else
					BONDLINK[linki][0]=res->block_i[0];
				BONDLINK[linki][1]=res->block_i[res->blockN-1];
				linki++;
			}
		}
	}
	return linki;
}
int getlink7(PROTEIN *pro, int BONDLINK[][2])
{
	int linki=0;
	int ir;
	int ib;
	int iss;
	int ss1,ss2;
	RES *r1,*r2;

	for(ir=0;ir<pro->resN;ir++){
		for(ib=2;ib<pro->r[ir].blockN;ib++){
			BONDLINK[linki][0]=pro->r[ir].block_i[ib-1];
			BONDLINK[linki][1]=pro->r[ir].block_i[ib];
			linki++;
		}
	}
	/*printf("why %d %d\n",linki,pro->ssbN);*/
	for(iss=0;iss<pro->ssbN;iss++){
		ss1=pro->ssb[iss*2];
		ss2=pro->ssb[iss*2+1];
		/*printf("%d ssb %d %d %d\n",linki,iss+1,ss1,ss2);*/
		r1=pro->a[ss1].atres;
		r2=pro->a[ss2].atres;
		BONDLINK[linki][0]=r1->block_i[1];
		BONDLINK[linki][1]=r2->block_i[1];
		linki++;
	}
	sfree(pro->ssb);
	return linki;
}
#define BONDLINKK3 -12.0
#define BONDLINKK4 -12.0
void getmatrixvalue(int NN, int ln1, int FRAGLINK[][2], int ln2, int BONDLINK1[][2],int ln3, int BONDLINK2[][2],double *FRAGMATRIX)
{
	int i,j;
	int il;
	double n;

	for(i=0;i<NN;i++)
		for(j=0;j<NN;j++)
			FRAGMATRIX[i*NN+j]=0.0;
	for(il=0;il<ln1;il++){
		i=FRAGLINK[il][0];
		j=FRAGLINK[il][1];
		FRAGMATRIX[i*NN+j]=-1.0;
		FRAGMATRIX[j*NN+i]=-1.0;
	}
	for(il=0;il<ln2;il++){
		i=BONDLINK1[il][0];
		j=BONDLINK1[il][1];
		FRAGMATRIX[i*NN+j]=BONDLINKK3;
		FRAGMATRIX[j*NN+i]=BONDLINKK3;
	}
	for(il=0;il<ln3;il++){
		i=BONDLINK2[il][0];
		j=BONDLINK2[il][1];
		FRAGMATRIX[i*NN+j]=BONDLINKK4;
		FRAGMATRIX[j*NN+i]=BONDLINKK4;
	}
	for(i=0;i<NN;i++){
		n=0.0;
		for(j=0;j<NN;j++){
			n-=FRAGMATRIX[i*NN+j];
		}
		FRAGMATRIX[i*NN+i]=n;
		/*printf("%d %8.3f\n",i+1,n);*/
	}
}

void genmatrix(PROTEIN *pro, double *matrix)
{
	int ln1,ln2,ln3;
	int (*FRAGLINK)[2];
	int (*BONDLINK1)[2];
	int (*BONDLINK2)[2];

	snew(FRAGLINK,pro->blockN*5);
	ln1=getlink5(pro, FRAGLINK);
	snew(BONDLINK1,pro->blockN*2);
	ln2=getlink6(pro, BONDLINK1);
	snew(BONDLINK2,pro->blockN*2+pro->ssbN);
	ln3=getlink7(pro, BONDLINK2);
	getmatrixvalue(pro->blockN,ln1,FRAGLINK,ln2,BONDLINK1,ln3,BONDLINK2,matrix);
	sfree(FRAGLINK);
	sfree(BONDLINK1);
	sfree(BONDLINK2);
}

void strq(double a[],int n,double q[],double b[],double c[])
{ int i,j,k,u;
  double h,f,g,h2;
  for (i=0; i<=n-1; i++)
  for (j=0; j<=n-1; j++)
    { u=i*n+j; q[u]=a[u]; /*if (i==j) printf("%d  %8.3f\n",u,a[u]);*/}
  for (i=n-1; i>=1; i--)
    { h=0.0;
      if (i>1)
        for (k=0; k<=i-1; k++)
          { u=i*n+k; h=h+q[u]*q[u];}
      if (h+1.0==1.0)
        { c[i]=0.0;
          if (i==1) c[i]=q[i*n+i-1];
          b[i]=0.0;
        }
      else
        { c[i]=sqrt(h);
          u=i*n+i-1;
          if (q[u]>0.0) c[i]=-c[i];
          h=h-q[u]*c[i];
          q[u]=q[u]-c[i];
          f=0.0;
          for (j=0; j<=i-1; j++)
            { q[j*n+i]=q[i*n+j]/h;
              g=0.0;
              for (k=0; k<=j; k++)
                g=g+q[j*n+k]*q[i*n+k];
              if (j+1<=i-1)
                for (k=j+1; k<=i-1; k++)
                  g=g+q[k*n+j]*q[i*n+k];
              c[j]=g/h;
              f=f+g*q[j*n+i];
            }
          h2=f/(h+h);
          for (j=0; j<=i-1; j++)
            { f=q[i*n+j];
              g=c[j]-h2*f;
              c[j]=g;
              for (k=0; k<=j; k++)
                { u=j*n+k;
                  q[u]=q[u]-f*c[k]-g*q[i*n+k];
                }
            }
          b[i]=h;
        }
    }
  for (i=0; i<=n-2; i++) c[i]=c[i+1];
  c[n-1]=0.0;
  b[0]=0.0;
  for (i=0; i<=n-1; i++)
    { if ((b[i]!=0.0)&&(i-1>=0))
        for (j=0; j<=i-1; j++)
          { g=0.0;
            for (k=0; k<=i-1; k++)
              g=g+q[i*n+k]*q[k*n+j];
            for (k=0; k<=i-1; k++)
              { u=k*n+j;
                q[u]=q[u]-g*q[k*n+i];
              }
          }
      u=i*n+i;
      b[i]=q[u]; q[u]=1.0;
      if (i-1>=0)
        for (j=0; j<=i-1; j++)
          { q[i*n+j]=0.0; q[j*n+i]=0.0;}
    }
  return;
}

int sstq(int n,double b[],double c[],double q[],double eps,int l)
{ int i,j,k,m,it,u,v;
  double d,f,h,g,p,r,e,s;
  c[n-1]=0.0; d=0.0; f=0.0;
  for (j=0; j<=n-1; j++)
    { it=0;
      h=eps*(fabs(b[j])+fabs(c[j]));
      if (h>d) d=h;
      m=j;
      while ((m<=n-1)&&(fabs(c[m])>d)) m=m+1;
      if (m!=j)
        { do
            { if (it==l)
                { printf("fail\n");
                  return(-1);
                }
              it=it+1;
              g=b[j];
              p=(b[j+1]-g)/(2.0*c[j]);
              r=sqrt(p*p+1.0);
              if (p>=0.0) b[j]=c[j]/(p+r);
              else b[j]=c[j]/(p-r);
              h=g-b[j];
              for (i=j+1; i<=n-1; i++)
                b[i]=b[i]-h;
              f=f+h; p=b[m]; e=1.0; s=0.0;
              for (i=m-1; i>=j; i--)
                { g=e*c[i]; h=e*p;
                  if (fabs(p)>=fabs(c[i]))
                    { e=c[i]/p; r=sqrt(e*e+1.0);
                      c[i+1]=s*p*r; s=e/r; e=1.0/r;
                    }
                  else
		      { e=p/c[i]; r=sqrt(e*e+1.0);
                      c[i+1]=s*c[i]*r;
                      s=1.0/r; e=e/r;
                    }
                  p=e*b[i]-s*g;
                  b[i+1]=h+s*(e*g+s*b[i]);
                  for (k=0; k<=n-1; k++)
                    { u=k*n+i+1; v=u-1;
                      h=q[u]; q[u]=s*q[v]+e*h;
                      q[v]=e*q[v]-s*h;
                    }
                }
              c[j]=s*p; b[j]=e*p;
            }
          while (fabs(c[j])>d);
        }
      b[j]=b[j]+f;
    }
  for (i=0; i<=n-1; i++)
    { k=i; p=b[i];
      if (i+1<=n-1)
        { j=i+1;
          while ((j<=n-1)&&(b[j]<=p))
            { k=j; p=b[j]; j=j+1;}
        }
      if (k!=i)
        { b[k]=b[i]; b[i]=p;
          for (j=0; j<=n-1; j++)
            { u=j*n+i; v=j*n+k;
              p=q[u]; q[u]=q[v]; q[v]=p;
            }
        }
    }
  return(1);
}

int eig_sstq(double a[],int n,double v[],double eps,int jt)
{
  int k;
	int i;
	double *b,*c;
	snew(b,n);
	snew(c,n);
	strq(a,n,v,b,c);
	k=sstq(n,b,c,v,eps,jt);
	for(i=0;i<n*n;i++){
		a[i]=0.0;
	}
	for(i=0;i<n;i++){
		a[i*n+i]=b[i];
	}
	sfree(b);
	sfree(c);
	return k;
}
void swap(double *a,double *b)
{
	double t;
	t=*a;*a=*b;*b=t;
}
int sorteigen(double a[],int n,double v[])
{
	int i,j,k,l;

	for(i=0;i<n-1;i++){
		k=i;
		for(j=i+1;j<n;j++){
			if(a[k*n+k]>a[j*n+j])
				k=j;
		}
		if(k!=i){
			swap(&(a[k*n+k]),&(a[i*n+i]));
			for(l=0;l<n;l++)
				swap(&(v[l*n+i]),&(v[l*n+k]));
		}
	}
	/*for(i=0;i<n-1;i++) printf("%d  %8.3f\n",i+1,a[i*n+i]);*/

	return 1;
}
void fluctuation(double *value,double *vector, PROTEIN *pro)
{
	double fluc;
	int ir,fr;
	int ifr,jfr;
	int ia;
	int N=pro->blockN;

	for(ia=0;ia<pro->atmN;ia++)
		pro->a[ia].fluc=1.0;
	for(ir=0;ir<pro->resN;ir++){
		for(ifr=0;ifr<pro->r[ir].blockN;ifr++){
			fr=pro->r[ir].block_i[ifr];
			fluc=0.0;
			for(jfr=1;jfr<N;jfr++){
				if(value[jfr*N+jfr]>5e-3)
					fluc=fluc+vector[fr*N+jfr]*vector[fr*N+jfr]/value[jfr*N+jfr];
			}
			/*printf("%d  %8.3f\n",fr+1,fluc);*/
			for(ia=0;ia<pro->r[ir].block_at[ifr][0];ia++){
				if(pro->r[ir].ishet==SMOL) pro->a[pro->r[ir].block_at[ifr][1]+ia].fluc=fluc;
				else pro->a[pro->r[ir].block_at[ifr][ia+1]].fluc=fluc;
			}
		}
	}
}

#define EPS 1e-15
#define MAXITERATE	 2000000

void calfluctuation(double *matrix, PROTEIN *pro)
{
	double *vector;
	snew(vector, pro->blockN*pro->blockN);
	eig_sstq(matrix,pro->blockN,vector,EPS,MAXITERATE);
	sorteigen(matrix,pro->blockN,vector);
	fluctuation(matrix,vector,pro);
	sfree(vector);
}
