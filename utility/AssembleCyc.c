#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct{
	int start;
	int N;	
	int CA;
	int C;
	int O;
	int sstart;
	int end;
	char nam[10];
}RES;
typedef struct{
	int an;
	float xyz[100][3];
	char nam[100][4];
	int resn;
	RES res[5];
}FRAGMENT;

typedef struct{
	float N[3];	
	float CA[3];
	float C[3];
	float O[3];
	float H[3];
	char resnam[10];
	int satn;
	float side[10][3];
	char nam[10][4];
}CYCRES;

#define MAXHBN  15
typedef struct{
	char linkname[2][100];
	int typ[2];
	float linkscore[2];
	int resn;
	int len[4];
	int hbn;
	int index[MAXHBN][2];
 	float hbdis[MAXHBN];
	CYCRES res[15];
}CYCPEP;

/*****************GEOMETRY*************************/
void getdx(float x1[3],float x2[3],float dx[3])
{
	dx[0]=x1[0]-x2[0];
	dx[1]=x1[1]-x2[1];
	dx[2]=x1[2]-x2[2];
}
float get_len(float x[3])
{
	float r;

	r=x[0]*x[0]+x[1]*x[1]+x[2]*x[2];

	return sqrt(r);
}

void direction( float a[3], float b[3], float c[3])
{
	int i;
	float r;
	
	getdx(a,b,c);
	r=get_len(c);
	for(i=0;i<3;i++) 
		c[i]=c[i]/r;
}

#define PI 3.14159265358979323846
#define TODEG(A)     ((A)*57.295779513)

float iprod(float a[3],float b[3])
{
	return (a[0]*b[0]+a[1]*b[1]+a[2]*b[2]);
}
float cos_angle(const float a[3],const float b[3])
{
  float   cos;
  int    m;
  double aa,bb,ip,ipa,ipb; 
  
  ip=ipa=ipb=0.0;
  for(m=0; m<3; m++) {
    aa   = a[m];
    bb   = b[m];
    ip  += aa*bb;
    ipa += aa*aa;
    ipb += bb*bb;
  }
  cos=ip/sqrt(ipa*ipb);
  if (cos > 1.0) 
    return  1.0; 
  if (cos <-1.0) 
    return -1.0;
  
  return cos;
}
float cal_angle(float x1[3], float x2[3], float x3[3])
{
	float x12[3], x32[3];
	float cosa;
	float a;

	getdx(x1,x2,x12);
	getdx(x3,x2,x32);
	cosa=cos_angle(x12,x32);
	a=acos(cosa);

	return TODEG(a);
}
void oprod(const float a[3],const float b[3],float c[3])
{
  c[0]=a[1]*b[2]-a[2]*b[1];
  c[1]=a[2]*b[0]-a[0]*b[2];
  c[2]=a[0]*b[1]-a[1]*b[0];
}

float cal_dih(float x1[3], float x2[3], float x3[3], float x4[3])
{
	float x12[3],x32[3],x34[3], m[3], n[3];
	float ipr,phi,cos_phi,sign;
	int i;

	getdx(x1,x2,x12);  
	getdx(x3,x2,x32);	
	getdx(x3,x4,x34);	

	oprod(x12,x32,m); 
	oprod(x32,x34,n);
  	cos_phi=cos_angle(m,n);
  	phi=acos(cos_phi);
	ipr=iprod(x12,n);
	sign=(ipr<0.0)?-1.0:1.0;
	phi=sign*phi; 

  	return TODEG(phi);
}

#define HBONDLEN 0.98
void calHNxyz(float H[3], float N[3], float CA[3], float C[3])
{
	float v1[3],v2[3],v3[3];
	int i;
	float len;

	direction(CA, N, v1);
     	direction(C, N, v2);
	for(i=0;i<3;i++){
		v3[i]=v1[i]+v2[i];
	}
	len=get_len(v3);
	for(i=0;i<3;i++){
		v3[i]=v3[i]/len;
	}
	for(i=0;i<3;i++) H[i]=N[i]-v3[i]*HBONDLEN;
}
void addHN(CYCPEP *p)
{
	int ir;

	calHNxyz(p->res[0].H, p->res[0].N, p->res[0].CA, p->res[p->resn-1].C);
	for(ir=1;ir<p->resn;ir++){
		calHNxyz(p->res[ir].H, p->res[ir].N, p->res[ir].CA, p->res[ir-1].C);
	}
}

float distance2(float x1[3], float x2[3])
{
        return (x1[0]-x2[0])*(x1[0]-x2[0])+(x1[1]-x2[1])*(x1[1]-x2[1])+(x1[2]-x2[2])*(x1[2]-x2[2]);
}


/***********************HYDROGEN BONDING*******************************************/
#define HBDISMAX  (2.5*2.5)
#define HBDISMIN  (1.6*1.6)

int dispair(CYCPEP *p, int pair[][2], float *hbdis)
{
	int ir,jr;
	float dis;
	int n=0;
	
	for(ir=0;ir<p->resn;ir++){
		for(jr=0;jr<p->resn;jr++){
			if(jr==ir) continue;
			dis=distance2(p->res[ir].H, p->res[jr].O);
			if(dis<HBDISMAX&&dis>HBDISMIN){
				pair[n][0]=ir;
				pair[n][1]=jr;
				hbdis[n]=sqrt(dis);
				n++;
			}
		}
	}

	return n;
}
int checkangle(CYCPEP *p, int n, int pair[][2], int *index)
{
	float ang1,ang2;
	float dih;
	int i;
	int ir,jr;
	int c=0;

	for(i=0;i<n;i++){
		ir=pair[i][0];
		jr=pair[i][1];
		ang1=cal_angle(p->res[ir].N, p->res[ir].H, p->res[jr].O);
		/*printf("ang1 %8.3f", ang1);*/
		if(ang1<90) continue;
		ang2=cal_angle(p->res[ir].H, p->res[jr].O, p->res[jr].C);
		/*printf("ang2 %8.3f", ang2);*/
		if(ang2<90) continue;
		dih=cal_dih(p->res[ir].H, p->res[jr].O, p->res[jr].C, p->res[jr].CA);
		/*printf("dih %8.3f", dih);*/
		if(dih>50&&dih<130) continue;
		index[c]=i;
		c++;
	}
		/*printf("\n");*/
	return c;	
}
/*************************************************************************/

void getxyz(char line[], float xyz[3])
{
	int i;
	for(i=0;i<3;i++)
		xyz[i]=atof(line+30+8*i);
}

#define LINELEN 81
void readfrag(char *pn, FRAGMENT *p)
{
	FILE *pf;
	char line[LINELEN];
	int n=0;
	int rn=0;
	
	if((pf=fopen(pn,"r"))==NULL){
		printf("ERROR: Can not open fragment structure file %s\n",pn);
		exit(0);
	}
	
	while(fgets(line, LINELEN,pf)){
		if(!strncmp(line,"ATOM",4)){
			if(rn==0||strncmp(line+17,p->res[rn-1].nam,9)){
				strncpy(p->res[rn].nam,line+17,9);
				p->res[rn].nam[9]='\0';
				p->res[rn].start=n;
				p->res[rn].sstart=n+4;
				if(rn!=0) p->res[rn-1].end=n-1;
				rn++;
			}
			getxyz(line,p->xyz[n]);
			strncpy(p->nam[n],line+13,3);
			p->nam[n][3]='\0';
			if(!strncmp(line+13,"N  ",3)){
				p->res[rn-1].N=n;
			}
			else if(!strncmp(line+13,"CA ",3)){
				p->res[rn-1].CA=n;
			}
			else if(!strncmp(line+13,"C  ",3)){
				p->res[rn-1].C=n;
			}
			else if(!strncmp(line+13,"O  ",3)){
				p->res[rn-1].O=n;
			}
			else if(!strncmp(line+13,"CB ",3)){
				p->res[rn-1].sstart=n;
			}
			n++;
		}
	}
	fclose(pf);
	p->an=n;
	p->resn=rn;
	p->res[rn-1].end=n-1;
}
void getfname(char *line, char *linkname)
{
	int i,j;
	i=14;j=0;

	while(line[i]!='\n'&&line[i]!=' '){
		linkname[j]=line[i];
		i++;j++;
	}
	linkname[j]='\0';
}
#define LINELEN 81
#define BESTCLUSTERNUM 3
void readlink(char *pn, FRAGMENT p[2][BESTCLUSTERNUM], int typ[2][BESTCLUSTERNUM], char fragname[2][100], char linkname[2][BESTCLUSTERNUM][100], float fragscore[2], float linkscore[2][BESTCLUSTERNUM], int tertyp[2][BESTCLUSTERNUM][2],int clusterN[2])
{
	FILE *pf;
	char line[LINELEN];
	int n=0;
	int rn=0;
	int id;
	int lid=0;
	int hasgetfragname1=0;
	int hasgetfragname2=0;
	int i,j;
	
	if((pf=fopen(pn,"r"))==NULL){
		printf("ERROR: Can not open protein linking file %s\n",pn);
		exit(0);
	}
	while(fgets(line, LINELEN,pf)){
		if(!strncmp(line,"TITLE",5)){
			id=atoi(line+8)-1;
			lid=atoi(line+11)-1;
			typ[id][lid]==atoi(line+13);
			clusterN[id]=lid+1;
		}
		else if(!strncmp(line,"REMARK  dock1",13)){
			if(hasgetfragname1==0){
				i=26;j=0;
				while(line[i]!='\n'&&line[i]!=' '){
					fragname[0][j]=line[i];
					i++;
					j++;
				}
				fragname[0][j]='\0';
				fragscore[0]=atof(line+16);
				hasgetfragname1=1;
			}
			tertyp[id][lid][0]=atoi(line+13);
		}
		else if(!strncmp(line,"REMARK  dock2",13)){
			if(hasgetfragname2==0){
				i=26;j=0;
				while(line[i]!='\n'&&line[i]!=' '){
					fragname[1][j]=line[i];
					i++;
					j++;
				}
				fragname[1][j]='\0';
				fragscore[1]=atof(line+16);
				hasgetfragname2=1;
			}
			tertyp[id][lid][1]=atoi(line+13);
		}
		else if(!strncmp(line,"REMARK  link",12)){
			getfname(line,linkname[id][lid]);
		}
		else if(!strncmp(line,"REMARK  score",13)){
			linkscore[id][lid]=atof(line+14);
		}
		if(!strncmp(line,"ATOM",4)){
			if(rn==0||strncmp(line+17,p[id][lid].res[rn-1].nam,9)){
				strncpy(p[id][lid].res[rn].nam,line+17,9);
				p[id][lid].res[rn].nam[9]='\0';
				p[id][lid].res[rn].start=n;
				p[id][lid].res[rn].sstart=n+4;
				if(rn!=0) p[id][lid].res[rn-1].end=n-1;
				rn++;
			}
			getxyz(line,p[id][lid].xyz[n]);
			strncpy(p[id][lid].nam[n],line+13,3);
			p[id][lid].nam[n][3]='\0';
			if(!strncmp(line+13,"N  ",3)){
				p[id][lid].res[rn-1].N=n;
			}
			else if(!strncmp(line+13,"CA ",3)){
				p[id][lid].res[rn-1].CA=n;
			}
			else if(!strncmp(line+13,"C  ",3)){
				p[id][lid].res[rn-1].C=n;
			}
			else if(!strncmp(line+13,"O  ",3)){
				p[id][lid].res[rn-1].O=n;
			}
			else if(!strncmp(line+13,"CB ",3)){
				p[id][lid].res[rn-1].sstart=n;
			}
			n++;
		}
		if(!strncmp(line,"END",3)){
			p[id][lid].an=n;
			p[id][lid].resn=rn;
			p[id][lid].res[rn-1].end=n-1;
			/*printf("%d %d\n", n,rn);*/
			n=0;rn=0;
		}
	}
	fclose(pf);
}
void cpx(float a[3],float b[3])
{
	int i;
	for(i=0;i<3;i++) a[i]=b[i];
}
void assemble(FRAGMENT *p1, FRAGMENT *l2, FRAGMENT *p2, FRAGMENT *l1, int t1, int t2, int tertyp1[2], int tertyp2[2],char chid,CYCPEP *cyc)
{
	int atno=1,resn=0;
	int ir,ip;
	int sr,er;
	int sa,ea;
	int i,s,e;
	FRAGMENT *p;
      char resname[3];

	for(ip=0;ip<4;ip++){
		if(ip==0){ p=p1; sr=0;er=p1->resn;if(t1==1) sa=0; else sa=1; if(t2==1) ea=0; else ea=1;} 
		else  if(ip==1){ p=l2;if(t2==1){ sr=1;er=l2->resn-1;sa=0; ea=0;} else{sr=0;er=l2->resn;sa=2; ea=2;} } 
		else if(ip==2){ p=p2; sr=0;er=p2->resn;if(t2==1) sa=0; else sa=1; if(t1==1) ea=0; else ea=1;} 
		else if(ip==3){ p=l1; if(t1==1){ sr=1;er=l1->resn-1;sa=0; ea=0;} else{sr=0;er=l1->resn;sa=2; ea=2;}} 
		for(ir=sr;ir<er;ir++){
			if(!(sa==2&&ir==sr)){
				resname[0]=p->res[ir].nam[0];
				resname[1]=p->res[ir].nam[1];
				resname[2]=p->res[ir].nam[2];
			}
			if(ea==2&&ir==er-1){
				if(p==l2){
					resname[0]=p2->res[0].nam[0];
					resname[1]=p2->res[0].nam[1];
					resname[2]=p2->res[0].nam[2];
				}
				else if(p==l1){
					resname[0]=p1->res[0].nam[0];
					resname[1]=p1->res[0].nam[1];
					resname[2]=p1->res[0].nam[2];
					resn=0;
				}
			}
			sprintf(cyc->res[resn].resnam,"%c%c%c %c%4d",resname[0],resname[1],resname[2],chid,resn+1);
			cyc->res[resn].satn=0;
 			for(i=p->res[ir].start;i<=p->res[ir].end;i++){
				if(ir==sr&&sa!=0&&i==p->res[ir].N) continue;
				if(ir==sr&&sa==2&&i==p->res[ir].CA) continue;

				if(ir==er-1){
					if(ea==1){
						if(i==p->res[ir].C||i==p->res[ir].O) continue;
					}
					if(ea==2){
						/*printf("%s\n",p->res[ir].nam);*/
						if(i>=p->res[ir].CA) break;
					}
				}
				if(sa==2&&ir==sr){
					if(i!=p->res[ir].C&&i!=p->res[ir].O) continue;
				}
				if((p==p1)&&ir==sr&&tertyp1[0]==0&&i>=p->res[ir].sstart) break;
				if((p==p2)&&ir==sr&&tertyp2[0]==0&&i>=p->res[ir].sstart) break;
				if((p==p2)&&ir==er-1&&tertyp1[1]==0&&i>=p->res[ir].sstart) break;
				if((p==p1)&&ir==er-1&&tertyp2[1]==0&&i>=p->res[ir].sstart) break;

				if(!strcmp(p->nam[i],"N  ")) cpx(cyc->res[resn].N, p->xyz[i]);
				else if(!strcmp(p->nam[i],"CA ")) cpx(cyc->res[resn].CA, p->xyz[i]);
				else if(!strcmp(p->nam[i],"C  ")) cpx(cyc->res[resn].C, p->xyz[i]);
				else if(!strcmp(p->nam[i],"O  ")) cpx(cyc->res[resn].O, p->xyz[i]);
				
				else{
					cpx(cyc->res[resn].side[cyc->res[resn].satn], p->xyz[i]);
					strcpy(cyc->res[resn].nam[cyc->res[resn].satn], p->nam[i]); 
					cyc->res[resn].nam[cyc->res[resn].satn][3]='\0';
					cyc->res[resn].satn++;
				}
				
				atno++;
			}
			resn++;
			if(ir==er-1&&ea!=0){
				resn--; 
			}
		}
	}
}

void printpepstructure(FILE *pf,int id,char fragname[2][100], float dockscore[2], CYCPEP *cyc)
{
	int atno=1,resn=1;
	int ir,ia;
	int ih;

	fprintf(pf,"TITLE %2d\n",id+1);
	fprintf(pf,"REMARK  dock1 %s\n",fragname[0]);
	fprintf(pf,"REMARK  dock2 %s\n",fragname[1]);
	fprintf(pf,"REMARK  link1 %s %d\n",cyc->linkname[0], cyc->typ[0]);
	fprintf(pf,"REMARK  link2 %s %d\n",cyc->linkname[1], cyc->typ[1]);
	fprintf(pf,"REMARK  score %8.3f %8.3f%8.3f%8.3f%8.3f\n",dockscore[0]+dockscore[1]+cyc->linkscore[0]+cyc->linkscore[1],dockscore[0],dockscore[1],cyc->linkscore[0],cyc->linkscore[1]);
	fprintf(pf,"REMARK  resn %2d %3d%3d%3d%3d\n",cyc->resn,cyc->len[0],cyc->len[1],cyc->len[2],cyc->len[3]);
	fprintf(pf,"REMARK  backbone HB No. %2d\n",cyc->hbn);
	for(ih=0;ih<cyc->hbn;ih++){
		fprintf(pf,"REMARK  HB %2d [%s]--[%s] %8.3f\n",ih+1,cyc->res[cyc->index[ih][0]].resnam,cyc->res[cyc->index[ih][1]].resnam, cyc->hbdis[ih]);
	}
	for(ir=0;ir<cyc->resn;ir++){
		fprintf(pf,"ATOM %6d  N   %s    %8.3f%8.3f%8.3f\n",atno,cyc->res[ir].resnam,cyc->res[ir].N[0],cyc->res[ir].N[1],cyc->res[ir].N[2]); atno++;
		fprintf(pf,"ATOM %6d  CA  %s    %8.3f%8.3f%8.3f\n",atno,cyc->res[ir].resnam,cyc->res[ir].CA[0],cyc->res[ir].CA[1],cyc->res[ir].CA[2]); atno++;
		fprintf(pf,"ATOM %6d  C   %s    %8.3f%8.3f%8.3f\n",atno,cyc->res[ir].resnam,cyc->res[ir].C[0],cyc->res[ir].C[1],cyc->res[ir].C[2]); atno++;
		fprintf(pf,"ATOM %6d  O   %s    %8.3f%8.3f%8.3f\n",atno,cyc->res[ir].resnam,cyc->res[ir].O[0],cyc->res[ir].O[1],cyc->res[ir].O[2]); atno++;
		for(ia=0;ia<cyc->res[ir].satn;ia++){
			fprintf(pf,"ATOM %6d  %s %s    %8.3f%8.3f%8.3f\n",atno,cyc->res[ir].nam[ia],cyc->res[ir].resnam,cyc->res[ir].side[ia][0],cyc->res[ir].side[ia][1],cyc->res[ir].side[ia][2]); atno++;
		}
	}

	fprintf(pf,"END\n");
}
int checkcollision(FRAGMENT *p1,FRAGMENT *p2)
{
        int ia,ja;
        float dis;

        for(ia=0;ia<p1->an;ia++){
                for(ja=0;ja<p2->an;ja++){
                        dis=distance2(p1->xyz[ia], p2->xyz[ja]);
                        if(dis<2.5*2.5){
                                return 0;
                        }
                }
        }
        return 1;
}

#define MAXMODEL 25
void calHBond(CYCPEP *cyc)
{
	int pair[MAXHBN][2];
	int hbn,ih;
 	float hbdis[MAXHBN];
	int index[MAXHBN];

	addHN(cyc);
	hbn=dispair(cyc, pair, hbdis);
	hbn=checkangle(cyc, hbn, pair, index);
	cyc->hbn=hbn;
	for(ih=0;ih<hbn;ih++){
		cyc->index[ih][0]=pair[index[ih]][0];
		cyc->index[ih][1]=pair[index[ih]][1];
		cyc->hbdis[ih]=hbdis[index[ih]];
	}
}

int main(int argc, char *argv[])
{
	FRAGMENT p[2][BESTCLUSTERNUM],p1,p2;
	CYCPEP cyc[MAXMODEL];
	char fragname[2][100];
	char linkname[2][BESTCLUSTERNUM][100];
	int typ[2][BESTCLUSTERNUM];
	float fragscore[2];
	float linkscore[2][BESTCLUSTERNUM];
	int tertyp[2][BESTCLUSTERNUM][2];
	int clusterN[2];
	int il,jl;
	int col;
	int maxlinkres=atoi(argv[3]);
	int modn=0,imod;
	FILE *pf;

	readlink(argv[1], p,typ,fragname,linkname,fragscore,linkscore,tertyp,clusterN);
	printf("%s %s %s\n",argv[1],fragname[0],fragname[1]);
	readfrag(fragname[0], &p1);
	readfrag(fragname[1], &p2);
	for(il=0;il<clusterN[0];il++) for(jl=0;jl<clusterN[1];jl++){
		if(modn>=MAXMODEL){
			printf("WARNING model number > MAXMODEL %d\n",MAXMODEL);
			break;
		}
		if(p[0][il].resn+p[1][jl].resn>maxlinkres) continue;
		col=checkcollision(&(p[1][jl]),&(p[0][il]));
		if(col==0) continue;
		assemble(&p1, &(p[1][jl]),&p2,&(p[0][il]),typ[0][il],typ[1][jl],tertyp[0][il],tertyp[1][jl],argv[4][0],cyc+modn);
		strcpy(cyc[modn].linkname[0],linkname[0][il]); strcpy(cyc[modn].linkname[1],linkname[1][jl]);
		cyc[modn].typ[0]=typ[0][il]; cyc[modn].typ[1]=typ[1][jl];
		cyc[modn].linkscore[0]=linkscore[0][il]; cyc[modn].linkscore[1]=linkscore[1][jl];
		cyc[modn].len[0]=p1.resn; cyc[modn].len[1]=p[1][jl].resn;cyc[modn].len[2]=p2.resn; cyc[modn].len[3]=p[0][il].resn;
		cyc[modn].resn=cyc[modn].len[0]+cyc[modn].len[1]+cyc[modn].len[2]+cyc[modn].len[3]-4;
		modn++;
	}
	if((pf=fopen(argv[2],"w"))==NULL){
		printf("ERROR: Can not create the file %s\n",argv[2]);
		exit(0);
	}
	for(imod=0;imod<modn;imod++){
		calHBond(cyc+imod);
		printpepstructure(pf, imod,fragname,fragscore, cyc+imod);
	}
	fclose(pf);
}

