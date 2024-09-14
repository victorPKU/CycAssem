#ifndef _geometry_h
#define _geometry_h

void cpnx(float xi[][3],float xo[][3], int n);
void cpx(float xi[3],float xo[3]);
float distance(float x1[3], float x2[3]);
float distance2(float x1[3], float x2[3]);
float cal_angle(float x1[3], float x2[3], float x3[3]);
float cal_dih(float x1[3], float x2[3], float x3[3], float x4[3]);
void get_center(int n, float x[][3], float center[3]);
void direction( float a[3], float b[3], float c[3]);
void translate(int n, float x0[][3],float x[][3],float v[3]);
void irtranslate(int n, float x0[][3],float x[][3],float v[3]);

#define PI 3.14159265358979323846
#define FOURPI		(4.*PI)
#define TORAD(A)     ((A)*0.017453293)
#define TODEG(A)     ((A)*57.295779513)

#endif

