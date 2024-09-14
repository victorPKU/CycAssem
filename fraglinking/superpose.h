#ifndef _superpose_h
#define _superpose_h

void do_rot(int an,float x[][3],float xt[][3],float R[3][3]);
void calc_fit_R(int an,float xp[][3],float x[][3],float R[3][3]);
float calcrmsd(int an,float xp[][3],float x[][3]);
float fit_terminal(float tptermxyz[4][3],float lptatx[4][3],float trans[3], float rotR[3][3]);

#endif

