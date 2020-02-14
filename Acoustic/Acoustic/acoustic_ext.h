//#define N 100

double Wmns[2], Wpls[2];
double dWmns[2], dWpls[2];
double *P_c, *P_n, *U_c, *U_n;
double *dP_c, *dP_n, *dU_c, *dU_n;
double *C, *Ro; // Material


void createArrays(void);

void freeArrays(void);

void init(void);

void singleStep(void);

void result(void);

void toInvariants(int ind, double Z_l, double Z_r);
void toDInvariants(int ind, double Z_l, double Z_r);

void singleStepCIR(void);

void interpolatedCoefs(double *u, double *du, int dir, int ind, double *coefs);

double interpolatedValue(double *coefs, int derivate);

double interpolatedValueCIR(double *u, int dir, int ind);

void fillStencilValues(int ind);

void debugPrint(void);

void debugGnuplot(void);

void debugGnuplotFileName(char *filename1, char *filename2);