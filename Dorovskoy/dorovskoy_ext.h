//#define N 100

double W1[2], W2[2], W3[2], W4[2];
double dW1[2], dW2[2], dW3[2], dW4[2];
double *V_s_c, *V_f_c, *H_c, *P_c;
double *V_s_n, *V_f_n, *H_n, *P_n;
double *dV_s_c, *dV_f_c, *dH_c, *dP_c;
double *dV_s_n, *dV_f_n, *dH_n, *dP_n;

double mu = 2.2815E9;
double alpha = 2507.905873;
double K = 2.642306867E9;
double rho_s = 1350.0;
double rho_f = 100.0;
double rho_0 = rho_s + rho_f;

//double phi = 0.1;
//double Cp1 = 2000;
//double Cp2 = 450;
//double Cs = 1300;


void createArrays(void);

void freeArrays(void);

void init(void);

void singleStep(void);

void toInvariants(int ind);
void toDInvariants(int ind);

void interpolatedCoefs(double *u, double *du, double c, int dir, int ind, double *coefs);

double interpolatedValue(double *coefs, int derivate);

void fillStencilValues(int ind);

void debugPrint(void);

void debugGnuplot(void);

void debugGnuplotFileName(char *filename1);

void singleStepCIR(void);

double interpolatedValueCIR(double *u, double c, int dir, int ind);