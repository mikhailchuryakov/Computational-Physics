//#define N 100

double **P_c, **P_n, **U_c, **U_n, **V_c, **V_n;
double **Px_c, **Px_n, **Ux_c, **Ux_n, **Vx_c, **Vx_n;
double **Py_c, **Py_n, **Uy_c, **Uy_n, **Vy_c, **Vy_n;
double **Pxy_c, **Pxy_n, **Uxy_c, **Uxy_n, **Vxy_c, **Vxy_n;
double **C, **Ro; // Material

double func(double x, double y);

double dfunc(double x, double y);

void createArrays(void);

void freeArrays(void);

void init(void);

void singleStep(void);

void omega(double **P_c, double **U_c, double *W, int iy, int ix, int axis);

void domega(double **dP_c, double **dU_c, double *dW, int iy, int ix, int axis);

void reconstruct(double *W, double *dW, double *Wn, double *dWn, int iy, int ix);

void backward(double **Pn, double **Un, double **dPn, double **dUn, double *Wn, double *dWn, int iy, int ix);

void update(void);

void singleStepCIR(void);

void interpolatedCoefs(double *u, double *du, int dir, int iy, int ix, double *coefs);

double interpolatedValue(double *coefs, int derivate);

double interpolatedValueCIR(double *u, int dir, int iy, int ix);

void fillStencilValues();

void debugGnuplotFileName(char *filename1, char *filename2);