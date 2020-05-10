//#define N 100

#define ZERO 0
#define FIRST 1
#define AXIS_X 0
#define AXIS_Y 1
#define LEFT 0
#define RIGHT 1
#define eps 1E-6

#define grid (200)

extern int Nx, Ny;

// template
#define S 1
// Courant
#define k 0.4
#define LENGTHx 2.0
#define LENGTHy 2.0
#define hx (LENGTHx / (Nx))
#define hy (LENGTHy / (Ny))
#define taux (k * hx / 1.0)
#define tauy (k * hy / 1.0)

#define yidxL (iy == 0 ? (Ny - 1) : (iy - 1))
#define yidxR (iy == (Ny - 1) ? 0 : (iy + 1))
#define xidxL (ix == 0 ? (Nx - 1) : (ix - 1))
#define xidxR (ix == (Nx - 1) ? 0 : (ix + 1))

extern double **P_c, **P_n;
extern double **Px_c, **Px_n;
extern double **Py_c, **Py_n;
extern double **Pxy_c, **Pxy_n;
extern double **C; // Material

double func(double x, double y);

double dfunc(double x, double y);

void createArrays(void);

void freeArrays(void);

void init(void);

void singleStepCIP(void);

void omega(double **P_c, double *W, int iy, int ix, int axis);

void domega(double **dP_c, double *dW, int iy, int ix, int axis);

void reconstruct(double *W, double *dW, double *Wn, double *dWn, int iy, int ix);

void backward(double **Pn, double **dPn, double *Wn, double *dWn, int iy, int ix);

void update(void);

void singleStepCIR(void);

void interpolatedCoefs(double *u, double *du, int dir, int iy, int ix, double *coefs);

double interpolatedValue(double *coefs, int derivate);

void debugGnuplotFileName(char *filename1);
