#define _USE_MATH_DEFINES

#include <stdio.h>
#include <cmath>
#include <stdlib.h>
#include <omp.h>

#include "acoustic2d_ext.h"

#define grid (100)

int Nx = grid, Ny = grid;

#define LENGTHx 2.0
#define LENGTHy 2.0

#define k 0.4
#define hx (LENGTHx / (Nx))
#define hy (LENGTHy / (Ny))
#define taux (k * hx / 1.0)
#define tauy (k * hy / 1.0)

// template
#define S 1
//int stencil_l[S + 1] = { -1, 0 };
//int stencil_r[S + 1] = { 0, 1 };

//double Wpls[2], Wmns[2], dWpls[2], dWmns[2];
// Wn[2],dWn[2] = {W_left, W_right}
//double Wn[2], dWn[2];

#define LEFT 0
#define RIGHT 1
//

#define M (2.0 * (Nx * 0.4 * hx / 2.0) / C[0][0] / taux)
//#define M 0.0


#define eps 1E-6
double func(double x, double y) {
	
	if (-0.4 - x <= eps && 0.1 - x >= eps){
		return 1.0 * pow(sin(2.0 * M_PI *(x - 0.6)), 4.0);
	}
	else {
		return 0;
	}
	
	//return exp(-100 * (x)*(x)) * exp(-100 * (y)*(y));
	//return exp(-100 * (x + 0.5)*(x + 0.5));
	//return 0.0;
}

double func_x(double x, double y) {
	
	if (-0.4 - x <= eps && 0.1 - x >= eps){
		return 8.0 * M_PI * pow(sin(2.0 * M_PI *(x - 0.6)), 3.0) * cos(2.0 * M_PI *(x - 0.6));
	}
	else {
		return 0;
	}
	
	//return -200 * (x) * exp(-100 * (x)*(x)) * exp(-100 * (y)*(y));
	//return 0.0;
}

double func_y(double x, double y) {
	return -200 * (y) * exp(-100 * (x)*(x)) * exp(-100 * (y)*(y));
}

double func_xy(double x, double y) {
	return 40000 * (x) * (y) * exp(-100 * (x)*(x)) * exp(-100 * (y)*(y));
}

int main() {
	char str[20], strCIR[20];
	double time;
	int number_of_processes = 0;
	int del = 25;

	printf("start\n");
	for (Nx = grid, Ny = grid; Nx <= 400 + 1; Nx = Nx * 2, Ny = Ny * 2) {
#pragma omp parallel
		{
			int step = 0;
#pragma omp single
			{
				time = omp_get_wtime();
				createArrays();
				printf("Arrays created\n");
				init();
				printf("initialized\n");
				printf("%d\t%lf\t%f\t%lf\t%lf\n", Nx, hx, 1.0 / hx, taux, M);
				number_of_processes = omp_get_num_threads();
				//printf("%d\t", number_of_processes);
			}
			for (step = 0; step < M; step++) {
#pragma omp single
				{
					if (step % 100 == 0) printf("%d\n", step);
					/*
					if ((step % del) == 0) {
						sprintf(str, "gif/P%d.dat", step / del + 1);
						sprintf(strCIR, "gif/U%d.dat", step / del + 1);
						debugGnuplotFileName(str, strCIR);
					}
					*/
				}
				singleStep();
			}
		}
		
		sprintf(str, "gif/P%d.dat", (int)M / del + 1);
		sprintf(strCIR, "gif/U%d.dat", (int)M / del + 1);
		debugGnuplotFileName(str, strCIR);
		
		sprintf(str, "P%d.dat", Nx);
		sprintf(strCIR, "U%d.dat", Nx);
		debugGnuplotFileName(str, strCIR);

		freeArrays();
		time = (omp_get_wtime() - time);
		printf("Concurrent OpenMP calculation with %d threads took %f seconds.\n", number_of_processes, time);
	}
	return 0;
}

void init(void) {
	int ix, iy;
	for (iy = 0; iy < Ny; iy++) {
		for (ix = 0; ix < Nx; ix++) {
			double x = -(LENGTHx / 2.0) + ix * hx;
			double y = -(LENGTHy / 2.0) + iy * hy;
			
			C[iy][ix] = 1.0;
			Ro[iy][ix] = 1.0;

			// nodes
			P_c[iy][ix] = func(x, y);
			U_c[iy][ix] = P_c[iy][ix] / (C[iy][ix] * Ro[iy][ix]);
			V_c[iy][ix] = P_c[iy][ix] / (C[iy][ix] * Ro[iy][ix]);
			V_n[iy][ix] = V_c[iy][ix];

			Px_c[iy][ix] = func_x(x, y);
			Ux_c[iy][ix] = 0.0;// Px_c[iy][ix] / (C[iy][ix] * Ro[iy][ix]);
			Vx_c[iy][ix] = 0.0;//Px_c[iy][ix] / (C[iy][ix] * Ro[iy][ix]);
			Vx_n[iy][ix] = Vx_c[iy][ix];

			Py_c[iy][ix] = 0.0;//func_y(x, y);
			Uy_c[iy][ix] = 0.0;//Py_c[iy][ix] / (C[iy][ix] * Ro[iy][ix]);
			Vy_c[iy][ix] = 0.0;//Py_c[iy][ix] / (C[iy][ix] * Ro[iy][ix]);
			Vy_n[iy][ix] = Vy_c[iy][ix];

			Pxy_c[iy][ix] = 0.0;//func_xy(x, y);
			Uxy_c[iy][ix] = 0.0;//Pxy_c[iy][ix] / (C[iy][ix] * Ro[iy][ix]);
			Vxy_c[iy][ix] = 0.0;//Pxy_c[iy][ix] / (C[iy][ix] * Ro[iy][ix]);
			Vxy_n[iy][ix] = Vxy_c[iy][ix];

			//printf("%lf\t%lf\t%lf\n", x, dP_c[ind], dU_c[ind]);
//			if (fabs(x - (0.5)) < eps){
	//			printf("%lf\t%lf\t%lf\n", x, C[iy][ix], C[iy][ix]);
		//	}
		}
	}
}

#define ZERO 0
#define FIRST 1
#define AXIS_X 0
#define AXIS_Y 1
void singleStep(void) {
	int ix, iy;
	int start = omp_get_thread_num();
	int stride = omp_get_num_threads();
	double *W, *dW;
	double *Wn, *dWn;
	// W[4] = {Wpls[0], Wpls[1], Wmns[0], Wmns[1]}
	W = (double *)malloc(4 * sizeof(double));
	// dW[4] = {dWpls[0], dWpls[1], dWmns[0], dWmns[1]}
	dW = (double *)malloc(4 * sizeof(double));

	Wn = (double *)malloc(2 * sizeof(double));
	dWn = (double *)malloc(2 * sizeof(double));

	// X-step
#pragma omp barrier
	for (iy = start; iy < Ny; iy += stride) {
		for (ix = 0; ix < Nx; ix++) {
			// update *W
			omega(P_c, U_c, W, iy, ix, AXIS_X);
			// update *dW
			domega(Px_c, Ux_c, dW, iy, ix, AXIS_X);
			// update $Wn, $dWn
			reconstruct(W, dW, Wn, dWn, iy, ix);
			// update $Pn, $Un, $dPn, $dUn
			backward(P_n, U_n, Px_n, Ux_n, Wn, dWn, iy, ix);

			omega(Py_c, Uy_c, W, iy, ix, AXIS_X);
			domega(Pxy_c, Uxy_c, dW, iy, ix, AXIS_X);
			reconstruct(W, dW, Wn, dWn, iy, ix);
			backward(Py_n, Uy_n, Pxy_n, Uxy_n, Wn, dWn, iy, ix);
		}	
	}

	update();
	//printf("%lf\t", Pxy_n[20][20]);
	/*
	// Y-step
#pragma omp barrier
	for (ix = start; ix < Nx; ix += stride) {
		for (iy = 0; iy < Ny; iy++) {
			omega(P_c, V_c, W, iy, ix, AXIS_Y);
			domega(Px_c, Vx_c, dW, iy, ix, AXIS_Y);
			reconstruct(W, dW, Wn, dWn, iy, ix);
			backward(P_n, V_n, Px_n, Vx_n, Wn, dWn, iy, ix);
			

			omega(Py_c, Vy_c, W, iy, ix, AXIS_Y);
			domega(Pxy_c, Vxy_c, dW, iy, ix, AXIS_Y);
			reconstruct(W, dW, Wn, dWn, iy, ix);
			backward(Py_n, Vy_n, Pxy_n, Vxy_n, Wn, dWn, iy, ix);

			//printf("%lf\t%lf\t%lf\t%lf\n", P_n[iy][ix], V_n[iy][ix], Px_n[iy][ix], Vx_n[iy][ix]);
			//printf("%lf\t%lf\t%lf\t%lf\n", Py_n[iy][ix], Vy_n[iy][ix], Pxy_n[iy][ix], Vxy_n[iy][ix]);
			//printf("%lf\t%lf\t%lf\t%lf\n", P_n[iy][ix], U_n[iy][ix], Px_n[iy][ix], Ux_n[iy][ix]);
			//printf("%lf\t%lf\t%lf\t%lf\n", Py_n[iy][ix], Uy_n[iy][ix], Pxy_n[iy][ix], Uxy_n[iy][ix]);
		}
	}

	update();
	*/
	free(W);
	free(dW);
	free(Wn);
	free(dWn);
	//printf("%lf\n", Pxy_n[20][20]);
}

#define yidxL (iy == 0 ? (Ny - 1) : (iy - 1))
#define yidxR (iy == (Ny - 1) ? 0 : (iy + 1))
#define xidxL (ix == 0 ? (Nx - 1) : (ix - 1))
#define xidxR (ix == (Nx - 1) ? 0 : (ix + 1))

void omega(double **P_c, double **U_c, double *W, int iy, int ix, int axis){
	double Z_l, Z_r;
	Z_l = C[iy][ix] * Ro[iy][ix];
	Z_r = C[iy][ix] * Ro[iy][ix];
	if (axis == AXIS_X){
		W[0] = (P_c[iy][xidxL] + Z_l * U_c[iy][xidxL]) / (2 * Z_l);
		W[1] = (P_c[iy][ix] + Z_l * U_c[iy][ix]) / (2 * Z_l);

		W[2] = (-P_c[iy][ix] + Z_r * U_c[iy][ix]) / (2 * Z_r);
		W[3] = (-P_c[iy][xidxR] + Z_r * U_c[iy][xidxR]) / (2 * Z_r);
		//printf("%lf\t%lf\n", Wmns[0], Wmns[1]);
	}
	else {
		W[0] = (P_c[yidxL][ix] + Z_l * U_c[yidxL][ix]) / (2 * Z_l);
		W[1] = (P_c[iy][ix] + Z_l * U_c[iy][ix]) / (2 * Z_l);

		W[2] = (-P_c[iy][ix] + Z_r * U_c[iy][ix]) / (2 * Z_r);
		W[3] = (-P_c[yidxR][ix] + Z_r * U_c[yidxR][ix]) / (2 * Z_r);
		//printf("%lf\t%lf\n", Wmns[0], Wmns[1]);
	}
}

void domega(double **dP_c, double **dU_c, double *dW, int iy, int ix, int axis){
	double Z_l, Z_r;
	Z_l = C[iy][ix] * Ro[iy][ix];
	Z_r = C[iy][ix] * Ro[iy][ix];
	if (axis == AXIS_X){
		dW[0] = (dP_c[iy][xidxL] + Z_l * dU_c[iy][xidxL]) / (2 * Z_l);
		dW[1] = (dP_c[iy][ix] + Z_l * dU_c[iy][ix]) / (2 * Z_l);

		dW[2] = (-dP_c[iy][ix] + Z_r * dU_c[iy][ix]) / (2 * Z_r);
		dW[3] = (-dP_c[iy][xidxR] + Z_r * dU_c[iy][xidxR]) / (2 * Z_r);
		//printf("%lf\t%lf\n", Wmns[0], Wmns[1]);
	}
	else {
		dW[0] = (dP_c[yidxL][ix] + Z_l * dU_c[yidxL][ix]) / (2 * Z_l);
		dW[1] = (dP_c[iy][ix] + Z_l * dU_c[iy][ix]) / (2 * Z_l);

		dW[2] = (-dP_c[iy][ix] + Z_r * dU_c[iy][ix]) / (2 * Z_r);
		dW[3] = (-dP_c[yidxR][ix] + Z_r * dU_c[yidxR][ix]) / (2 * Z_r);
		//printf("%lf\t%lf\n", dWmns[0], dWmns[1]);
	}
}

void reconstruct(double *W, double *dW, double *Wn, double *dWn, int iy, int ix){
	// coefs[5] = { x, a, b, c, d }
	double coefs1[5], coefs2[5];
	double stencil_values_w1[S + 1];
	double stencil_values_w2[S + 1];
	double stencil_Dvalues_w1[S + 1];
	double stencil_Dvalues_w2[S + 1];

	unsigned int j;
	for (j = 0; j < S + 1; j++) {
		stencil_values_w2[j] = W[j];
		stencil_Dvalues_w2[j] = dW[j];

		stencil_values_w1[j] = W[j + S + 1];
		stencil_Dvalues_w1[j] = dW[j + S + 1];
	}

	interpolatedCoefs(stencil_values_w2, stencil_Dvalues_w2, LEFT, iy, ix, coefs2);
	interpolatedCoefs(stencil_values_w1, stencil_Dvalues_w1, RIGHT, iy, ix, coefs1);

	Wn[0] = interpolatedValue(coefs2, ZERO); // W_pls
	Wn[1] = interpolatedValue(coefs1, ZERO); // W_mns

	dWn[0] = interpolatedValue(coefs2, FIRST); // dW_pls
	dWn[1] = interpolatedValue(coefs1, FIRST); // dW_mns
}

void backward(double **Pn, double **Un, double **dPn, double **dUn, double *Wn, double *dWn, int iy, int ix) {
	double Z_l, Z_r;
	Z_l = C[iy][ix] * Ro[iy][ix];
	Z_r = C[iy][ix] * Ro[iy][ix];

	Pn[iy][ix] = 2.0 * Z_r * Z_l * (Wn[0] - Wn[1]) / (Z_l + Z_r);
	Un[iy][ix] = 2.0 * (Z_l * Wn[0] + Z_r * Wn[1]) / (Z_l + Z_r);
	//printf("%lf\t%lf\n", Wn[0], Wn[1]);

	dPn[iy][ix] = 2.0 * Z_r * Z_l * (dWn[0] - dWn[1]) / (Z_l + Z_r);
	dUn[iy][ix] = 2.0 * (Z_l * dWn[0] + Z_r * dWn[1]) / (Z_l + Z_r);

	//dPn[iy][ix] = 2.0 * (C[iy][ix] * Z_l * dWn[0] - C[iy][ix] * Z_r * dWn[1]) / (Z_l + Z_r);
	//dUn[iy][ix] = 2.0 * (C[iy][ix] * Z_l * Z_r * dWn[0] + C[iy][ix] * Z_l * Z_r * dWn[1]) / (Z_l + Z_r);
}

void interpolatedCoefs(double *u, double *du, int dir, int iy, int ix, double *coefs) {
	// direction
	if (dir == LEFT) {
		coefs[1] = (du[1] + du[0]) / hx / hx - 2.0 * (u[1] - u[0]) / hx / hx / hx;
		coefs[2] = (2.0 * du[1] + du[0]) / hx - 3.0 * (u[1] - u[0]) / hx / hx;
		coefs[3] = du[1];
		coefs[4] = u[1];
		coefs[0] = -C[iy][ix] * taux;
	}
	else {
		coefs[1] = (du[1] + du[0]) / hx / hx - 2.0 * (u[1] - u[0]) / hx / hx / hx;
		coefs[2] = -(2.0 * du[0] + du[1]) / hx + 3.0 * (u[1] - u[0]) / hx / hx;
		coefs[3] = du[0];
		coefs[4] = u[0];
		coefs[0] = C[iy][ix] * taux;
	}
}

void update() {
	int ix, iy;
	int start = omp_get_thread_num();
	int stride = omp_get_num_threads();

	for (iy = start; iy < Ny; iy += stride) {
		for (ix = 0; ix < Nx; ix++) {
			P_c[iy][ix] = P_n[iy][ix];
			U_c[iy][ix] = U_n[iy][ix];
			V_c[iy][ix] = V_n[iy][ix];

			Px_c[iy][ix] = Px_n[iy][ix];
			Ux_c[iy][ix] = Ux_n[iy][ix];
			Vx_c[iy][ix] = Vx_n[iy][ix];

			Py_c[iy][ix] = Py_n[iy][ix];
			Uy_c[iy][ix] = Uy_n[iy][ix];
			Vy_c[iy][ix] = Vy_n[iy][ix];

			Pxy_c[iy][ix] = Pxy_n[iy][ix];
			Uxy_c[iy][ix] = Uxy_n[iy][ix];
			Vxy_c[iy][ix] = Vxy_n[iy][ix];
		}
	}
}

double interpolatedValue(double *coefs, int derivate) {
	double a = coefs[1], b = coefs[2], c = coefs[3], d = coefs[4], x = coefs[0];
	double val;
	if (derivate == FIRST) {
		val = 3.0 * a * x * x + 2.0 * b * x + c;
	}
	else {
		val = a * x * x * x + b * x * x + c * x + d;
	}
	return val;
}

void debugGnuplotFileName(char *filename1, char *filename2){
	FILE *pFileP, *pFileU;
	pFileP = fopen(filename1, "w");
	//pFileU = fopen(filename2, "w");
	int ix, iy;
	for (iy = 0; iy < Ny; iy++) {
		for (ix = 0; ix < Nx; ix++) {
			//if (abs(P_c[iy][ix]) > 1E-5) 
				fprintf(pFileP, "%.12lf\t%.12lf\t%.12lf\n", -(LENGTHx / 2.0) + hx * ix, -(LENGTHx / 2.0) + hy * iy, P_c[iy][ix]);
		}
	}
	fclose(pFileP);
//	fclose(pFileU);
}

void createArrays(){
	int i;

	C = (double **)malloc(Ny * sizeof(double *));
	Ro = (double **)malloc(Ny * sizeof(double *));

	P_c = (double **)malloc(Ny * sizeof(double *));
	P_n = (double **)malloc(Ny * sizeof(double *));
	U_c = (double **)malloc(Ny * sizeof(double *));
	U_n = (double **)malloc(Ny * sizeof(double *));
	V_c = (double **)malloc(Ny * sizeof(double *));
	V_n = (double **)malloc(Ny * sizeof(double *));

	Px_c = (double **)malloc(Ny * sizeof(double *));
	Px_n = (double **)malloc(Ny * sizeof(double *));
	Ux_c = (double **)malloc(Ny * sizeof(double *));
	Ux_n = (double **)malloc(Ny * sizeof(double *));
	Vx_c = (double **)malloc(Ny * sizeof(double *));
	Vx_n = (double **)malloc(Ny * sizeof(double *));

	Py_c = (double **)malloc(Ny * sizeof(double *));
	Py_n = (double **)malloc(Ny * sizeof(double *));
	Uy_c = (double **)malloc(Ny * sizeof(double *));
	Uy_n = (double **)malloc(Ny * sizeof(double *));
	Vy_c = (double **)malloc(Ny * sizeof(double *));
	Vy_n = (double **)malloc(Ny * sizeof(double *));

	Pxy_c = (double **)malloc(Ny * sizeof(double *));
	Pxy_n = (double **)malloc(Ny * sizeof(double *));
	Uxy_c = (double **)malloc(Ny * sizeof(double *));
	Uxy_n = (double **)malloc(Ny * sizeof(double *));
	Vxy_c = (double **)malloc(Ny * sizeof(double *));
	Vxy_n = (double **)malloc(Ny * sizeof(double *));

	for (i = 0; i < Ny; i++){
		C[i] = (double *)malloc(Nx * sizeof(double));
		Ro[i] = (double *)malloc(Nx * sizeof(double));

		P_c[i] = (double *)malloc(Nx * sizeof(double));
		P_n[i] = (double *)malloc(Nx * sizeof(double));
		U_c[i] = (double *)malloc(Nx * sizeof(double));
		U_n[i] = (double *)malloc(Nx * sizeof(double));
		V_c[i] = (double *)malloc(Nx * sizeof(double));
		V_n[i] = (double *)malloc(Nx * sizeof(double));

		Px_c[i] = (double *)malloc(Nx * sizeof(double));
		Px_n[i] = (double *)malloc(Nx * sizeof(double));
		Ux_c[i] = (double *)malloc(Nx * sizeof(double));
		Ux_n[i] = (double *)malloc(Nx * sizeof(double));
		Vx_c[i] = (double *)malloc(Nx * sizeof(double));
		Vx_n[i] = (double *)malloc(Nx * sizeof(double));

		Py_c[i] = (double *)malloc(Nx * sizeof(double));
		Py_n[i] = (double *)malloc(Nx * sizeof(double));
		Uy_c[i] = (double *)malloc(Nx * sizeof(double));
		Uy_n[i] = (double *)malloc(Nx * sizeof(double));
		Vy_c[i] = (double *)malloc(Nx * sizeof(double));
		Vy_n[i] = (double *)malloc(Nx * sizeof(double));

		Pxy_c[i] = (double *)malloc(Nx * sizeof(double));
		Pxy_n[i] = (double *)malloc(Nx * sizeof(double));
		Uxy_c[i] = (double *)malloc(Nx * sizeof(double));
		Uxy_n[i] = (double *)malloc(Nx * sizeof(double));
		Vxy_c[i] = (double *)malloc(Nx * sizeof(double));
		Vxy_n[i] = (double *)malloc(Nx * sizeof(double));
	}
	
}

void freeArrays(){
	int i;
	for (i = 0; i < Ny; i++) {
		free(C[i]);
		free(Ro[i]);

		free(Px_c[i]);
		free(Px_n[i]);
		free(Ux_c[i]);
		free(Ux_n[i]);
		free(Vx_c[i]);
		free(Vx_n[i]);

		free(Py_c[i]);
		free(Py_n[i]);
		free(Uy_c[i]);
		free(Uy_n[i]);
		free(Vy_c[i]);
		free(Vy_n[i]);

		free(Pxy_c[i]);
		free(Pxy_n[i]);
		free(Uxy_c[i]);
		free(Uxy_n[i]);
		free(Vxy_c[i]);
		free(Vxy_n[i]);
	}

	free(Px_c);
	free(Px_n);
	free(Ux_c);
	free(Ux_n);
	free(Vx_c);
	free(Vx_n);

	free(Py_c);
	free(Py_n);
	free(Uy_c);
	free(Uy_n);
	free(Vy_c);
	free(Vy_n);

	free(Pxy_c);
	free(Pxy_n);
	free(Uxy_c);
	free(Uxy_n);
	free(Vxy_c);
	free(Vxy_n);
}



// CIR
void singleStepCIR(void) {
	int iy, ix;
	int start = omp_get_thread_num();
	int stride = omp_get_num_threads();
	double stencil_values_w1[S + 1];
	double stencil_values_w2[S + 1];

	double *W = (double *)malloc(4 * sizeof(double));

	// X-step
	for (iy = start; iy < Ny; iy += stride) {
		for (ix = 0; ix < Nx; ix++) {
			double Wpls_n, Wmns_n;
			double Z_l, Z_r;
			unsigned int j;

			Z_l = C[iy][ix] * Ro[iy][ix];
			Z_r = C[iy][ix] * Ro[iy][ix];

			omega(P_c, U_c, W, iy, ix, AXIS_X);

			for (j = 0; j < S + 1; j++) {
				stencil_values_w2[j] = W[j];

				stencil_values_w1[j] = W[j + S + 1];
			}

			Wpls_n = interpolatedValueCIR(stencil_values_w2, LEFT, iy, ix);
			Wmns_n = interpolatedValueCIR(stencil_values_w1, RIGHT, iy, ix);

			P_n[iy][ix] = 2.0 * Z_r * Z_l * (Wpls_n - Wmns_n) / (Z_l + Z_r);
			U_n[iy][ix] = 2.0 * (Z_l * Wpls_n + Z_r * Wmns_n) / (Z_l + Z_r);
		}
	}

#pragma omp barrier

#pragma omp single
	{
		for (iy = 0; iy < Ny; iy++) {
			for (ix = 0; ix < Nx; ix++) {
				P_c[iy][ix] = P_n[iy][ix];
				U_c[iy][ix] = U_n[iy][ix];
				V_c[iy][ix] = V_n[iy][ix];
			}
		}
	}

	// Y-step
	for(ix = start; ix < Nx; ix += stride) {
		for (iy = 0; iy < Ny; iy++) {
			double Wpls_n, Wmns_n;
			double Z_l, Z_r;
			unsigned int j;

			Z_l = C[iy][ix] * Ro[iy][ix];
			Z_r = C[iy][ix] * Ro[iy][ix];

			omega(P_c, V_c, W, iy, ix, AXIS_Y);

			for (j = 0; j < S + 1; j++) {
				stencil_values_w2[j] = W[j];

				stencil_values_w1[j] = W[j + S + 1];
			}

			Wpls_n = interpolatedValueCIR(stencil_values_w2, LEFT, iy, ix);
			Wmns_n = interpolatedValueCIR(stencil_values_w1, RIGHT, iy, ix);

			P_n[iy][ix] = 2.0 * Z_r * Z_l * (Wpls_n - Wmns_n) / (Z_l + Z_r);
			V_n[iy][ix] = 2.0 * (Z_l * Wpls_n + Z_r * Wmns_n) / (Z_l + Z_r);
			//printf("%lf\t%lf\n", P_n[iy][ix], U_n[iy][ix]);
		}
	}

#pragma omp barrier

#pragma omp single
	{
		for (iy = 0; iy < Ny; iy++) {
			for (ix = 0; ix < Nx; ix++) {
				P_c[iy][ix] = P_n[iy][ix];
				U_c[iy][ix] = U_n[iy][ix];
				V_c[iy][ix] = V_n[iy][ix];
			}
		}
	}

	free(W);
}

double interpolatedValueCIR(double *u, int dir, int iy, int ix) {
	double a, b, x;
	if (dir == LEFT){
		a = (u[1] - u[0]) / hx;
		b = u[1];
		x = -C[iy][ix] * taux;
	}
	else {
		a = (u[1] - u[0]) / hx;
		b = u[0];
		x = C[iy][ix] *  taux;
	}
	double val = a * x + b;
	return val;
}

