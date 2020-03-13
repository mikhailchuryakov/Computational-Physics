#define _USE_MATH_DEFINES

#include <stdio.h>
#include <cmath>
#include <stdlib.h>
#include <omp.h>

#include "acoustic2d_ext.h"

int Nx = grid, Ny = grid;

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

	free(W);
	free(dW);
	free(Wn);
	free(dWn);
	//printf("%lf\n", Pxy_n[20][20]);
}

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