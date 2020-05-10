#define _USE_MATH_DEFINES

#include <stdio.h>
#include <cmath>
#include <stdlib.h>
#include <omp.h>

#include "advection2d_ext.h"

int Nx = grid, Ny = grid;

void singleStepCIP(void) {
	int ix, iy;
	int start = omp_get_thread_num();
	int stride = omp_get_num_threads();
	double *W, *dW;
	double Wn, dWn;
	double *coefs2;
	// W[4] = {Wpls[0], Wpls[1]}
	W = (double *)malloc(2 * sizeof(double));
	// dW[4] = {dWpls[0], dWpls[1]}
	dW = (double *)malloc(2 * sizeof(double));
	coefs2 = (double *)malloc(5 * sizeof(double));

	// X-step
#pragma omp barrier
	for (iy = start; iy < Ny; iy += stride) {
		for (ix = 0; ix < Nx; ix++) {
			
			W[0] = P_c[iy][xidxL];
			W[1] = P_c[iy][ix];
			dW[0] = Px_c[iy][xidxL];
			dW[1] = Px_c[iy][ix];
			interpolatedCoefs(W, dW, LEFT, iy, ix, coefs2);
			Wn = interpolatedValue(coefs2, ZERO); // W_pls
			dWn = interpolatedValue(coefs2, FIRST); // dW_pls
			P_n[iy][ix] = Wn;
			Px_n[iy][ix] = dWn;

			W[0] = Py_c[iy][xidxL];
			W[1] = Py_c[iy][ix];
			dW[0] = Pxy_c[iy][xidxL];
			dW[1] = Pxy_c[iy][ix];
			interpolatedCoefs(W, dW, LEFT, iy, ix, coefs2);
			Wn = interpolatedValue(coefs2, ZERO); // W_pls
			dWn = interpolatedValue(coefs2, FIRST); // dW_pls
			Py_n[iy][ix] = Wn;
			Pxy_n[iy][ix] = dWn;
		}
	}
	update();

	// Y-step
#pragma omp barrier
	for (ix = start; ix < Nx; ix += stride) {
		for (iy = 0; iy < Ny; iy++) {
			W[0] = P_c[yidxL][ix];
			W[1] = P_c[iy][ix];
			dW[0] = Px_c[yidxL][ix];
			dW[1] = Px_c[iy][ix];
			interpolatedCoefs(W, dW, LEFT, iy, ix, coefs2);
			Wn = interpolatedValue(coefs2, ZERO); // W_pls
			dWn = interpolatedValue(coefs2, FIRST); // dW_pls
			P_n[iy][ix] = Wn;
			Px_n[iy][ix] = dWn;

			W[0] = Py_c[yidxL][ix];
			W[1] = Py_c[iy][ix];
			dW[0] = Pxy_c[yidxL][ix];
			dW[1] = Pxy_c[iy][ix];
			interpolatedCoefs(W, dW, LEFT, iy, ix, coefs2);
			Wn = interpolatedValue(coefs2, ZERO); // W_pls
			dWn = interpolatedValue(coefs2, FIRST); // dW_pls
			Py_n[iy][ix] = Wn;
			Pxy_n[iy][ix] = dWn;

		}
	}

	update();

	free(W);
	free(dW);
	free(coefs2);
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
		printf("loose");
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

void update() {
	int ix, iy;
	int start = omp_get_thread_num();
	int stride = omp_get_num_threads();

	for (iy = start; iy < Ny; iy += stride) {
		for (ix = 0; ix < Nx; ix++) {
			P_c[iy][ix] = P_n[iy][ix];
			Px_c[iy][ix] = Px_n[iy][ix];
			Py_c[iy][ix] = Py_n[iy][ix];
			Pxy_c[iy][ix] = Pxy_n[iy][ix];
		}
	}
}