#define _USE_MATH_DEFINES

#include <stdio.h>
#include <cmath>
#include <stdlib.h>
#include <omp.h>

#include "acoustic2d_ext.h"


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
	for (ix = start; ix < Nx; ix += stride) {
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
		x = C[iy][ix] * taux;
	}
	double val = a * x + b;
	return val;
}