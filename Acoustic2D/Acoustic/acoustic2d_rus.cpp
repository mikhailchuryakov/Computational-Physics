#define _USE_MATH_DEFINES

#include <stdio.h>
#include <cmath>
#include <stdlib.h>
#include <omp.h>

#include "acoustic2d_ext.h"

void singleStepRus(void) {
	int iy, ix;
	int start = omp_get_thread_num();
	int stride = omp_get_num_threads();
	double stencil_values_w1[S + 3];
	double stencil_values_w2[S + 3];

	double *W = (double *)malloc(8 * sizeof(double));

	// X-step
	for (iy = start; iy < Ny; iy += stride) {
		for (ix = 0; ix < Nx; ix++) {
			double Wpls_n, Wmns_n;
			double Z_l, Z_r;
			unsigned int j;

			Z_l = C[iy][ix] * Ro[iy][ix];
			Z_r = C[iy][ix] * Ro[iy][ix];

			omega_Rus(P_c, U_c, W, iy, ix, AXIS_X);

			for (j = 0; j < S + 3; j++) {
				stencil_values_w2[j] = W[j];

				stencil_values_w1[j] = W[j + S + 3];
			}

			Wpls_n = interpolatedValueRus(stencil_values_w2, LEFT, iy, ix);
			Wmns_n = interpolatedValueRus(stencil_values_w1, RIGHT, iy, ix);

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

			omega_Rus(P_c, V_c, W, iy, ix, AXIS_Y);

			for (j = 0; j < S + 3; j++) {
				stencil_values_w2[j] = W[j];

				stencil_values_w1[j] = W[j + S + 3];
			}

			Wpls_n = interpolatedValueRus(stencil_values_w2, LEFT, iy, ix);
			Wmns_n = interpolatedValueRus(stencil_values_w1, RIGHT, iy, ix);

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

#define yidxL2 (iy == 0 ? (Ny - 2) : (iy == 1 ? (Ny - 1) : (iy - 2)))
#define yidxR2 (iy == (Ny - 1) ? 1 : (iy == (Nx - 2) ? 0 : (iy + 2)))
#define xidxL2 (ix == 0 ? (Nx - 2) : (ix == 1 ? (Nx - 1) : (ix - 2)))
#define xidxR2 (ix == (Nx - 1) ? 1 : (ix == (Nx - 2) ? 0 : (ix + 2)))
void omega_Rus(double **P_c, double **U_c, double *W, int iy, int ix, int axis){
	double Z_l, Z_r;
	Z_l = C[iy][ix] * Ro[iy][ix];
	Z_r = C[iy][ix] * Ro[iy][ix];
	if (axis == AXIS_X){
		W[0] = (P_c[iy][xidxL2] + Z_l * U_c[iy][xidxL2]) / (2 * Z_l);
		W[1] = (P_c[iy][xidxL] + Z_l * U_c[iy][xidxL]) / (2 * Z_l);
		W[2] = (P_c[iy][ix] + Z_l * U_c[iy][ix]) / (2 * Z_l);
		W[3] = (P_c[iy][xidxR] + Z_l * U_c[iy][xidxR]) / (2 * Z_l);

		W[4] = (-P_c[iy][xidxL] + Z_r * U_c[iy][xidxL]) / (2 * Z_r);
		W[5] = (-P_c[iy][ix] + Z_r * U_c[iy][ix]) / (2 * Z_r);
		W[6] = (-P_c[iy][xidxR] + Z_r * U_c[iy][xidxR]) / (2 * Z_r);
		W[7] = (-P_c[iy][xidxR2] + Z_r * U_c[iy][xidxR2]) / (2 * Z_r);
		//printf("%lf\t%lf\n", Wmns[0], Wmns[1]);
	}
	else {
		W[0] = (P_c[yidxL2][ix] + Z_l * U_c[yidxL2][ix]) / (2 * Z_l);
		W[1] = (P_c[yidxL][ix] + Z_l * U_c[yidxL][ix]) / (2 * Z_l);
		W[2] = (P_c[iy][ix] + Z_l * U_c[iy][ix]) / (2 * Z_l);
		W[3] = (P_c[yidxR][ix] + Z_l * U_c[yidxR][ix]) / (2 * Z_l);

		W[4] = (-P_c[yidxL][ix] + Z_r * U_c[yidxL][ix]) / (2 * Z_r);
		W[5] = (-P_c[iy][ix] + Z_r * U_c[iy][ix]) / (2 * Z_r);
		W[6] = (-P_c[yidxR][ix] + Z_r * U_c[yidxR][ix]) / (2 * Z_r);
		W[7] = (-P_c[yidxR2][ix] + Z_r * U_c[yidxR2][ix]) / (2 * Z_r);
		//printf("%lf\t%lf\n", Wmns[0], Wmns[1]);
	}
}

double interpolatedValueRus(double *u, int dir, int iy, int ix) {
	double a, b, c, d, x;
	//double val;
	if (dir == LEFT){
		a = (-u[0] + 3.0 * u[1] - 3.0 * u[2] + u[3]) / hx / hx / hx / 6.0;
		b = (u[1] - 2 * u[2] + u[3]) / hx / hx / 2.0;
		c = (u[0] - 6.0 * u[1] + 3.0 * u[2] + 2.0 * u[3]) / hx / 6.0;
		d = u[2];
		x = -C[iy][ix] * taux;
		
		//val = u[2] + C[iy][ix] * (u[1] + u[3]) / 2.0 + C[iy][ix] * C[iy][ix] * (u[1] + u[3] - 2.0*u[2]) / 2.0 + C[iy][ix] * (C[iy][ix] * C[iy][ix] - 1)*(u[0] - 3.0*u[1] + 3.0*u[2] - u[3])/6.0;
	}
	else {
		
		a = (-u[0] + 3.0 * u[1] - 3.0 * u[2] + u[3]) / hx / hx / hx / 6.0;
		b = (u[0] - 2 * u[1] + u[2]) / hx / hx / 2.0;
		c = -(2.0*u[0] + 3.0 * u[1] - 6.0 * u[2] + u[3]) / hx / 6.0;
		d = u[1];
		x = C[iy][ix] * taux;
		
		//val = u[1] + C[iy][ix] * (u[2] + u[0]) / 2.0 + C[iy][ix] * C[iy][ix] * (u[2] + u[0] - 2.0*u[1]) / 2.0 + C[iy][ix] * (C[iy][ix] * C[iy][ix] - 1)*(u[3] - 3.0*u[2] + 3.0*u[1] - u[0])/6.0;
	}
	double val = a * x * x * x + b * x * x + c * x + d;
	return val;
}