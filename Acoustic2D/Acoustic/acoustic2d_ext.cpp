#define _USE_MATH_DEFINES

#include <stdio.h>
#include <cmath>
#include <stdlib.h>
#include <omp.h>

#include "acoustic2d_ext.h"


#define M (2.0 * (Nx * 0.4 * hx / 2.0) / C[0][0] / taux)

double **P_c, **P_n, **U_c, **U_n, **V_c, **V_n;
double **Px_c, **Px_n, **Ux_c, **Ux_n, **Vx_c, **Vx_n;
double **Py_c, **Py_n, **Uy_c, **Uy_n, **Vy_c, **Vy_n;
double **Pxy_c, **Pxy_n, **Uxy_c, **Uxy_n, **Vxy_c, **Vxy_n;
double **C, **Ro; // Material

int main() {
	char str[20], strCIR[20];
	double time;
	int number_of_processes = 0;
	int del = 25;

	printf("start\n");
	for (Nx = grid, Ny = grid; Nx <= 200 + 1; Nx = Nx * 2, Ny = Ny * 2) {
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
					
					if ((step % del) == 0) {
						sprintf(str, "gif/P%d.dat", step / del + 1);
						sprintf(strCIR, "gif/U%d.dat", step / del + 1);
						debugGnuplotFileName(str, strCIR);
					}
					
				}
				singleStepRus();
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
	return -200 * (y)* exp(-100 * (x)*(x)) * exp(-100 * (y)*(y));
}

double func_xy(double x, double y) {
	return 40000 * (x)* (y)* exp(-100 * (x)*(x)) * exp(-100 * (y)*(y));
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
			V_c[iy][ix] = 0.0;
			V_n[iy][ix] = V_c[iy][ix];

			Px_c[iy][ix] = func_x(x, y);
			Ux_c[iy][ix] = Px_c[iy][ix] / (C[iy][ix] * Ro[iy][ix]);
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

