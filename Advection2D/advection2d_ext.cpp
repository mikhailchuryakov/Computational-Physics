#define _USE_MATH_DEFINES

#include <stdio.h>
#include <cmath>
#include <stdlib.h>
#include <omp.h>

#include "advection2d_ext.h"


#define M (2.0 * (Nx * 0.4 * hx / 2.0) / C[0][0] / taux)

double **P_c, **P_n;
double **Px_c, **Px_n;
double **Py_c, **Py_n;
double **Pxy_c, **Pxy_n;
double **C; // Material

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
				printf("%d\n", number_of_processes);
			}
			for (step = 0; step < M; step++) {
#pragma omp single
				{
					if (step % 100 == 0) printf("%d\n", step);
					
					if ((step % del) == 0) {
						sprintf(str, "gif/P%d.dat", step / del + 1);
						debugGnuplotFileName(str);
					}
					
				}
				singleStepCIP();
			}
		}
		
		sprintf(str, "gif/P%d.dat", (int)M / del + 1);
		
		sprintf(str, "P%d.dat", Nx);
		debugGnuplotFileName(str);

		freeArrays();
		time = (omp_get_wtime() - time);
		printf("Concurrent OpenMP calculation with %d threads took %f seconds.\n", number_of_processes, time);
	}
	return 0;
}

double func(double x, double y) {
	
	if (-0.4 - x <= eps && 0.1 - x >= eps && -0.4 - y <= eps && 0.1 - y >= eps){
		return 1.0 * pow(sin(2.0 * M_PI *(x - 0.6)), 4.0);
	}
	else {
		return 0;
	}
	
	//return exp(-400 * (x)*(x)) * exp(-400 * (y)*(y));
}

double func_x(double x, double y) {
	
	if (-0.4 - x <= eps && 0.1 - x >= eps && -0.4 - y <= eps && 0.1 - y >= eps){
		return 8.0 * M_PI * pow(sin(2.0 * M_PI *(x - 0.6)), 3.0) * cos(2.0 * M_PI *(x - 0.6));
	}
	else {
		return 0;
	}
	
	//return -800 * (x) * exp(-400 * (x)*(x)) * exp(-400 * (y)*(y));
}

double func_y(double x, double y) {
	return -800 * (y)* exp(-400 * (x)*(x)) * exp(-400 * (y)*(y));
}

double func_xy(double x, double y) {
	return 640000 * (x)* (y)* exp(-400 * (x)*(x)) * exp(-400 * (y)*(y));
}

void init(void) {
	int ix, iy;
	for (iy = 0; iy < Ny; iy++) {
		for (ix = 0; ix < Nx; ix++) {
			double x = -(LENGTHx / 2.0) + ix * hx;
			double y = -(LENGTHy / 2.0) + iy * hy;
			
			C[iy][ix] = 1.0;

			// nodes
			P_c[iy][ix] = func(x, y);
			Px_c[iy][ix] = func_x(x, y);
			Py_c[iy][ix] = 0.0;//func_y(x, y);
			Pxy_c[iy][ix] = 0.0;//func_xy(x, y);

		}
	}
}


void debugGnuplotFileName(char *filename1){
	FILE *pFileP;
	pFileP = fopen(filename1, "w");
	int ix, iy;
	for (iy = 0; iy < Ny; iy++) {
		for (ix = 0; ix < Nx; ix++) {
			if (abs(P_c[iy][ix]) > 1E-5) 
				fprintf(pFileP, "%.12lf\t%.12lf\t%.12lf\n", -(LENGTHx / 2.0) + hx * ix, -(LENGTHx / 2.0) + hy * iy, P_c[iy][ix]);
		}
	}
	fclose(pFileP);
}

void createArrays(){
	int i;

	C = (double **)malloc(Ny * sizeof(double *));

	P_c = (double **)malloc(Ny * sizeof(double *));
	P_n = (double **)malloc(Ny * sizeof(double *));
	
	Px_c = (double **)malloc(Ny * sizeof(double *));
	Px_n = (double **)malloc(Ny * sizeof(double *));
	
	Py_c = (double **)malloc(Ny * sizeof(double *));
	Py_n = (double **)malloc(Ny * sizeof(double *));
	
	Pxy_c = (double **)malloc(Ny * sizeof(double *));
	Pxy_n = (double **)malloc(Ny * sizeof(double *));

	for (i = 0; i < Ny; i++){
		C[i] = (double *)malloc(Nx * sizeof(double));

		P_c[i] = (double *)malloc(Nx * sizeof(double));
		P_n[i] = (double *)malloc(Nx * sizeof(double));

		Px_c[i] = (double *)malloc(Nx * sizeof(double));
		Px_n[i] = (double *)malloc(Nx * sizeof(double));

		Py_c[i] = (double *)malloc(Nx * sizeof(double));
		Py_n[i] = (double *)malloc(Nx * sizeof(double));

		Pxy_c[i] = (double *)malloc(Nx * sizeof(double));
		Pxy_n[i] = (double *)malloc(Nx * sizeof(double));
	}
	
}

void freeArrays(){
	int i;
	for (i = 0; i < Ny; i++) {
		free(C[i]);

		free(Px_c[i]);
		free(Px_n[i]);

		free(Py_c[i]);
		free(Py_n[i]);

		free(Pxy_c[i]);
		free(Pxy_n[i]);
	}
	free(C);

	free(Px_c);
	free(Px_n);

	free(Py_c);
	free(Py_n);

	free(Pxy_c);
	free(Pxy_n);
}

