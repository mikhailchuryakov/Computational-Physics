#define _USE_MATH_DEFINES

#include <stdio.h>
#include <stdlib.h>
#include <cmath>

#include "norma2d.h"

int grid = 100;
int Nx = grid, Ny = grid;

#define h (2.0 / (Nx))

void createArrays() {
	int i;

	W = (double **)malloc(Ny * sizeof(double *));
	A = (double **)malloc(Ny * sizeof(double *));

	for (i = 0; i < Ny; i++){
		W[i] = (double *)malloc(Nx * sizeof(double));
		A[i] = (double *)malloc(Nx * sizeof(double));
	}
	Hx = (double *)malloc(Nx * sizeof(double));
	Hy = (double *)malloc(Ny * sizeof(double));
}

void freeArrays() {
	int i;
	for (i = 0; i < Ny; i++) {
		free(W[i]);
		free(A[i]);
	}
	free(W);
	free(Hx);
	free(Hy);
	free(A);
}


int main() {
	char str[20];
	FILE *pFile;
	pFile = fopen("result.txt", "a");

	for (Nx = grid, Ny = grid; Nx <= 400 + 1; Nx = Nx * 2, Ny = Ny * 2) {
		createArrays();

		sprintf(str, "P%d.dat", (int) Nx);
		printf("%s\n", str);

		initializeArrays(str);
		solveAnalitic();

		fprintAnalitic();

		fprintf(pFile, "%d\t%g\t%g\n", (int) Nx, normL1(), normLINF());
		freeArrays();
	}

	fclose(pFile);

	return 0;
}

double normL1(){
	double val = 0;
	int iy, ix;

	for (iy = 0; iy < Ny; iy++) {
		for (ix = 0; ix < Nx; ix++) {
			val += fabs(W[iy][ix] - A[iy][ix]) * h * h;
		}
	}

	return val;
}

double normLINF(){
	double max = 0;
	int iy, ix;
	double val;

	for (iy = 0; iy < Ny; iy++) {
		for (ix = 0; ix < Nx; ix++) {
			//printf("%d\t%lf-%lf\t%lf\n", ind, W[ind], A[ind], W[ind] - A[ind]);
			if ((val = fabs(W[iy][ix] - A[iy][ix])) > max) {
				printf("%d\t%d\t%lf-%lf\t%lf\n", iy, ix, W[iy][ix], A[iy][ix], W[iy][ix] - A[iy][ix]);
				max = val;
			}
		}
	}
	return max;
}

#define eps 1E-6
void solveAnalitic(void){
	int iy, ix;
	double x;

	for (iy = 0; iy < Ny; iy++) {
		for (ix = 0; ix < Nx; ix++) {
			
			if (0.4 - Hx[ix] <= eps && 0.9 - Hx[ix] >= eps){
				A[iy][ix] = 1.0 * pow(sin(2.0 * M_PI *(Hx[ix] + 0.6)), 4.0);
			}
			else{
				A[iy][ix] = 0.0;
			}

			/*
			if (-0.8 - H[ind] <= eps && -0.3 - H[ind] >= eps){
			A[ind] = -1.0 / 5.0 * pow(sin(2.0 * M_PI *(H[ind] - 0.2)), 4.0);
			//A[ind] = 0;
			}
			else if (0.15 - H[ind] <= eps && 0.4 - H[ind] >= eps) {
			//printf("%lf\n", H[ind]);
			A[ind] = 2.0 / 2.5 * pow(sin(4 * M_PI *(H[ind] - 1.15)), 4.0);
			//A[ind] = 1.0 * pow(sin(4 * M_PI *(H[ind] - 1.15)), 4.0);
			//A[ind] = 0;
			}
			else {
			A[ind] = 0;
			}
			*/
			/*
			//A[ind] = -1.0 / 5.0 * exp(-100 * (x + 0.5)*(x + 0.5)) + 2.0 / 2.5 * exp(-400 * (x - 0.25)*(x - 0.25));
			if ((x >= -0.5 && x <= -0.2) || (x >= 0.2 && x <= 0.5)) {
				A[ind] = 0.6;
			}
			else if (x > -0.2 && x < 0.2) {
				A[ind] = 3.0;
			}
			else {
				A[ind] = 1.0;
			}
			*/
		}
	}
}

void initializeArrays(char* file1){
	int iy, ix;
	FILE *first;
	first = fopen(file1, "r");

	if (first != NULL) {
		for (iy = 0; iy < Ny; iy++) {
			for (ix = 0; ix < Nx; ix++) {
				char ahx[30];
				char ahy[30];
				char f[30];
				fscanf(first, "%s\t%s\t%s", ahx, ahy, f);
				if (iy == 0) Hx[ix] = atof(ahx);
				if (ix == 0) Hy[iy] = atof(ahy);
				W[iy][ix] = atof(f);
				/*
				if (-0.8 - H[ind] <= eps && -0.3 - H[ind] >= eps){
				W[ind] = atof(f);
				//W[ind] = 0;
				}
				else if (0.15 - H[ind] <= eps && 0.4 - H[ind] >= eps) {
				//W[ind] = atof(f);
				W[ind] = 0;
				}
				else {
				W[ind] = 0;
				}
				*/
				//printf("%.9f\t%.9f\t%.9f\n", Hx[ix], Hy[iy], W[iy][ix]);
			}
		}
	}
	else {
		printf("ERROR\n");
		exit(-1);
	}

	fclose(first);
}


void fprintAnalitic(){
	int iy, ix;
	FILE *pFile;
	pFile = fopen("analitic.dat", "w");
	for (iy = 0; iy < Ny; iy++) {
		for (ix = 0; ix < Nx; ix++) {
			fprintf(pFile, "%.9f\t%.9f\t%.9f\n", -1.0 + h * ix, -1.0 + h * iy, A[iy][ix]);
		}
	}
	fclose(pFile);
}