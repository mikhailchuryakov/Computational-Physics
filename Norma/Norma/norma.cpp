#define _USE_MATH_DEFINES

#include <stdio.h>
#include <stdlib.h>
#include <cmath>

#include "norma.h"

double N;

#define h (2.0 / (N))

void createArrays() {
	W = (double *) malloc(N * sizeof(double));
	H = (double *) malloc(N * sizeof(double));
	A = (double *) malloc(N * sizeof(double));
}

void freeArrays() {
	free(W);
	free(H);
	free(A);
}


int main() {
	char str[20];
	FILE *pFile;
	pFile = fopen("result.txt", "a");

	for (N = 100; N <= 3201; N = N * 2) {
		createArrays();

		sprintf(str, "P%d.dat", (int) N);
		printf("%s\n", str);

		initializeArrays(str);
		solveAnalitic();

		fprintAnalitic();

		fprintf(pFile, "%d\t%g\t%g\n", (int) N, normL1(), normLINF());
		freeArrays();
	}

	fclose(pFile);

	return 0;
}

double normL1(){
	double val = 0;
	int ind;

	for (ind = 0; ind < N; ind++){
		val += fabs(W[ind] - A[ind]) * h;
	}

	return val;
}

double normLINF(){
	double max = 0;
	int ind;
	double val;

	for (ind = 0; ind < N; ind++) {
		//printf("%d\t%lf-%lf\t%lf\n", ind, W[ind], A[ind], W[ind] - A[ind]);
		if ((val = fabs(W[ind] - A[ind])) > max) {
			printf("%d\t%lf-%lf\t%lf\n", ind, W[ind], A[ind], W[ind]-A[ind]);
			max = val;
		}
	}
	return max;
}

#define eps 1E-6
void solveAnalitic(void){
	int ind;
	double x;

	for (ind = 0; ind < N; ind++) {
		/*
		//if (-0.7 - H[ind] <= eps && -0.2 - H[ind] >= eps){
			//A[ind] = 1.0 * pow(sin(2.0 * M_PI *(H[ind] - 0.3)), 4.0);
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
		x = -1.0 + ind * h;
		//A[ind] = -1.0 / 5.0 * exp(-100 * (x + 0.5)*(x + 0.5)) + 2.0 / 2.5 * exp(-400 * (x - 0.25)*(x - 0.25));
		/*
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
		A[ind] = 0.5088899240798759 *exp(-400.0 * pow(-0.5 + x, 2.0)) - 0.008889568928092085 * exp(-400.0 * pow(-0.1125 + x, 2.0)) -
			0.008889562612751877 * exp(-400.0 * pow(0.1125 + x, 2.0)) + 0.5088892074609678 * exp(-400.0 * pow(0.5 + x, 2.0));
	}
}


/*void solveAnalitic(void){
	int ind;

	for (ind = 0; ind < N; ind++) {
		if ((ind >= N / 4) && (ind <= 3 * N / 4)) {
			A[ind] = 1.0 / (2.0 * 1);
		}
		else {
			A[ind] = 0.0;
		}
	}
}*/

void initializeArrays(char* file1){
	int ind = 0;
	FILE *first;
	first = fopen(file1, "r");

	if (first != NULL) {
		for (ind = 0; ind < N; ind++) {
			char ah[30];
			char f[30];
			fscanf(first, "%s\t%s", ah, f);
			H[ind] = atof(ah);
			W[ind] = atof(f);
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
		}
	}
	else {
		printf("ERROR\n");
		exit(-1);
	}

	fclose(first);
}


void fprintAnalitic(){
	FILE *pFile;
	pFile = fopen("analitic.dat", "w");
	unsigned int ind;
	for (ind = 0; ind < N; ind++){
		fprintf(pFile, "%.9f\t%.9f\n", -1.0 + h * ind, A[ind]);
	}
	fclose(pFile);
}