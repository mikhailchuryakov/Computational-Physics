#define _USE_MATH_DEFINES

#include <stdio.h>
#include <cmath>
#include <stdlib.h>

#include "acoustic_ext.h"

int N;

#define LENGTH 2.0

#define k 0.2
#define h (LENGTH / (N))
#define tau (k * h / 1.0)

// template
#define S 1
//int stencil_l[S + 1] = { -1, 0 };
//int stencil_r[S + 1] = { 0, 1 };
double stencil_values_w1[S + 1];
double stencil_values_w2[S + 1];
double stencil_Dvalues_w1[S + 1];
double stencil_Dvalues_w2[S + 1];
#define LEFT 0
#define RIGHT 1
//

#define M (((N * 0.45 * h / 2.0) / C[0] + (N * 0.5 * h / 2.0) / C[N - 1]) / tau)
//#define M 0.0


#define eps 1E-6
double func(double x) {
	/*
	if (-0.7 - x <= eps && -0.2 - x >= eps){
		return 1.0 * pow(sin(2.0 * M_PI *(x - 0.3)), 4.0);
	}
	else {
		return 0;
	}
	*/
	//return exp(-100 * (x + 0.5)*(x + 0.5));
	return 0.5 * exp(-400 * (x)*(x));
}

double dfunc(double x) {
	/*
	if (-0.7 - x <= eps && -0.2 - x >= eps){
		return 8.0 * M_PI * pow(sin(2.0 * M_PI *(x - 0.3)), 3.0) * cos(2.0 * M_PI *(x - 0.3));
	}
	else {
		return 0;
	}
	*/
	//return -200 * (x + 0.5) * exp(-100 * (x + 0.5)*(x + 0.5));
	return -400 * (x) * exp(-400 * (x)*(x));
}

int main() {
	unsigned int step, del;
	char str[20], strCIR[20];
	del = 50;

	for (N = 800; N <= 801; N = N * 2) {
		createArrays();
		init();

		printf("%d\t%lf\t%f\t%lf\t%lf\n", N, h, 1.0 / h, tau, M);
		for (step = 0; step < M; step++) {
			
			if ((step % del) == 0) {
				sprintf(str, "gif/P%d.dat", step / del + 1);
				sprintf(strCIR, "gif/U%d.dat", step / del + 1);
				debugGnuplotFileName(str, strCIR);
			}
			
			singleStep();
		}
		
		sprintf(str, "gif/P%d.dat", (int)M / del + 1);
		sprintf(strCIR, "gif/U%d.dat", (int)M / del + 1);
		debugGnuplotFileName(str, strCIR);
		
		sprintf(str, "P%d.dat", N);
		sprintf(strCIR, "U%d.dat", N);
		debugGnuplotFileName(str, strCIR);

		freeArrays();
	}
	return 0;
}

void init(void) {
	int ind;
	for (ind = 0; ind < N; ind++) {
		double x = -(LENGTH / 2.0) + ind * h;
		/*
		// substance
		if (x <= 0.0 - 0.5*h) {
			C[ind] = 1.0;
			Ro[ind] = 1.0;
		}
		else {
			C[ind] = 0.5;
			Ro[ind] = 1.0;
		}
		*/
		/*
		if (x >= -0.2 && x < 0.2 - h / 2) {
			C[ind] = 2.0;
			Ro[ind] = 1.0;
		}
		else {
			C[ind] = 1.0;
			Ro[ind] = 1.0;
		}
		*/
		
		if ((x >= -0.5 && x < -0.2 - h / 2) || (x >= 0.2 && x < 0.5 - h / 2)) {
			C[ind] = 0.6;
			Ro[ind] = 1.0;
		}
		else if (x > -0.2 && x < 0.2 - h / 2) {
			C[ind] = 3.0;
			Ro[ind] = 1.0;
		}
		else {
			C[ind] = 1.0;
			Ro[ind] = 1.0;
		}
		

		// nodes
		P_c[ind] = func(x);
		U_c[ind] = P_c[ind] / (C[ind] * Ro[ind]);

		dP_c[ind] = dfunc(x) / Ro[ind == 0 ? N - 1 : ind - 1];
		dU_c[ind] = dP_c[ind] * C[ind] * Ro[ind];
		//dP_c[ind] = dfunc(x);
		//dU_c[ind] = dfunc(x);
		//printf("%lf\t%lf\t%lf\n", x, dP_c[ind], dU_c[ind]);
		if (fabs(x -(0.5)) < eps){
			printf("%lf\t%lf\t%lf\n", x, C[ind == 0 ? N - 1 : ind - 1], C[ind]);
		}
	}
}

#define ZERO 0
#define FIRST 1
void singleStep(void) {
	int ind;

	for (ind = 0; ind < N; ind++) {
		// coefs[5] = { x, a, b, c, d }
		double coefs1[5], coefs2[5];
		double Wpls_n, Wmns_n, dWpls_n, dWmns_n;
		double Z_l, Z_r;
		double Kl, Kr;

		if (ind == 0){
			Z_l = C[N - 1] * Ro[N - 1];
		}
		else {
			Z_l = C[ind - 1] * Ro[ind - 1];
		}
		Z_r = C[ind] * Ro[ind];

		toInvariants(ind, Z_l, Z_r);
		toDInvariants(ind, Z_l, Z_r);

		fillStencilValues(ind);

		interpolatedCoefs(stencil_values_w1, stencil_Dvalues_w1, RIGHT, ind, coefs1);
		interpolatedCoefs(stencil_values_w2, stencil_Dvalues_w2, LEFT, ind, coefs2);

		Wmns_n = interpolatedValue(coefs1, ZERO);
		dWmns_n = interpolatedValue(coefs1, FIRST);

		Wpls_n = interpolatedValue(coefs2, ZERO);
		dWpls_n = interpolatedValue(coefs2, FIRST);


		P_n[ind] = 2.0 * Z_r * Z_l * (Wpls_n - Wmns_n) / (Z_l + Z_r);
		U_n[ind] = 2.0 * (Z_l * Wpls_n + Z_r * Wmns_n) / (Z_l + Z_r);

		Kl = C[ind == 0 ? N - 1 : ind - 1] * Z_l;
		Kr = C[ind] * Z_r;

		/*
		if (ind == 0){
			dP_n[ind] = (Z_l * Z_l * dWpls_n - C[ind] * Ro[N - 1] * Z_r * dWmns_n) / (Z_l + Z_r) / Ro[N - 1] + (C[N - 1] * Ro[ind] * Z_l * dWpls_n - Z_r * Z_r * dWmns_n) / (Z_l + Z_r) / Ro[ind];
			dU_n[ind] = (Z_l * Z_r * dWpls_n + C[ind] * Ro[N - 1] * Z_r * dWmns_n) / (Z_l*Z_l + Z_r*Z_l) * Kl + (C[N - 1] * Ro[ind] * Z_l * dWpls_n + Z_l * Z_r * dWmns_n) / (Z_r*Z_l + Z_r*Z_r) * Kr;
		}
		else {
			dP_n[ind] = (Z_l * Z_l * dWpls_n - C[ind] * Ro[ind - 1] * Z_r * dWmns_n) / (Z_l + Z_r) / Ro[ind - 1] + (C[ind - 1] * Ro[ind] * Z_l * dWpls_n - Z_r * Z_r * dWmns_n) / (Z_l + Z_r) / Ro[ind];
			dU_n[ind] = (Z_l * Z_r * dWpls_n + C[ind] * Ro[ind - 1] * Z_r * dWmns_n) / (Z_l*Z_l + Z_r*Z_l) * Kl + (C[ind - 1] * Ro[ind] * Z_l * dWpls_n + Z_l * Z_r * dWmns_n) / (Z_r*Z_l + Z_r*Z_r) * Kr;
		}
		*/

		if (ind == 0){
			dP_n[ind] = 2.0 * (C[N - 1] * Z_l * dWpls_n - C[ind] * Z_r * dWmns_n) / (Z_l + Z_r);
			dU_n[ind] = 2.0 * (C[N - 1] * Z_l * Z_r * dWpls_n + C[ind] * Z_l * Z_r * dWmns_n) / (Z_l + Z_r);
		}
		else {
			dP_n[ind] = 2.0 * (C[ind - 1] * Z_l * dWpls_n - C[ind] * Z_r * dWmns_n) / (Z_l + Z_r);
			dU_n[ind] = 2.0 * (C[ind - 1] * Z_l * Z_r * dWpls_n + C[ind] * Z_l * Z_r * dWmns_n) / (Z_l + Z_r);
		}
		
	}

	for (ind = 0; ind < N; ind++) {
		P_c[ind] = P_n[ind];
		dP_c[ind] = dP_n[ind];
		U_c[ind] = U_n[ind];
		dU_c[ind] = dU_n[ind];
	}
}

void toInvariants(int ind, double Z_l, double Z_r){
	if (ind == 0) {
		Wpls[0] = (P_c[N - 1] + Z_l * U_c[N - 1]) / (2 * Z_l);
	}
	else{
		Wpls[0] = (P_c[ind - 1] + Z_l * U_c[ind - 1]) / (2 * Z_l);
	}
	Wpls[1] = (P_c[ind] + Z_l * U_c[ind]) / (2 * Z_l);
	
	Wmns[0] = (-P_c[ind] + Z_r * U_c[ind]) / (2 * Z_r);
	if (ind == N - 1) {
		Wmns[1] = (-P_c[0] + Z_r * U_c[0]) / (2 * Z_r);
	}
	else {
		Wmns[1] = (-P_c[ind + 1] + Z_r * U_c[ind + 1]) / (2 * Z_r);
	}
}


void toDInvariants(int ind, double Z_l, double Z_r){
	double Kl, Kr;
	Kl = C[ind == 0 ? N - 1 : ind - 1] * C[ind == 0 ? N - 1 : ind - 1] * Ro[ind == 0 ? N - 1 : ind - 1];
	Kr = C[ind] * C[ind] * Ro[ind];

	if (ind == 0) {
		dWpls[0] = (dP_c[N - 1] * Ro[N - 1] + Z_l * dU_c[N - 1] / Kl) / (2.0 * Z_l);
	}
	else {
		dWpls[0] = (dP_c[ind - 1] * Ro[ind - 1] + Z_l * dU_c[ind - 1] / Kl) / (2.0 * Z_l);
	}
	dWpls[1] = (dP_c[ind] * Ro[ind == 0 ? N - 1 : ind - 1] + Z_l * dU_c[ind] / Kl) / (2.0 * Z_l);

	dWmns[0] = (-dP_c[ind] * Ro[ind] + Z_r * dU_c[ind] / Kr) / (2.0 * Z_r);
	if (ind == N - 1) {
		dWmns[1] = (-dP_c[0] * Ro[ind] + Z_r * dU_c[0] / Kr) / (2.0 * Z_r);
	}
	else{
		dWmns[1] = (-dP_c[ind + 1] * Ro[ind] + Z_r * dU_c[ind + 1] / Kr) / (2.0 * Z_r);
	}
}

/*
void toDInvariants(int ind, double Z_l, double Z_r){
	double Kl, Kr;
	Kl = C[ind == 0 ? N - 1 : ind - 1] * C[ind == 0 ? N - 1 : ind - 1] * Ro[ind == 0 ? N - 1 : ind - 1];
	Kr = C[ind] * C[ind] * Ro[ind];

	if (ind == 0) {
		dWpls[0] = (dP_c[N - 1] + Z_l * dU_c[N - 1]) / (2.0 * Z_l);
	}
	else {
		dWpls[0] = (dP_c[ind - 1] + Z_l * dU_c[ind - 1]) / (2.0 * Z_l);
	}
	dWpls[1] = (dP_c[ind] + Z_l * dU_c[ind]) / (2.0 * Z_l);

	dWmns[0] = (-dP_c[ind] + Z_r * dU_c[ind]) / (2.0 * Z_r);
	if (ind == N - 1) {
		dWmns[1] = (-dP_c[0] + Z_r * dU_c[0]) / (2.0 * Z_r);
	}
	else{
		dWmns[1] = (-dP_c[ind + 1] + Z_r * dU_c[ind + 1]) / (2.0 * Z_r);
	}
}
*/
void fillStencilValues(int ind) {
	unsigned int j;
	for (j = 0; j < S + 1; j++) {
		stencil_values_w2[j] = Wpls[j];
		stencil_Dvalues_w2[j] = dWpls[j];

		stencil_values_w1[j] = Wmns[j];
		stencil_Dvalues_w1[j] = dWmns[j];
	}
}

void interpolatedCoefs(double *u, double *du, int dir, int ind, double *coefs) {

	if (dir == LEFT) {
		coefs[1] = (du[1] + du[0]) / h / h - 2.0 * (u[1] - u[0]) / h / h / h;
		coefs[2] = (2.0 * du[1] + du[0]) / h - 3.0 * (u[1] - u[0]) / h / h;
		coefs[3] = du[1];
		coefs[4] = u[1];
		coefs[0] = -C[ind == 0 ? N - 1 : ind - 1] * tau;
	}
	else {
		coefs[1] = (du[1] + du[0]) / h / h - 2.0 * (u[1] - u[0]) / h / h / h;
		coefs[2] = -(2.0 * du[0] + du[1]) / h + 3.0 * (u[1] - u[0]) / h / h;
		coefs[3] = du[0];
		coefs[4] = u[0];
		coefs[0] = C[ind] * tau;
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
	pFileU = fopen(filename2, "w");
	int ind;
	for (ind = 0; ind < N; ind++){
		fprintf(pFileP, "%.12lf\t%.12lf\n", -(LENGTH / 2.0) + h * ind, P_c[ind]);
	}
	for (ind = 0; ind < N; ind++){
		fprintf(pFileU, "%.12lf\t%.12f\n", -(LENGTH / 2.0) + h * ind, U_c[ind]);
	}
	fclose(pFileP);
	fclose(pFileU);
}

void createArrays(){
	C = (double *)malloc(N * sizeof(double));
	Ro = (double *)malloc(N * sizeof(double));
	P_c = (double *)malloc(N * sizeof(double));
	P_n = (double *)malloc(N * sizeof(double));
	U_c = (double *)malloc(N * sizeof(double));
	U_n = (double *)malloc(N * sizeof(double));
	dP_c = (double *)malloc(N * sizeof(double));
	dP_n = (double *)malloc(N * sizeof(double));
	dU_c = (double *)malloc(N * sizeof(double));
	dU_n = (double *)malloc(N * sizeof(double));
}

void freeArrays(){
	free(C);
	free(Ro);
	free(P_c);
	free(U_c);
	free(P_n);
	free(U_n);
	free(dP_c);
	free(dU_c);
	free(dP_n);
	free(dU_n);
}



// CIR
void singleStepCIR(void) {
	int ind = 0;

	for (ind = 0; ind < N; ind++) {
		double Wpls_n, Wmns_n;
		double Z_l, Z_r;

		if (ind == 0){
			Z_l = C[N - 1] * Ro[N - 1];
		}
		else {
			Z_l = C[ind - 1] * Ro[ind - 1];
		}
		Z_r = C[ind] * Ro[ind];

		toInvariants(ind, Z_l, Z_r);
		
		fillStencilValues(ind);

		Wmns_n = interpolatedValueCIR(stencil_values_w1, RIGHT, ind);
		Wpls_n = interpolatedValueCIR(stencil_values_w2, LEFT, ind);

		P_n[ind] = 2.0 * Z_r * Z_l * (Wpls_n - Wmns_n) / (Z_l + Z_r);
		U_n[ind] = 2.0 * (Z_l * Wpls_n + Z_r * Wmns_n) / (Z_l + Z_r);
	}

	for (ind = 0; ind < N; ind++){
		P_c[ind] = P_n[ind];
		U_c[ind] = U_n[ind];
	}
}

double interpolatedValueCIR(double *u, int dir, int ind) {
	double a, b, x;
	if (dir == LEFT){
		a = (u[1] - u[0]) / h;
		b = u[1];
		x = -C[ind == 0 ? N - 1 : ind - 1] * tau;
	}
	else {
		a = (u[1] - u[0]) / h;
		b = u[0];
		x = C[ind] *  tau;
	}
	double val = a * x + b;
	return val;
}

