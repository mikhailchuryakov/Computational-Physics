#define _USE_MATH_DEFINES

#include <stdio.h>
#include <cmath>
#include <stdlib.h>

#include "dorovskoy_ext.h"

int N;

double tau = 3.0 * K * rho_0 * (rho_f - rho_s)*rho_s + rho_0*rho_0*rho_s * (4.0 * mu + 3.0 * alpha*rho_s*(rho_f + rho_s));
double psi = rho_0*rho_0*rho_s*rho_s * (-12.0 * rho_f*rho_s *(-3.0 * K * K + alpha * rho_0 * (4.0 * mu * rho_0 + 3.0 * K * (rho_f + rho_s))) +
	pow(3.0 * K * (rho_f - rho_s) + rho_0 * (4.0 * mu + 3.0 * alpha * rho_s * (rho_f + rho_s)), 2.0));
double betta = rho_0 * rho_0 * rho_s * rho_s * (9.0 * K * K * (rho_f + rho_s) * (rho_f + rho_s) +
	rho_0 * rho_0 * (9.0 * alpha * alpha * rho_f * rho_f * rho_s * rho_s +
	6.0 * alpha * rho_f * rho_s * (-4.0 * mu + 3.0 * alpha * rho_s * rho_s) +
	pow(4.0 * mu + 3.0 * alpha * rho_s * rho_s, 2.0)) -
	6.0 * K * rho_0 * (4.0 * mu * rho_s + 3.0 * alpha * rho_f * rho_f * rho_s +
	3.0 * alpha * rho_s * rho_s * rho_s + rho_f * (-4.0 * mu + 6.0 * alpha * rho_s * rho_s)));
double gamma = tau;
double omega = 3.0 * K * rho_0 * (rho_f + rho_s)*rho_s + rho_0*rho_0*rho_s * (4.0 * mu - 3.0 * alpha*rho_s*(rho_f + rho_s));
double delta = 6 * K * K * rho_s - 3.0 * K * alpha * rho_0 * rho_s * (rho_f + 3.0 * rho_s) +
	alpha * rho_0 * rho_0 * rho_s * (-4.0 * mu + 3.0 * alpha * rho_s * (rho_f + rho_s));
double epsilon = 6.0 * rho_0 * rho_0 * rho_s * (-K + alpha * rho_0 * rho_s);
double phi = sqrt(6.0) * (K - alpha * rho_0 * rho_s) / rho_0;

double lambda2 = sqrt(tau - sqrt(psi)) / (sqrt(6.0) * rho_0 * rho_s);
double lambda1 = -lambda2;
double lambda4 = sqrt(tau + sqrt(psi)) / (sqrt(6.0) * rho_0 * rho_s);
double lambda3 = -lambda4;

double a1 = sqrt(-sqrt(betta) + gamma) * phi / (-alpha * sqrt(betta) + delta);
double a2 = sqrt(sqrt(betta) + gamma) * phi / (alpha * sqrt(betta) + delta);
double b1 = sqrt(6.0) * rho_s / (sqrt(-sqrt(betta) + gamma));
double b2 = sqrt(6.0) * rho_s / (sqrt(sqrt(betta) + gamma));
double c1 = (-sqrt(betta) + omega) / epsilon;
double c2 = (sqrt(betta) + omega) / epsilon;

double a = (1.0 + omega / sqrt(betta)) / 4.0;
double am = (1.0 - omega / sqrt(betta)) / 4.0;
double b = epsilon / (4.0 * sqrt(betta));
double c = sqrt(gamma - sqrt(betta)) * (delta * delta - alpha*alpha*betta) / (4.0*sqrt(betta)*phi*(alpha*gamma - delta));
double cm = sqrt(gamma + sqrt(betta)) * (delta * delta - alpha*alpha*betta) / (4.0*sqrt(betta)*phi*(alpha*gamma - delta));
double d = sqrt(gamma - sqrt(betta)) * (gamma + sqrt(betta)) * (alpha * sqrt(betta) - delta) / (4.0 * sqrt(6.0)*sqrt(betta)*rho_s*(alpha*gamma - delta));
double dm = sqrt(gamma + sqrt(betta)) * (sqrt(betta) - gamma) * (alpha * sqrt(betta) + delta) / (4.0 * sqrt(6.0)*sqrt(betta)*rho_s*(alpha*gamma - delta));

#define LENGTH 4000.0

#define k 0.2
#define h (LENGTH / (N))
#define tauT (k * h / 2000.0)

// template
#define S 1
//int stencil_l[S + 1] = { -1, 0 };
//int stencil_r[S + 1] = { 0, 1 };
double stencil_values_w1[S + 1];
double stencil_values_w2[S + 1];
double stencil_values_w3[S + 1];
double stencil_values_w4[S + 1];
double stencil_Dvalues_w1[S + 1];
double stencil_Dvalues_w2[S + 1];
double stencil_Dvalues_w3[S + 1];
double stencil_Dvalues_w4[S + 1];
#define LEFT 0
#define RIGHT 1
//

#define M (((N * 0.45 * h / 2.0) / 450.0 + (N * 0.5 * h / 2.0) / 2000.0) / tauT)
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
	return 1 * exp(-400.0 * (x)*(x));
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
	return -800.0 * x * exp(-400.0 * (x)*(x));
}

int main() {
	unsigned int step, del;
	char str[20];
	del = 25;

	printf("it started\n");
	printf("%lf\n", tau);
	printf("%lf\n", sqrt(psi));
	printf("%lf\n", betta);
	printf("%lf\n", gamma);
	printf("%lf\n", omega);
	printf("%lf\n", delta);
	printf("%lf\n", epsilon);
	printf("%lf\n", phi);
	printf("\n");
	printf("%lf\n", lambda1);
	printf("%lf\n", lambda2);
	printf("%lf\n", lambda3);
	printf("%lf\n", lambda4);
	printf("\n");
	printf("%.12lf\n", a1);
	printf("%.12lf\n", a2);
	printf("%.12lf\n", b1);
	printf("%.12lf\n", b2);
	printf("%lf\n", c1);
	printf("%lf\n", c2);
	printf("\n");
	printf("%lf\n", a);
	printf("%lf\n", am);
	printf("%lf\n", b);
	printf("%lf\n", d);
	printf("%lf\n", dm);
	printf("%lf\n", c);
	printf("%lf\n", cm);
	printf("\n");

	for (N = 800; N <= 801; N = N * 2) {
		createArrays();
		init();

		printf("%d\t%lf\t%f\t%lf\t%lf\n", N, h, 1.0 / h, tauT, M);
		for (step = 0; step < M; step++) {
			
			if ((step % del) == 0) {
				sprintf(str, "gif/P%d.dat", step / del + 1);
				debugGnuplotFileName(str);
			}
			
			singleStep();
		}
		
		sprintf(str, "gif/P%d.dat", (int)M / del + 1);
		debugGnuplotFileName(str);
		
		sprintf(str, "P%d.dat", N);
		debugGnuplotFileName(str);

		freeArrays();
	}
	return 0;
}

void init(void) {
	int ind;
	for (ind = 0; ind < N; ind++) {
		double x = -(LENGTH / 2.0) + ind * h;

		V_s_c[ind] = func(x);
		V_f_c[ind] = func(x);
		H_c[ind] = func(x);
		P_c[ind] = func(x);

		dV_s_c[ind] = dfunc(x);
		dV_f_c[ind] = dfunc(x);
		dH_c[ind] = dfunc(x);
		dP_c[ind] = dfunc(x);

	}
}

#define ZERO 0
#define FIRST 1
void singleStep(void) {
	int ind;

	for (ind = 0; ind < N; ind++) {
		// coefs[5] = { x, a, b, c, d }
		double coefs1[5], coefs2[5], coefs3[5], coefs4[5];
		double W1n, W2n, W3n, W4n;
		double dW1n, dW2n, dW3n, dW4n;

		toInvariants(ind);
		toDInvariants(ind);

		fillStencilValues(ind);

		interpolatedCoefs(stencil_values_w1, stencil_Dvalues_w1, lambda1, LEFT, ind, coefs1);
		interpolatedCoefs(stencil_values_w2, stencil_Dvalues_w2, lambda2, RIGHT, ind, coefs2);
		interpolatedCoefs(stencil_values_w3, stencil_Dvalues_w3, lambda3, LEFT, ind, coefs3);
		interpolatedCoefs(stencil_values_w4, stencil_Dvalues_w4, lambda4, RIGHT, ind, coefs4);

		W1n = interpolatedValue(coefs1, ZERO);
		dW1n = interpolatedValue(coefs1, FIRST);
		W2n = interpolatedValue(coefs2, ZERO);
		dW2n = interpolatedValue(coefs2, FIRST);
		W3n = interpolatedValue(coefs3, ZERO);
		dW3n = interpolatedValue(coefs3, FIRST);
		W4n = interpolatedValue(coefs4, ZERO);
		dW4n = interpolatedValue(coefs4, FIRST);

		V_s_n[ind] = a1 * W1n - a1 * W2n + a2 * W3n - a2 * W4n;
		V_f_n[ind] = -b1 * W1n + b1 * W2n - b2 * W3n + b2 * W4n;
		H_n[ind] = c1 * W1n + c1 * W2n + c2 * W3n + c2 * W4n;
		P_n[ind] = W1n + W2n + W3n + W4n;
		
		dV_s_n[ind] = a1 * dW1n - a1 * dW2n + a2 * dW3n - a2 * dW4n;
		dV_f_n[ind] = -b1 * dW1n + b1 * dW2n - b2 * dW3n + b2 * dW4n;
		dH_n[ind] = c1 * dW1n + c1 * dW2n + c2 * dW3n + c2 * dW4n;
		dP_n[ind] = dW1n + dW2n + dW3n + dW4n;
	}

	for (ind = 0; ind < N; ind++) {
		V_s_c[ind] = V_s_n[ind];
		V_f_c[ind] = V_f_n[ind];
		H_c[ind] = H_n[ind];
		P_c[ind] = P_n[ind];

		dV_s_c[ind] = dV_s_n[ind];
		dV_f_c[ind] = dV_f_n[ind];
		dH_c[ind] = dH_n[ind];
		dP_c[ind] = dP_n[ind];
	}
}

#define lind (ind ==  0 ? (N - 1) : (ind - 1))
#define rind (ind ==  N - 1 ? 0 : (ind + 1))
void toInvariants(int ind){
	W1[0] = a * P_c[lind] - b * H_c[lind] + c * V_s_c[lind] - d * V_f_c[lind];
	W3[0] = am * P_c[lind] + b * H_c[lind] - cm * V_s_c[lind] + dm * V_f_c[lind];

	W1[1] = a * P_c[ind] - b * H_c[ind] + c * V_s_c[ind] - d * V_f_c[ind];
	W3[1] = am * P_c[ind] + b * H_c[ind] - cm * V_s_c[ind] + dm * V_f_c[ind];
	W2[0] = a * P_c[ind] - b * H_c[ind] - c * V_s_c[ind] + d * V_f_c[ind];
	W4[0] = am * P_c[ind] + b * H_c[ind] + cm * V_s_c[ind] - dm * V_f_c[ind];
		
	W2[1] = a * P_c[rind] - b * H_c[rind] - c * V_s_c[rind] + d * V_f_c[rind];
	W4[1] = am * P_c[rind] + b * H_c[rind] + cm * V_s_c[rind] - dm * V_f_c[rind];
}


void toDInvariants(int ind){
	dW1[0] = a * dP_c[lind] - b * dH_c[lind] + c * dV_s_c[lind] - d * dV_f_c[lind];
	dW3[0] = am * dP_c[lind] + b * dH_c[lind] - cm * dV_s_c[lind] + dm * dV_f_c[lind];

	dW1[1] = a * dP_c[ind] - b * dH_c[ind] + c * dV_s_c[ind] - d * dV_f_c[ind];
	dW3[1] = am * dP_c[ind] + b * dH_c[ind] - cm * dV_s_c[ind] + dm * dV_f_c[ind];
	dW2[0] = a * dP_c[ind] - b * dH_c[ind] - c * dV_s_c[ind] + d * dV_f_c[ind];
	dW4[0] = am * dP_c[ind] + b * dH_c[ind] + cm * dV_s_c[ind] - dm * dV_f_c[ind];

	dW2[1] = a * dP_c[rind] - b * dH_c[rind] - c * dV_s_c[rind] + d * dV_f_c[rind];
	dW4[1] = am * dP_c[rind] + b * dH_c[rind] + cm * dV_s_c[rind] - dm * dV_f_c[rind];
}


void fillStencilValues(int ind) {
	unsigned int j;
	for (j = 0; j < S + 1; j++) {
		stencil_values_w1[j] = W1[j];
		stencil_Dvalues_w1[j] = dW1[j];
		stencil_values_w2[j] = W2[j];
		stencil_Dvalues_w2[j] = dW2[j];
		stencil_values_w3[j] = W3[j];
		stencil_Dvalues_w3[j] = dW3[j];
		stencil_values_w4[j] = W4[j];
		stencil_Dvalues_w4[j] = dW4[j];
	}
}

void interpolatedCoefs(double *u, double *du, double c, int dir, int ind, double *coefs) {

	if (dir == LEFT) {
		coefs[1] = (du[1] + du[0]) / h / h - 2.0 * (u[1] - u[0]) / h / h / h;
		coefs[2] = (2.0 * du[1] + du[0]) / h - 3.0 * (u[1] - u[0]) / h / h;
		coefs[3] = du[1];
		coefs[4] = u[1];
		coefs[0] = c * tauT;
	}
	else {
		coefs[1] = (du[1] + du[0]) / h / h - 2.0 * (u[1] - u[0]) / h / h / h;
		coefs[2] = -(2.0 * du[0] + du[1]) / h + 3.0 * (u[1] - u[0]) / h / h;
		coefs[3] = du[0];
		coefs[4] = u[0];
		coefs[0] = c * tauT;
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

// CIR
void singleStepCIR(void) {
	int ind = 0;

	for (ind = 0; ind < N; ind++) {
		double W1n, W2n, W3n, W4n;

		toInvariants(ind);
		//printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", W1[0], W1[1], W2[0], W2[1], W3[0], W3[1]);
		fillStencilValues(ind);

		W1n = interpolatedValueCIR(stencil_values_w1, lambda1, LEFT, ind);
		W2n = interpolatedValueCIR(stencil_values_w2, lambda2, RIGHT, ind);
		W3n = interpolatedValueCIR(stencil_values_w3, lambda3, LEFT, ind);
		W4n = interpolatedValueCIR(stencil_values_w4, lambda4, RIGHT, ind);

		V_s_n[ind] = a1 * W1n - a1 * W2n + a2 * W3n - a2 * W4n;
		V_f_n[ind] = -b1 * W1n + b1 * W2n - b2 * W3n + b2 * W4n;
		H_n[ind] = c1 * W1n + c1 * W2n + c2 * W3n + c2 * W4n;
		P_n[ind] = W1n + W2n + W3n + W4n;
	}

	for (ind = 0; ind < N; ind++){
		V_s_c[ind] = V_s_n[ind];
		V_f_c[ind] = V_f_n[ind];
		H_c[ind] = H_n[ind];
		P_c[ind] = P_n[ind];
	}
}

double interpolatedValueCIR(double *u, double c, int dir, int ind) {
	double a, b, x;
	if (dir == LEFT){
		a = (u[1] - u[0]) / h;
		b = u[1];
		x = c * tauT;
	}
	else {
		a = (u[1] - u[0]) / h;
		b = u[0];
		x = c * tauT;
	}
	double val = a * x + b;
	return val;
}

void debugGnuplotFileName(char *filename1){
	FILE *pFileP;
	pFileP = fopen(filename1, "w");
	int ind;
	for (ind = 0; ind < N; ind++){
		fprintf(pFileP, "%.12lf\t%.12lf\n", -(LENGTH / 2.0) + h * ind, V_s_c[ind]);
	}
	fclose(pFileP);
}

void createArrays(){
	V_s_c = (double *)malloc(N * sizeof(double));
	V_f_c = (double *)malloc(N * sizeof(double));
	H_c = (double *)malloc(N * sizeof(double));
	P_c = (double *)malloc(N * sizeof(double));
	dV_s_c = (double *)malloc(N * sizeof(double));
	dV_f_c = (double *)malloc(N * sizeof(double));
	dH_c = (double *)malloc(N * sizeof(double));
	dP_c = (double *)malloc(N * sizeof(double));

	V_s_n = (double *)malloc(N * sizeof(double));
	V_f_n = (double *)malloc(N * sizeof(double));
	H_n = (double *)malloc(N * sizeof(double));
	P_n = (double *)malloc(N * sizeof(double));
	dV_s_n = (double *)malloc(N * sizeof(double));
	dV_f_n = (double *)malloc(N * sizeof(double));
	dH_n = (double *)malloc(N * sizeof(double));
	dP_n = (double *)malloc(N * sizeof(double));
}

void freeArrays(){
	free(V_s_c);
	free(V_f_c);
	free(H_c);
	free(P_c);
	free(V_s_n);
	free(V_f_n);
	free(H_n);
	free(P_n);
	free(dV_s_c);
	free(dV_f_c);
	free(dH_c);
	free(dP_c);
	free(dV_s_n);
	free(dV_f_n);
	free(dH_n);
	free(dP_n);
}