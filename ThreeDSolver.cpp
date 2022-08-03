#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <random>
#include <algorithm>
#include <omp.h>
#include "ThreeDMG.h"

int main(int argc, char **argv) {
	int n = 31;
	int n2 = n * n;
	int n3 = n2 * n;
	int NNZ = 0;
	int maxlayers = 0;
	int maxiters = 1000;
	int sweeps = 2;
	int cycles = 2;
	double EPS = 10e-8;
	double dx = 1.0;
	std::vector<double> A;
	std::vector<int> JA;
	std::vector<int> IA{0};
	std::vector<double> R(n3, 0.0);
	int current_n = n;
	for (int layer = 0; layer <= maxlayers; layer++) {
		for (int i = 0; i < (int)pow(current_n, 3); i++) {
			if (i >= (int)pow(current_n, 2)) {
				A.push_back(-1 / dx / dx);
				JA.push_back(i - pow(current_n, 2));
				NNZ++;
			}			
			if (i % (int)pow(current_n, 2) >= current_n) {
				A.push_back(-1 / dx / dx);
				JA.push_back(i - current_n);
				NNZ++;
			}
			if (i % current_n != 0) {
				A.push_back(-1 / dx / dx);
				JA.push_back(i - 1);
				NNZ++;
			}
			A.push_back(6 / dx / dx);
			JA.push_back(i);
			NNZ++;
			if (i % current_n != (current_n - 1)) {
				A.push_back(-1 / dx / dx);
				JA.push_back(i + 1);
				NNZ++;
			}
			if (i % (int)pow(current_n, 2) < ((int)pow(current_n, 2) - current_n)) {
				A.push_back(-1 / dx / dx);
				JA.push_back(i + current_n);
				NNZ++;
			}
			if (i < (int)pow(current_n, 3) - (int)pow(current_n, 2)) {
				A.push_back(-1 / dx / dx);
				JA.push_back(i + pow(current_n, 2));
				NNZ++;
			}
			IA.push_back(NNZ);
		}
		
		if (current_n % 2 == 0 || current_n < 2) {
			maxlayers = layer;
			break;
		}
		current_n /= 2;
		dx *= 2;
	}
	std::vector<double> B(n3,1.0);
	/*std::random_device device_random_;
	std::default_random_engine generator_(device_random_());
	std::normal_distribution<> distribution(0.0,.0001);
	for (int i=0; i<n3; i++) {
		B.push_back(distribution(generator_));
	}*/

	std::vector<double> PHI(n3, 0.0);
	for (int i=0; i < 100; i++)	PHI = ThreeDMG(0,maxlayers,B,PHI,A,IA,JA,n,sweeps,0);

	double resid = 0.0; 
	for (int i = 0; i < n3; i++) {
		double val = B[i];
		for (int j=IA[i]; j<IA[i+1]; j++) {
			val -= PHI[JA[j]] * A[j];
		}
		resid += val;
	}

	std::cout << resid / (n3) << std::endl;

	return 0;
}