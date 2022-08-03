#include <iostream>
#include <cmath>
#include <vector>
#include <omp.h>

namespace ThreeD {
	std::vector<double> restriction(std::vector<double> R, int n);
	std::vector<double> interpolation(std::vector<double> R2, int n2);
}
std::vector<double> ThreeD::restriction(std::vector<double> R, int n) {
	int n2 = n / 2;
	std::vector<double> R2(n2*n2*n2,0.0);
	#pragma omp parallel for
	for (int i = 0; i < n2; i++) {
		for (int j = 0; j < n2; j++) {
			for (int k = 0; k < n2; k++) {
				// Middle layer
				R2[i*n2*n2+j*n2+k] += R[n * n * (2 * i + 1) + n * (2 * j + 1) + (2 * k + 1)] * 8;
				R2[i*n2*n2+j*n2+k] += R[n * n * (2 * i + 1) + n * (2 * j    ) + (2 * k + 1)] * 4;
				R2[i*n2*n2+j*n2+k] += R[n * n * (2 * i + 1) + n * (2 * j + 2) + (2 * k + 1)] * 4;
				R2[i*n2*n2+j*n2+k] += R[n * n * (2 * i + 1) + n * (2 * j + 1) + (2 * k    )] * 4;
				R2[i*n2*n2+j*n2+k] += R[n * n * (2 * i + 1) + n * (2 * j + 1) + (2 * k + 2)] * 4;
				R2[i*n2*n2+j*n2+k] += R[n * n * (2 * i + 1) + n * (2 * j    ) + (2 * k    )] * 2;
				R2[i*n2*n2+j*n2+k] += R[n * n * (2 * i + 1) + n * (2 * j + 2) + (2 * k    )] * 2;
				R2[i*n2*n2+j*n2+k] += R[n * n * (2 * i + 1) + n * (2 * j    ) + (2 * k + 2)] * 2;
				R2[i*n2*n2+j*n2+k] += R[n * n * (2 * i + 1) + n * (2 * j + 2) + (2 * k + 2)] * 2;
				// Bottom layer
				R2[i*n2*n2+j*n2+k] += R[n * n * (2 * i    ) + n * (2 * j + 1) + (2 * k + 1)] * 4;
				R2[i*n2*n2+j*n2+k] += R[n * n * (2 * i    ) + n * (2 * j    ) + (2 * k + 1)] * 2;
				R2[i*n2*n2+j*n2+k] += R[n * n * (2 * i    ) + n * (2 * j + 2) + (2 * k + 1)] * 2;
				R2[i*n2*n2+j*n2+k] += R[n * n * (2 * i    ) + n * (2 * j + 1) + (2 * k    )] * 2;
				R2[i*n2*n2+j*n2+k] += R[n * n * (2 * i    ) + n * (2 * j + 1) + (2 * k + 2)] * 2;
				R2[i*n2*n2+j*n2+k] += R[n * n * (2 * i    ) + n * (2 * j    ) + (2 * k    )] * 1;
				R2[i*n2*n2+j*n2+k] += R[n * n * (2 * i    ) + n * (2 * j + 2) + (2 * k    )] * 1;
				R2[i*n2*n2+j*n2+k] += R[n * n * (2 * i    ) + n * (2 * j    ) + (2 * k + 2)] * 1;
				R2[i*n2*n2+j*n2+k] += R[n * n * (2 * i    ) + n * (2 * j + 2) + (2 * k + 2)] * 1;
				// Top layer 
				R2[i*n2*n2+j*n2+k] += R[n * n * (2 * i + 2) + n * (2 * j + 1) + (2 * k + 1)] * 4;
				R2[i*n2*n2+j*n2+k] += R[n * n * (2 * i + 2) + n * (2 * j    ) + (2 * k + 1)] * 2;
				R2[i*n2*n2+j*n2+k] += R[n * n * (2 * i + 2) + n * (2 * j + 2) + (2 * k + 1)] * 2;
				R2[i*n2*n2+j*n2+k] += R[n * n * (2 * i + 2) + n * (2 * j + 1) + (2 * k    )] * 2;
				R2[i*n2*n2+j*n2+k] += R[n * n * (2 * i + 2) + n * (2 * j + 1) + (2 * k + 2)] * 2;
				R2[i*n2*n2+j*n2+k] += R[n * n * (2 * i + 2) + n * (2 * j    ) + (2 * k    )] * 1;
				R2[i*n2*n2+j*n2+k] += R[n * n * (2 * i + 2) + n * (2 * j + 2) + (2 * k    )] * 1;
				R2[i*n2*n2+j*n2+k] += R[n * n * (2 * i + 2) + n * (2 * j    ) + (2 * k + 2)] * 1;
				R2[i*n2*n2+j*n2+k] += R[n * n * (2 * i + 2) + n * (2 * j + 2) + (2 * k + 2)] * 1;
				R2[i*n2*n2+j*n2+k] /= 64;
			}
		}
	}

	return R2;

}

std::vector<double> ThreeD::interpolation(std::vector<double> R2, int n2) {
	int n = n2 * 2 + 1;
	std::vector<double> R(n*n*n,0.0);
	for (int i = 0; i < n2; i++) {
		for (int j = 0; j < n2; j++) {
			for (int k = 0; k < n2; k++) {
				// Middle layer
				R[n * n * (2 * i + 1) + n * (2 * j + 1) + (2 * k + 1)] += R2[i*n2*n2+j*n2+k]    ;
				R[n * n * (2 * i + 1) + n * (2 * j    ) + (2 * k + 1)] += R2[i*n2*n2+j*n2+k] / 2;
				R[n * n * (2 * i + 1) + n * (2 * j + 2) + (2 * k + 1)] += R2[i*n2*n2+j*n2+k] / 2;
				R[n * n * (2 * i + 1) + n * (2 * j + 1) + (2 * k    )] += R2[i*n2*n2+j*n2+k] / 2;
				R[n * n * (2 * i + 1) + n * (2 * j + 1) + (2 * k + 2)] += R2[i*n2*n2+j*n2+k] / 2;
				R[n * n * (2 * i + 1) + n * (2 * j    ) + (2 * k    )] += R2[i*n2*n2+j*n2+k] / 4;
				R[n * n * (2 * i + 1) + n * (2 * j + 2) + (2 * k    )] += R2[i*n2*n2+j*n2+k] / 4;
				R[n * n * (2 * i + 1) + n * (2 * j    ) + (2 * k + 2)] += R2[i*n2*n2+j*n2+k] / 4;
				R[n * n * (2 * i + 1) + n * (2 * j + 2) + (2 * k + 2)] += R2[i*n2*n2+j*n2+k] / 4;
				// Bot layer
				R[n * n * (2 * i    ) + n * (2 * j + 1) + (2 * k + 1)] += R2[i*n2*n2+j*n2+k] / 2;
				R[n * n * (2 * i    ) + n * (2 * j    ) + (2 * k + 1)] += R2[i*n2*n2+j*n2+k] / 4;
				R[n * n * (2 * i    ) + n * (2 * j + 2) + (2 * k + 1)] += R2[i*n2*n2+j*n2+k] / 4;
				R[n * n * (2 * i    ) + n * (2 * j + 1) + (2 * k    )] += R2[i*n2*n2+j*n2+k] / 4;
				R[n * n * (2 * i    ) + n * (2 * j + 1) + (2 * k + 2)] += R2[i*n2*n2+j*n2+k] / 4;
				R[n * n * (2 * i    ) + n * (2 * j    ) + (2 * k    )] += R2[i*n2*n2+j*n2+k] / 8;
				R[n * n * (2 * i    ) + n * (2 * j + 2) + (2 * k    )] += R2[i*n2*n2+j*n2+k] / 8;
				R[n * n * (2 * i    ) + n * (2 * j    ) + (2 * k + 2)] += R2[i*n2*n2+j*n2+k] / 8;
				R[n * n * (2 * i    ) + n * (2 * j + 2) + (2 * k + 2)] += R2[i*n2*n2+j*n2+k] / 8;
				// Top layer 
				R[n * n * (2 * i + 2) + n * (2 * j + 1) + (2 * k + 1)] += R2[i*n2*n2+j*n2+k] / 2;
				R[n * n * (2 * i + 2) + n * (2 * j    ) + (2 * k + 1)] += R2[i*n2*n2+j*n2+k] / 4;
				R[n * n * (2 * i + 2) + n * (2 * j + 2) + (2 * k + 1)] += R2[i*n2*n2+j*n2+k] / 4;
				R[n * n * (2 * i + 2) + n * (2 * j + 1) + (2 * k    )] += R2[i*n2*n2+j*n2+k] / 4;
				R[n * n * (2 * i + 2) + n * (2 * j + 1) + (2 * k + 2)] += R2[i*n2*n2+j*n2+k] / 4;
				R[n * n * (2 * i + 2) + n * (2 * j    ) + (2 * k    )] += R2[i*n2*n2+j*n2+k] / 8;
				R[n * n * (2 * i + 2) + n * (2 * j + 2) + (2 * k    )] += R2[i*n2*n2+j*n2+k] / 8;
				R[n * n * (2 * i + 2) + n * (2 * j    ) + (2 * k + 2)] += R2[i*n2*n2+j*n2+k] / 8;
				R[n * n * (2 * i + 2) + n * (2 * j + 2) + (2 * k + 2)] += R2[i*n2*n2+j*n2+k] / 8;
			}
		}
	}
	return R;
}