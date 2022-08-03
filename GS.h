#ifndef GS_H
#define GS_H

#include <iostream>
#include <cmath>
#include <vector>
#include <omp.h>

std::vector<double> Gauss_Seidel(std::vector<double> B, std::vector<double> PHI,
	std::vector<double> A,std::vector<int> IA, 
	std::vector<int> JA, int ntot, int sweeps, int sum) {
	
	double D = A[IA[sum]], val;
	int j,k;
	for (int iters=0; iters<sweeps; iters++) {
		#pragma omp parallel 
		{
		#pragma omp for private(k,val) 
			for (j = 0; j < ntot; j++) {
				val = B[j];
				for (k = IA[sum+j]; k < IA[sum+j+1]; k++) {
					val -= PHI[JA[k]] * A[k]; 
				}
				PHI[j] += val / D;
			}
			

		}
	}

	return PHI;

}

#endif