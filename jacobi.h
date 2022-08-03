#ifndef JACOBI_H
#define JACOBI_H

#include <iostream>
#include <cmath>
#include <vector>
#include <omp.h>
#include "ThreeDMG.h"
extern long long total;

std::vector<double> Jacobi(std::vector<double> B, std::vector<double> PHI,
	std::vector<double> A,std::vector<int> IA, 
	std::vector<int> JA, int ntot, int sweeps, int sum) {
	double D = A[IA[sum]], val;	
	int k,j;	
	std::vector<double> midPHI(ntot);
	for (int iters=0; iters<sweeps; iters++) {	
		midPHI = PHI;
	#pragma omp parallel 
		{	
			#pragma omp for private(k,val,j)
			for (j = 0; j < ntot; j++) {
				val = B[j];
				for (k = IA[sum+j]; k < IA[sum+j+1]; k++) {
					val -= midPHI[JA[k]] * A[k]; 
				}
				PHI[j] += val / D;
			}
		}
		
	}
	
	return PHI;
}

#endif