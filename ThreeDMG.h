#ifndef THREEDMG_H
#define THREEDMG_H
#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>
#include "GS.h"
#include "ThreeDMapping.h"
#include "jacobi.h"
extern long long total;


std::vector<double> ThreeDMG(int level, int maxlevel, std::vector<double>B,std::vector<double> PHI,
	std::vector<double> A, std::vector<int> IA, std::vector<int> JA, int n, int sweeps, int sum) {
	int i, n3 = n*n*n;
	if (level == maxlevel) {
		PHI = Gauss_Seidel(B,PHI,A,IA,JA,n3,sweeps*2,sum);
		return PHI;
	}
	else {
		// Pre Smoothing
		PHI = Gauss_Seidel(B,PHI,A,IA,JA,n3,sweeps,sum);
		std::vector<double> R(n3);

		// Calculating Residual
		double val;
		#pragma omp parallel for private(val)
		for (i=0; i < n3; i++) {
			val = B[i];
			for (int j=IA[i+sum]; j<IA[i+sum+1]; j++) {
				val -= PHI[JA[j]] * A[j];
			}
			R[i] = val;
		}

		// Restricting Residual
		std::vector<double> B2;
		B2 = ThreeD::restriction(R,n);
		
		//Solving on Coarse Grid
		int nnew = n / 2;
		std::vector<double> PHI2(nnew*nnew*nnew, 0.0);
		PHI2 = ThreeDMG(level+1, maxlevel, B2, PHI2, A, IA, JA, nnew, sweeps, sum+n3);

		// Interpolating
		std::vector<double> PHI_PROL;
		PHI_PROL = ThreeD::interpolation(PHI2, nnew);
		// Adding to solution
		#pragma omp parallel for
		for (i=0; i<n3; i++) {
			PHI[i] += PHI_PROL[i];
		}

		// Post Smoothing
		PHI = Gauss_Seidel(B,PHI,A,IA,JA,n3,sweeps,sum);
		return PHI;
	}
}

#endif