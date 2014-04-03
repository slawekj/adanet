/*
 *   This is an implementation of ADANET algorithm for Gene Regulatory Network
 *   inference from mRNA expression data, in form of a Matlab package.
 *   Copyright (C) 2014  Janusz Slawek
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program, see LICENSE.
 */

#include "mex.h"
#include "decision_stump.h"

/*
 * Input arguments
 */

#define X_IN prhs[0]
#define Y_IN prhs[1]
#define M_IN prhs[2]

/*
 * Output arguments
 */

#define W_OUT plhs[0]

/*
 *
 * prhs An array of right-hand input arguments.
 * plhs An array of left-hand output arguments.
 * nrhs The number of right-hand arguments, or the size of the prhs array.
 * nlhs The number of left-hand arguments, or the size of the plhs array.
 *
 */

void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
	mwSize n,p,m;
	double *x, *y, *w;

	/* Check for proper number of arguments */
	
	if (nrhs != 3) {
		mexErrMsgTxt("Three arguments required.");
	} else if (nlhs > 1) {
		mexErrMsgTxt("Too many output arguments.");
	}

	/* Check for dimensions of W */
	n = mxGetM(X_IN);
	p = mxGetN(X_IN);
	x = mxGetPr(X_IN);
	y = mxGetPr(Y_IN);
	m = mxGetN(M_IN);

	W_OUT = mxCreateDoubleMatrix(1,p, mxREAL);
	w = mxGetPr(W_OUT);

	/* Calculate W */
	train_decision_stump(n,p,x,y,m,w);

	/* return */
	return;
}

