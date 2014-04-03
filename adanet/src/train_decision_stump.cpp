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

#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cfloat>
#include "decision_stump.h"

inline int compare(const void *a, const void *b) {
	double result = **(double **) a - **(double **) b;
	if (result > 0)
		return 1;
	if (result < 0)
		return -1;
	return 0;
}

void train_decision_stump(int N, int P, double const *x,
		double const *y_zero_one, int M, double *w) {
	/*
	 * reset random seed
	 */
	srand(time(NULL));

	/*
	 * allocate memory of arrays ...
	 */

	const double **x_to_sort = (const double **) calloc(N,
			sizeof(const double *));
	int *x_sorted_index = (int *) calloc(N * P, sizeof(int));
	double *D = (double *) calloc(N, sizeof(double));
	double *y_neg_pos = (double *) calloc(N, sizeof(double));
	double *y_predicted = (double *) calloc(N, sizeof(double));

	/*
	 * ... and scalars ...
	 */
	int col;
	int row;
	int next_row;
	int sumN;
	int sumP;
	double total_w;
	int iteration;
	double best_error;
	double local_error_01;
	double local_error_10;
	const double *origin;
	double error_01_l;
	double error_01_r;
	int s;
	double local_threshold;
	int best_column;
	double best_sgn;
	double best_threshold;
	double alpha;

	/*
	 * initialize y_neg_pos
	 */
	for (row = 0; row < N; row++) {
		y_neg_pos[row] = 2 * y_zero_one[row] - 1.0;
	}

	/*
	 * get indices of sorted columns of x
	 *
	 */

	for (col = 0; col < P; col++) {
		for (row = 0; row < N; row++) {
			x_to_sort[row] = x + col * N + row;
		}
		origin = x_to_sort[0];
		qsort(x_to_sort, N, sizeof(double *), compare);
		for (row = 0; row < N; row++) {
			x_sorted_index[col * N + row] = x_to_sort[row] - origin;
		}
	}

	/*
	 * initialize D
	 */
	sumN = 0;
	sumP = 0;
	for (row = 0; row < N; row++) {
		if (y_zero_one[row] > 0.0) {
			sumP++;
		} else {
			sumN++;
		}
	}

	/* check if we have both classes */
	if (!sumN || !sumP) {
		return;
	}

	for (row = 0; row < N; row++) {
		if (y_zero_one[row] > 0.0) {
			D[row] = 1.0 * sumN;
		} else {
			D[row] = 1.0 * sumP;
		}
	}
	total_w = 0.0;
	for (row = 0; row < N; row++) {
		total_w += D[row];
	}

	/* normalize D */
	for (row = 0; row < N; row++) {
		D[row] /= total_w;
	}

	/* initialize weights of features to 0 */
	for (col = 0; col < P; col++) {
		w[col] = 0.0;
	}

	for (iteration = 0; iteration < M; iteration++) {
		/* total_w is always 1.0 because it gets normalized */
		total_w = 1.0;

		best_error = 1.0;
		for (col = 0; col < P; col++) {
			/*
			 * initialize variables
			 */
			error_01_l = 0.0;
			error_01_r = 0.0;
			for (row = 0; row < N; row++) {
				if (y_zero_one[row] == 0.0) {
					error_01_r += D[row];
				}
			}
			/*
			 * go through all the thresholds
			 */
			for (s = 0; s < N - 1; s++) {
				row = x_sorted_index[col * N + s];
				next_row = x_sorted_index[col * N + s + 1];

				error_01_l += y_zero_one[row] * D[row];
				error_01_r -= (1.0 - y_zero_one[row]) * D[row];
				local_error_01 = error_01_l + error_01_r;
				local_error_10 = 1.0 - local_error_01;

				/*
				 * test if it's a valid split point
				 * OK: 1 2 | 3 4
				 * WRONG: 1 1 | 1 3
				 *
				 */
				if (x[col * N + row] < x[col * N + next_row]) {
					/*
					 * it's a valid split
					 */
					local_threshold = (x[col * N + row] + x[col * N + next_row])
							/ 2.0;

					if (local_error_01 < best_error) {
						best_error = local_error_01;
						best_sgn = 1.0;
						best_column = col;
						best_threshold = local_threshold;
					}
					if (local_error_10 < best_error) {
						best_error = local_error_10;
						best_sgn = -1.0;
						best_column = col;
						best_threshold = local_threshold;
					}
				}
			}
		}
		if (best_error == 0.0) {
			/*
			 * Special case, error equal to zero
			 * no need to look any further
			 * break the loop, iterations done
			 *
			 */
			w[best_column] += 3.0;
			break;
		} else if (best_error < 0.5) {
			/*
			 * Find class prediction
			 */
			for (row = 0; row < N; row++) {
				if (x[best_column * N + row] < best_threshold) {
					y_predicted[row] = best_sgn * -1.0;
				} else {
					y_predicted[row] = best_sgn * 1.0;
				}
			}

			/*
			 * Find alpha
			 */
			alpha = 0.5 * log((1 - best_error) / best_error);

			/*
			 * Update D
			 */
			for (row = 0; row < N; row++) {
				D[row] *= exp(-alpha * y_neg_pos[row] * y_predicted[row]);
			}
			total_w = 0.0;
			for (row = 0; row < N; row++) {
				total_w += D[row];
			}

			/* normalize D */
			for (row = 0; row < N; row++) {
				D[row] /= total_w;
			}

			/*
			 * Update weight
			 */
			w[best_column] += alpha;
		}
	}

	/*
	 * Scale w
	 */
	total_w = 0.0;
	for (col = 0; col < P; col++) {
		total_w += w[col];
	}
	if (total_w > 0.0) {
		for (col = 0; col < P; col++) {
			w[col] /= total_w;
		}
	}

	/* deallocate all variables allocated dynamically */
	free(y_predicted);
	free(y_neg_pos);
	free(D);
	free(x_to_sort);
	free(x_sorted_index);

	return;
}
