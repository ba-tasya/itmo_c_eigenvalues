#include "return_codes.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double EPS = 1e-9;

int is_zero(double x)
{
	return (fabs(x) < EPS);
}

void rotate_right(double *matrix, double c, double s, size_t x, size_t y, size_t n)
{
	double tempx = 0;
	double tempy = 0;
	for (size_t i = 0; i < n; ++i)
	{
		tempx = matrix[i * n + x];
		tempy = matrix[i * n + y];
		matrix[i * n + x] = c * tempx + s * tempy;
		matrix[i * n + y] = -s * tempx + c * tempy;
	}
}

void rotate_left(double *matrix, double c, double s, size_t x, size_t y, size_t n)
{
	double tempx = 0;
	double tempy = 0;
	for (size_t j = 0; j < n; ++j)
	{
		tempx = matrix[x * n + j];
		tempy = matrix[y * n + j];
		matrix[x * n + j] = c * tempx + s * tempy;
		matrix[y * n + j] = -s * tempx + c * tempy;
	}
}

void calc_cs(double a, double b, double *c, double *s)
{
	if (a == 0 && b == 0)
	{
		*c = 1;
		*s = 0;
		return;
	}
	double tmp = 1.0f / (sqrt(a * a + b * b));
	*c = tmp * a;
	*s = tmp * b;
}

void hessen(double *matrix, size_t n)
{
	if (n == 1)
	{
		return;
	}
	double c, s;
	for (size_t column = 0; column < n - 2; ++column)
	{
		for (size_t row = column + 2; row < n; ++row)
		{
			calc_cs(matrix[(column + 1) * n + column], matrix[row * n + column], &c, &s);
			rotate_left(matrix, c, s, column + 1, row, n);
			rotate_right(matrix, c, s, column + 1, row, n);
		}
	}
}

void givesn(double *matrixA, double *rotates, size_t n)
{
	double c, s;
	for (size_t column = 0; column < n - 1; ++column)
	{
		calc_cs(matrixA[column * n + column], matrixA[(column + 1) * n + column], &c, &s);
		rotate_left(matrixA, c, s, column, column + 1, n);
		rotates[column * 2] = c;
		rotates[column * 2 + 1] = s;
	}
	for (size_t column = 0; column < n - 1; ++column)
	{
		rotate_right(matrixA, rotates[column * 2], rotates[column * 2 + 1], column, column + 1, n);
	}
}

int check_matrix(double *matrix, size_t n)
{
	for (size_t i = 1; i < n; ++i)
	{
		if (!is_zero(matrix[i * n + (i - 1)]) && i != n - 1 && !is_zero(matrix[(i + 1) * n + i]))
		{
			return 0;
		}
	}
	return 1;
}

void eigvals(double *matrixA, double *rotates, size_t n)
{
	hessen(matrixA, n);
	while (!check_matrix(matrixA, n))
	{
		givesn(matrixA, rotates, n);
	}
}

int main(int argc, char **argv)
{
	int result = SUCCESS;
	FILE *in = NULL, *out = NULL;
	double *matrix = NULL, *rotates = NULL;
	if (argc != 3)
	{
		fprintf(stderr, "Expected two arguments: input file name and output file name.");
		result = ERROR_PARAMETER_INVALID;
		goto CleanUp;
	}
	in = fopen(argv[1], "r");
	if (!in)
	{
		fprintf(stderr, "Failed to open input file");
		result = ERROR_CANNOT_OPEN_FILE;
		goto CleanUp;
	}
	size_t n;
	int ret = fscanf(in, "%zu", &n);
	if (ret == EOF || ret == 0)
	{
		fprintf(stderr, "Invalid data");
		result = ERROR_DATA_INVALID;
		goto CleanUp;
	}
	matrix = malloc(sizeof(double) * n * n);
	if (!matrix)
	{
		fprintf(stderr, "Failed to allocate memory");
		result = ERROR_OUT_OF_MEMORY;
		goto CleanUp;
	}
	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < n; j++)
		{
			int ret = fscanf(in, "%lf", &matrix[i * n + j]);
			if (ret == EOF || ret == 0)
			{
				fprintf(stderr, "Invalid data");
				result = ERROR_DATA_INVALID;
				goto CleanUp;
			}
			if (matrix[i * n + j] != 0 && EPS > fabs(matrix[i * n + j]) / 1e9)
			{
				EPS = fabs(matrix[i * n + j]) / 1e9;
			}
		}
	}
	out = fopen(argv[2], "w");
	if (!out)
	{
		fprintf(stderr, "Failed to open output file");
		result = ERROR_CANNOT_OPEN_FILE;
		goto CleanUp;
	}
	rotates = malloc(sizeof(double) * (n - 1) * 2);
	if (!rotates)
	{
		fprintf(stderr, "Failed to allocate memory");
		result = ERROR_OUT_OF_MEMORY;
		goto CleanUp;
	}
	eigvals(matrix, rotates, n);
	for (size_t i = 0; i < n; ++i)
	{
		if (i < n - 1 && !is_zero(matrix[(i + 1) * n + i]))
		{
			double a = matrix[i * n + i], b = matrix[i * n + (i + 1)], c = matrix[(i + 1) * n + i],
				   d = matrix[(i + 1) * n + (i + 1)];
			double D = (a * a - 2 * a * d + d * d + 4 * b * c);

			if (D >= 0)
			{
				fprintf(out, "%g\n", (a + d + sqrt(D)) / 2);
				fprintf(out, "%g\n", (a + d - sqrt(D)) / 2);
			}
			else
			{
				fprintf(out, "%g +%gi\n", (a + d) / 2, sqrt(-D) / 2);
				fprintf(out, "%g %gi\n", (a + d) / 2, -sqrt(-D) / 2);
			}
			++i;
		}
		else
		{
			fprintf(out, "%g\n", matrix[i * n + i]);
		}
	}
CleanUp:
	free(matrix);
	free(rotates);
	if (in)
		fclose(in);
	if (out)
		fclose(out);
	return result;
}