#include "makespl.h"
#include "piv_ge_solver.h"

#include <stdio.h>
#include <stdlib.h>
#include <float.h>

/* UWAGA: liczbę używanych f. bazowych można ustawić przez wartość
          zmiennej środowiskowej APPROX_BASE_SIZE
*/

/*
 * Funkcje bazowe: n - liczba funkcji a,b - granice przedzialu aproksymacji i
 * - numer funkcji x - wspolrzedna dla ktorej obliczana jest wartosc funkcji
 */
double
fi( int i, double x )
{
	int j;
	double f, pf, tmp;
	pf = 1;
	f = 2 * x;

	if( i == 0 )
		return pf;

	if( i == 1 )
		return f;

	for( j = 1; j < i; j++ )
	{
		tmp = f;
		f = f * 2 * x;
		f -= 2 * j * pf;
		pf = tmp;
	}

	return f;
}

/* Pierwsza pochodna fi */
double
dfi( int i, double x )
{
	double df;

	if( i == 0 )
		return 0;

	df = fi( i-1, x );
	df = df * 2 * i;
	return df;
}

/* Druga pochodna fi */
double
d2fi( int i, double x )
{
	double d2f;

	if( i == 0 )
		return 0;

	d2f = dfi( i-1, x );
	d2f = d2f * 2 * i;
	return d2f;
}

/* Trzecia pochodna fi */
double
d3fi( int i, double x )
{
	double d3f;

	if( i == 0 )
		return 0;

	d3f = d2fi( i-1, x );
	d3f = d3f * 2 * i;
	return d3f;
}


void
make_spl(points_t * pts, spline_t * spl)
{

	matrix_t       *eqs= NULL;
	double         *x = pts->x;
	double         *y = pts->y;
	double		a = x[0];
	double		b = x[pts->n - 1];
	int		i, j, k;
	int		nb = pts->n - 3 > 10 ? 10 : pts->n - 3;
	char *nbEnv= getenv( "APPROX_BASE_SIZE" );

	if( nbEnv != NULL && atoi( nbEnv ) > 0 )
		nb = atoi( nbEnv );

	eqs = make_matrix(nb, nb + 1);

	for (j = 0; j < nb; j++)
	{
		for (i = 0; i < nb; i++)
			for (k = 0; k < pts->n; k++)
				add_to_entry_matrix(eqs, j, i, fi(i, x[k]) * fi(j, x[k]));

		for (k = 0; k < pts->n; k++)
			add_to_entry_matrix(eqs, j, nb, y[k] * fi(j, x[k]));
	}


	if (piv_ge_solver(eqs))
	{
		spl->n = 0;
		return;
	}

	if (alloc_spl(spl, nb) == 0)
	{
		for (i = 0; i < spl->n; i++)
		{
			double xx = spl->x[i] = a + i*(b-a)/(spl->n-1);
			xx+= 10.0*DBL_EPSILON;  // zabezpieczenie przed ulokowaniem punktu w poprzednim przedziale
			spl->f[i] = 0;
			spl->f1[i] = 0;
			spl->f2[i] = 0;
			spl->f3[i] = 0;
			for (k = 0; k < nb; k++) {
				double		ck = get_entry_matrix(eqs, k, nb);
				spl->f[i]  += ck * fi  (k, xx);
				spl->f1[i] += ck * dfi (k, xx);
				spl->f2[i] += ck * d2fi(k, xx);
				spl->f3[i] += ck * d3fi(k, xx);
			}
		}
	}
}
