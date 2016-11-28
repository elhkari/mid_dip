//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define OUTFILE "dip1out.txt"

#define a 10.
//#define b 10.
#define L 1.

#define NUM 301
#define NUM2 601
#define c 5.e-15
#define drho ((a)/(NUM - c))
#define dz ((a)/(NUM - c))
/*
#define lNUM 0

#define S (0.)
#define sdrho ( S/(lNUM - 0.5))
#define sdz ( S/(lNUM - 0.5))

#define drho ((a-S)/(NUM - lNUM))
#define dz ((a-S)/(NUM - lNUM))
*/
double rho[NUM], z[NUM2], psi[NUM][NUM2], vel[NUM][NUM2], V[NUM][NUM2];

double E, num, denum, maxgE, E_kin, E_pot, C;
double gE[NUM][NUM2];

void prepfuncs()
{
	int i, j;
	for (i = 0; i < NUM; i++)
	{
		for (j = 0; j < NUM2; j++)
		{
			V[i][j] = 1./(sqrt(rho[i] * rho[i] + z[j] * z[j]));
			psi[i][j] = exp(-sqrt(rho[i] * rho[i] + z[j] * z[j]));
			vel[i][j] = 0;
		}
	}
}


void calc_E()
{
	double sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9, sum10, sum11, sum12;
	int i, j;

	sum1 = sum2 = sum3 = sum4 = sum5 = sum6 = sum7 = sum8 = sum9 = sum10 = sum11 = sum12 = 0;

	for (i = 0; i < NUM - 1; i++)
	{
		sum9 += psi[i][0] * psi[i][0] * rho[i];
		sum10 += psi[i + 1][0] * psi[i + 1][0] * rho[i + 1];
		
		for (j = 0; j < NUM2 - 1; j++)
		{
			sum1 += (psi[i + 1][j] - psi[i][j])*(psi[i + 1][j] - psi[i][j])*(rho[i + 1] + rho[i]);
			sum2 += (psi[i + 1][j + 1] - psi[i][j + 1])*(psi[i + 1][j + 1] - psi[i][j + 1])*(rho[i + 1] + rho[i]);

			sum3 += (psi[i][j + 1] - psi[i][j])*(psi[i][j + 1] - psi[i][j])*rho[i];
			sum4 += (psi[i + 1][j + 1] - psi[i + 1][j])*(psi[i + 1][j + 1] - psi[i + 1][j])*rho[i + 1];

			sum5 += psi[i][j] * psi[i][j] * rho[i] * V[i][j];
			sum6 += psi[i][j + 1] * psi[i][j + 1] * rho[i] * V[i][j+1];
			sum7 += psi[i + 1][j] * psi[i + 1][j] * rho[i + 1] * V[i+1][j];
			sum8 += psi[i + 1][j + 1] * psi[i + 1][j + 1] * rho[i + 1] * V[i+1][j+1];

			sum11 += (psi[i][j] * psi[i][j] + psi[i][j + 1] * psi[i][j + 1]) * rho[i];
			sum12 += (psi[i + 1][j] * psi[i + 1][j] + psi[i + 1][j + 1] * psi[i + 1][j + 1]) * rho[i + 1];
		}
	}


	sum1 = dz*sum1/(8.*drho);
	sum2 = dz*sum2/(8.*drho);
	
	sum3 = drho*sum3 /(4.* dz);
	sum4 = drho*sum4 /(4.* dz);
	
	sum5 = drho*dz*sum5 / 4.;
	sum6 = drho*dz*sum6 / 4.;
	sum7 = drho*dz*sum7 / 4.;
	sum8 = drho*dz*sum8 / 4.;

	sum9 = drho*L*sum9 / 4.;
	sum10 = drho*L*sum10 / 4.;

	sum11 = drho*dz*sum11 / 4.;
	sum12 = drho*dz*sum12 / 4.;

	//num = sum1 + sum2 + sum3 + sum4 - sum5 - sum6 - sum7 - sum8 + sum9 + sum10;
	E_kin = sum1 + sum2 + sum3 + sum4;
	E_pot = sum5 + sum6 + sum7 + sum8;
	C = sum9 + sum10;
	num = E_kin - E_pot + sum9 + sum10;
	denum = sum11 + sum12;

	E = num / denum;
}

void calc_gE()
{
	double dv, dw;
	int i, j;
////////////////////////////////////////////////////////////////////////////////////////
	i = 0;
//
	j = 0;

	dv = (rho[i + 1] + rho[i])*(psi[i + 1][j] - psi[i][j])*(-1.)*dz / (4.*drho) +
		rho[i] * (psi[i][j + 1] - psi[i][j])*(-1.)*drho / (2.*dz) -
		psi[i][j] * V[i][j] * rho[i] * drho*dz / 2.  +
		L*psi[i][j]*rho[i]*drho / 2.;

	dw = psi[i][j] * rho[i] * drho*dz / 2.;

	gE[i][j] = (dw*num - dv*denum) / (denum*denum);
	//
	j = NUM2 - 1;

	dv = (rho[i + 1] + rho[i])*(psi[i + 1][j] - psi[i][j])*(-1.)*dz / (4.*drho) +
		rho[i] * (psi[i][j] - psi[i][j-1])*drho / (2.*dz) -
		psi[i][j] * V[i][j] * rho[i] * drho*dz / 2.;

	dw = psi[i][j] * rho[i] * drho*dz / 2.;

	gE[i][j] = (dw*num - dv*denum) / (denum*denum);
	//

	for (j = 1; j < NUM2 - 1; j++)
	{
		dv = (rho[i+1] + rho[i])*(psi[i+1][j] - psi[i][j])*(-1.)*dz / (2.*drho) +
			rho[i]*(psi[i][j+1] - psi[i][j])*(-1.)*drho / (2.*dz) +
			rho[i] * (psi[i][j] - psi[i][j-1])*drho / (2.*dz) -
				 psi[i][j]*V[i][j]*rho[i]*drho*dz;

			dw = psi[i][j]*rho[i]*drho*dz;

			gE[i][j] = (dw*num - dv*denum) / (denum*denum);
	}
////////////
	i = NUM - 1;
	//
	j = 0;

	dv = (rho[i] + rho[i-1])*(psi[i][j] - psi[i-1][j])*dz / (4.*drho) +
		rho[i] * (psi[i][j + 1] - psi[i][j])*(-1.)*drho / (2.*dz) -
		psi[i][j] * V[i][j] * rho[i] * drho*dz / 2. +
		L*psi[i][j] * rho[i] * drho / 2.;

	dw = psi[i][j] * rho[i] * drho*dz / 2.;

	gE[i][j] = (dw*num - dv*denum) / (denum*denum);
	//
	j = NUM2 - 1;

	dv = (rho[i] + rho[i-1])*(psi[i][j] - psi[i-1][j])*dz / (4.*drho) +
		rho[i] * (psi[i][j] - psi[i][j - 1])*drho / (2.*dz) -
		psi[i][j] * V[i][j] * rho[i] * drho*dz / 2.;

	dw = psi[i][j] * rho[i] * drho*dz / 2.;

	gE[i][j] = (dw*num - dv*denum) / (denum*denum);
	//

	for (j = 1; j < NUM2 - 1; j++)
	{
		dv = (rho[i] + rho[i-1])*(psi[i][j] - psi[i-1][j])*dz / (2.*drho) +
			rho[i] * (psi[i][j + 1] - psi[i][j])*(-1.)*drho / (2.*dz) +
			rho[i] * (psi[i][j] - psi[i][j-1])*drho / (2.*dz) -
			psi[i][j] * V[i][j] * rho[i] * drho*dz;

		dw = psi[i][j] * rho[i] * drho*dz;

		gE[i][j] = (dw*num - dv*denum) / (denum*denum);
	}
//////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////
	j = 0;
	
	for (i = 1; i < NUM - 1; i++)
	{
		dv = (rho[i + 1] + rho[i])*(psi[i + 1][j] - psi[i][j])*(-1.)*dz / (4. * drho) +
			(rho[i] + rho[i - 1])*(psi[i][j] - psi[i - 1][j])*dz / (4. * drho) +
			rho[i] * (psi[i][j + 1] - psi[i][j])*(-1.)*drho / dz -
			psi[i][j] * V[i][j] * rho[i] * drho*dz + 
			L*psi[i][j]*rho[i]*drho;

		dw = psi[i][j] * rho[i] * drho*dz;

		gE[i][j] = (dw*num - dv*denum) / (denum*denum);
	}
	////////////
	j = NUM2 - 1;
	
	for (i = 1; i < NUM - 1; i++)
	{
		dv = (rho[i + 1] + rho[i])*(psi[i + 1][j] - psi[i][j])*(-1.)*dz / (4. * drho) +
			(rho[i] + rho[i - 1])*(psi[i][j] - psi[i - 1][j])*dz / (4. * drho) +
			rho[i] * (psi[i][j] - psi[i][j-1])*drho / dz -
			psi[i][j] * V[i][j] * rho[i] * drho*dz;

		dw = psi[i][j] * rho[i] * drho*dz;

		gE[i][j] = (dw*num - dv*denum) / (denum*denum);
	}
	//////////////////////////////////////////////////////////////////////////////

	for (i = 1; i < NUM - 1; i++)
	{
		for (j = 1; j < NUM2 - 1; j++)
		{
			dv = (rho[i+1] + rho[i])*(psi[i+1][j] - psi[i][j])*(-1.)*dz / (2*drho) + 
				 (rho[i] + rho[i-1])*(psi[i][j] - psi[i-1][j])*dz / (2 * drho) + 
				 rho[i]*(psi[i][j+1] - psi[i][j])*(-1.)*drho / dz +
				 rho[i] * (psi[i][j] - psi[i][j-1])*drho / dz -
				 2*psi[i][j]*V[i][j]*rho[i]*drho*dz;

			dw = 2.*psi[i][j]*rho[i]*drho*dz;

			gE[i][j] = (dw*num - dv*denum) / (denum*denum);
		}
	}
}

//нормировка ВФ на 1
void normto1()
{
	int i, j;
	for (i = 0; i < NUM; i++)
	{
		for (j = 0; j < NUM2; j++)
		{
			psi[i][j] = psi[i][j] / sqrt(denum);
		}
	}
	calc_E();
}

int main()
{
	int i, j, k, l, step, sec;
	double G, dt, F, v, n1, n2, n3;
	time_t tstart = time(NULL);
	
/*	
	for (i = 0; i < NUM; i++)
	{
		rho[i] = c*drho + drho*i;
	}

	for (j = 0; j < NUM2; j++)
	{
		z[j] = -a + j*dz;
	}
	printf("%.14le	%.14le\n", rho[0], rho[300]);
	printf("%.14le	%.14le	%.14le	%.14le	%.14le	%.14le\n", z[0], z[280], z[299], z[300], z[320], z[600]);
*/	

	
	for (i = 0; i < NUM; i++)
	{
		rho[i] = c*drho + drho*i*i*i/((NUM - 1.)*(NUM - 1.));
	}

	k = 0;
	for (j = NUM2/2; j < NUM2; j++)
	{
		z[j] = c*dz + dz*k*k*k / ((NUM - 1.)*(NUM - 1.));
		k++;
	}
	k = 0;
	for (j = (NUM2/2 - 1); j >= 0; j--)
	{
		z[j] = -c*dz - dz*k*k*k / ((NUM - 1.)*(NUM - 1.));
		k++;
	}
	printf("%.14le	%.14le\n", rho[0], rho[300]);
	printf("%.14le	%.14le	%.14le	%.14le	%.14le	%.14le\n", z[0], z[280], z[299], z[300], z[320], z[600]);
	
////////////////////////////////////////////////////////////////////////////////////////////////////////////
	prepfuncs();

	n1 = n2 = 0;
	step = 0;
	G = 0.995;
	dt = 8e-2;

	calc_E();
	calc_gE();

	maxgE = gE[0][0];
	for (i = 0; i < NUM; i++)
	{
		for (j = 0; j < NUM2; j++)
		{
			if (maxgE < fabs(gE[i][j]))
				maxgE = fabs(gE[i][j]);
			n1 = i; n2 = j;
		}
	}

	n3 = a + 1.;
	for (j = 0; j < NUM2; j++)
	{
		if (n3 > fabs(z[j]))
			n3 = fabs(z[j]);
	}
	printf("%.14le	%.14le	%.14le	%.14le	%.14le	%.14le\n", E, E_kin, E_pot, maxgE, rho[0], n3);



	while (maxgE > 1e-10)
	{

		for (k = 0; k < NUM; k++)
		{
			for (l = 0; l < NUM2; l++)
			{
				vel[k][l] = G*(vel[k][l] + gE[k][l] * dt);
				psi[k][l] = psi[k][l] + vel[k][l] * dt;
			}
		}

		calc_E();
		calc_gE();

		maxgE = gE[0][0];
		for (i = 0; i < NUM; i++)
		{
			for (j = 0; j < NUM2; j++)
			{
				if (maxgE < fabs(gE[i][j]))
				{
					maxgE = fabs(gE[i][j]);
					n1 = i; n2 = j;
				}
			}
		}

		if (step == 5000)
		{
			printf("%.14le	%.14le	%.14le	%.14le	%.14le	%.14le\n", E, E_kin, E_pot, maxgE, n1, n2);
			step = 0;
		}

		step++;
	}

	printf("%.14le	%.14le	%.14le	%.14le	%.14le	%.14le\n", E, E_kin, E_pot, maxgE, num, denum);
	normto1();
	printf("%.14le	%.14le	%.14le	%.14le	%.14le	%.14le\n", E, E_kin, E_pot, maxgE, num, denum);

	time_t tstop = time(NULL);
	sec = tstop - tstart;
	printf("spent_time = %ld\n", sec);
	FILE *f;
	f = fopen(OUTFILE, "w");
//	for (i = 0; i < NUM; i++)
		//fprintf(f, "%.14le	%.14le\n", rho[i], z[i]);

	return 0;
}




