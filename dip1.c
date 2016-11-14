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

#define drho (a/(NUM - 0.5))
#define dz (a/(NUM - 0.5))

double rho[NUM], z[NUM2], psi[NUM][NUM2], vel[NUM][NUM2], V[NUM][NUM2];

double E, num, denum, maxgE;
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


	sum1 = dz*sum1/(4.*drho);
	sum2 = dz*sum2/(4.*drho);
	
	sum3 = drho*sum3 /(4.* dz);
	sum4 = drho*sum4 /(4.* dz);
	
	sum5 = drho*dz*sum5 / 4.;
	sum6 = drho*dz*sum6 / 4.;
	sum7 = drho*dz*sum7 / 4.;
	sum8 = drho*dz*sum8 / 4.;

	sum9 = drho*L*sum9 / 4.;
	sum10 = drho*L*sum10 / 4.;

	sum11 = drho*dz*sum11 / 8.;
	sum12 = drho*dz*sum12 / 8.;

	num = sum1 + sum2 + sum3 + sum4 - sum5 - sum6 - sum7 - sum8 + sum9 + sum10;
	denum = sum11 + sum12;

	E = num / denum;
}

void calc_gE()
{
	double dv, dw;
	int i, j;

	for (i = 1; i < NUM - 1; i++)
	{
		for (j = 1; j < NUM2 - 1; j++)
		{
			dv = (rho[i+1] + rho[i])*(psi[i+1][j] - psi[i][j])*(-1.)*dz / (2*drho) + 
				 (rho[i] + rho[i-1])*(psi[i][j] - psi[i-1][j])*dz / (2 * drho) + 
				 rho[i]*(psi[i][j+1] - psi[i][j])*(-1.)*drho / dz +
				 rho[i] * (psi[i][j] - psi[i][j-1])*drho / dz -
				 2*psi[i][j]*V[i][j]*rho[i]*drho*dz;

			dw = psi[i][j]*rho[i]*drho*dz;

			gE[i][j] = (dw*num - dv*denum) / (denum*denum);
		}
	}
}


int main()
{
	int i, j, k, l, step, sec;
	double G, dt, cosT, F, v;
	time_t tstart = time(NULL);

	for (i = 0; i < NUM; i++)
	{
		rho[i] = 0.5*drho + drho*i;
	}

	for (j = 0; j < NUM2; j++)
	{
		z[j] = -a + j*dz;
	}

	prepfuncs();

	step = 0;
	G = 0.995;
	dt = 8e-3;

	calc_E();
	calc_gE();

	maxgE = gE[1][1];
	for (i = 1; i < NUM - 1; i++)
	{
		for (j = 1; j < NUM2 - 1; j++)
		{
			if (maxgE < fabs(gE[i][j]))
				maxgE = fabs(gE[i][j]);
		}
	}

	printf("%.14le		%.14le\n", E, maxgE);



	while (maxgE > 1e-10)
	{

		for (k = 1; k < NUM - 1; k++)
		{
			for (l = 1; l < NUM2 - 1; l++)
			{
				vel[k][l] = G*(vel[k][l] + gE[k][l] * dt);
				psi[k][l] = psi[k][l] + vel[k][l] * dt;
			}
		}

		calc_E();
		calc_gE();

		maxgE = gE[1][1];
		for (i = 1; i < NUM - 1; i++)
		{
			for (j = 1; j < NUM2 - 1; j++)
			{
				if (maxgE < fabs(gE[i][j]))
					maxgE = fabs(gE[i][j]);
			}
		}

		if (step == 5000)
		{
			printf("%.14le		%.14le\n", E, maxgE);
			step = 0;
		}

		step++;
	}

	time_t tstop = time(NULL);
	sec = tstop - tstart;
	printf("spent_time = %ld\n", sec);
	FILE *f;
	f = fopen(OUTFILE, "w");
//	for (i = 0; i < NUM; i++)
		//fprintf(f, "%.14le	%.14le\n", rho[i], z[i]);

	return 0;
}




