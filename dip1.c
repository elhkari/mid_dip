//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define OUTFILE "nt2out.txt"

#define a 10.
//#define b 10.
#define L 1.
#define NUM 301
#define NUM2 601

#define drho (a/NUM)
#define dz (a/NUM)

double rho[NUM], z[NUM2], psi[NUM][NUM2], V[NUM][NUM2], 
dpsi_rho[NUM][NUM], dpsi_z[NUM][NUM], dR_rho[NUM][NUM], dR_z[NUM][NUM];

double E, num, denum, maxgE, maxgE_rho, maxgE_z;
double gE_rho[NUM], gE_z[NUM];

void prepfuncs()
{
	int i, j;
	for (i = 0; i < NUM; i++)
	{
		for (j = 0; j < NUM2; j++)
		{
			V[i][j] = 1./(sqrt(rho[i] * rho[i] + z[j] * z[j]));
			psi[i][j] = exp(-R[i][j]);
			dpsi_rho[i][j] = psi[i][j] * (-rho[i] / R[i][j]);
			dpsi_z[i][j] = psi[i][j] * (-z[j] / R[i][j]);
			dR_rho[i][j] = rho[i] / R[i][j];
			dR_z[i][j] = z[j] / R[i][j];
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
		
		for (j = 0; j < NUM - 1; j++)
		{
			sum1 += (psi[i + 1][j] - psi[i][j])*(psi[i + 1][j] - psi[i][j])*(rho[i + 1] + rho[i]);
			sum2 += (psi[i + 1][j + 1] - psi[i][j + 1])*(psi[i + 1][j + 1] - psi[i][j + 1])*(rho[i + 1] + rho[i]);

			sum3 += (psi[i][j + 1] - psi[i][j])*(psi[i][j + 1] - psi[i][j])*rho[i];
			sum4 += (psi[i + 1][j + 1] - psi[i + 1][j])*(psi[i + 1][j + 1] - psi[i + 1][j])*rho[i + 1];

			sum5 += psi[i][j] * psi[i][j] * rho[i] / R[i][j];
			sum6 += psi[i][j + 1] * psi[i][j + 1] * rho[i] / R[i][j+1];
			sum7 += psi[i + 1][j] * psi[i + 1][j] * rho[i + 1] / R[i+1][j];
			sum8 += psi[i + 1][j + 1] * psi[i + 1][j + 1] * rho[i + 1] / R[i+1][j+1];

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

void calc_gE_rho()
{
	double dv, dw, sum1, sum2, sum3, sum4, sum5;
	double sum1a[NUM], sum3a[NUM], sum5a[NUM];
	int i, j;

	

	for (i = 1; i < NUM - 1; i++)
	{
		sum1 = sum2 = sum3 = sum4 = sum5 = 0;

		for (j = 1; j < NUM-1; j++)
		{
			/*
			sum1a[j] = 2. * (psi[i + 1][j] - psi[i][j])*(rho[i + 1] + rho[i])*(-dpsi_rho[i][j]) +
				(psi[i + 1][j] - psi[i][j])*(psi[i + 1][j] - psi[i][j]) +
				2. * (psi[i][j] - psi[i - 1][j])*dpsi_rho[i][j] +
				(psi[i][j] - psi[i - 1][j])*(psi[i][j] - psi[i - 1][j]);

			sum2 += 2. * (psi[i][j + 1] - psi[i][j])*(dpsi_rho[i][j + 1] - dpsi_rho[i][j])*rho[i] +
				(psi[i][j + 1] - psi[i][j])*(psi[i][j + 1] - psi[i][j]);

			sum3a[j] = 2. * psi[i][j] * rho[i] * dpsi_rho[i][j] / R[i][j] -
				psi[i][j] * psi[i][j] * rho[i] * dR_rho[i][j] / R[i][j] / R[i][j] +
				psi[i][j] * psi[i][j] / R[i][j];

			sum5a[j] = 2. * psi[i][j] * rho[i] * dpsi_rho[i][j] + psi[i][j]*psi[i][j];
			*/

			sum1 += 2. * (psi[i + 1][j] - psi[i][j])*(rho[i + 1] + rho[i])*(-dpsi_rho[i][j]) +
				(psi[i + 1][j] - psi[i][j])*(psi[i + 1][j] - psi[i][j]) +
				2. * (psi[i][j] - psi[i - 1][j])*dpsi_rho[i][j] +
				(psi[i][j] - psi[i - 1][j])*(psi[i][j] - psi[i - 1][j]);

			sum2 += 2. * (psi[i][j + 1] - psi[i][j])*(dpsi_rho[i][j + 1] - dpsi_rho[i][j])*rho[i] +
				(psi[i][j + 1] - psi[i][j])*(psi[i][j + 1] - psi[i][j]);

			sum3 += 2. * psi[i][j] * rho[i] * dpsi_rho[i][j] / R[i][j] -
				psi[i][j] * psi[i][j] * rho[i] * dR_rho[i][j] / R[i][j] / R[i][j] +
				psi[i][j] * psi[i][j] / R[i][j];

			sum5 += 2. * psi[i][j] * rho[i] * dpsi_rho[i][j] + psi[i][j] * psi[i][j];
		}
		sum4 += 2. * psi[i][0] * dpsi_rho[i][0] + psi[i][0] * psi[i][0];

		sum1 = 2. * sum1;
		sum3 = 2. * sum3;
		sum5 = 2. * sum5;

		/*
		for (j = 1; j < NUM - 1; j++)
		{
			
			sum1 += 2. * sum1a[j];
			sum3 += 2. * sum3a[j];
			sum5 += 2. * sum5a[j];
		}

		j = NUM - 1;
		sum1a[j] = 2. * (psi[i + 1][j] - psi[i][j])*(rho[i + 1] + rho[i])*(-dpsi_rho[i][j]) +
			(psi[i + 1][j] - psi[i][j])*(psi[i + 1][j] - psi[i][j]) +
			2. * (psi[i][j] - psi[i - 1][j])*dpsi_rho[i][j] +
			(psi[i][j] - psi[i - 1][j])*(psi[i][j] - psi[i - 1][j]);

		sum3a[j] = 2. * psi[i][j] * rho[i] * dpsi_rho[i][j] / R[i][j] -
			psi[i][j] * psi[i][j] * rho[i] * dR_rho[i][j] / R[i][j] / R[i][j] +
			psi[i][j] * psi[i][j] / R[i][j];

		sum5a[j] = 2. * psi[i][j] * rho[i] * dpsi_rho[i][j] + psi[i][j] * psi[i][j];
		

		sum1 += sum1a[0] + sum1a[NUM - 1];
		sum3 += sum3a[0] + sum3a[NUM - 1];
		sum5 += sum5a[0] + sum5a[NUM - 1];
		*/

		sum1 = sum1*dz / 8. / drho;
		sum2 = sum2*drho / 2. / dz;
		sum3 = sum3*drho*dz / 2.;
		sum4 = sum4*drho*L / 2.;
		sum5 = sum5*drho*dz / 4.;

		dv = sum1 + sum2 - sum3 + sum4;
		dw = sum5;

		gE_rho[i] = (dw*num - dv*denum) / (denum*denum);
	}
}

void calc_gE_z()
{
	double dv, dw, sum1, sum2, sum3, sum4;
	double sum2a[NUM], sum3a[NUM], sum4a[NUM];
	int i, j;

	

	for (j = 1; j < NUM - 1; j++)
	{
		sum1 = sum2 = sum3 = sum4 = 0;
		for (i = 1; i < NUM - 1; i++)
		{
			sum1 += (psi[i + 1][j] - psi[i][j])*(dpsi_z[i + 1][j] - dpsi_z[i][j])*(rho[i + 1] + rho[i]);

			sum2 += rho[i] * ((psi[i][j + 1] - psi[i][j])*(-dpsi_z[i][j]) +
				(psi[i + 1][j] - psi[i][j])*(psi[i + 1][j] - psi[i][j]) +
				(psi[i][j] - psi[i][j - 1])*dpsi_z[i][j]);

			sum3 += rho[i] * (2. * psi[i][j] * dpsi_z[i][j] / R[i][j] -
				psi[i][j] * psi[i][j] * dR_z[i][j] / R[i][j] / R[i][j]);

			sum4 += psi[i][j] * rho[i] * dpsi_z[i][j];

			/*
			sum1 += (psi[i + 1][j] - psi[i][j])*(dpsi_z[i + 1][j] - dpsi_z[i][j])*(rho[i + 1] + rho[i]);
				
			sum2a[i] = rho[i] * ((psi[i][j + 1] - psi[i][j])*(-dpsi_z[i][j]) +
				(psi[i + 1][j] - psi[i][j])*(psi[i + 1][j] - psi[i][j]) +
				(psi[i][j] - psi[i][j - 1])*dpsi_z[i][j]);

			sum3a[i] = rho[i]*(2. * psi[i][j] * dpsi_z[i][j] / R[i][j] -
				psi[i][j] * psi[i][j] * dR_z[i][j] / R[i][j] / R[i][j]);

			sum4a[i] = psi[i][j] * rho[i] * dpsi_z[i][j];
			*/
		}
		sum2 = 2. * sum2;
		sum3 = 2. * sum3;
		sum4 = 2. * sum4;
		/*
		for (i = 1; i < NUM - 1; i++)
		{
			sum2 += 2. * sum2a[j];
			sum3 += 2. * sum3a[j];
			sum4 += 2. * sum4a[j];
		}

		i = NUM - 1;
		sum2a[i] = rho[i] * ((psi[i][j + 1] - psi[i][j])*(-dpsi_z[i][j]) +
			(psi[i + 1][j] - psi[i][j])*(psi[i + 1][j] - psi[i][j]) +
			(psi[i][j] - psi[i][j - 1])*dpsi_z[i][j]);

		sum3a[i] = rho[i] * (2. * psi[i][j] * dpsi_z[i][j] / R[i][j] -
			psi[i][j] * psi[i][j] * dR_z[i][j] / R[i][j] / R[i][j]);

		sum4a[i] = psi[i][j] * rho[i] * dpsi_z[i][j];

		sum2 += sum2a[0] + sum2a[NUM - 1];
		sum3 += sum3a[0] + sum3a[NUM - 1];
		sum4 += sum4a[0] + sum4a[NUM - 1];
		*/

		sum1 = sum1*dz / 2. / drho;
		sum2 = sum2*drho / 2. / dz;
		sum3 = sum3*drho*dz / 2.;
		sum4 = sum4*drho*dz / 2.;

		dv = sum1 + sum2 - sum3;
		dw = sum4;

		gE_z[i] = (dw*num - dv*denum) / (denum*denum);
	}
}

int main()
{
	int i, j, k, step, sec;
	double G, dt, cosT, F, v;
	time_t tstart = time(NULL);

	for (i = 0; i < NUM; i++)
	{
		rho[i] = -b +  drho*i;
		vel_rho[i] = 0.;
	}

	for (j = 0; j < NUM; j++)
	{
		z[j] = -a + j*dz;
		vel_z[i] = 0.;
	}

	G = 0.995;
	dt = 8e-3;

	calc_E();
	calc_gE_rho();
	calc_gE_z();

	maxgE_rho = gE_rho[1];
	for (i = 2; i < NUM - 1; i++)
	{
		if (maxgE_rho < fabs(gE_rho[i]))
			maxgE_rho = fabs(gE_rho[i]);
	}

	maxgE_z = gE_z[1];
	for (j = 2; j < NUM - 1; j++)
	{
		if (maxgE_z < fabs(gE_z[j]))
			maxgE_z = fabs(gE_z[j]);
	}

	if (maxgE_rho < maxgE_z)
	{
		maxgE = maxgE_z;
	}
	else
	{
		maxgE = maxgE_rho;
	}

	printf("%.14le		%.14le\n", E, maxgE);



	while (maxgE > 1e-10)
	{

		for (k = 1; k < NUM - 1; k++)
		{
			vel_rho[k] = G*(vel_rho[k] + gE_rho[k] * dt);
			rho[k] = rho[k] + vel_rho[k] * dt;
		}

		for (k = 1; k < NUM - 1; k++)
		{
			vel_z[k] = G*(vel_z[k] + gE_z[k] * dt);
			z[k] = z[k] + vel_z[k] * dt;
		}

		calc_E();
		calc_gE_rho();
		calc_gE_z();

		maxgE_rho = gE_rho[1];
		for (i = 2; i < NUM - 1; i++)
		{
			if (maxgE_rho < fabs(gE_rho[i]))
				maxgE_rho = fabs(gE_rho[i]);
		}

		maxgE_z = gE_z[1];
		for (j = 2; j < NUM - 1; j++)
		{
			if (maxgE_z < fabs(gE_z[j]))
				maxgE_z = fabs(gE_z[j]);
		}

		if (maxgE_rho < maxgE_z)
		{
			maxgE = maxgE_z;
		}
		else
		{
			maxgE = maxgE_rho;
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
	for (i = 0; i < NUM; i++)
		fprintf(f, "%.14le	%.14le\n", rho[i], z[i]);

	return 0;
}




