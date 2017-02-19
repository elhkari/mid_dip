#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define OUTFILE "dip1out.txt"

#define a 10.
#define as 0.5 //длина отрезка с более мелким шагом
#define L 1.

#define NUM 301
#define NUM2 (2*NUM - 1)
#define NUMs 120 //кол-во точек с более мелким шагом

#define c 5.e-02 //коэф. определяет расстояние от 0 до ближайшей точки

#define hrhos (as/NUMs)
#define hzs (as/NUMs)

#define hrho ((a - as)/(NUM - NUMs))
#define hz ((a - as)/(NUM - NUMs))

#define G 0.995
#define dt 8.e-2

double rho[NUM], z[NUM2], drho[NUM], dz[NUM2], psi[NUM][NUM2], vel[NUM][NUM2], V[NUM][NUM2];

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

	for (i = 0; i < NUM-1; i++)
	{
		drho[i] = (rho[i + 1] - rho[i]);
	}
	for (j = 0; j < NUM2-1; j++)
	{
		dz[j] = (z[j + 1] - z[j]);
	}
}


void calc_E()
{
	double sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9, sum10, sum11, sum12;
	
	int i, j;

	sum1 = sum2 = sum3 = sum4 = sum5 = sum6 = sum7 = sum8 = sum9 = sum10 = sum11 = sum12 = 0;

	for (i = 0; i < NUM - 1; i++)
	{
		sum9 += psi[i][0] * psi[i][0] * rho[i]*drho[i];
		sum10 += psi[i + 1][0] * psi[i + 1][0] * rho[i + 1]*drho[i];
		
		for (j = 0; j < NUM2 - 1; j++)
		{
			sum1 += (psi[i + 1][j] - psi[i][j])*(psi[i + 1][j] - psi[i][j])*(rho[i + 1] + rho[i])*dz[j]/drho[i];
			sum2 += (psi[i + 1][j + 1] - psi[i][j + 1])*(psi[i + 1][j + 1] - psi[i][j + 1])*(rho[i + 1] + rho[i])*dz[j] / drho[i];

			sum3 += (psi[i][j + 1] - psi[i][j])*(psi[i][j + 1] - psi[i][j])*rho[i]*drho[i]/dz[j];
			sum4 += (psi[i + 1][j + 1] - psi[i + 1][j])*(psi[i + 1][j + 1] - psi[i + 1][j])*rho[i + 1] * drho[i] / dz[j];

			sum5 += psi[i][j] * psi[i][j] * rho[i] * V[i][j]*drho[i]*dz[j];
			sum6 += psi[i][j + 1] * psi[i][j + 1] * rho[i] * V[i][j+1] * drho[i]*dz[j];
			sum7 += psi[i + 1][j] * psi[i + 1][j] * rho[i + 1] * V[i+1][j] * drho[i]*dz[j];
			sum8 += psi[i + 1][j + 1] * psi[i + 1][j + 1] * rho[i + 1] * V[i+1][j+1] * drho[i]*dz[j];

			sum11 += (psi[i][j] * psi[i][j] + psi[i][j + 1] * psi[i][j + 1]) * rho[i] * drho[i]*dz[j];
			sum12 += (psi[i + 1][j] * psi[i + 1][j] + psi[i + 1][j + 1] * psi[i + 1][j + 1]) * rho[i + 1] * drho[i]*dz[j];
		}
	}


	sum1 = sum1/8.;
	sum2 = sum2/8.;
	
	sum3 = sum3 /4.;
	sum4 = sum4 /4.;
	
	sum5 = sum5 / 4.;
	sum6 = sum6 / 4.;
	sum7 = sum7 / 4.;
	sum8 = sum8 / 4.;

	sum9 = L*sum9 / 4.;
	sum10 = L*sum10 / 4.;

	sum11 = sum11 / 4.;
	sum12 = sum12 / 4.;

	E_kin = sum1 + sum2 + sum3 + sum4;
	E_pot = sum5 + sum6 + sum7 + sum8;
	C = sum9 + sum10;
	num = E_kin - E_pot + C;
	denum = sum11 + sum12;

	E = num / denum;
}

void clear_gE()
{
	int i, j;
	for (i = 0; i < NUM; i++)
	{
		for (j = 0; j < NUM2; j++)
		{
			gE[i][j] = 0;
		}
	}
}

void calc_gE()
{
	double tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12;
	double dv1, dv2, dv3, dv4, dv5, dv6, dv7, dv8;
	double dw1, dw2, dw3, dw4;
	int i, j;

	clear_gE();
	
	//ДАЛЬШЕ НАДО ПОМНИТЬ, ЧТО МЫ ЗАПИСЫВАЕМ (-1)*ГРАДИЕНТ
	//gE = (dw*num - dv*denum) / (denum*denum);

	for (i = 0; i < NUM-1; i++)
	{
		gE[i][0] -= L*psi[i][0] * rho[i] * drho[i] / (2.*denum);
		gE[i+1][0] -= L*psi[i+1][0] * rho[i+1] * drho[i] / (2.*denum);

		for (j = 0; j < NUM2-1; j++)
		{
			//ПРОИЗВОДНАЯ ЧИСЛИТЕЛЯ

			//КИНЕТИЧЕСКИЕ СЛАГАЕМЫЕ
			tmp1 = (rho[i + 1] + rho[i])*(psi[i + 1][j] - psi[i][j])*dz[j] / (4. * drho[i]);
			// - (i,j) //+ (i+1,j) 
			dv1 = -tmp1 / denum;
			gE[i][j] -= dv1;
			gE[i + 1][j] += dv1;
			///////////////////////////////////////////////
			
			tmp2 = (rho[i + 1] + rho[i])*(psi[i + 1][j + 1] - psi[i][j + 1])*dz[j] / (4. * drho[i]); 
			// - (i,j+1) //+(i+1,j+1)
			dv2 = -tmp2 / denum;
			gE[i][j+1] -= dv2;
			gE[i + 1][j+1] += dv2;
			///////////////////////////////////////////////
			
			tmp3 = rho[i] * (psi[i][j + 1] - psi[i][j])*drho[i] / (2. * dz[j]);
			//-(i,j) //+(i,j+1)
			dv3 = -tmp3 / denum;
			gE[i][j] -= dv3;
			gE[i][j + 1] += dv3;

			/////////////////////////////////////////////

			tmp4 = rho[i+1] * (psi[i+1][j + 1] - psi[i+1][j])*drho[i] / (2. * dz[j]);
			//-(i+1,j) //+(i+1,j+1)
			dv4 = -tmp4 / denum;
			gE[i+1][j] -= dv4;
			gE[i+1][j + 1] += dv4;

			//ПОТЕНЦИАЛЬНЫЕ СЛАГАЕМЫЕ
			tmp5 = psi[i][j] * V[i][j] * rho[i] * drho[i] * dz[j] / 2.;
			dv5 = tmp5 / denum;
			gE[i][j] += dv5;

			tmp6 = psi[i][j+1] * V[i][j+1] * rho[i] * drho[i] * dz[j] / 2.;
			dv6 = tmp6 / denum;
			gE[i][j+1] += dv6;

			tmp7 = psi[i+1][j] * V[i+1][j] * rho[i+1] * drho[i] * dz[j] / 2.;
			dv7 = tmp7 / denum;
			gE[i+1][j] += dv7;

			tmp8 = psi[i+1][j+1] * V[i+1][j+1] * rho[i+1] * drho[i] * dz[j] / 2.;
			dv8 = tmp8 / denum;
			gE[i+1][j+1] += dv8;

			//ПРОИЗВОДНАЯ ЗНАМЕНАТЕЛЯ
			tmp9 = psi[i][j] * rho[i] * drho[i] * dz[j] / 2.;
			dw1 = tmp9*num / (denum*denum);
			gE[i][j] += dw1;
			
			tmp10 = psi[i][j+1] * rho[i] * drho[i] * dz[j] / 2.;
			dw2 = tmp10*num / (denum*denum);
			gE[i][j+1] += dw2;
			
			tmp11 = psi[i+1][j] * rho[i+1] * drho[i] * dz[j] / 2.;
			dw3 = tmp11*num / (denum*denum);
			gE[i+1][j] += dw3;

			tmp12 = psi[i+1][j+1] * rho[i+1] * drho[i] * dz[j] / 2.;
			dw4 = tmp12*num / (denum*denum);
			gE[i+1][j+1] += dw4;
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
	//double G, dt, F, v, n1, n2, n3;
	double F, v, n1, n2, n3;
	time_t tstart = time(NULL);
	
	/*
	for (i = 0; i < NUM; i++)
	{
		rho[i] = c*hrho + hrho*i;
	}

	for (j = 0; j < NUM2; j++)
	{
		z[j] = -a + j*hz;
	}
	printf("%.14le	%.14le\n", rho[0], rho[300]);
	printf("%.14le	%.14le	%.14le	%.14le	%.14le	%.14le\n", z[0], z[280], z[299], z[300], z[320], z[600]);
	*/

	
	for (i = 0; i < NUMs; i++)
	{
		rho[i] = c*hrhos + hrhos*i;
	}
	for (i = NUMs; i < NUM; i++)
	{
		rho[i] = rho[i-1] + hrho;
	}

	k = 1;
	z[NUM2 / 2] = c*hzs;
	for (j = (NUM2/2 + 1); j < (NUM2/2 + NUMs); j++)
	{
		z[j] = c*hzs + hzs*k;
		k++;
	}
	for (j = (NUM2 / 2 + NUMs); j < NUM2; j++)
	{
		z[j] = z[j-1] + hz;
	}
	k = 2;
	for (j = (NUM2/2 - 1); j >= 0; j--)
	{
		z[j] = -z[j+k];
		k= k+2;
	}

	/*
	for (j = 0; j < (NUM2 / 2 - NUMs); j++)
	{
	z[j] = -a + hz*j;
	k++;
	}
	for (j = (NUM2 / 2 - NUMs); j < (NUM2 / 2 + NUMs); j++)
	{
	z[j] = z[j - 1] + hzs;
	}
	for (j = (NUM2 / 2 + NUMs); j < NUM2; j++)
	{
	z[j] = z[j - 1] + hz;
	}
	*/

	printf("%.14le	%.14le\n", rho[0], rho[300]);
	printf("%.14le	%.14le	%.14le	%.14le	%.14le	%.14le\n\n", z[NUM2 / 2 -1],
		z[NUM2 / 2], z[NUM2 / 2 + 1], z[NUM2 / 2 + NUMs -1], z[NUM2 / 2 + NUMs], z[NUM2 / 2 + NUMs + 1]);
	

////////////////////////////////////////////////////////////////////////////////////////////////////////////
	prepfuncs();

	n1 = n2 = 0;
	step = 0;
	//G = 0.995;
	//dt = 8e-2;

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

	normto1();
	printf("%.14lf	%.14lf	%.14lf	%.14le\n", E, E_kin, E_pot, maxgE);


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

		if (step == 1000)
		{
			normto1();
			printf("%.14lf	%.14lf	%.14lf	%.14le	%.14lf	%.14lf\n", E, E_kin, E_pot, maxgE, n1, n2);
			step = 0;
		}

		step++;
	}

	printf("%.14lf	%.14lf	%.14lf	%.14le	%.14lf	%.14lf\n", E, E_kin, E_pot, maxgE, num, denum);
	normto1();
	printf("%.14lf	%.14lf	%.14lf	%.14le	%.14lf	%.14lf\n", E, E_kin, E_pot, maxgE, num, denum);

	time_t tstop = time(NULL);
	sec = tstop - tstart;
	printf("spent_time = %ld\n\a", sec);
/*	FILE *f;
	f = fopen(OUTFILE, "w");
	for (i = 0; i < NUM; i++)
		fprintf(f, "%.14le	%.14le\n", rho[i], z[i]);
*/
	return 0;
}