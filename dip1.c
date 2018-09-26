#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define OUTFILE "dip1out.txt"

#define a 10.
//#define z0 10.	//расстояние от 0 до плоскости
//#define as  3.0*1.0	//длина приращения последнего шага
#define L (0.)

#define NUM 100
#define NUM2 (2*NUM)

#define h1  (a/((NUM-1)+((NUM-2)*(NUM-1)*hrho)/2.))  // расстояние от нуля до первой точки
#define c 0.699

#define x 8. //на сколько поднимаем плоскость
#define M ((int)(x/hz))
//#define h2  (z0/(NUM+((NUM-1)*NUM*hrho)/2.))  // расстояние от нуля до первой точки (в сторону плоскости) 

#define hrho (as/(NUM-2))
#define hz (a/(NUM-1+c))
//#define hz0 as/(NUM-2)

#define G 0.995
#if NUM == 100
	#define dt 2.e-1
	#define as  (3.0*1.0)
#elif NUM == 150
	#define dt 2.e-1
	#define as  (3.0*1.5)
#elif NUM == 200
	#define dt 2.e-1
	#define as  (3.0*2.0)
#elif NUM == 300
	#define dt 2.e-1
	#define as  (3.0*3.0)
#endif

double rho[NUM], z[NUM2], drho[NUM], dz[NUM2], psi[NUM][NUM2], vel[NUM][NUM2], V[NUM][NUM2];

double E, num, denum, maxgE, E_kin, E_pot, C;
double gE[NUM][NUM2];

void prepfuncs()
{
	int i, j;
	for (i = 0; i < NUM; i++)
	{
		for (j = M; j < NUM2; j++)
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
	for (j = M; j < NUM2-1; j++)
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
		sum9 += psi[i][M] * psi[i][M] * rho[i]*drho[i];
		sum10 += psi[i + 1][M] * psi[i + 1][M] * rho[i + 1]*drho[i];
		
		for (j = M; j < NUM2 - 1; j++)
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
		gE[i][M] -= L*psi[i][M] * rho[i] * drho[i] / (2.*denum);
		gE[i+1][M] -= L*psi[i+1][M] * rho[i+1] * drho[i] / (2.*denum);

		for (j = M; j < NUM2-1; j++)
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
		for (j = M; j < NUM2; j++)
		{
			psi[i][j] = psi[i][j] / sqrt(denum);
		}
	}
	calc_E();
}

int main()
{
	int i, j, k, l, step, sec;
	double F, v, n1, n2, n3;
	time_t tstart = time(NULL);
	
	rho[0] = 0;
	for (i = 1; i < NUM; i++)
	{
		rho[i] = rho[i-1] + h1*(1+(i-1)*hrho);
	}
	
	k = 1;
	z[NUM2 / 2] = c*hz;
	for (j = NUM2 / 2 + 1; j < NUM2; j++)
	{
		z[j] = z[j-1] + hz;
		k++;
	}

	l = 1;
	z[NUM2 / 2 - 1] = -c*hz;
	for (j = (NUM2 / 2 - 2); j >= M; j--)
	{
		z[j] = z[j+1]-hz;
		l++;
	}
	printf("M = %d	hz = %.14le	h1 = %.14le	dh = %.14le\n", M, hz, h1, hrho);
	printf("%.14le	%.14le\n", rho[0], rho[NUM - 1]);
	printf("%.14le	%.14le	%.14le	%.14le	%.14le	%.14le\n\n", z[M],
		z[NUM2 / 2 -1], z[NUM2 / 2], z[NUM2 / 2 + 1], z[NUM2 -2], z[NUM2-1]);
	

////////////////////////////////////////////////////////////////////////////////////////////////////////////
	prepfuncs();

	n1 = n2 = 0;
	step = 0;

	calc_E();
	calc_gE();

	maxgE = gE[0][M];
	for (i = 0; i < NUM; i++)
	{
		for (j = M; j < NUM2; j++) //пров
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
			for (l = M; l < NUM2; l++)
			{
				vel[k][l] = G*(vel[k][l] + gE[k][l] * dt);
				psi[k][l] = psi[k][l] + vel[k][l] * dt;
			}
		}

		calc_E();
		calc_gE();

		maxgE = gE[0][M];
		for (i = 0; i < NUM; i++)
		{
			for (j = M; j < NUM2; j++)
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
