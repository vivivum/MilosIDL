#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "idl_export.h"		/* IDL external definitions */

#define NTERMS 10

#define PI 3.141592653
#define LIMITE_INFERIOR_PRECISION_SVD pow(2.0,-39)
#define LIMITE_INFERIOR_PRECISION_TRIG pow(2.0,-39)
#define LIMITE_INFERIOR_PRECISION_SINCOS pow(2.0,-39)
#define LIMITA_PRECISION_fxp(number) ( fabs(number) > LIMITE_INFERIOR_PRECISION_SVD ? number : 0 )
#define LIMITA_PRECISION_trig(number) ( fabs(number) > LIMITE_INFERIOR_PRECISION_TRIG ? number : 0 )
#define LIMITA_PRECISION_SINCOS(number) ( fabs(number) > LIMITE_INFERIOR_PRECISION_TRIG ? number : 0 )

void ROTACION(double angulo,double in1_1,double in1_2,double in2_1,double in2_2,double* out1_1,double* out1_2,double* out2_1,double* out2_2){

	double seno, coseno;

	seno=LIMITA_PRECISION_SINCOS(sin(angulo));
	coseno=LIMITA_PRECISION_SINCOS(cos(angulo));

	//rotacion vector columna 1
	(*out1_1)=LIMITA_PRECISION_fxp(in1_1*coseno-in2_1*seno);
	(*out2_1)=LIMITA_PRECISION_fxp(in1_1*seno  +in2_1*coseno);

	//rotacion vector columna 2
	(*out1_2)=LIMITA_PRECISION_fxp(in1_2*coseno-in2_2*seno);
	(*out2_2)=LIMITA_PRECISION_fxp(in1_2*seno  +in2_2*coseno);

}

double calculaAngulo(double in1_1,double in1_2,double in2_2)
{

	double angulo,in1_2b,dif,arc;

	if(in1_1==in2_2){  //Si los elementos de la diagonal son iguales rotamos PI/4.
		angulo=LIMITA_PRECISION_fxp(PI/4);
	}
	else{
		in1_2b=LIMITA_PRECISION_fxp(2*in1_2);
		dif=LIMITA_PRECISION_fxp(in2_2-in1_1);

		arc=in1_2b/dif;

		angulo=LIMITA_PRECISION_trig(atan(arc))/2;

//		printf("\nAngulo: %e \n",angulo);

		if(angulo>(PI/4)){
			angulo=LIMITA_PRECISION_fxp(PI/4);
		}

		if(angulo<(-PI/4)){
			angulo=LIMITA_PRECISION_fxp(-PI/4);
		}

	}

	return angulo;
}

double PD(double in1_1,double in1_2,double in2_2,double *out1_1,double *out1_2,double *out2_2){

	double angulo,angulo1;
	double kk;

	angulo=calculaAngulo(in1_1,in1_2,in2_2);

	ROTACION(angulo,in1_1,in1_2,in1_2,in2_2,out1_1,out1_2,&kk,out2_2);
	ROTACION(angulo,*out1_1,kk,*out1_2,*out2_2,out1_1,out1_2,&kk,out2_2);

	//forzamos los no diagonles a cero
	*out1_2=0;

	return angulo;
}

void PND(double angulofila,double angulocolumna,double in1_1,double in1_2,double in2_1,double in2_2,double *out1_1,double *out1_2,double *out2_1,double *out2_2){

	double aux;

	ROTACION(angulofila,in1_1,in1_2,in2_1,in2_2,out1_1,out1_2,out2_1,out2_2);
	ROTACION(angulocolumna,*out1_1,*out2_1,*out1_2,*out2_2,out1_1,out1_2,out2_1,out2_2);

	aux=*out1_2;
	*out1_2=*out2_1;
	*out2_1=aux;

}


void PV(double angulo,double in1_1,double in1_2,double in2_1,double in2_2,double* out1_1,double* out1_2,double* out2_1,double* out2_2){

	double aux,kk1,kk2,kk3,kk4;

	ROTACION(angulo,in1_1,in2_1,in1_2,in2_2,&kk1,&kk2,&kk3,&kk4);

	*out1_1 = kk1;
	*out1_2 = kk3;
	*out2_1 = kk2;
	*out2_2 = kk4;
}


void InicializarMatrizIdentidad(double A[NTERMS][NTERMS],int tamanio){
	int i=0,j;

	for(i=0;i<tamanio;i++){
		for(j=0;j<tamanio;j++)
			A[i][j] = 0;
	}

	for(i=0;i<tamanio;i++){
			A[i][i] = 1;
	}
}


void InicializarMatrizCeros(double A[NTERMS][NTERMS],int tamanio){
	int i=0,j;

	for(i=0;i<tamanio;i++){
		for(j=0;j<tamanio;j++)
			A[i][j]=0;
	}

}

int multmatrix(double *a,int naf,int nac, double *b,int nbf,int nbc,double *result,int *fil,int *col){

    int i,j,k;
    double sum;

	if(nac==nbf){
		(*fil)=naf;
		(*col)=nbc;

//		free(*result);
//		*result=calloc((*fil)*(*col),sizeof(double));
//		printf("a ver ..\n");

		for ( i = 0; i < naf; i++)
		    for ( j = 0; j < nbc; j++){
				sum=0;
				for ( k = 0;  k < nbf; k++){
//					printf("i: %d,j:%d,k=%d .. a[%d][%d]  .. b[%d][%d]\n",i,j,k,i,k,k,j);
					sum += a[i*nac+k] * b[k*nbc+j];
				}
//				printf("Sum\n");
				result[(*col)*i+j] = sum;

      		}

		return 1;
	}
	return 0;

}

void svdcordic(double a[NTERMS*NTERMS], double w[NTERMS], double v[NTERMS][NTERMS],int max_iter)
{

	int m = NTERMS;
	int n = NTERMS;
  int tamanio = NTERMS;
	int ind_f,ind_c,i,j;
	int iter;

	double angulos[(int)NTERMS/2];
	double* a_out=calloc(NTERMS*NTERMS,sizeof(double));
	double autovectores[NTERMS][NTERMS];
	double autovectores_out[NTERMS][NTERMS];
	double kk1,kk2,kk3,kk4;

	InicializarMatrizIdentidad(autovectores,NTERMS);
	InicializarMatrizCeros(autovectores_out,NTERMS);

	int posFila[NTERMS][NTERMS]={0,0,0,0,0,0,0,0,0,0,
											2,2,2,2,2,2,2,2,2,2,
											4,4,4,4,4,4,4,4,4,4,
											1,1,1,1,1,1,1,1,1,1,
											6,6,6,6,6,6,6,6,6,6,
											3,3,3,3,3,3,3,3,3,3,
											8,8,8,8,8,8,8,8,8,8,
											5,5,5,5,5,5,5,5,5,5,
											9,9,9,9,9,9,9,9,9,9,
											7,7,7,7,7,7,7,7,7,7};

	int posCol[NTERMS][NTERMS]={0,2,4,1,6,3,8,5,9,7,
										  0,2,4,1,6,3,8,5,9,7,
  										0,2,4,1,6,3,8,5,9,7,
										  0,2,4,1,6,3,8,5,9,7,
										  0,2,4,1,6,3,8,5,9,7,
										  0,2,4,1,6,3,8,5,9,7,
										  0,2,4,1,6,3,8,5,9,7,
										  0,2,4,1,6,3,8,5,9,7,
										  0,2,4,1,6,3,8,5,9,7,
										  0,2,4,1,6,3,8,5,9,7};

	for(iter=0;iter<max_iter;iter++){

		angulos[0]=PD(*(a+tamanio*0+0),*(a+tamanio*0+1),*(a+tamanio*1+1),
						(a_out+tamanio*0+0),(a_out+tamanio*0+2),(a_out+tamanio*2+2));

		angulos[1]=PD(*(a+tamanio*2+2),*(a+tamanio*2+3),*(a+tamanio*3+3),
						(a_out+tamanio*4+4),(a_out+tamanio*1+4),(a_out+tamanio*1+1));

		angulos[2]=PD(*(a+tamanio*4+4),*(a+tamanio*4+5),*(a+tamanio*5+5),
						(a_out+tamanio*6+6),(a_out+tamanio*3+6),(a_out+tamanio*3+3));


		angulos[3]=PD(*(a+tamanio*6+6),*(a+tamanio*6+7),*(a+tamanio*7+7),
						(a_out+tamanio*8+8),(a_out+tamanio*5+8),(a_out+tamanio*5+5));


		angulos[4]=PD(*(a+tamanio*8+8),*(a+tamanio*8+9),*(a+tamanio*9+9),
						(a_out+tamanio*9+9),(a_out+tamanio*7+9),(a_out+tamanio*7+7));

		//PND for filas
		//fila 0
		PND(angulos[0],angulos[1],*(a+tamanio*0+2),*(a+tamanio*0+3),*(a+tamanio*1+2),*(a+tamanio*1+3),
						(a_out+tamanio*0+4),(a_out+tamanio*0+1),(a_out+tamanio*2+4),(a_out+tamanio*1+2));

		PND(angulos[0],angulos[2],*(a+tamanio*0+4),*(a+tamanio*0+5),*(a+tamanio*1+4),*(a+tamanio*1+5),
						(a_out+tamanio*0+6),(a_out+tamanio*0+3),(a_out+tamanio*2+6),(a_out+tamanio*2+3));

		PND(angulos[0],angulos[3],*(a+tamanio*0+6),*(a+tamanio*0+7),*(a+tamanio*1+6),*(a+tamanio*1+7),
						(a_out+tamanio*0+8),(a_out+tamanio*0+5),(a_out+tamanio*2+8),(a_out+tamanio*2+5));

		PND(angulos[0],angulos[4],*(a+tamanio*0+8),*(a+tamanio*0+9),*(a+tamanio*1+8),*(a+tamanio*1+9),
						(a_out+tamanio*0+9),(a_out+tamanio*0+7),(a_out+tamanio*2+9),(a_out+tamanio*2+7));

		//fila 1
		PND(angulos[1],angulos[2],*(a+tamanio*2+4),*(a+tamanio*2+5),*(a+tamanio*3+4),*(a+tamanio*3+5),
						(a_out+tamanio*4+6),(a_out+tamanio*3+4),(a_out+tamanio*1+6),(a_out+tamanio*1+3));

		PND(angulos[1],angulos[3],*(a+tamanio*2+6),*(a+tamanio*2+7),*(a+tamanio*3+6),*(a+tamanio*3+7),
						(a_out+tamanio*4+8),(a_out+tamanio*4+5),(a_out+tamanio*1+8),(a_out+tamanio*1+5));

		PND(angulos[1],angulos[4],*(a+tamanio*2+8),*(a+tamanio*2+9),*(a+tamanio*3+8),*(a+tamanio*3+9),
						(a_out+tamanio*4+9),(a_out+tamanio*4+7),(a_out+tamanio*1+9),(a_out+tamanio*1+7));

		//fila 2
		PND(angulos[2],angulos[3],*(a+tamanio*4+6),*(a+tamanio*4+7),*(a+tamanio*5+6),*(a+tamanio*5+7),
						(a_out+tamanio*6+8),(a_out+tamanio*5+6),(a_out+tamanio*3+8),(a_out+tamanio*3+5));

		PND(angulos[2],angulos[4],*(a+tamanio*4+8),*(a+tamanio*4+9),*(a+tamanio*5+8),*(a+tamanio*5+9),
						(a_out+tamanio*6+9),(a_out+tamanio*6+7),(a_out+tamanio*3+9),(a_out+tamanio*3+7));

		//fila 3
		PND(angulos[3],angulos[4],*(a+tamanio*6+8),*(a+tamanio*6+9),*(a+tamanio*7+8),*(a+tamanio*7+9),
						(a_out+tamanio*8+9),(a_out+tamanio*7+8),(a_out+tamanio*5+9),(a_out+tamanio*5+7));

		//-----------------------------------------------------------------------------

		//PV por columnas
		int ind_angulo;
		ind_angulo=0;
		for(ind_c=0;ind_c<NTERMS;ind_c=ind_c+2){
			for(ind_f=0;ind_f<NTERMS;ind_f=ind_f+2){

				PV(angulos[ind_angulo],autovectores[ind_f][ind_c],autovectores[ind_f][ind_c+1],
											 autovectores[ind_f+1][ind_c],autovectores[ind_f+1][ind_c+1],
													&kk1,&kk2,&kk3,&kk4);

								autovectores_out[posFila[ind_f][ind_c]][posCol[ind_f][ind_c]]=kk1;
								autovectores_out[posFila[ind_f][ind_c+1]][posCol[ind_f][ind_c+1]]=kk2;

								autovectores_out[posFila[ind_f+1][ind_c]][posCol[ind_f+1][ind_c]]=kk3;
								autovectores_out[posFila[ind_f+1][ind_c+1]][posCol[ind_f+1][ind_c+1]]=kk4;
			}
			ind_angulo++;
		}

		//copiar a_out en a
		for(ind_f=0;ind_f<NTERMS;ind_f++){
			for(ind_c=0;ind_c<NTERMS;ind_c++){
				a[ind_f*NTERMS+ind_c] = a_out[ind_f*NTERMS+ind_c];
			}
		}

		//copiar autovectores_aux en autovectores
		for(ind_f=0;ind_f<NTERMS;ind_f++){
			for(ind_c=0;ind_c<NTERMS;ind_c++){
				autovectores[ind_f][ind_c] = autovectores_out[ind_f][ind_c];
			}
		}

	}  //end for max_iter

	//copiar autovectores en v
	for(ind_f=0;ind_f<NTERMS;ind_f++){
		for(ind_c=0;ind_c<NTERMS;ind_c++){
			v[ind_f][ind_c] = autovectores[ind_f][ind_c];
		}
	}

	for(ind_c=0;ind_c<NTERMS;ind_c++){
			w[ind_c]=a[ind_c*NTERMS+ind_c];
	}

	free(a_out);

}


int mil_svd(double *double1_var, double *double2_var, short *short_var,
			double *double3_var)
{

	//	double1_var = h, double2_var = beta, double3_var = delta

int i,j;
double epsilon,top,factor,maximo,minimo;
static double h_svd[NTERMS*NTERMS];
static double v2[NTERMS][NTERMS],w2[NTERMS],v[NTERMS*NTERMS],w[NTERMS];
static double vaux[NTERMS*NTERMS],waux[NTERMS];
static double aux[NTERMS];
static double aux2[NTERMS];
static double beta[NTERMS];
int max_iter;
int aux_nf,aux_nc;

epsilon= 1e-25;
top=1.0;
factor=0.0;
maximo=0.0;
minimo=1000000000.0;
max_iter = 18;

//	BUSCO EL MAXIMO DE EL VECTOR DE ENTRADA double1_var = h
for(j=0;j<NTERMS*NTERMS;j++){
		if(fabs(double1_var[j])>maximo){
			maximo=fabs(double1_var[j]);
		}
	}
	factor=maximo;

	// replico double1_var (h) en h_svd
	for(j=0;j<NTERMS*NTERMS;j++){
		h_svd[j]=double1_var[j];
		}
		// replico double1_var (h) en h_svd
		for(j=0;j<NTERMS;j++){
			beta[j]=double2_var[j];
			}
//NORMALIZACION DE LA MATRIX

for(j=0;j<NTERMS*NTERMS;j++){
		h_svd[j]=h_svd[j]/factor;
}

//SVD

	svdcordic(h_svd,w2,v2,max_iter);

	for(i=0;i<NTERMS-1;i++){
		for(j=0;j<NTERMS-1;j++){
			v[i*NTERMS+j]=v2[i][j];
		}
	}

	for(j=0;j<NTERMS-1;j++){
		w[j]=w2[j]*factor;
	}

//REDUNDANTE
	for(j=0;j<NTERMS*NTERMS;j++){
			vaux[j]=v[j];
	}

	for(j=0;j<NTERMS;j++){
			waux[j]=w[j];
	}


multmatrix(beta,1,NTERMS,vaux,NTERMS,NTERMS,aux2,&aux_nf,&aux_nc);

for(i=0;i<NTERMS;i++){
	aux2[i]= aux2[i]*((fabs(waux[i]) > epsilon) ? (1/waux[i]) : 0.0);
}

multmatrix(vaux,NTERMS,NTERMS,aux2,NTERMS,1,aux,&aux_nf,&aux_nc);

	for(j=0;j<NTERMS*NTERMS;j++){
		double1_var[j] = h_svd[j];
		}

		*short_var = (short *) max_iter;

		for(j=0;j<NTERMS;j++){
				double2_var[j] = beta[j]*2.0;//*factor;
		}

		for(j=0;j<NTERMS;j++){
				double3_var[j] = aux[j];//waux[j];
		}

  return 1;

}

int mil_svd_c(int argc, void* argv[])

{
  if(argc != 4) return 0;

  return mil_svd((double *) argv[0], (double *) argv[1],
			     (short *) argv[2], (double *) argv[3]);
}
