#include "tipos.h"



int select_indice(int *score)
{
    int cont=1,f=0,prob_sel;
    
	prob_sel=getrandom(0,score[0]);
	while(f==0)
	{
		if(prob_sel<=score[cont])	f=1;
		else	prob_sel -= score[cont++];
	}
	return cont;
}


int vert_diferente(int *orden,int guia,int **elite,int dim)
{               
/* Devuelve el primer vertice con posiciones diferentes.
   devuelve un -1 si todos son iguales */
	
	int j,cont=0;

	j=getrandom(1,dim);	
	while(cont<=dim && orden[j]==elite[guia][j])
	{
		cont++;
		j++;
		if(j==dim+1) j=1;
	}
	if( orden[j]==elite[guia][j])
		return -1;
	else
		return j;
}


int compara_con_elite(int dim,int *orden,int **e_pos,int a)
{                         
/* Devuelve la suma de las diferencias de posiciones entre las solucion elite a y
   la solucion de orden. Si devuelve un 0 esque son la misma. */
   	
	int i,dif=0;
	
	for(i=1;i<=dim;i++)
		dif += abs(e_pos[a][orden[i]]-i);
	return dif;
}                                                                     
                                                                     



void guarda_elite(int dim,int *orden,int **elite,int **e_pos,int pos)
{
	int i;
	
	for(i=1;i<=dim;i++)
	{
		elite[pos][i]=orden[i];
		e_pos[pos][orden[i]]=i;
	}
}	

void reverse(int *orden,int dim)
{
	int i,*aux;
	
	aux = reserva_vector_int(dim);
	
	for(i=1;i<=dim;i++) aux[i]=orden[i];
	
	for(i=1;i<=dim;i++) orden[i]=aux[dim+1-i];
	
	free(aux);
}   

void orden_aleatorio(int *orden,int dim)
{
	int j,a,i,pos;
	
	for(i=1;i<=dim;i++) orden[i]=0;
	
	for(j=1;j<=dim;j++)
	{
		/* posicion de j (de las libres) */
		a=getrandom(1,dim+1-j); 
		pos=1;
		for(i=1 ; i<a ; i++)
		{
			while(orden[pos]!=0)
				pos++;
		     pos++;
		}
		while(orden[pos]!=0)
				pos++;
		orden[pos]=j;
	}
		
}
	
	

long calcula_coste(long **matriz,int *orden,long dim)
{
	long coste=0;                                      
	long i,j;
	
	for(i=1   ; i< dim ; i++)
	for(j=i+1 ; j<=dim ; j++)
		coste += (long) matriz[orden[i]][orden[j]];
	
	return coste;
}
	
	
long **reserva_matriz(int dim)
{
	long i,**aux;

    aux=(long**)calloc(dim+1,sizeof(long*));
    if(!aux) abortar("Reserva matriz");

	for(i=1;i<=dim;i++)
		aux[i] = reserva_vector_long(dim);  
	
	return aux;
}  
    
long *reserva_vector_long(int dim)
{
	long *aux;

	aux=(long*)calloc(dim+2,sizeof(long));
	if(!aux) abortar("Reserva vector long");
    
	return aux;
} 

int *reserva_vector_int(int dim)
{
	int *aux;

	aux=(int*)calloc(dim+2,sizeof(int));
	if(!aux) abortar("Reserva vector int");
    
	return aux;
} 

float *reserva_vector_float(int dim)
{
	float *aux;

	aux=(float*)calloc(dim+2,sizeof(float));
	if(!aux) abortar("Reserva vector float");
    
	return aux;
} 


void abortar(char *texto)
{
	printf("\n%s",texto);
	exit(6);
}
