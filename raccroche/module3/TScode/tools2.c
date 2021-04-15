#include "tipos.h"

long **lee_new_random(char *fichero,int *dim)
{
	FILE *puntfile;
	int i,j;
	long **matriz;

	if((puntfile = fopen(fichero,"r"))==NULL)
	    	abortar("\nData file not found");
	fscanf(puntfile,"%d\n",dim);
	matriz = reserva_matriz(*dim);
	for(i=1 ; i<=*dim ; i++)
	{
		for(j=1 ; j<=*dim ; j++)
	    	fscanf(puntfile,"%ld",&matriz[i][j]);
	  	fscanf(puntfile,"\n");
	}
	fclose(puntfile);

	return matriz;
}

void actualiza_frec(int *frec,int dim,int t)
{
	int s;

    if(frec[t]==0)
    {
       	for(s=1;s<=dim;s++)
       	{
       		frec[s] +=100;
       		frec[0] +=100;
       	}
    }
    frec[t]--;
    frec[0]--;
}



void calcula_frec(int *frec,int dim)
{
	int max,i;

	max=0;
	for(i=1 ; i<=dim ;i++)
			if(frec[i]>max)		max=frec[i];

    frec[0]=0;
	for(i=1 ; i<=dim ; i++)
	{
			frec[i] = max - frec[i];
			frec[0] += frec[i];
	}
}


long path_relinking(int *orden,int *orden_inv,int dim,int elite_num,
				    long **dif,long coste,int **e_pos)
{

	int j,s,i,guia,min,dife,*duplica;
    int real_moves=0;
    long mejor_coste,d_coste;

    mejor_coste = coste;
    duplica 	= reserva_vector_int(dim);

	for(i=1;i<=dim;i++) /* We go through all the indices */
	{
		j = orden[i];

		/* Select as the new position of j, the closest between the guides */

		min=10000;
		for(s=1;s<=elite_num;s++)
		{
			dife=abs(i-e_pos[s][j]);
			if(dife<min) {
		  		min=dife;
		  		guia=s;
			}
	    }
	    if(min>0)
	    {
	    	real_moves++;
			coste += insert(i,e_pos[guia][j],orden,dif);
		    inserta2(i,e_pos[guia][j],orden,orden_inv);
		    if(coste>mejor_coste) mejor_coste = coste;

			if(real_moves % 4 ==0)
			{
		    	/* Apply Greedy to an intermediate solution*/
				 d_coste=coste;
		         for(i=1;i<=dim;i++) duplica[i]=orden[i];
		         d_coste = greedy(duplica,dif,dim,d_coste);
		         if(d_coste>mejor_coste) mejor_coste=d_coste;
            }
        }
    }
    if(real_moves >0)
	{
		 /* Apply Greedy to the latest solution */
			d_coste=coste;
		    for(i=1;i<=dim;i++) duplica[i]=orden[i];
		    d_coste = greedy(duplica,dif,dim,d_coste);
		    if(d_coste>mejor_coste) mejor_coste=d_coste;
    }
    free(duplica);
    return mejor_coste;
}




void reemplaza_elite(long *mvg,int *mig,int dim,int *orden,int **elite,
				     int **e_pos,long *c_elite,long coste,int elite_num)
{
 /* Replace a guide by the optimal initial location if it is different from the ones stored  */

	long minimo;
	int i;

	guarda_elite(dim,orden,elite,e_pos,*mig);
	c_elite[*mig]=coste;

    minimo=1000000000;
	for(i=1;i<=elite_num;i++)
		if(c_elite[i]<minimo)
		{
			minimo=c_elite[i];
			*mig=i;
			*mvg=c_elite[i];
		}
}



void inicia_elites(int dim,int *orden,int *elite_num,int **e_pos,
				   int **elite,long *c_elite,int *mig,long *mvg,
				   long coste,int numeroelites)
{
		int i,find;
		long minimo;

		if(*elite_num==0)
		{
			guarda_elite(dim,orden,elite,e_pos,++(*elite_num));
			c_elite[*elite_num]=coste;
		}
		else
		{   /* If it is different from the stored ones we store it*/

			find=0;
			for(i=1;i<=*elite_num;i++)
				if(compara_con_elite(dim,orden,e_pos,i)==0)
						find=1;
			if(find==0)
			{
				guarda_elite(dim,orden,elite,e_pos,++(*elite_num));
				c_elite[*elite_num]=coste;
			}

			/* When we put the last guide, we calculate the lowest value*/

			if(*elite_num==numeroelites)
			{
				minimo=1000000000;
				for(i=1;i<=*elite_num;i++)
					if(c_elite[i]<minimo)
					{
						minimo=c_elite[i];
						*mig=i;
						*mvg=c_elite[i];
					}
			}
	}
}



long optimo_local(int *orden_l,int *orden_invl,int dim,long coste,long **dif)
{
	/* Calculate the best local from the best_local (order_l)   mejor_local (orden_l)*/

	int i;
	long old_coste;

   	do {
    		  old_coste = coste;
		      for(i=1;i<=dim;i++)
		        	  	coste += improving(orden_l,dif,dim,i,orden_invl);
	} while(coste>old_coste);

    return coste;
}




void actual_a_bestlocal(int *orden,int *orden_inv,int *orden_l,
						int *orden_invl,int dim)
{
	int s;

	for(s=1;s<=dim;s++)
	{
		orden_l[s]=orden[s];
		orden_invl[s]=orden_inv[s];
	}
}


void bestlocal_a_actual(int *orden,int *orden_inv,int *orden_l,
						int *orden_invl,int dim)
{
	int s;

	for(s=1;s<=dim;s++)
	{
		orden[s]=orden_l[s];
		orden_inv[s]=orden_invl[s];
	}
}




long **lee_reinelt(char *fichero,int *dim)
{
	FILE *puntfile;
	int modulo,i,j;
	long **matriz;


		if((puntfile = fopen(fichero,"r"))==NULL)
	    	abortar("\nData file not found");
	    while(fgetc(puntfile)!='\n')
	       ; /* lee comentario */

			fscanf(puntfile,"%d\n",dim);

		matriz = reserva_matriz(*dim);

		/* Dos formatos distintos de datos:
		   Columnas de 10 en 10 para matrices 50x50 y 60x60
		   Columnas de  6 en  6	para matrices 44x44 y 56x56 */

		if(*dim % 10 == 0)
		  modulo =10;
		else
		  modulo = 6;
	    for(i=1 ; i<=*dim ; i++)
		for(j=1 ; j<=*dim ; j++)
	  	{
		 if(i*j % modulo == 1)
		  		fscanf(puntfile,"\n");
	    	fscanf(puntfile,"%ld",&matriz[i][j]);
	  	}
	    fclose(puntfile);

	    return matriz;
}


void relllena_aleat(long **matriz,int dim,int r,int semilla)
{
/** Matrix array with random values **/

		int i,j;

		srand( semilla );
		for(i=1 ; i<=dim ; i++)
		for(j=1 ; j<=dim ; j++)
			if(i!=j)
				matriz[i][j] = getrandom(0,r);
}


void inicia(long **matriz,long **dif,long *score,int *score_int,int dim)
{
	int i,j;
	long menor,suma;

	for(i=1 ; i<=dim ; i++)
	for(j=1 ; j<=dim ; j++)
	{
		dif[i][j] = matriz[i][j] - matriz[j][i];
		if(i!=j)
		{
			score[i] += matriz[i][j];
			score[j] += matriz[i][j];
		}
	}


	/*** Inicializacion  Scores **********/

	menor=2000000000L;
	suma=0;
	for(i=1 ; i<=dim ;i++)
	{
		if(score[i]==0) score[i]=1;
		if(score[i]<menor && score[i]!=1)
			menor = score[i];
		suma += score[i];
	}

    if(suma<32500 || menor==0) menor=1;

	if(menor<=0) abortar("Hay un score negativo");

	for(i=1;i<=dim;i++)
	{
		 score[i] /= menor;
	     score[0] += score[i];
	}
	while (score[0]>32500)
	{
	  	menor=30000;
	  	score[0]=0;
	  	for(i=1;i<=dim;i++)
	  		if(score[i]<menor && score[i]>1)
	  			menor=score[i];
	  	for(i=1;i<=dim;i++)
	  	{
	  		score[i] /= menor;
	  		if(score[i]<1) score[i]=1;
	  		score[0] += score[i];
	  	}
	 }


	for(i=0 ; i<=dim ; i++)
		score_int[i] = (int) score[i];
}

int **matriz_enteros(int filas,int columnas)
{
	int **elite,i;

	elite=(int**)calloc(filas+1,sizeof(int*));
    if(!elite) abortar("Reserva elite");

	for(i=1;i<=filas;i++)
		elite[i] = reserva_vector_int(columnas);

	return elite;
}
