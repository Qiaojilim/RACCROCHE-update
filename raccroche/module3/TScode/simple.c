#include "tipos.h"

void becker(long **matriz,int *orden,int dim)
{
	long fila,columna;
	double *q,max;
	int i,s,*aux,best_i,j;
	
	q=(double*)calloc(dim+2,sizeof(double));
	if(!q) abortar("Problemas en Becker"); 
	
	aux = reserva_vector_int(dim+1);
	    
	for(s=1 ; s<=dim ; s++)
	{   
		for(i=1;i<=dim;i++)	q[i]=0;
		i=1;
		for(i=1;i<=dim;i++)
		{
			if(aux[i]!=0)	continue;
			{   
				fila=columna=0;  
				for(j=1;j<=dim;j++)
				{   
					if(aux[j]!=0)	continue;           
					fila += matriz[i][j];
					columna += matriz[j][i];
				} 
				if(columna==0) q[i]=0;
				else
					q[i] = (double)fila / (double)columna;
			}
        }
		
		/* Encontrar el maximo */
		
		max=-10000;    
    	for(i=1 ; i<=dim; i++)
    	{       
   			if(aux[i]==0 && q[i]>max)
    		{
    			max=q[i];
    			best_i=i;
    		}
    	}
    	orden[s]=best_i;
    	aux[best_i]=1;
    }
	free(aux);
}	