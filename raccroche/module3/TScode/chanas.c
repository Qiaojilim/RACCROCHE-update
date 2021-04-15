#include "tipos.h"

long chanas_koby(long **matriz,int *orden,int dim, long coste, long *iter)
{
	long acum,valor,best_valor;
	int parcial,ultimo,pos,best_pos,iterar,primer_empate;  
	
	*iter = -1;
		
	/* (SORT(REVERSE))SORT(P) */ 
	
	do{	    
	   acum = 0;
	   for(parcial=1;parcial<dim;parcial++) 
	    {
	     best_valor = -10000;
	     ultimo=parcial+1; 
	     
	     for(pos=1;pos<=ultimo;pos++)
	      {
	       valor = inserta_parcial(matriz,ultimo,pos,orden,parcial);  
	      
	       if(best_valor < valor)
	        {
	         best_pos = pos;
	         best_valor = valor;
	        }     
	      } 
	      inserta(ultimo,best_pos,orden);   
	      acum += best_valor;
	      
	    } 
	    
	   iterar = 0;
	   if( acum != coste)
	    {
	     coste = acum; 
	     primer_empate = 1;
	     iterar = 1; 
	    } 
	   else if(primer_empate)
	    {
	     primer_empate = 0;
	     reverse(orden,dim);
	     iterar = 1;
	    } 
	   (*iter)++;   
	    
	   }while(iterar);  
	   
	   
    
    return acum;
	
}

long inserta_parcial(long **matriz, int ultimo, int pos, int *orden, int parcial)
{ 
 /* Evalua el coste de insertar orden[ultimo] en la posicion pos de orden */
 /* No cambia orden */
 
 long valor = 0;
 int i;
 
 for(i=1;i<pos;i++)  
   valor += matriz[orden[i]][orden[ultimo]];
 for(i=pos;i<=parcial;i++)  
   valor += matriz[orden[ultimo]][orden[i]];     
 
 return valor;
}	