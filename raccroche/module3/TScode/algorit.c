#include "tipos.h"


long greedy_lento(int *orden,long **dif,int dim,long coste)
{ 
/* Estudia todos los elementos para ver aquel que produciria un aumento
   mayor, y lo realiza. Despues vuelve a estudiar todos los elementos...
*/   
   
	long c,old_c;
	long best_valor,valor;
	int best_j,best_pos,j,pos;        

	c=coste;
	
	do { 
		
	    old_c = c;
	    best_valor=0;
		for(j=1 ; j<=dim ; j++)
		{
			valor = mejor_insert(j,orden,dif,dim,&pos); 
			if(valor > best_valor)	
			{
				best_valor = valor;
				best_pos   = pos;
				best_j     = j;
			}	
			
		}  
		if(best_valor>0)
		{
			c += (long) best_valor;
			inserta(best_j,best_pos,orden);
		}		
		
		
	} while(c > old_c);		
	
	return c;		                  
}


	                  
long greedy(int *orden,long **dif,int dim,long coste)
{
	long c,old_c,valor;
	int j,pos;
	        
	c=coste;
	
	do {
	    old_c = c;
		for(j=1 ; j<=dim ; j++)
		{
			valor = mejor_insert(j,orden,dif,dim,&pos);
			if (valor>0)
			{
				inserta(j,pos,orden);
				c += (long) valor;
			}
		}
	} while(c > old_c);		
	
	return c;		                  
}	                  
	                  
long mejor_insert(int j,int *orden,long **dif,int dim,int *pos)
{
/* Dado el elemento de la posicion j de orden, se busca la mejor
   posicion en la que insertarlo. Se devuelve el update del coste 
   y la mejor posicion en pos. Puede que la mejor posicion sea 
   empeorar ya que no considera el dejarlo donde esta.
*/

	long best_c,c=0;
	int k,best_pos;
		
	if(j>2)  /* Inicializacion */
	{
		best_c   = dif[orden[j]][orden[j-1]]; 
		best_pos = j-1;         
	}
	else
	{
	 	best_c   = dif[orden[j+1]][orden[j]]; 
		best_pos = j+1;   
	}

	/* Insertarlo en una posicion anterior */
	for(k=j-1;k>=1;k--)
	{
		c += dif[orden[j]][orden[k]];                  
		if(c > best_c)
		{
			best_c = c;
	        best_pos = k;
	    }
	}
	
	/* Insertarlo en una posterior */
	c=0; 	                  
	for(k=j+1;k<=dim;k++)
	{
		c += dif[orden[k]][orden[j]];   	
	    if(c > best_c)
	    {
	    	best_c = c;
	        best_pos = k;
	    }
	}  
	
	*pos = best_pos;
	return best_c;
}


void inserta2(int j,int pos,int *orden,int *orden_inv)
{
/** Performs the insertion of j after order */
// insert j at position pos

	int k,a;
		
	if ( pos > j)
	{   
	    a=orden[j];       
		for(k=j; k<pos ; k++) 
		{
			orden[k]=orden[k+1];
			orden_inv[orden[k+1]]=k;
		}
		orden[pos]=a;
		orden_inv[a]=pos;
	}		
	else
	{	
		a=orden[j];  
	    for(k=j;k>pos;k--) 
	    {
	    	orden[k]=orden[k-1];            
	    	orden_inv[orden[k-1]]=k;
	    }
	    orden[pos]=a;
	    orden_inv[a]=pos;
	}
}

void inserta(int j,int pos,int *orden)
{
/** Realiza la insercion de j en pos de orden */

	int k,a;
		
	if ( pos > j)
	{   
	    a=orden[j];       
		for(k=j; k<pos ; k++) 
			orden[k]=orden[k+1];
		orden[pos]=a;
	}		
	else
	{	
		a=orden[j];  
	    for(k=j;k>pos;k--) 
	    	orden[k]=orden[k-1];            
	    orden[pos]=a;
	}
}		
		

long greedy_switch(int *orden,long **dif,int dim,long coste)
{
	long c,old_c,valor,valor2;
	int j;
	        
	c=coste;
	
	do {
	    old_c = c;
	    
	    /** Para el primero **/
	    valor = insert(1,2,orden,dif);
	    if (valor>0) {  
	    	inserta(1,2,orden);
	    	c += (long) valor;    }

		for(j=2 ; j<dim ; j++)
		{    
			valor  = insert(j,j+1,orden,dif);
			valor2 = insert(j,j-1,orden,dif);
			if (valor>valor2 && valor>0) {
				inserta(j,j+1,orden);
				c += (long) valor;    }
			else if(valor2>valor && valor2>0) {
				inserta(j,j-1,orden);
				c += (long) valor2;    }
		}  
		
		/** Para el ultimo **/
		valor = insert(dim,dim-1,orden,dif);
	    if (valor>0) {  inserta(dim,dim-1,orden);
	    				c += (long) valor;    }
	} while(c > old_c);		
	
	return c;
}		     
		                  

long greedy_switch_lento(int *orden,long **dif,int dim,long coste)
{
	long c,old_c,valor,best_valor,valor2;
	int j,best_j,best_pos;
	        
	c=coste;
	
	do {
	    old_c = c;
	    best_valor=0;
	    
	    /** Para el primero **/
	    valor = insert(1,2,orden,dif);
	    if (valor>0) {
				best_valor = valor;
				best_pos   = 2;
				best_j     = 1;
		}		

		for(j=2 ; j<dim ; j++)
		{    
			valor  = insert(j,j+1,orden,dif);
			valor2 = insert(j,j-1,orden,dif);
			if (valor>valor2 && valor>best_valor) {
					best_valor = valor;
					best_pos   = j+1;
					best_j     = j;  }
			else if(valor2>valor && valor2>best_valor)
			{
				best_valor = valor2;
				best_pos   = j-1;
				best_j     = j;
			}		
		}  
		
		/** Para el ultimo **/
		valor = insert(dim,dim-1,orden,dif);
	    if (valor>0) {
				best_valor = valor;
				best_pos   = dim-1;
				best_j     = dim;
		}		
		
		if(best_valor>0) {  
			inserta(best_j,best_pos,orden);
	    	c += (long) best_valor;
	    }     
	    
	} while(c > old_c);		
	
	return c;
}		     






long insert(int j,int i,int *orden,long **dif)
{
/* Evaluate the insertion of the element of the position j of order, in the position i.
   Returns the cost difference */

	long c=0;	                  
	int k;
	                  
	if(i==j) return 0;	                                 
	                               
	if(i<j)
		for(k=i;k<j;k++)
			c += dif[orden[j]][orden[k]];                  
	else
		for(k=j+1;k<=i;k++)
			c += dif[orden[k]][orden[j]];   	
	return c;
}		        
