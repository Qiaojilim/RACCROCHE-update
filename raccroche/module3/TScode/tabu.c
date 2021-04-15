#include "tipos.h"

#define BLOCK 5

long improving(int *orden,long **dif,int dim,int j,int *orden_inv)
{
	long valor;
	int pos;
	        
	valor = mejor_insert(j,orden,dif,dim,&pos);
	if (valor>0)
	{
		inserta2(j,pos,orden,orden_inv);
		return valor;		                  
	}
	
	return 0;
}
   
   
long sure_insert(int *orden,long **dif,int dim)
{
	long valor;
	int j,pos;
	
	j=getrandom(1,dim);        
	valor = mejor_insert(j,orden,dif,dim,&pos);
	inserta(j,pos,orden);
	
	return valor;		                  
}
	

int indice(int dim)
{
	static int i=1;
	int a,t,j;

	a=getrandom( 1,BLOCK);

	for(t=1 ; t<=BLOCK ;t++)
	{
		if(a==t)
			j=i;
		if(i==dim)
			i=1;
		else
			i++;
	}
	return j;
}  


long sure_switch(int *orden,long **dif,int dim)
{ 
	int j;               
	long valor,valor2;
	
	j = getrandom(1,dim);

	if(j==1)
	{	
		valor = insert(1,2,orden,dif);
	   	inserta(1,2,orden);
	}
    else if( j>1 && j<dim)
    {
		valor  = insert(j,j+1,orden,dif);
		valor2 = insert(j,j-1,orden,dif);
		if (valor > valor2 )  
				inserta(j,j+1,orden);
		else  {
				inserta(j,j-1,orden);
				valor = valor2;    
		}   
	}
	else
	{
		valor = insert(dim,dim-1,orden,dif);
	    inserta(dim,dim-1,orden);
	}
	return valor;
}