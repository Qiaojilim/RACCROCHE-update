/****  	TABU SEARCH PROCEDURE FOR THE LINEAR ORDERING PROBLEM
 * 		NEEDS 2 ARGUMENTS:  FILE_NAME_(NORMAL_FORM) NUMBER_GLOBAL_ITERATIONS  (>= 10)
 * 		RETURNS 2 NUMBERS:  VALUE OF BEST SOLUTION FOUND, TIME ELAPSED IN SECONDS
 */

#include "tipos.h"
#define NUM_ELITE 4    //memorize 4 elite solutions

int main(int argc,char **argv)
{
	long **matriz,**dif,bcoste;
	int i,j,t,*orden,*orden_invl,*orden_inv,*orden_l, *orden_final;
	long coste, mejor_coste,c_init,mejor_local,*score,*tabu;
	int *score_int,*frec,pos,mig,pr; /* Index of the guide with worse value*/
	long change_total,iter,change,tenure,iter_inten;
	int **elite,**e_pos,elite_num=0,find;
	long mvg,*c_elite,global;
	float *medias,posi_media;
	time_t t_tabu,f_tabu, total_t;
	int dim,relink,fin_inten,fin_diver;

	matriz=lee_new_random(argv[1],&dim);
	global=atol(argv[2]);

	/************   Initializations   *************************/
	
	orden 			= reserva_vector_int(dim);
    orden_inv 		= reserva_vector_int(dim);
    orden_l		 	= reserva_vector_int(dim);
    orden_invl		= reserva_vector_int(dim);
    score 			= reserva_vector_long(dim);
    score_int 		= reserva_vector_int(dim);
	frec      		= reserva_vector_int(dim);
    tabu 			= reserva_vector_long(dim);
    dif 			= reserva_matriz(dim);			// this is a 2d matrix (dim * dim), index from 1 to dim
	c_elite 		= reserva_vector_long(dim);
    elite 			= matriz_enteros(NUM_ELITE,dim); //elite solutions: 2D matrix (4 rows, dim columns), index of rows 1,2,3,4
    e_pos 			= matriz_enteros(NUM_ELITE,dim); // 2d matrix
    medias			= reserva_vector_float(dim) ;
	

    inicia(matriz,dif,score,score_int,dim);    // given matriz, dif is the skew symmetric matrix of differences

	for(i=1 ; i<=dim ; i++)
		orden_l[i]=orden_invl[i]=orden_inv[i]=orden[i]=i;
		
	coste = mejor_coste = c_init = calcula_coste(matriz,orden,dim);
	
	//printf("\nThe sum of superdianonal elements of the input matrix: %ld\nDim of the matrix is %d\n", coste, dim);  //289493
	//printf("Initial order is \n");	
 	//for(i=1 ; i<=dim ; i++)
 	//	printf("%d ", orden[i]);   // order from 1 to 1015 indices


	tenure		= 2*(long)sqrt(dim);     // the moved object becomes tabu-active for TabuTenure iterations, and therefore it cannot be selected for insertions 
	fin_inten	= 10*dim; 	// dim       end of intensification phase
	fin_diver	= 10*dim; 	// dim/2     end of diversification phase
	pr      	= 4;		// Full Algorithm

	relink = 0;
    srand(1471);

	/************   Algoritmo  Tabu ********************************/

	t_tabu = clock();
	change_total=iter=0;

	

	while(change_total<global)  //global iterations
	{
		iter++;
		for(j=1;j<=dim;j++) {
			frec[j]=0;
			tabu[j]=-1000;
			medias[j]=orden_inv[j];
		}
		change=iter_inten=0;
		mejor_local=coste;
		actual_a_bestlocal(orden,orden_inv,orden_l,orden_invl,dim); //orden_l[s]=orden[s]; orden_invl[s]=orden_inv[s];
		
		///*************
		//printf("\n Intensification \nchange_total=%ld\n",change_total);

		while(change < fin_inten)
		{
		    /* Intensification based on better insertion */

			j = select_indice(score_int);  /* indice original */
			if(iter_inten - tabu[j] > tenure )
			{
				i = orden_inv[j];tabu[j] = ++iter_inten;frec[j]++;
				coste += mejor_insert(i,orden,dif,dim,&pos);
				inserta2(i,pos,orden,orden_inv);
				if(coste>mejor_local)
				{
					actual_a_bestlocal(orden,orden_inv,orden_l,orden_invl,dim);
					mejor_local=coste;
					///*************
					orden_final = orden;
					change=0;

				}
				else	change++;

				medias[j] += pos;

		    }
		    if(change == fin_inten ) {
		    	mejor_local=optimo_local(orden_l,orden_invl,dim,mejor_local,dif);
		        bestlocal_a_actual(orden,orden_inv,orden_l,orden_invl,dim);
		        coste = mejor_local;
		        ///*************
				orden_final = orden;

		    }
		}

		if(mejor_local > mejor_coste) {    //best local score greater than best score
			 	mejor_coste=mejor_local;
			 	change_total=0;
			 	///*************
				//printf("\n*** %ldvs%ld", mejor_local, mejor_coste);
				//printf("\n*** orden_final \n");
				for(i=1 ; i<=dim ; i++)
 					printf("%d\n", orden_final[i]);   // order from 1 to 1015 indices


		}
		else 	change_total++;


	/************* Path Relinking ***********************************/

		if(pr>1 && elite_num < NUM_ELITE)
		  inicia_elites(dim,orden,&elite_num,e_pos,elite,c_elite,&mig,&mvg,coste,NUM_ELITE);

		if( (pr==2 || pr==4) && elite_num >= NUM_ELITE)
		{
		        find=0;
				for(i=1;i<=elite_num;i++)
					if(compara_con_elite(dim,orden,e_pos,i)==0) find=1;

				if(find==0)
				{
					bcoste=path_relinking(orden_l,orden_invl,dim,elite_num,dif,coste,e_pos);
				  	if(bcoste>mejor_coste)
				  	{
		  				mejor_coste = bcoste;
		  				relink++;
		  			}
		  		}
				if(mejor_local > mvg )
		   			reemplaza_elite(&mvg,&mig,dim,orden,elite,e_pos,c_elite,coste,elite_num);
		}


	/*************  Long Term Diversification *************************/


		if(pr>2 && change_total > global/2)
		{
			///*************
			//printf("\n Diversification - IF branch \n");
		   for(i=1;i<=dim;i++)
		  		medias[i] /= (1+frec[i]);

		  for(i=1;i<=dim;i++)
		  {
		  		posi_media=0;
		  		for(j=1;j<=elite_num;j++)
		  			posi_media += (float) e_pos[j][i];
		  		posi_media /= elite_num;
		  		medias[i] = (medias[i]+posi_media)/2;
   		  }

   		  coste = mejor_local;
   		  for(j=1;j<=dim;j++)
   		  {
   		  		i = orden_inv[j];
				pos = dim - (int) medias[j];
				if(pos<1)   pos=1;
				if(pos>dim) pos=dim;
				coste += insert(i,pos,orden,dif); /* Evaluate the insertion of the element of the position i of order, in the position pos.
   													Returns the cost difference */
				inserta2(i,pos,orden,orden_inv); //Performs the insertion of i at pos
				if(coste > mejor_coste)
				{
					mejor_coste=coste;
					///*************
					//printf("%ld  ", mejor_coste);
					
				}
   		  }
   		}
   		else
   		{
			/********** Diversification of the usual Tabu ****************/
			///*************
			//printf("\n Diversification - ELSE branch \n");

			calcula_frec(frec,dim);

			for(j=0 ; j< fin_diver ;j++)
			{
				t = select_indice(frec);
				i = orden_inv[t];
				coste += mejor_insert(i,orden,dif,dim,&pos);
				inserta2(i,pos,orden,orden_inv);
				actualiza_frec(frec,dim,t);
				if(coste > mejor_coste)
				{
					mejor_coste=coste;
					///*************
					//printf("%ld  ", mejor_coste);
				}
     		}
        }
    } //end while

    f_tabu = clock();

	total_t = (double)(f_tabu - t_tabu) / CLOCKS_PER_SEC;
//printf("\n%s %ld %4.2f",argv[1],mejor_coste,((float)f_tabu - t_tabu)/CLK_TCK );
//printf("\nMatrix is %s; value of best solution %ld; the number of seconds used by the CPU %4.2ld  \n\n",argv[1],mejor_coste,total_t);

    return 0;
}





