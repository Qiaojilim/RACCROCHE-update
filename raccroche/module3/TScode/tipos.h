void escribe_result(char **Nombre);


void actualiza_frec(int *frec,int dim,int t);
void calcula_frec(int *frec,int dim);
long path_relinking(int *orden,int *orden_inv,int dim,int elite_num,
				    long **dif,long coste,int **e_pos);
void reemplaza_elite(long *mvg,int *mig,int dim,int *orden,int **elite,
				     int **e_pos,long *c_elite,long coste,int elite_num);
void inicia_elites(int dim,int *orden,int *elite_num,int **e_pos,
				   int **elite,long *c_elite,int *mig,long *mvg,
				   long coste,int numeroelites);
long optimo_local(int *orden_l,int *orden_invl,int dim,long coste,long **dif);
void actual_a_bestlocal(int *orden,int *orden_inv,int *orden_l,
						int *orden_invl,int dim);
void bestlocal_a_actual(int *orden,int *orden_inv,int *orden_l,
						int *orden_invl,int dim);
long **lee_reinelt(char *fichero,int *dim);
void lee_sgb(long **matriz,int dim,int semilla);
void relllena_aleat(long **matriz,int dim,int r,int semilla);
void inicia(long **matriz,long **dif,long *score,int *score_int,int dim);
int **matriz_enteros(int filas,int columnas);						
long **lee_new_random(char *fichero,int *dim);						
						



long chanas_koby(long **matriz,int *orden,int dim, long coste, long *iter);
long inserta_parcial(long **matriz, int ultimo, int pos, int *orden, int parcial);

long mejor_insert(int j,int *orden,long **dif,int dim,int *pos);
long greedy(int *orden,long **dif,int dim,long coste);  
long greedy_lento(int *orden,long **dif,int dim,long coste);
long insert(int j,int i,int *orden,long **dif);  
void inserta(int j,int pos,int *orden);
void inserta2(int j,int pos,int *orden,int *orden_inv);
long greedy_switch(int *orden,long **dif,int dim,long coste);
long greedy_switch_lento(int *orden,long **dif,int dim,long coste);
 
int select_indice(int *score);     
int vert_diferente(int *orden,int guia,int **elite,int dim);
int compara_con_elite(int dim,int *orden,int **e_pos,int a); 
void guarda_elite(int dim,int *orden,int **elite,int **e_pos,int pos);     
long **reserva_matriz(int dim);
long *reserva_vector_long(int dim); 
float *reserva_vector_float(int dim);
int *reserva_vector_int(int dim); 
void abortar(char *texto);
long calcula_coste(long **matriz,int *orden,long dim); 
void reverse(int *orden,int dim);
void orden_aleatorio(int *orden,int dim);                                                  
                                                  
long improving(int *orden,long **dif,int dim,int j,int *orden_inv);
long sure_insert(int *orden,long **dif,int dim);
int indice(int dim); 
long sure_switch(int *orden,long **dif,int dim);

void becker(long **matriz,int *orden,int dim);

#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>


#define getrandom( min, max ) ((rand() % (int)(((max)+1) - (min))) + (min))
