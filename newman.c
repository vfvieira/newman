#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "petscmat.h"
#include "petscvec.h"
#include "petscis.h"
#include "network.h"
#include "newman.h"

PetscScalar *unity_array;
int nbranches;

void GetEigBSparse(Mat *net, Vec *k_sub, struct NetData ndata, PetscInt *nodes_comm1, PetscInt *nodes_comm2, IS ids, int *id_comm1, int *id_comm2, Vec *eig_result){
	
	//struct timeval t_start, t_int, t_end;
	//gettimeofday(&t_start, NULL);
	
	//double t_diff;
	int i;
	
	// Variáveis "globais"
	PetscInt edges = ndata.edges;
	PetscScalar weight = ndata.weight;
	
	PetscInt display = ndata.display;
		
	PetscInt nodes_subnet;
	ISGetSize(ids,&nodes_subnet);
	
	const PetscInt *ids_subnet;
	PetscMalloc(nodes_subnet*sizeof(PetscInt),&ids_subnet);
	ISGetIndices(ids, &ids_subnet); // nao consegui destruir... (destruida com ids?)
	
	Vec eig_init, diff, diff2; // Vetores para iteração no Power Method
	PetscScalar *x0; // Vetores usados como referência por eig_init e eig_result
	PetscInt pos;
	
	int iter, max_iter;
	max_iter = 2*nodes_subnet;// * nodes_subnet; // Número máximo de iterações do Power Method
	PetscReal lambda, dominant, lambda_previous, lambda_current; // Base para normalização
	
	PetscReal error = INFINITY;
	PetscReal epsilon = 0.000001;
	
	Vec unity;
	VecCreateSeqWithArray(PETSC_COMM_WORLD, 1, nodes_subnet, unity_array, &unity); // ok free
	
	pos = 0;
	
	Vec ax0; // Guarda primeira parcela da iteração
	PetscScalar kTx;

	// Construindo termo para normalização de B
	Vec diagB,diagBx0;
	VecCreateSeq(PETSC_COMM_WORLD,nodes_subnet,&diagB); // ok free
	VecCreateSeq(PETSC_COMM_WORLD,nodes_subnet,&diagBx0); // ok free
	
	MatMult(*net,unity,diagB); // Calculando os graus dos nós
		
	PetscScalar sumK;
	VecSum(*k_sub,&sumK);
	
	sumK = -(sumK/(2 * weight));
		
	VecAXPY(diagB,sumK,*k_sub); // diagB eh o terceiro termo
	
	//===============================================================================
	// Vou encontrar o autovetor dominante iterando somente sobre estruturas esparsas.
	// Assim, pretendo economizar memória e processamento.
	//===============================================================================
	//Determinando vetor inicial
	//Tentar vetores diferentes...
	PetscMalloc(nodes_subnet*sizeof(PetscScalar),&x0); // ok free
	for (i=0; i< nodes_subnet; i++){x0[i] = (i % 2);}
	
	VecCreateSeqWithArray(PETSC_COMM_WORLD, 1, nodes_subnet, x0, &eig_init); // ok free
	VecCopy(eig_init,*eig_result);
	
	
	VecCreateSeq(PETSC_COMM_WORLD,nodes_subnet,&diff); // ok free
	VecCreateSeq(PETSC_COMM_WORLD,nodes_subnet,&diff2); // ok free
	VecCreateSeq(PETSC_COMM_WORLD,nodes_subnet,&ax0); // ok free
	
	// Iteração do Power Method sobre estruturas esparsas
	iter = 0;
	while(iter < max_iter && error > epsilon){
				
		VecCopy(*eig_result,eig_init); // eig_init <- eig_result
		
		MatMult(*net,eig_init,ax0);
		// Primeira parcela está em ax0

		VecDot(*k_sub,eig_init,&kTx);
		// VecScale altera o vetor, por isso estou usando eig_result como auxiliar.
		VecCopy(*k_sub,*eig_result);
		VecScale(*eig_result, (PetscScalar) ((kTx)/(2 * weight)) ); 
		// Agora eu estou com (k (kT x0)) / (2 m) (segunda parcela) em eig_result
	
		VecAXPBY(*eig_result,(PetscScalar)1,(PetscScalar)-1,ax0);
		
		// Termo para normalização de B
		VecPointwiseMult(diagBx0,diagB,eig_init);
		VecAXPBY(*eig_result,(PetscScalar)-1,(PetscScalar)1,diagBx0); // Trabalhando com a terceira parcela
		// Agora o resultado está em eig_result (eig_result = B x0)
		
		VecNorm(*eig_result,NORM_INFINITY,&lambda);
		VecScale(*eig_result,(PetscScalar) 1/lambda);
			
		// Computando o resíduo
		VecCopy(eig_init,diff);
		VecCopy(*eig_result,diff2);
		VecAbs(diff);
		VecAbs(diff2);
		VecAXPBY(diff,1,-1,diff2);
		VecNorm(diff,NORM_2,&error);
		
		iter++;
	}
	VecGetValues(*eig_result,1,&pos,&lambda_current);
	VecGetValues(eig_init,1,&pos,&lambda_previous);
	
	if(lambda_previous * lambda_current < 0){
		dominant = -1;
		lambda = -lambda;
	}
	else{
		dominant = 1;
	}
	
	// eig_result agora convergiu para autovetor dominante (em magnitude).
	if(display >= 2){
		PetscPrintf(PETSC_COMM_WORLD,"\n \t Largest eigenvalue-> Domin.: ");
		PetscPrintf(PETSC_COMM_WORLD,"%f ",lambda);
	}
	//gettimeofday(&t_int, NULL);
    //t_diff = t_int.tv_sec - t_start.tv_sec + (double)(t_int.tv_usec - t_start.tv_usec)/CLOCKS_PER_SEC;
    //printf("\n\t\t ***** Dominant: %.2f segundos... (%d nós) \n", t_diff,nodes_subnet);
    
	if(dominant < 0){
		if(display >= 2){
			PetscPrintf(PETSC_COMM_WORLD,"(Negative. Shifting...)");	
		}
		Vec shift;
		VecCreateSeq(PETSC_COMM_WORLD,nodes_subnet,&shift);// ok free
		
		PetscReal new_lambda;
		// Inicializando vetor inicial
		for(i=0; i<nodes_subnet; i++){	x0[i] = (i % 2);}
		VecCopy(eig_init,*eig_result);
		
		iter = 0;
		error = INFINITY;
		// Loop para achar maior algébrico
		while(iter < max_iter && error > epsilon){
			VecCopy(*eig_result,eig_init); // eig_init <- eig_result
			VecCopy(eig_init,shift);
			MatMult(*net,eig_init,ax0);
			// Primeira parcela em ax0
	
			VecDot(*k_sub,eig_init,&kTx);
			VecCopy(*k_sub,*eig_result); // VecScale altera o vetor, por isso estou usando eig_result como auxiliar.
			VecScale(*eig_result, (PetscScalar) ((kTx)/(2 * weight)) ); 
			// Agora eu estou com (k (kT x0)) / (2 m) (segunda parcela) em eig_result

			VecAXPBY(*eig_result,(PetscScalar)1,(PetscScalar)-1,ax0);
			// Agora o resultado está em eig_result (eig_result = B x0)

			// Aplicando o shift para pular do maior em magnitude para o maior algébrico
			VecScale(shift,(PetscScalar) lambda/2);
			VecAXPBY(*eig_result,(PetscScalar)-1,(PetscScalar)1,shift);

			// Termo de normalização de B (terceira parcela)
			VecPointwiseMult(diagBx0,diagB,eig_init);
			VecAXPBY(*eig_result,(PetscScalar)-1,(PetscScalar)1,diagBx0);
			
			VecNorm(*eig_result,NORM_INFINITY,&new_lambda);
			VecScale(*eig_result,(PetscScalar) 1/new_lambda);
			
			// Computando o resíduo
			VecCopy(eig_init,diff);
			VecCopy(*eig_result,diff2);
			VecAbs(diff);
			VecAbs(diff2);
			VecAXPBY(diff,1,-1,diff2);
			VecNorm(diff,NORM_2,&error);

			iter++;
		}
		new_lambda = new_lambda + (lambda/2);
		
		if(display >= 2){
			PetscPrintf(PETSC_COMM_WORLD," Algeb.: %f",new_lambda);
			PetscPrintf(PETSC_COMM_WORLD,"\n");
		}
		VecDestroy(&shift);
	}
	
	(*nodes_comm1) = 0;
	(*nodes_comm2) = 0;
	
	PetscScalar *x1;
	PetscInt *id_x1;
	PetscMalloc(nodes_subnet*sizeof(PetscScalar),&x1); // ok free
	PetscMalloc(nodes_subnet*sizeof(PetscInt),&id_x1); // ok free
		
	for(i = 0 ; i < nodes_subnet ; i++){
		id_x1[i] = i;
	}
		
	if(nodes_subnet > 0){
		VecGetValues(*eig_result,nodes_subnet,id_x1,x1);
	}
	
	// Associando comunidades e contabilizando
	for(i = 0 ; i < nodes_subnet ; i++){
		// Olhando valores pelo x1
		if(x1[i] >= 0){
			x1[i] = 1;
			id_comm1[*nodes_comm1] = ids_subnet[i];
			(*nodes_comm1)++;
		}
		else{
			x1[i] = -1;
			id_comm2[*nodes_comm2] = ids_subnet[i];
			(*nodes_comm2)++;
		}
	}
	
	VecSetValues(*eig_result,nodes_subnet,id_x1,x1,INSERT_VALUES);
	
	// Calculando a quantidade de vértices em cada sub-comunidade
	PetscScalar commsdiff;
	VecDot(*eig_result,unity,&commsdiff);
		
	VecDestroy(&diagB);
	VecDestroy(&diagBx0);
	VecDestroy(&eig_init);
	VecDestroy(&diff);
	VecDestroy(&diff2);
	VecDestroy(&ax0);
	VecDestroy(&unity);
	PetscFree(x1);
	PetscFree(id_x1);
	PetscFree(x0);
	//PetscFree(ids_subnet);
		
	//gettimeofday(&t_end, NULL);
    //t_diff = t_end.tv_sec - t_int.tv_sec + (double)(t_end.tv_usec - t_int.tv_usec)/CLOCKS_PER_SEC;
}

void GetDeltaQ(Mat *net, Vec *k_sub, struct NetData ndata, Vec *eig_result, IS ids, PetscScalar *delta_q){

	// Variáveis "globais"
	PetscInt edges = ndata.edges;
	PetscScalar weight = ndata.weight;
	
	PetscInt nodes_subnet;
	ISGetSize(ids,&nodes_subnet);
	
	Vec unity;
	VecCreateSeqWithArray(PETSC_COMM_WORLD, 1, nodes_subnet, unity_array, &unity); // ok free
	
	Vec k_sub_copy;
	VecDuplicate(*k_sub,&k_sub_copy); // ok free
	VecCopy(*k_sub,k_sub_copy);
	
	Vec diagB,diagBy;
	// Termo para normalização de B.
	VecCreateSeq(PETSC_COMM_WORLD,nodes_subnet,&diagB); // ok free
	VecCreateSeq(PETSC_COMM_WORLD,nodes_subnet,&diagBy); // ok free
	MatMult(*net,unity,diagB);  // Pegando as somas das linhas de net
	
	PetscScalar sumK;
	VecSum(k_sub_copy,&sumK);
	sumK = -(sumK/(2 * weight));
	VecAXPY(diagB,sumK,k_sub_copy); // diagB eh o terceiro termo

	Vec ay;
	VecCreateSeq(PETSC_COMM_WORLD,nodes_subnet,&ay);  // ok free
	PetscScalar kTy;
	
	VecDot(k_sub_copy,*eig_result,&kTy);
	VecScale(k_sub_copy, (PetscScalar) ((kTy)/(2 * weight)) );
	// Segunda parcela (k (kTy/2E)) em k
	
	MatMult(*net,*eig_result,ay);
	// Primeira parcela em ay (A y)
	
	//VecAXPBY(Vec y,PetscScalar alpha,PetscScalar beta,Vec x)
	//y = alpha x + beta y
	VecAXPBY(ay,(PetscScalar)-1,(PetscScalar)1,k_sub_copy);
	
	VecPointwiseMult(diagBy,*eig_result,diagB);
	VecAXPBY(ay,(PetscScalar)-1,(PetscScalar)1,diagBy);
	
	// Agora o resultado está em ay (ay = B y)
	VecDot(*eig_result,ay,delta_q);
	*delta_q = (*delta_q) / (4 * weight);
	
	VecDestroy(&diagB);
	VecDestroy(&diagBy);
	VecDestroy(&ay);
	VecDestroy(&unity);
	VecDestroy(&k_sub_copy);
}

void FineTuningImproved(Mat *net, Vec *k, struct NetData ndata, PetscInt *nodes_comm1, PetscInt *nodes_comm2, IS ids, int *id_comm1, int *id_comm2, PetscScalar *delta_q, Vec *eig_result){
		
	// Essa versão vai ser serial. Ou seja, só estou pensando com 1 processador e esse vai fazer todas as mudanças.
	// Posteriormente eu tenho que fazer:
	//	Cada processador avalia o ganho da movimentação do seu pedaço e vê qual é a melhor movimentação local.
	//	Depois todos enviam para o rank 0 e ele vê qual é o melhor global, informa pra todo mundo qual é a nova modularidade e o processador responsável pelo vértice realiza a alteração no seu pedaço local.
	
	int variation = ndata.variation;
	
	PetscInt nodes_subnet;
	ISGetSize(ids,&nodes_subnet);
	
	const PetscInt *ids_subnet;
	
	Vec best_eig_result;
		
	VecAssemblyBegin(*eig_result);
	VecAssemblyEnd(*eig_result);
	
	VecDuplicate(*eig_result,&best_eig_result); // ok free
	VecCopy(*eig_result,best_eig_result);
	
	PetscScalar edges = ndata.edges;
	PetscScalar weight = ndata.weight;
	
	PetscInt *moved;
	PetscMalloc((nodes_subnet)*sizeof(PetscInt),&moved); //ok free
	
	int i, count;
	PetscInt best_delta_q_id;
	PetscScalar new_delta_q, best_delta_q, best_delta_q_id_value;
		
	for(i = 0 ; i < nodes_subnet ; i++){
		moved[i] = 0;
	}
	
	PetscMalloc(nodes_subnet*sizeof(PetscInt),&ids_subnet); //ok free
	ISGetIndices(ids, &ids_subnet);
	
	PetscScalar *eig_result_array;
	PetscMalloc(nodes_subnet*sizeof(PetscScalar),&eig_result_array); //ok free
	for(i = 0 ; i < nodes_subnet ; i++){
		eig_result_array[i] = 0;
	}
		
	PetscInt *id_local;
	PetscMalloc(nodes_subnet*sizeof(PetscInt),&id_local); //ok free
	for(i = 0 ; i < nodes_subnet ; i++){
		id_local[i] = i;
	}
	
	VecGetValues(best_eig_result,nodes_subnet,id_local,eig_result_array);
	
	PetscScalar *zero_array;
	PetscMalloc(nodes_subnet*sizeof(PetscScalar),&zero_array);
	for(i = 0 ; i < nodes_subnet ; i++){
		zero_array[i] = 0;
	}
	
	Vec zero;
	VecCreateSeqWithArray(PETSC_COMM_WORLD, 1,  nodes_subnet, zero_array, &zero);
	
	PetscScalar *k_sub_array;
	PetscMalloc(nodes_subnet*sizeof(PetscScalar),&k_sub_array);
	for(i = 0 ; i < nodes_subnet ; i++){
		k_sub_array[i] = 0;
	}
		
	PetscScalar *Bx_array;
	PetscMalloc(nodes_subnet*sizeof(PetscScalar),&Bx_array); //ok free
	
	// ******************************************
	// Algumas operações preliminares
	// (Estas operações são importantes para o desenvolvimento da ideia mas são comuns ao algoritmo todo, então 
	// serão feitas uma só vez)
	// ******************************************
	
	// Pegando o sub-vetor de k correspondente aos índices da comunidade corrente
	Vec k_sub;//, k_aux;
	VecDuplicate(*k,&k_sub); // ok free
	
	VecCopy(*k,k_sub);
	
	VecGetValues(k_sub,nodes_subnet,id_local,k_sub_array);
	
	Vec k_sub_copy;
	VecDuplicate(k_sub,&k_sub_copy); // ok free
	
	VecCopy(k_sub,k_sub_copy);
	
	// Mexi nessa parte porque estava dando pau. Na versão original, diagA é obtida com matgetdiagonal.
	// Mas como a diagonal de A é sempre zero, eu amarrei isso.
	Vec diagA;
	//VecDuplicate(k_sub,&diagA); // ok free //comentei pra ver se funciona. pode ser que seja importante.
	//MatGetDiagonal(*net,diagA);  //comentei pra ver se funciona. pode ser que seja importante.s
	VecCreateSeqWithArray(PETSC_COMM_WORLD, 1,  nodes_subnet, zero_array, &diagA);
	
	//VecView(diagA,0);
	//exit(0);
	
	Vec Bx;
	VecDuplicate(k_sub,&Bx); // ok free
	
	const PetscScalar *B_col_array;
	PetscMalloc(nodes_subnet*sizeof(const PetscScalar),&B_col_array);
	
	const PetscInt *cols;
	PetscMalloc(nodes_subnet*sizeof(const PetscInt),&cols);
	PetscInt ncols;
	
	Vec B_col;
	VecDuplicate(k_sub,&B_col);
	
	Vec Bx_nova;
	VecDuplicate(k_sub,&Bx_nova);
	
	PetscScalar kTy;
	
	Vec diagAx;
	VecDuplicate(k_sub,&diagAx); // ok free
	
	// ******************************************
	// Fim das operações preliminares
	// ******************************************
	
	PetscScalar si;
	
	count = 0;
		
	// ******************************************
	// Faço o cálculo de Bx uma vez só. 
	// Nas próximas eu faço só uma atualização de Bx...
	// ******************************************
	MatMult(*net,best_eig_result,Bx);
	// Primeira parcela: Vector Bx está guardando Ax

	VecPointwiseMult(diagAx,best_eig_result,diagA);
	// Segunda parcela: Vector diagAx está guardando Diag(A)x

	VecDot(k_sub,best_eig_result,&kTy);

	VecCopy(k_sub,k_sub_copy);
	VecScale(k_sub_copy, (PetscScalar) ((kTy)/(2.0 * weight)) );
	//Terceira parcela: Vector k_sub_copy está guardando (kkT x)/(2m)

	// Já vou adiantando a soma das 3 primeiras parcelas para poder aproveitar as estruturas de dados
	VecAXPBY(Bx,(PetscScalar)-1.0,(PetscScalar)1.0,diagAx);
	VecAXPBY(Bx,(PetscScalar)-1.0,(PetscScalar)1.0,k_sub_copy);

	// Agora Bx está guardando Ax - Diag(A)x - (kkT x)/(2m)

	VecCopy(k_sub,k_sub_copy);
	VecPointwiseMult(k_sub_copy, k_sub_copy,k_sub_copy);
	VecPointwiseMult(k_sub_copy,best_eig_result,k_sub_copy);
	VecScale(k_sub_copy, (PetscScalar) ((1.0)/(2.0 * weight)) );
	// Quarta parcela: Vector k_sub_copy está guradando Diag(kkT) x / (2m)

	VecAXPBY(Bx,(PetscScalar)1.0,(PetscScalar)1.0,k_sub_copy);
	// Agora Bx está guardando Ax - Diag(A)x - (kkT x)/(2m) + Diag(kkT) x / (2m)
	VecAssemblyBegin(Bx);
	VecAssemblyEnd(Bx);

	VecGetValues(Bx,nodes_subnet,id_local,Bx_array);
	// ******************************************
	// **** FIM DO CÁLCULO DE Bx ****
	// ******************************************
	
	int new_best_eig = 0;
	int max_count;
	
	if(variation == 1){
		max_count = nodes_subnet;
	}
	if(variation == 2){
		max_count = nodes_subnet / 10;
	}
	if(variation == 3){
		max_count = nodes_subnet / 5;
	}
	
	PetscBool iguais;

	while(count < max_count){
						
		// ******************************************
		// Fase inicial (calculando B * eig_result)
		// ******************************************																				
		
		//VecAXPBY(Vec y,PetscScalar alpha,PetscScalar beta,Vec x)
		//y = alpha x + beta y
		
		// Preciso calcular Bx = Ax - Diag(A)x - (kkT x)/(2m) + Diag(kkT) x / (2m)
		if(new_best_eig == 1){
			
			/*		
			// **********************************************
			// Início da ideia nova
			// **********************************************
			VecSetValue(zero,best_delta_q_id,1,INSERT_VALUES);
			VecAssemblyBegin(zero);
			VecAssemblyEnd(zero);
			MatMult(*net,zero,B_col);
			
			
			//MatGetRowUpperTriangular(*net);
			//MatGetRow(*net,best_delta_q_id,&nodes_subnet,&id_local,&B_col_array);
			//MatRestoreRow(*net,best_delta_q_id,&ncols,&cols,&B_col_array);
			//MatRestoreRowUpperTriangular(*net);
			
			//VecSetValues(B_col,nodes_subnet,cols,B_col_array,INSERT_VALUES);
			//VecAssemblyBegin(B_col); 
			//VecAssemblyEnd(B_col);
			
			
			VecCopy(k_sub,k_sub_copy);
			VecScale(k_sub_copy, (PetscScalar) (k_sub_array[best_delta_q_id]/(2.0 * weight)) );
			
			VecAXPBY(B_col,-1,1,k_sub_copy);
			VecSetValue(B_col,best_delta_q_id,0,INSERT_VALUES);
			VecAssemblyBegin(B_col);
			VecAssemblyEnd(B_col);
			
			
			VecAXPBY(Bx,-2 * best_delta_q_id_value,1,B_col);
			VecGetValues(Bx,nodes_subnet,id_local,Bx_array);
			
			VecSetValue(zero,best_delta_q_id,0,INSERT_VALUES);
			
			VecCopy(Bx,Bx_nova);
			
							
			// **********************************************
			// Fim da ideia nova
			// **********************************************
			*/
			
			
			
			// **********************************************
			// Início da ideia velha
			// **********************************************
			MatMult(*net,best_eig_result,Bx);
			// Primeira parcela: Vector Bx está guardando Ax
		
			VecPointwiseMult(diagAx,best_eig_result,diagA);
			// Segunda parcela: Vector diagAx está guardando Diag(A)x
		
			VecDot(k_sub,best_eig_result,&kTy);

			VecCopy(k_sub,k_sub_copy);
			VecScale(k_sub_copy, (PetscScalar) ((kTy)/(2.0 * weight)) );
			//Terceira parcela: Vector k_sub_copy está guardando (kkT x)/(2m)
		
			// Já vou adiantando a soma das 3 primeiras parcelas para poder aproveitar as estruturas de dados
			VecAXPBY(Bx,(PetscScalar)-1.0,(PetscScalar)1.0,diagAx);
			VecAXPBY(Bx,(PetscScalar)-1.0,(PetscScalar)1.0,k_sub_copy);
		
			// Agora Bx está guardando Ax - Diag(A)x - (kkT x)/(2m)
		
			VecCopy(k_sub,k_sub_copy);
			VecPointwiseMult(k_sub_copy, k_sub_copy,k_sub_copy);
			VecPointwiseMult(k_sub_copy,best_eig_result,k_sub_copy);
			VecScale(k_sub_copy, (PetscScalar) ((1.0)/(2.0 * weight)) );
			// Quarta parcela: Vector k_sub_copy está guradando Diag(kkT) x / (2m)
		
			VecAXPBY(Bx,(PetscScalar)1.0,(PetscScalar)1.0,k_sub_copy);
			// Agora Bx está guardando Ax - Diag(A)x - (kkT x)/(2m) + Diag(kkT) x / (2m)
			VecAssemblyBegin(Bx);
			VecAssemblyEnd(Bx);
			
			VecGetValues(Bx,nodes_subnet,id_local,Bx_array);
			// **********************************************
			// Fim da ideia velha
			// **********************************************
			
			
			/*VecEqual(Bx,Bx_nova,&iguais);
			
			
			if(iguais == PETSC_TRUE){
				printf("\niguais...");
			}
			else{
				printf("\n\t\t\tDIFERENTES!!\n");
				
				
				VecView(Bx,0);
				VecView(Bx_nova,0);
				
				exit(0);
			}
			*/
			
			
			
			
			new_best_eig = 0;
		}
		
		// ******************************************
		// Fim da fase inicial (calculando B * eig_result)
		//  ******************************************
		best_delta_q = -1; // menor valor possível
		best_delta_q_id = 0;
		best_delta_q_id_value = 1.0;
		for(i = 0 ; i < nodes_subnet ; i++){
			//current = ids_subnet[i];
					
			si = eig_result_array[i];
			
			if(moved[i] == 0){
				new_delta_q = 0;
										
				new_delta_q = -1.0 * si/(weight) * Bx_array[i];
								
				if(new_delta_q > best_delta_q){
					best_delta_q = new_delta_q;
					best_delta_q_id = i;
					best_delta_q_id_value = si;
					//printf("\nnovo best_delta_q_id_value = %f\n", best_delta_q_id_value);
				}
			}
		}
		
		//printf("\n Best deltaQ: %f (index %d)\n", best_delta_q, best_delta_q_id);
		//printf("\n Best deltaQ id value: %f (index %d)\n", best_delta_q_id_value, best_delta_q_id);
				
		moved[best_delta_q_id] = 1;
		
		if(best_delta_q_id_value == 1.0){
			eig_result_array[best_delta_q_id] = -1.0;
			VecSetValue(*eig_result,best_delta_q_id,-1.0,INSERT_VALUES);
		}
		else if(best_delta_q_id_value == -1.0){
			eig_result_array[best_delta_q_id] = 1.0;
			VecSetValue(*eig_result,best_delta_q_id,1.0,INSERT_VALUES);
		}
		else{
			printf("\nFUDEU!! (best_delta_q_id_value = %f)\n", best_delta_q_id_value);
			VecView(*eig_result,0);
			printf("\n array: ");
			for(i = 0 ; i < nodes_subnet ; i++){
				printf("%f ",eig_result_array[i]);
			}
			exit(0);
		}
		VecAssemblyBegin(*eig_result);
		VecAssemblyEnd(*eig_result);
		
		if( (*delta_q) + best_delta_q > (*delta_q) ){
			(*delta_q) = (*delta_q) + best_delta_q;
			VecCopy(*eig_result,best_eig_result);
			new_best_eig = 1;
		}
		
		if(variation == 4){
			if(*delta_q < 0){
				count = nodes_subnet;
			}
		}
		count++;
	}
	
	VecCopy(best_eig_result,*eig_result);
	
	//printf("\n \t Overall Best deltaQ: %f \n", *delta_q);
		
	PetscScalar *x1;
	PetscInt *id_x1;
	PetscMalloc(nodes_subnet*sizeof(PetscScalar),&x1); //ok free
	PetscMalloc(nodes_subnet*sizeof(PetscInt),&id_x1); //ok free
	
		
	VecGetValues(*eig_result,nodes_subnet,id_local,eig_result_array);
	(*nodes_comm1) = 0;
	(*nodes_comm2) = 0;
	// Associando comunidades e contabilizando
	for(i = 0 ; i < nodes_subnet ; i++){
		
		// Olhando valores pelo x1
		if(eig_result_array[i] == 1){
			id_comm1[*nodes_comm1] = ids_subnet[i];
			(*nodes_comm1)++;
		}
		else{
			id_comm2[*nodes_comm2] = ids_subnet[i];
			(*nodes_comm2)++;
		}
	}
	
	for(i = (*nodes_comm1) ; i < nodes_subnet ; i++){
		id_comm1[i] = 0;
	}
	for(i = (*nodes_comm2) ; i < nodes_subnet ; i++){
		id_comm2[i] = 0;
	}
	
	VecDestroy(&best_eig_result);
	PetscFree(moved);
		
	//PetscFree(ids_subnet);
	PetscFree(eig_result_array);
	PetscFree(id_local);
	PetscFree(Bx_array);
	VecDestroy(&best_eig_result);
	VecDestroy(&k_sub);
	VecDestroy(&k_sub_copy);
	VecDestroy(&diagA);
	VecDestroy(&Bx);
	VecDestroy(&diagAx);
	PetscFree(x1);
	PetscFree(id_x1);
}

void SplitCommunities(Mat *net, Vec *k_complete, struct NetData ndata, struct branch_node *bnode, IS ids, PetscInt *comms, PetscInt *ncomms, PetscScalar *q, PetscScalar *delta_q){
		
	// Variáveis "globais"
	PetscInt display = ndata.display;
	
	PetscInt variation = ndata.variation;
	
	PetscInt nodes_subnet;
	ISGetSize(ids,&nodes_subnet);
		
	
	if(display >= 2){
		PetscPrintf(PETSC_COMM_WORLD,"\n \t Try to split community with %d nodes: ", nodes_subnet);
	}
	
	Vec eig_result;
	VecCreateSeq(PETSC_COMM_WORLD, nodes_subnet, &eig_result); // ok free
	
	Mat subnet;
	MatGetSubMatrix(*net,ids,ids,MAT_INITIAL_MATRIX,&subnet); // ok free
	
	
	//PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT);
	//MatView(subnet,PETSC_VIEWER_STDOUT_WORLD);
	//VecView(k_complete,0);
	
	
	// Pegando o sub-vetor de k correspondente aos índices da comunidade corrente
	Vec k_sub, k_aux;
	// Cria k_aux do mesmo tipo do k_complete
	VecDuplicate(*k_complete,&k_aux); // ok free
	
	// Copia k_complete. Necessário pois o VecGetSubVector() mantém referência para o original.
	VecCopy(*k_complete,k_aux);  
	VecGetSubVector(k_aux,ids,&k_sub); // ok free
		
	// Dando o split na rede em duas comunidades
	PetscInt *id_comm1, *id_comm2, nodes_comm1, nodes_comm2;
	PetscMalloc(nodes_subnet*sizeof(PetscInt),&id_comm1); // ok free
	PetscMalloc(nodes_subnet*sizeof(PetscInt),&id_comm2); // ok free
		
	GetEigBSparse(&subnet,&k_sub,ndata,&nodes_comm1,&nodes_comm2, ids, id_comm1, id_comm2, &eig_result);
	GetDeltaQ(&subnet,&k_sub,ndata,&eig_result,ids,delta_q);
	
	// Aplicando o ajuste fino com Kernighan-Lin e a ideia de melhora da Sun et al
	if(variation == 1 || variation == 2 || variation == 3 || variation == 4){
		FineTuningImproved(&subnet, &k_sub, ndata, &nodes_comm1, &nodes_comm2, ids, id_comm1, id_comm2, delta_q, &eig_result);
	}
		
	// Não vou mais precisar dessas estruturas
	MatDestroy(&subnet);
	VecDestroy(&k_sub);
	VecDestroy(&k_aux);
		
	if(display >= 2){
		PetscPrintf(PETSC_COMM_WORLD,"\n == Modularity gain in splitting: %f.",*delta_q);
	}
	
	bnode->delta_q = (*delta_q);
		
	if(*delta_q > 0.0000001 && nodes_comm1 > 1 && nodes_comm2 > 1){
				
		// Partindo pras sub-redes
		
		if(display >= 2){
			PetscPrintf(PETSC_COMM_WORLD," Modularity gain is positive. OK to split. \n");
		}
		// Criando conjunto de índices para cada comunidade a ser quebrada.
		if(display >= 2){
			PetscPrintf(PETSC_COMM_WORLD,"\n == Community 1: %d node(s)",nodes_comm1);
		}
		IS ids_comm1;
		
		ISCreateGeneral(PETSC_COMM_WORLD,nodes_comm1,id_comm1,PETSC_OWN_POINTER,&ids_comm1); // ok free
		if(display >= 3){
			ISView(ids_comm1,PETSC_VIEWER_STDOUT_WORLD);
		}
		if(display >= 2){
			PetscPrintf(PETSC_COMM_WORLD,"\n == Community 2: %d node(s)",nodes_comm2);
		}

		IS ids_comm2;
		ISCreateGeneral(PETSC_COMM_WORLD,nodes_comm2,id_comm2,PETSC_OWN_POINTER,&ids_comm2); // ok free
		if(display >= 3){
			ISView(ids_comm2,PETSC_VIEWER_STDOUT_WORLD);
		}
		*q = *q + *delta_q; // Atualizando delta Q
			
		if(nodes_comm1 > 1){
			struct branch_node *bnode_left=NULL;
			bnode_left = (struct branch_node*) malloc(sizeof(struct branch_node));
			bnode_left->left = NULL;
			bnode_left->right = NULL;
			bnode_left->parent = bnode;
			bnode_left->delta_q = -1;
			bnode_left->branch_nodes = nodes_comm1;
			bnode_left->index = nbranches;
			nbranches++;
			bnode->left = bnode_left;
			
			/*
			left_index = 2 * current_index + 1;
			GenerateCommunityFile(bnode_left.branch_nodes, bnode_left.branch_index, ndata, id_comm1);
			*/
			
			// Tentar quebrar comunidade c1
			SplitCommunities(net,k_complete,ndata,bnode_left,ids_comm1,comms,ncomms,q,delta_q);
		}
		if(nodes_comm2 > 1){
			struct branch_node *bnode_right = NULL;
			bnode_right = (struct branch_node*) malloc(sizeof(struct branch_node));
			bnode_right->left = NULL;
			bnode_right->right = NULL;
			bnode_right->parent = bnode;
			bnode_right->delta_q = -1;
			bnode_right->branch_nodes = nodes_comm2;
			bnode_right->index = nbranches;
			nbranches++;
			bnode->right = bnode_right;
			
			/*
			right_index = 2 * current_index + 2;
			GenerateCommunityFile(bnode_right.branch_nodes, bnode_right.branch_index, ndata, id_comm2);
			*/
					
			// Tentar quebrar comunidade c2
			SplitCommunities(net,k_complete,ndata,bnode_right,ids_comm2,comms,ncomms,q,delta_q);
		}
	
		ISDestroy(&ids_comm1);
		ISDestroy(&ids_comm2);
	}
	else{
		if(display >= 2){
			PetscPrintf(PETSC_COMM_WORLD," Modularity gain is negative. No good to split. Updating communities... \n");
			PetscPrintf(PETSC_COMM_WORLD," ====================================================================== \n");
		}
		*ncomms = *ncomms + 1;
		if(display >= 1){
			PetscPrintf(PETSC_COMM_WORLD,"\t \t ======> Found community with %d nodes!\n", nodes_subnet);
		}
		const PetscInt *ids_subnet;
		PetscMalloc(nodes_subnet*sizeof(PetscInt),&ids_subnet);
		ISGetIndices(ids, &ids_subnet);
		
		int i;
		for(i = 0 ; i < nodes_subnet ; i++){
			comms[ids_subnet[i]] = *ncomms;
		}
	}
	VecDestroy(&eig_result);
}

void Newman(Mat *net, struct NetData ndata, FILE *out_file){
	
	// Variáveis "globais"
	PetscInt nodes = ndata.nodes;
	PetscInt display = ndata.display;
	
	PetscScalar *k_complete_array;
	PetscMalloc(nodes*sizeof(PetscScalar),&k_complete_array); // ok free
	int i;
	for(i = 0 ; i < nodes ; i++){
		k_complete_array[i] = 0;
	}
	Vec k_complete;
	VecCreateSeqWithArray(PETSC_COMM_WORLD,1,ndata.nodes,k_complete_array,&k_complete); // ok free
	
	//PetscInt nodes_subnet;
	PetscInt nodes_subnet;
	PetscInt *id_comm1, *id_comm2, *id_subnet;
	PetscScalar q, delta_q;
		
	q = 0; // Modularidade da divisão de comunidades
		
	// Alocando espaço para índices usados como entrada
	PetscMalloc(nodes*sizeof(PetscInt),&id_subnet); // ok free

	nodes_subnet = nodes;
	for(i = 0 ; i < nodes_subnet ; i++){
		id_subnet[i] = i;
	}
		
	nbranches = 0;
	
	//GenerateCommunityFile(nodes, bnode.branch_index, ndata, id_subnet);
		
	struct branch_node *bnode = NULL;
	bnode = (struct branch_node*) malloc(sizeof(struct branch_node));
	bnode->left = NULL;
	bnode->right = NULL;
	bnode->parent = NULL;
	bnode->delta_q = -1;
	bnode->branch_nodes = nodes;
	bnode->index = nbranches;
	nbranches++;
	
		
	IS ids;
	ISCreateGeneral(PETSC_COMM_WORLD,nodes_subnet,id_subnet,PETSC_OWN_POINTER,&ids); // ok free
		
	Mat subnet;
	MatGetSubMatrix(*net,ids,ids,MAT_INITIAL_MATRIX,&subnet); // ok free
	
	// Alocando espaço para índices das comunidades divididas
	PetscMalloc(nodes_subnet*sizeof(PetscInt),&id_comm1); // ok free
	PetscMalloc(nodes_subnet*sizeof(PetscInt),&id_comm2); // ok free
	
	if(display >= 2){
		PetscPrintf(PETSC_COMM_WORLD,"---> Splitting network in communities...");
	}
	
	unity_array = (PetscScalar*) malloc(sizeof(PetscScalar) * nodes); // ok free
	for(i = 0 ; i < nodes ; i++){
		unity_array[i] = 1;
	}
	
	Vec unity;
	VecCreateSeqWithArray(PETSC_COMM_WORLD, 1, nodes, unity_array, &unity); // ok free
	
	MatMult(subnet,unity,k_complete); // Calculando os graus
	
	PetscScalar weight;
	VecSum(k_complete,&weight);
	ndata.weight = weight / 2;
	
	
	
	PetscInt *comms;
	PetscMalloc(nodes*sizeof(PetscInt),&comms); // ok free
	for(i = 0 ; i < nodes ; i++){
		comms[i] = 0;
	}
	PetscInt ncomms = 0;
				
	SplitCommunities(&subnet,&k_complete,ndata,bnode,ids,comms,&ncomms,&q,&delta_q);
	
	if(ndata.hierarchy == 1){	
		PrintBranchTree(bnode);
		GenerateDendrogramScript(ncomms,nbranches,q,ndata,bnode);
	}
	
	
	// Isso tem que apagar? //BreadthFirstTree(nbranches,bnode);
	
	
	
	char comms_filename[100];
	sprintf(comms_filename,"out/%s_comms.txt",ndata.prefix_out);
	
	
	FILE *comms_file;
	comms_file = fopen(comms_filename, "w");
	
	if(comms_file == NULL){
		PetscPrintf(PETSC_COMM_WORLD,"Warning: Error writing file in %s.\n", comms_filename);
	}
	
	for(i = 0 ; i < nodes ; i++){
		fprintf(comms_file,"%d\n", comms[i]);
	}
	
	
	
	
	
	if(display >= 1){
		PetscPrintf(PETSC_COMM_WORLD,"\nFinal Modularity: %f \n",q);
		PetscPrintf(PETSC_COMM_WORLD,"Number of communities found: %d \n",ncomms);
	}
	fclose(comms_file);
	
	
	fprintf(out_file,"Final Modularity: %f \n",q);
	fprintf(out_file,"Number of communities found: %d \n",ncomms);
	
	
	
	
	
	//PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT);
	//MatView(*net,PETSC_VIEWER_STDOUT_WORLD);
	
	
	ISDestroy(&ids);
	VecDestroy(&k_complete);
	free(unity_array);
	PetscFree(k_complete_array);
	PetscFree(comms);
	VecDestroy(&unity);
	MatDestroy(&subnet);
}
