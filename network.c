#include "network.h"
#include "petscpc.h"
#include <unistd.h>
#include <sys/syscall.h>
#include <errno.h>

void PrintBranchTree(struct branch_node *bnode){

	if(bnode->left != NULL){PrintBranchTree(bnode->left);}
	if(bnode->right != NULL){PrintBranchTree(bnode->right);}
	
	if(bnode->delta_q < 0){
		bnode->delta_q = 0;
	}
	PetscScalar left_delta_q = 0;
	PetscScalar right_delta_q = 0;
	if(bnode->left != NULL){
		left_delta_q = (bnode->left)->delta_q;
	}
	if(bnode->right != NULL){
		right_delta_q = (bnode->right)->delta_q;
	}
	
	if(left_delta_q > right_delta_q){
		//bnode->delta_q = bnode->delta_q + left_delta_q;
		bnode->delta_q = bnode->delta_q + left_delta_q + right_delta_q;
	}
	else{
		//bnode->delta_q = bnode->delta_q + right_delta_q;
		bnode->delta_q = bnode->delta_q + left_delta_q + right_delta_q;
	}
	
	//printf("\n ---> (%d) delta_q = %f (%d nodes)", bnode->index, bnode->delta_q, bnode->branch_nodes);
}

void BreadthFirstTree(int nbranches, struct branch_node *bnode, int **merge, double *merge_heigth){
	
	// **************************************************************
	//	Início da busca em largura na árvore
	// **************************************************************
	struct branch_node *top = NULL;
	struct branch_node *tail = NULL;
	struct branch_node *node = NULL;
		
	int i;
	PetscScalar **branch_matrix;
	PetscMalloc((nbranches)*sizeof(PetscScalar*),&branch_matrix); // Free OK
	for(i = 0 ; i < nbranches ; i++){
		PetscMalloc((5)*sizeof(PetscScalar),&(branch_matrix[i])); // Free OK
	}
	
	int bcount = nbranches - 1;
			
	// Empilhando a raiz
	top = bnode;
	//printf("enfileirou %d\n\n",bnode->index);
	tail = bnode;
	top->queue_next = NULL;
	tail->queue_next = NULL;
		
	while(top != NULL){
		
		// Desempilhou o topo
		node = top;
		top = top->queue_next;
		if(top == NULL){tail = NULL;}
				
		//printf("\n\t ---> Visitando... (%d) delta_q = %f (%d nodes)", node->index, node->delta_q, node->branch_nodes);
		
		// Busca em largura é feita de cima pra baixo
		// Gurdando valores para poder ler de baixo pra cima
		branch_matrix[bcount][0] = node->index;
		if(node->parent != NULL){
			branch_matrix[bcount][1] = (node->parent)->delta_q;
			//branch_matrix[bcount][1] = 666;
		}
		else{
			branch_matrix[bcount][1] = 1;
		}
		branch_matrix[bcount][2] = node->branch_nodes;
		if(node->parent != NULL){
			branch_matrix[bcount][3] = (node->parent)->index;
		}
		else{
			branch_matrix[bcount][3] = -1;
		}
		if(node->left == NULL && node->right == NULL){
			branch_matrix[bcount][4] = 1;
		}
		else{
			branch_matrix[bcount][4] = 0;
		}
		bcount--;
		
		if(node->left != NULL){
			//Empilhando o left
			if(top == NULL){
				//printf("enfileirou left %d\n\n",node->left->index);
				top = node->left;
				top->queue_next = NULL;
				tail = top;
			}
			else{
				//printf("enfileirou left %d\n\n",node->left->index);
				tail->queue_next = (node->left);
				tail = (node->left);
				tail->queue_next = NULL;
			}
		}
		
		if(node->right != NULL){
			//Empilhando o right
			if(top == NULL){
				//printf("enfileirou right %d\n\n",node->right->index);
				top = node->right;
				tail->queue_next = NULL;
				tail = top;
			}
			else{
				//printf("enfileirou right %d\n\n",node->right->index);
				tail->queue_next = (node->right);
				tail = (node->right);
				tail->queue_next = NULL;
			}
		}
	}
	
	// **************************************************************
	//	Fim da busca em largura na árvore
	// **************************************************************
	//printf("\n");
	//for(i = 0 ; i < nbranches ; i++){
	//	printf("%f\t%f\t%f\t%f\t%f\n", branch_matrix[i][0], branch_matrix[i][1], branch_matrix[i][2], branch_matrix[i][3], branch_matrix[i][4]);
	//}
	
	
	// **************************************************************
	//	Início da montagem do dendrograma
	// **************************************************************
	int merge_index = -1;
	int parent_index;
	int j;
	int line_index;
	int mcount = 0;
	int count_leaf = 1;
	
	for(i = 0 ; i < nbranches ; i++){

		if(branch_matrix[i][3] != -1){
			if(branch_matrix[i][4] == 1){ // Eh folha
				merge_index = -count_leaf;
				count_leaf++;
				//printf("\n\tAchou uma folha...");
			}
			else{
				for(j = 0 ; j < mcount ; j++){
					if(merge[j][2] == branch_matrix[i][0]){ // achou linha que tem como pai o índice atual
						merge_index = j+1;
					}
				}
				//printf("\n\tNao eh folha...");
			}
		
			
			parent_index = branch_matrix[i][3];
					
			line_index = -1;
			// Descubrindo onde colocar o índice
			for(j = 0 ; j < mcount ; j++){
				if(merge[j][2] == parent_index){// Achou alguma linha com o mesmo pai
					merge[j][1] = merge_index;
					line_index = j;
				}
			}
	
	
			if(line_index == -1){
				merge[mcount][0] = merge_index;
				merge[mcount][2] = parent_index;
				if(branch_matrix[i][3] != -1){
					merge_heigth[mcount] = branch_matrix[i][1];
				}
				mcount++;
			}
		
		
		}
	}
	
	// **************************************************************
	//	Fim da montagem do dendrograma
	// **************************************************************
	
	for(i = 0 ; i < nbranches ; i++){
		PetscFree(branch_matrix[i]);
	}
	PetscFree(branch_matrix);
}

/*void DepthFirstTree(struct branch_node *bnode, int *count_leaf, int *mcount, int **merge, double *merge_heigth){

	if(bnode->left != NULL){DepthFirstTree(bnode->left,count_leaf,mcount,merge,merge_heigth);}
	if(bnode->right != NULL){DepthFirstTree(bnode->right,count_leaf,mcount,merge,merge_heigth);}
	
	printf("\n ---> (%d) delta_q = %f (%d nodes)", bnode->index, bnode->delta_q, bnode->branch_nodes);
	
	int merge_index = -1;
	int parent_index;
	int j;
	int line_index;
	
	if(bnode->parent != NULL){
		if(bnode->left == NULL && bnode->right == NULL){ // eh folha
			
			merge_index = -(*count_leaf);
			(*count_leaf)++;
			
			printf("\n\tAchou uma folha...");
		}
		else{ // nao eh folha
			for(j = 0 ; j < (*mcount) ; j++){
				if(merge[j][2] == bnode->index){ // achou linha que tem como pai o índice atual
					merge_index = j+1;
				}
			}
			printf("\n\tNao eh folha...");
		}
-hierarchy 0
	
		if(bnode->parent != NULL){
			parent_index = bnode->parent->index;
		}
		else{
			parent_index = -1;
		}
		
			
		line_index = -1;
		// Descubrindo onde colocar o índice
		for(j = 0 ; j < (*mcount) ; j++){
			if(merge[j][2] == parent_index){// Achou alguma linha com o mesmo pai
				merge[j][1] = merge_index;
				line_index = j;
			}
		}
	
	
		if(line_index == -1){
	
			merge[(*mcount)][0] = merge_index;
			merge[(*mcount)][2] = parent_index;
			if(bnode->parent != NULL){
				merge_heigth[(*mcount)] = bnode->parent->delta_q;
			}
			(*mcount)++;
		}
	}
	
}
*/


void GenerateDendrogramScript(int ncomms, int nbranches, PetscScalar q, struct NetData ndata, struct branch_node *bnode){
	
	int i;
	char r_script_filename[100];
	sprintf(r_script_filename,"out/%s_dendrogram.R",ndata.prefix_out);
	
	int **merge;
	PetscMalloc((ncomms-1)*sizeof(PetscInt*),&merge);
	for(i = 0 ; i < ncomms-1 ; i++){
		PetscMalloc((3)*sizeof(PetscInt),&(merge[i]));
	}
	
	double *merge_heigth;
	PetscMalloc((ncomms-1)*sizeof(PetscScalar),&merge_heigth);
		
	BreadthFirstTree(nbranches,bnode,merge,merge_heigth);
	
	
	/*printf("Merge:\n");
	for(i = 0 ; i < ncomms-1 ; i++){
		printf("%d\t%d\t%d\t(heigth = %f)\n",merge[i][0],merge[i][1],merge[i][2],merge_heigth[i]);
	}
	printf("\n");
	*/
	
	
	
	
	/*
	printf("\n\nMerge:\n");
	for(i = 0 ; i < ncomms-1 ; i++){
		printf("%d \t %d \t %f", merge[i][0], merge[i][1], merge_heigth[i]);
		printf("\n");
	}
	*/
			
	
	/*
	int mcount = 0;
	int count_leaf = 1;
	int merge_index = 0;
	int line_index = 0;
	
	int current_parent = 0;
	PetscInt nodes = ndata.nodes;
	BranchTreeToDendrogram(bnode,&count_leaf,&mcount,merge,merge_heigth);
	*/
			
	FILE *r_script_file; //Input network file
	r_script_file = fopen(r_script_filename, "w");
	
	if(r_script_file == NULL){
		PetscPrintf(PETSC_COMM_WORLD,"Warning: Error while reading network from %s.\n", r_script_filename);
	}
	
	fprintf(r_script_file,"a <- list()\n");
	fprintf(r_script_file,"a$merge <- matrix(c(");
	for(i = 0 ; i < ncomms-1 ; i++){
		fprintf(r_script_file,"%d,%d", merge[i][0], merge[i][1]);
		if(i != ncomms-2){
			fprintf(r_script_file,",");
			fprintf(r_script_file,"\n");	
		}
		
	}
	fprintf(r_script_file,"), nc=2, byrow=TRUE )\n");
	
	fprintf(r_script_file,"a$height <-  c(");
	for(i = 0 ; i < ncomms-1 ; i++){
		fprintf(r_script_file,"%f", merge_heigth[i]);
		if(i != ncomms-2){
			fprintf(r_script_file,",");
		}
	}
	fprintf(r_script_file,")\n");
	
	fprintf(r_script_file,"a$order <- 1:%d\n",ncomms);
	fprintf(r_script_file,"a$labels <- c(");
	for(i = 0 ; i < ncomms ; i++){
		fprintf(r_script_file,"'%d'",i+1);
		if(i != ncomms-1){
			fprintf(r_script_file,",");
		}
	}
	fprintf(r_script_file,")\n");
	fprintf(r_script_file,"class(a) <- \"hclust\"\n");
	fprintf(r_script_file,"hc <- as.dendrogram(a)\n");
	
	fprintf(r_script_file,"png(\"%s.png\",width=8,height=6,units=\"in\",res=800)\n",r_script_filename);
	fprintf(r_script_file,"par(mar=c(2,2,2,2),xaxs = \"i\",yaxs = \"i\",cex.axis=1,cex.lab=1)\n");
	
	fprintf(r_script_file,"op <- par(	cex = 1, lty = 1, lwd = 1.8, pty = \"m\" )\n");
	
	fprintf(r_script_file,"plot(hc, edgePar = list(pch = c(1,NA), cex = 5, lab.cex = 1))\n");
	
	fprintf(r_script_file,"dev.off()\n");
	
	fclose(r_script_file);
	
	for(i = 0 ; i < ncomms-1 ; i++){
		free(merge[i]);
	}
	
	free(merge_heigth);
	
	//#$#$#$#$#$#$#$#$
	
	 int rc;
     system("R CMD BATCH -file /home/vinicius/Dropbox/codigos_dsc/newman_seq/out/qqqq_dendrogram.R");
     	
	

}





void GenerateCommunityFile(int nodes_comm, int branch_index, struct NetData ndata, PetscInt *id_comm){
	char comm_filename[100];
	
	sprintf(comm_filename,"out/%s_b%d.txt",ndata.prefix_out,branch_index);
	
	FILE *comm_file;
	comm_file = fopen(comm_filename,"w");
	
	int i;
	
	for(i = 0 ; i < nodes_comm ; i++){
		fprintf(comm_file,"%d\n",id_comm[i]);
	}
	
	fclose(comm_file);
	
}



void ReadNetworkToCSR(FILE *f, struct NetData ndata, Mat *net, PetscInt *rows, PetscInt *cols, PetscScalar *values){

	// Variáveis "globais"
	PetscInt nodes = ndata.nodes;
	PetscInt edges = ndata.edges;
	PetscInt starting_node = ndata.starting_node;
	PetscInt display = ndata.display;
	
	// Variáveis para leitura dos dados no arquivo
	char interaction[150]; // Ficar de olho na quantidade de caracteres que podem aparecer em cada linha
	int count, r_count, current_from, last_from, gap;
	char *from_node_tok, *to_node_tok, *value_tok, *token;
	PetscInt from_node, to_node;
	PetscScalar value;
		
	// Variáveis da estrutura CSR
	//PetscInt *rows, *cols;
	//PetscScalar *values;
	
	if(display >= 1){
		PetscPrintf(PETSC_COMM_WORLD, "Reading network \n");
	}

	// Estruturas que vao guardar as strings lidas dos arquivos
	// Ficar de olho na quantidade de caracteres que podem aparecer
	from_node_tok = (char*) malloc(sizeof(char)*100); //destruida
    to_node_tok = (char*) malloc(sizeof(char)*100); //destruida
	value_tok = (char*) malloc(sizeof(char)*100); //destruida
	
	PetscMalloc( (nodes+1) * sizeof(PetscInt),&rows);
	PetscMalloc(edges * sizeof(PetscInt),&cols);
	PetscMalloc(edges * sizeof(PetscScalar),&values);
	
	last_from = -1;
	count = 0;
	r_count = 0;
	current_from = 0;
	while(fgets(interaction, sizeof interaction, f) != NULL){
				
		token = strtok(interaction, "\t"); // Leu o from
		strcpy(from_node_tok, (const char *) token);
		from_node = (PetscInt) atoi(from_node_tok);
		
		token = strtok(NULL, "\t"); // Leu o to
		strcpy(to_node_tok, (const char *) token);
		to_node = (PetscInt) atoi(to_node_tok);
		
		token = strtok(NULL, "\t"); // Leu o valor
		strcpy(value_tok, (const char *) token);
		value = (const PetscScalar) atof(value_tok);
		
		cols[count] = (const PetscInt) to_node - starting_node;
     	values[count] = value;
     		
     	last_from = current_from;
		current_from = from_node;
		    
	    gap = current_from - last_from;
	        
	    // Corrigindo o problema de nós sem arestas
	    while(gap != 0){
			rows[r_count] = (const PetscInt) count;
	    	r_count++;
	    	gap--;
	    }
	    count++;
	}
	
	// Colocando valor na última posição (ptr para o próximo)
	rows[r_count] = (const PetscInt) count;
	
	// Corrigindo o problema dos últimos nós não terem arestas
	int i;
	for(i = r_count ; i < (nodes+1) ; i++){
		rows[i] = count;
	}
	
	if(display >= 4){
		PetscPrintf(PETSC_COMM_WORLD,"rows:\t");
		for(i = 0 ; i < (nodes + 1) ; i++){
			printf("%d ", rows[i]);
		}
		printf("\n");
	
		PetscPrintf(PETSC_COMM_WORLD,"cols:\t");
		for(i = 0 ; i < edges ; i++){
			printf("%d ", cols[i]);
		}
		printf("\n");
		
		PetscPrintf(PETSC_COMM_WORLD,"values:\t");
		for(i = 0 ; i < edges ; i++){
			printf("%f ", values[i]);
		}
		printf("\n");
	
		PetscPrintf(PETSC_COMM_WORLD,"FIM DA EXIBICAO\n");
	}
	
	MatCreateSeqSBAIJWithArrays(PETSC_COMM_WORLD, 1, nodes, nodes, rows, cols, values, net); //destruida
	
	//Visualizacao
	if(display >= 4){
		printf("\n\nA:\n");
		//PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_DENSE);
		PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT);
		MatView(*net,PETSC_VIEWER_STDOUT_WORLD);
	}
	
	if(display >= 1){
		PetscPrintf(PETSC_COMM_WORLD, "CSR structure created.\n");
	}
				
	// Liberando a memória
	free(from_node_tok);
	free(to_node_tok);
	free(value_tok);
	//PetscFree(rows);
	//PetscFree(cols);
	//PetscFree(values);
}
