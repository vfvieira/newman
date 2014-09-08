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

void parse_args(int argc, char **argv, int *nodes, int *edges, int *starting_node, int *display, int *variation, int *hierarchy, char *filename_in, char *prefix_out){
		
	*nodes = -1;
	*edges = -1;
	*starting_node = -1;
	*display = -1;
	sprintf(prefix_out,"%s","output");
	*variation = 0;
	*hierarchy = 0;
		
	if(argc < 2){
		printf("\nError: Missing fiename.");
		exit(0);
	}
	if(argc % 2 == 0){
		printf("\nWrong usage.\n");
		exit(0);
	}
	
	PetscOptionsGetInt(PETSC_NULL,"-n",nodes,PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL,"-e",edges,PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL,"-s",starting_node,PETSC_NULL);
	PetscOptionsGetString(PETSC_NULL,"-i",filename_in,100,PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL,"-d",display,PETSC_NULL);
	PetscOptionsGetString(PETSC_NULL,"-o",prefix_out,100,PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL,"-v",variation,PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL,"-hierarchy",hierarchy,PETSC_NULL);

	/*int i;
	for(i = 1 ; i < argc ; i++){
		if(argv[i][0] == '-'){
			switch(argv[i][1]){
				case 'n': // setando número de nós
					*nodes = atoi(argv[i+1]);
					i++;
					break;
				case 'e':
					*edges = atoi(argv[i+1]);
					i++;
					break;
				case 's':
					*starting_node = atoi(argv[i+1]);
					i++;
					break;
				case 'i':
					sprintf(filename_in,"%s",argv[i+1]);
					i++;
					break;
				case 'd':
					*display = atoi(argv[i+1]);
					i++;
					break;
				case 'o':
					sprintf(prefix_out,"%s",argv[i+1]);
					i++;
					break;
				case 'v':
					*variation = atoi(argv[i+1]);
					i++;
					break; 
				default:
					printf("\nUnknown option -%s\n",argv[i+1]);
					i++;
					break;
			}
		}
	}
	*/
	
	if(*display == -1){*display = 1;}
	if(*starting_node == -1){*starting_node = 1;}
}

int main(int argc, char *argv[]){
	
	PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
	FILE *net_file; //Input network file
	struct NetData ndata;
	Mat net;
	
	char filename_in[100];
	char prefix_out[100];
	
	int nodes, edges, starting_node, display, variation, hierarchy;
	parse_args(argc,argv,&nodes,&edges,&starting_node,&display,&variation,&hierarchy,filename_in,prefix_out);
	
	sprintf(ndata.prefix_out,"%s",prefix_out);

	ndata.nodes = nodes;
	ndata.edges = edges;
	ndata.starting_node = starting_node;
	ndata.display = display;
	ndata.nbranches = 0;
	ndata.variation = variation;
	ndata.hierarchy = hierarchy;
	
	PetscInt i;
	struct timeval start, end;
	gettimeofday(&start, NULL);
	
	float diff;
	
	net_file = fopen(filename_in, "r");
	if(ndata.display >= 1){
		PetscPrintf(PETSC_COMM_WORLD,"Reading network from file %s...\n", filename_in);
	}
	if (net_file == NULL){
		PetscPrintf(PETSC_COMM_WORLD,"Error while reading network from %s.\n", filename_in);
		return 0;
	}
			
	if(ndata.display >= 1){
		PetscPrintf(PETSC_COMM_WORLD, "Network successfully read: %d nodes and %d edges. ", ndata.nodes, ndata.edges);
		PetscPrintf(PETSC_COMM_WORLD, "Here we go!\n");
	 	PetscPrintf(PETSC_COMM_WORLD,"Generating CSR matrices...\n");
   	}
   	
	PetscScalar *k_complete_array;
	PetscMalloc((nodes)*sizeof(PetscScalar),&k_complete_array); //destruida
	
	for(i = 0 ; i < nodes ; i++){
		k_complete_array[i] = 0;
	}
	Vec k_complete;
	VecCreateSeqWithArray(PETSC_COMM_WORLD,1,nodes,k_complete_array,&k_complete); //destruida
	
	// Variáveis da estrutura CSR
	PetscInt *rows, *cols;
	PetscScalar *values;
	
	ReadNetworkToCSR(net_file, ndata, &net, rows, cols, values);
	
	char out_filename[100];
	sprintf(out_filename,"out/%s_out.txt",ndata.prefix_out);
	FILE *out_file;
	out_file = fopen(out_filename, "w");
	
	Newman(&net,ndata,out_file);
	
	fflush(stdout);
	fflush(stderr);
	
   	fclose(net_file);
   	//PetscFree(rows);
	//PetscFree(cols);
	//PetscFree(values);
    PetscFree(k_complete_array);
    VecDestroy(&k_complete);
    MatDestroy(&net);
    
    gettimeofday(&end, NULL);
    diff = end.tv_sec - start.tv_sec + (double)(end.tv_usec - start.tv_usec)/CLOCKS_PER_SEC;
    
   if(ndata.display >= 1){
    	PetscPrintf(PETSC_COMM_WORLD,"Done. Total time: %.2f seconds\n", diff);
	}
	fprintf(out_file,"Done. Total time: %.2f seconds\n", diff);
	
	fclose(out_file);
    PetscFinalize();
	return 0;
}
