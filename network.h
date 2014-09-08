#ifndef _NETWORK_H
#define _NETWORK_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "petscmat.h"
#include "petscvec.h"
#include "petscis.h"

struct NetData{
	PetscInt nodes, edges, starting_node;
	PetscInt display;
	PetscInt nbranches;
	char prefix_out[100];
	PetscInt variation;
	PetscInt hierarchy;
	PetscScalar weight;
};

struct branch_node{
	int index;
	int branch_nodes;
	int has_child;
	double delta_q;
	
	struct branch_node *left;
	struct branch_node *right;
	struct branch_node *parent;
	
	struct branch_node *queue_next;
};


//void ReadNetworkToCSR(FILE *f, struct NetData ndata, Mat *net);
void ReadNetworkToCSR(FILE *f, struct NetData ndata, Mat *net, PetscInt *rows, PetscInt *cols, PetscScalar *values);
void GenerateDendrogramScript(int ncomms, int nbranches, PetscScalar q, struct NetData ndata, struct branch_node *bnode);
void GenerateCommunityFile(int nodes_comm, int branch_index, struct NetData ndata, PetscInt *id_comm);
void PrintBranchTree(struct branch_node *bnode);
void DepthFirstTree(struct branch_node *bnode, int *count_leaf, int *mcount, int **merge, double *merge_heigth);

void BreadthFirstTree(int nbranches,struct branch_node *bnode, int **merge, double *merge_heigth);
void Enqueue(struct branch_node *bnode, struct branch_node *top);

#endif  /* _NETWORK_H */
