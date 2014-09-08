#include "network.h"
void GetEigBSparse(Mat *net, Vec *k, struct NetData ndata, PetscInt *nodes_comm1, PetscInt *nodes_comm2, IS ids, int *id_comm1, int *id_comm2, Vec *eig_result);
void GetDeltaQ(Mat *net, Vec *k, struct NetData ndata, Vec *eig_result, IS ids, PetscScalar *delta_q);
void FineTuningImproved(Mat *net, Vec *k, struct NetData ndata, PetscInt *nodes_comm1, PetscInt *nodes_comm2, IS ids, int *id_comm1, int *id_comm2, PetscScalar *delta_q, Vec *eig_result);
void SplitCommunities(Mat *net, Vec *k_complete, struct NetData ndata, struct branch_node *bnode, IS ids, PetscInt *comms, PetscInt *ncomms, PetscScalar *q, PetscScalar *delta_q);
void Newman(Mat *net, struct NetData ndata, FILE *out_file);
