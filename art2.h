#ifndef ART2_H
#define ART2_H
#include "data.h"

typedef struct network_t network_t;
typedef struct layer1_t layer1_t;
typedef struct layer2_t layer2_t;
typedef struct matrix_t matrix_t;
typedef struct param_t param_t;
typedef struct node_t node_t;
typedef struct cluster_t cluster_t;

struct network_t
{
    param_t *params;
    double *input;
    int *inputs_track;
    double *R;
    int winner;
    int resonance;
    int reset;
    int used;
    int cycle;
    matrix_t *buttom_up;
    matrix_t *top_down;
    layer1_t *f1;
    layer2_t *f2;
    int clusers_number;
    cluster_t *clusters;
};



struct layer1_t
{
    double *W;
    double *X;
    double *U;
    double *V;
    double *P;
    double *Q;
};

struct node_t
{
    double state;
    double output;
    int inhibited;
};

struct layer2_t
{
    node_t *nodes;
};

struct matrix_t
{
    int height;
    int width;
    double **data;
};

struct param_t
{
    int dim;
    int categories;
    double sigma;
    double  a;
    double b;
    double c;
    double d;
    double vigilance;
    double e;
    double buttom_up_init_value;
    double alpha;
    char *db_learn;
    char *db_test;
    int db_size;
    int nepoch;
    int niter;

};


enum{
    DB_LEARN = 0,
    DB_TEST,
    PARAM_ALPHA,
    PARAM_SIGMA,
    PARAM_A,
    PARAM_B,
    PARAM_C,
    PARAM_D,
    PARAM_E,
    PARAM_VIGIL,
    PARAM_NBEPOCH,
    PARAM_CAT,
    PARAM_DIM,
    PARAM_ITERATION,
    PARAM_DBSIZE,
    PARAM_END
};

struct cluster_t
{
    int prototype_vector;
    int nmembers;
    int *data_vectors;
};

param_t * init_params(const char *filename);
void free_params(param_t *params);
network_t *new_network (param_t *params);
void free_network (network_t *network);
void print_vector(double *vec, int dim, const char *label);
void print_matrix(matrix_t *matrix, const char *label);
void slow_learn_art2(network_t *art2, content_t *vectors);
void quick_learn_art2(network_t *art2, content_t *vectors);
void print_clusters(network_t *art2, content_t *data);
#endif 