#include "art2.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

static layer1_t *new_layer1(int dim)
{
    layer1_t *layer = (layer1_t *)malloc(sizeof *layer);
    layer->W = (double *)calloc(dim, sizeof *layer->W);
    layer->X = (double *)calloc(dim, sizeof *layer->X);
    layer->U = (double *)calloc(dim, sizeof *layer->U);
    layer->V = (double *)calloc(dim, sizeof *layer->V);
    layer->Q = (double *)calloc(dim, sizeof *layer->Q);
    layer->P = (double *)calloc(dim, sizeof *layer->P);
    return layer;
}

static void free_layer1(layer1_t *layer)
{
    free(layer->W);
    free(layer->X);
    free(layer->U);
    free(layer->V);
    free(layer->P);
    free(layer->Q);
    free(layer);
    layer = NULL;
}

static layer2_t *new_layer2(int node_dim)
{
    layer2_t *layer = (layer2_t *)malloc(sizeof *layer);
    layer->nodes = (node_t *)calloc(node_dim, sizeof *layer->nodes);
    return layer;
}

static void free_layer2(layer2_t *layer)
{
    free(layer->nodes);
    free(layer);
    layer = NULL;
}

static matrix_t *new_matrix(int h, int w)
{
    matrix_t *matrix = (matrix_t *)malloc(sizeof *matrix);
    matrix->height = h;
    matrix->width = w;
    double **mat = (double **)calloc(h, sizeof *mat);
    for (int iline = 0; iline < h; ++iline)
    {
        mat[iline] = (double *)calloc(w, sizeof *mat[iline]);
    }
    matrix->data = mat;
    return matrix;
}

static void free_matrix(matrix_t *matrix)
{
    for (int iline = 0; iline < matrix->height; ++iline)
    {
        free(matrix->data[iline]);
    }
    free(matrix->data);
    free(matrix);
    matrix = NULL;
}

void matrix_assign(matrix_t *matrix, double value)
{
    double **data = matrix->data;
    for (int iline = 0; iline < matrix->height; ++iline)
    {
        for (int icol = 0; icol < matrix->width; ++icol)
        {
            data[iline][icol] = value;
        }
    }
}

network_t *new_network(param_t *params)
{
    network_t *network = (network_t *)malloc(sizeof *network);
    param_t *p = (param_t *)malloc(sizeof *p);
    bcopy(params, p, sizeof(param_t));
    network->params = p;
    network->input = NULL;
    double *r = (double *)calloc(params->dim, sizeof *r);
    network->R = r;
    network->f1 = new_layer1(params->dim);
    network->f2 = new_layer2(params->categories);
    network->buttom_up = new_matrix(params->categories, params->dim);
    network->top_down = new_matrix(params->dim, params->categories);
    matrix_assign(network->buttom_up, params->buttom_up_init_value);
    network->inputs_track = (int *)calloc(params->db_size, sizeof *network->inputs_track);
    network->clusters = NULL;
    network->clusers_number = 0;
    network->used = 1;
    network->reset = 1;
    network->resonance = 0;
    network->cycle = 1;
    network->winner = -1;
    return network;
}

void free_network(network_t *network)
{
    free(network->params);
    free(network->R);
    free_layer1(network->f1);
    free_layer2(network->f2);
    free_matrix(network->buttom_up);
    free_matrix(network->top_down);
    if (network->clusters != NULL)
    {
        for (int icluster = 0; icluster < network->clusers_number; ++icluster)
        {
            free(network->clusters[icluster].data_vectors);
        }
    }
    free(network->clusters);
    free(network->inputs_track);
    free(network);
    network = NULL;
}

void assign_value(double *vector, int dim, double value)
{
    for (int ivec = 0; ivec < dim; ++ivec)
    {
        vector[ivec] = value;
    }
}

void copy_f1(layer1_t *dst, layer1_t *src, int dim)
{
    for (int ivec = 0; ivec < dim; ++ivec)
    {
        dst->W[ivec] = src->W[ivec];
        dst->X[ivec] = src->X[ivec];
        dst->U[ivec] = src->U[ivec];
        dst->V[ivec] = src->V[ivec];
        dst->P[ivec] = src->P[ivec];
        dst->Q[ivec] = src->Q[ivec];
    }
}

double norm(double *vec, int dim)
{
    double sum = 0.0;
    for (int ivec = 0; ivec < dim; ++ivec)
        sum += (vec[ivec]) * (vec[ivec]);
    return sqrt(sum);
}

double fabs(double v)
{
    return (v < 0) ? v * -1.0 : v;
}

double diff(double *v1, double *v2, int dim)
{
    double diffirence = 0.0;
    for (int ivec = 0; ivec < dim; ++ivec)
    {
        diffirence += fabs(v2[ivec] - v1[ivec]);
    }
    return diffirence / (double)dim;
}

double f(double v, double sigma)
{
    return ((v > sigma) ? v : 0.0);
}

void print_vector(double *vec, int dim, const char *label)
{
    printf("%s : ", label);
    for (int ivec = 0; ivec < dim; ++ivec)
        printf("%0.3f ", vec[ivec]);
    printf("\n");
}

void print_matrix(matrix_t *matrix, const char *label)
{
    printf("%s : \n", label);
    for (int ih = 0; ih < matrix->height; ++ih)
    {
        for (int iw = 0; iw < matrix->width; ++iw)
        {
            printf("%0.3f ", matrix->data[ih][iw]);
        }
        printf("\n");
    }
}

void reset_stm_f1(network_t *art2)
{
    int dim = art2->params->dim;
    layer1_t *f1 = art2->f1;
    for (int ivec = 0; ivec < dim; ++ivec)
    {
        f1->W[ivec] = f1->X[ivec] = f1->V[ivec] = f1->U[ivec] = f1->P[ivec] = f1->Q[ivec] = 0.0;
    }
}

void equation_w(network_t *art2)
{
    int dim = art2->params->dim;
    layer1_t *f1 = art2->f1;
    param_t *params = art2->params;
    for (int ivec = 0; ivec < dim; ++ivec)
        f1->W[ivec] = art2->input[ivec] + ((art2->winner > -1) ? (params->a * f1->U[ivec]) : 0.0);
}

void equation_x(network_t *art2)
{
    int dim = art2->params->dim;
    layer1_t *f1 = art2->f1;
    param_t *params = art2->params;
    double norm_w = norm(f1->W, dim);
    for (int ivec = 0; ivec < dim; ++ivec)
        f1->X[ivec] = f1->W[ivec] / (params->e + norm_w);
}

void equation_v(network_t *art2)
{
    int dim = art2->params->dim;
    layer1_t *f1 = art2->f1;
    param_t *params = art2->params;
    for (int ivec = 0; ivec < dim; ++ivec)
        f1->V[ivec] = f(f1->X[ivec], params->sigma) + params->b * f(f1->Q[ivec], params->sigma);
}

void equation_u(network_t *art2)
{
    int dim = art2->params->dim;
    layer1_t *f1 = art2->f1;
    param_t *params = art2->params;
    double norm_v = norm(f1->V, dim);
    for (int ivec = 0; ivec < dim; ++ivec)
        f1->U[ivec] = f1->V[ivec] / (params->e + norm_v);
}

void equation_p(network_t *art2)
{
    int dim = art2->params->dim;
    layer1_t *f1 = art2->f1;
    param_t *params = art2->params;
    for (int ivec = 0; ivec < dim; ++ivec)
        f1->P[ivec] = f1->U[ivec] + ((art2->winner > -1) ? (params->d * art2->buttom_up->data[art2->winner][ivec]) : 0.0);
}

void equation_q(network_t *art2)
{
    int dim = art2->params->dim;
    layer1_t *f1 = art2->f1;
    param_t *params = art2->params;
    double norm_p = norm(f1->P, dim);
    for (int ivec = 0; ivec < dim; ++ivec)
        f1->Q[ivec] = f1->P[ivec] / (params->e + norm_p);
}

void equation_r(network_t *art2)
{
    int dim = art2->params->dim;
    layer1_t *f1 = art2->f1;
    param_t *params = art2->params;
    double norm_p = norm(f1->P, dim);
    double norm_u = norm(f1->U, dim);
    for (int ivec = 0; ivec < dim; ++ivec)
        art2->R[ivec] = (f1->U[ivec] + params->c * f1->P[ivec]) / (norm_u + params->c * norm_p);
}

void print_layer_f1(layer1_t *f1, int dim)
{
    printf("Layer f1\n");
    print_vector(f1->W, dim, "W ");
    print_vector(f1->X, dim, "X ");
    print_vector(f1->V, dim, "V ");
    print_vector(f1->U, dim, "U ");
    print_vector(f1->P, dim, "P ");
    print_vector(f1->Q, dim, "Q ");
}

void print_layer_f2(network_t *art2)
{
    printf("Layer f2\n");
    node_t *nodes = art2->f2->nodes;
    int dim_nodes = art2->params->categories;
    for (int inode = 0; inode < dim_nodes; ++inode)
    {
        printf("%0.3f ", nodes[inode].state);
    }
    printf("\n");
}

void ltm_butom_up_learn_slow(network_t *art2)
{
    int dim = art2->params->dim;
    layer1_t *f1 = art2->f1;
    param_t *params = art2->params;
    double **mat_bu = art2->buttom_up->data;
    for (int ivec = 0; ivec < dim; ++ivec)
    {
        mat_bu[art2->winner][ivec] = params->alpha * params->d * f1->U[ivec] +
                                     (1.0 + params->alpha * params->d * (params->d - 1.0)) * mat_bu[art2->winner][ivec];
    }
}

void ltm_butom_up_learn_quick(network_t *art2)
{
    int dim = art2->params->dim;
    layer1_t *f1 = art2->f1;
    param_t *params = art2->params;
    double **mat_bu = art2->buttom_up->data;
    for (int ivec = 0; ivec < dim; ++ivec)
    {
        mat_bu[art2->winner][ivec] = f1->U[ivec] / (1.0 - params->d);
    }
}

void ltm_top_down_learn_slow(network_t *art2)
{
    int dim = art2->params->dim;
    layer1_t *f1 = art2->f1;
    param_t *params = art2->params;
    double **mat_td = art2->top_down->data;
    for (int ivec = 0; ivec < dim; ++ivec)
    {
        mat_td[ivec][art2->winner] = params->alpha * params->d * f1->U[ivec] +
                                     (1.0 + params->alpha * params->d * (params->d - 1.0)) * mat_td[ivec][art2->winner];
    }
}

void ltm_top_down_learn_quick(network_t *art2)
{
    int dim = art2->params->dim;
    layer1_t *f1 = art2->f1;
    param_t *params = art2->params;
    double **mat_td = art2->top_down->data;
    for (int ivec = 0; ivec < dim; ++ivec)
    {
        mat_td[ivec][art2->winner] = f1->U[ivec] / (1.0 - params->d);
    }
}

void reset_stm_f2(network_t *art2)
{
    node_t *nodes = art2->f2->nodes;
    int dim_nodes = art2->params->categories;
    for (int inode = 0; inode < dim_nodes; ++inode)
    {
        nodes[inode].state = 0.0;
        nodes[inode].output = 0.0;
        nodes[inode].inhibited = 0;
    }
}

void preprocess_stm_f1(network_t *art2)
{
    int dim = art2->params->dim;
    layer1_t *f1 = art2->f1;
    equation_w(art2);
    assign_value(f1->Q, dim, 0.0);
    assign_value(f1->U, dim, 0.0);
    equation_x(art2);
    assign_value(f1->P, dim, 0.0);
    equation_v(art2);
}

void process_stm_f1(network_t *art2)
{
    equation_u(art2);
    equation_w(art2);
    equation_p(art2);
    equation_x(art2);
    equation_q(art2);
    equation_v(art2);
}

void process_stm_f2(network_t *art2)
{
    int dim = art2->params->dim;
    layer1_t *f1 = art2->f1;
    param_t *params = art2->params;
    double **bu = art2->buttom_up->data;
    node_t *nodes = art2->f2->nodes;
    for (int inode = 0; inode < params->categories; ++inode)
    {
        nodes[inode].state = 0.0;
        if (nodes[inode].inhibited == 0)
        {
            for (int ivec = 0; ivec < dim; ++ivec)
            {
                nodes[inode].state += bu[inode][ivec] * f1->P[ivec];
            }
        }
    }
}

void stm_f2_wta(network_t *art2)
{
    double max = -INFINITY;
    param_t *params = art2->params;
    node_t *nodes = art2->f2->nodes;
    for (int inode = 0; inode < params->categories; ++inode)
    {
        if (max < nodes[inode].state && nodes[inode].inhibited == 0)
        {
            max = nodes[inode].state;
            art2->winner = inode;
        }
    }
    if (max <= 0.0)
        art2->winner = -1;
}

void orienting_subsystem(network_t *art2)
{
    equation_u(art2);
    equation_p(art2);
    equation_r(art2);
}

int urand(int min, int max)
{
    return (int)rand() % (max - min) + min - 1;
}

int *shuffle_vector(int *vec, int size)
{
    int rnd = 0;
    if (vec == NULL)
    {
        int *v = (int *)calloc(size, sizeof *v);
        vec = v;
        for (int idx = 0; idx < size; ++idx)
            vec[idx] = idx;
    }
    for (int idx = 0; idx < size; ++idx)
    {
        rnd = urand(idx, size);
        vec[idx] ^= vec[rnd];
        vec[rnd] ^= vec[idx];
        vec[idx] ^= vec[rnd];
    }
    return vec;
}

int get_clusters_number(int *vec, int dim)
{
    int max = -1000;
    for (int ivec = 0; ivec < dim; ++ivec)
    {
        if (max < vec[ivec])
        {
            max = vec[ivec];
        }
    }
    return max + 1;
}

int count_cluster_members(network_t *art2, int prototype)
{
    int *inputs_track = art2->inputs_track,
        size = art2->params->db_size, count = 0;
    for (int ivec = 0; ivec < size; ++ivec)
    {
        if (inputs_track[ivec] == prototype)
            ++count;
    }
    return count;
}

void get_cluster_members(int *inputs_track, int size, cluster_t *cluster)
{
    int idx = 0;
    for (int ivec = 0; ivec < size && idx < cluster->nmembers; ++ivec)
    {
        if (inputs_track[ivec] == cluster->prototype_vector)
        {
            cluster->data_vectors[idx] = ivec;
            ++idx;
        }
    }
}

void make_clusers(network_t *art2)
{
    art2->clusers_number = get_clusters_number(art2->inputs_track, art2->params->db_size);
    cluster_t *clusters = calloc(art2->clusers_number, sizeof(*clusters));
    for (int icluster = 0; icluster < art2->clusers_number; ++icluster)
    {
        clusters[icluster].prototype_vector = icluster;
        clusters[icluster].nmembers = count_cluster_members(art2, icluster);
        clusters[icluster].data_vectors = calloc(clusters[icluster].nmembers, sizeof *clusters[icluster].data_vectors);
        get_cluster_members(art2->inputs_track, art2->params->db_size, &clusters[icluster]);
    }
    art2->clusters = clusters;
}

void slow_learn_art2(network_t *art2, content_t *vectors)
{
    int dim = art2->params->dim;
    param_t *params = art2->params;
    int *rnd_access = NULL;
    node_t *nodes = art2->f2->nodes;
    for (int iepoch = 0; iepoch < params->nepoch; ++iepoch)
    {
        rnd_access = shuffle_vector(rnd_access, params->db_size);
        for (int idata = 0; idata < params->db_size; ++idata)
        {
            art2->input = vectors[rnd_access[idata]].v;
            preprocess_stm_f1(art2);
            process_stm_f1(art2);
            process_stm_f2(art2);
            do
            {
                stm_f2_wta(art2);
                orienting_subsystem(art2);
                double norm_r = norm(art2->R, dim);
                if (norm_r < params->vigilance - params->e)
                {
                    //printf("NR: %0.6f p: %.6f epoch: %d %d\n", norm_r, params->vigilance - params->e, iepoch, art2->winner);
                    nodes[art2->winner].inhibited = 1;
                    art2->winner = -1;
                    art2->reset = 1;
                }
                if (norm_r >= params->vigilance - params->e)
                {
                    art2->reset = 0;
                    equation_w(art2);
                    equation_x(art2);
                    equation_q(art2);
                    equation_v(art2);
                }
            } while (art2->reset == 1);
            if (art2->winner == -1)
            {
                reset_stm_f1(art2);
                reset_stm_f2(art2);
                continue;
            }
            for (int iter = 0; iter < params->niter; ++iter)
            {
                art2->inputs_track[rnd_access[idata]] = art2->winner;
                ltm_butom_up_learn_slow(art2);
                ltm_top_down_learn_slow(art2);
                process_stm_f1(art2);
            }
            reset_stm_f2(art2);
        }
    }
    free(rnd_access);
    make_clusers(art2);
}

void quick_learn_art2(network_t *art2, content_t *vectors)
{
    int dim = art2->params->dim;
    param_t *params = art2->params;
    int *rnd_access = NULL;
    node_t *nodes = art2->f2->nodes;
    for (int iepoch = 0; iepoch < params->nepoch; ++iepoch)
    {
        rnd_access = shuffle_vector(rnd_access, params->db_size);
        for (int idata = 0; idata < params->db_size; ++idata)
        {
            art2->input = vectors[rnd_access[idata]].v;
            preprocess_stm_f1(art2);
            process_stm_f1(art2);
            process_stm_f2(art2);
            do
            {
                stm_f2_wta(art2);
                orienting_subsystem(art2);
                double norm_r = norm(art2->R, dim);
                if (norm_r < params->vigilance - params->e)
                {
                    //printf("NR: %0.6f p: %.6f epoch: %d %d\n", norm_r, params->vigilance - params->e, iepoch, art2->winner);
                    nodes[art2->winner].inhibited = 1;
                    art2->winner = -1;
                    art2->reset = 1;
                }
                if (norm_r >= params->vigilance - params->e)
                {
                    art2->reset = 0;
                    equation_w(art2);
                    equation_x(art2);
                    equation_q(art2);
                    equation_v(art2);
                }
            } while (art2->reset == 1);
            if (art2->winner == -1)
            {
                reset_stm_f1(art2);
                reset_stm_f2(art2);
                continue;
            }
            art2->inputs_track[rnd_access[idata]] = art2->winner;
            ltm_butom_up_learn_quick(art2);
            ltm_top_down_learn_quick(art2);
            process_stm_f1(art2);
            reset_stm_f2(art2);
        }
    }
    free(rnd_access);
    make_clusers(art2);
}

void print_clusters(network_t *art2, content_t *data)
{
    cluster_t *clusters = art2->clusters;
    int dim = art2->params->dim;
    printf("Clusters number: %d\n", art2->clusers_number);

    for (int icluster = 0; icluster < art2->clusers_number; ++icluster)
    {
        printf("\n\n");
        print_vector(art2->buttom_up->data[clusters[icluster].prototype_vector], dim, "Prototype vector");
        printf("\n");
        for (int imember = 0; imember < clusters[icluster].nmembers; ++imember)
        {
            print_vector(data[clusters[icluster].data_vectors[imember]].v, dim, data[clusters[icluster].data_vectors[imember]].label);
        }
    }
}

