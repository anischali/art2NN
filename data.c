#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "data.h"
#include "art2.h"


void normlize_vector(double *vec, int dim, double norm_vec)
{
    for (int ivec = 0; ivec < dim; ++ivec)
        vec[ivec] /= norm_vec;
}



data_t *parse_database(const char *db_filename, param_t *params)
{
    FILE *database = fopen(db_filename, "r");
    data_t *data = (data_t *)malloc(sizeof *data);
    data->dim = params->dim;
    data->size = params->db_size;
    content_t *content = (content_t *)malloc(params->db_size * sizeof *content);
    size_t len = 0;
    char *line = NULL, *token = NULL, *save_ptr = NULL, *end_ptr = NULL;
    ssize_t nread;
    int count = 0, count_fields = 0;
    for (int icontent = 0; icontent < params->db_size; ++icontent)
    {
        content[icontent].v = (double *)calloc(params->dim, sizeof *content[icontent].v);
    }
    while ((nread = getline(&line, &len, database)) != -1)
    {
        token = line;
        count_fields = 0;
        double norm = 0.0;
        for (char *str = token; count_fields <= params->dim; str = NULL)
        {
            token = strtok_r(str, ",", &save_ptr);
            if (token == NULL || !(count < params->db_size))
                break;
            if (token != NULL && count_fields < params->dim)
            {
                content[count].v[count_fields] = strtod(token, &end_ptr);
                norm += content[count].v[count_fields] * content[count].v[count_fields];
                end_ptr = NULL;
            }
            else
            {
                //content[count].norm = sqrt(norm);
                //normlize_vector(content[count].v, params->dim, content[count].norm);
                token[strlen(token)-1] = '\0';
                content[count].label = strdup(token);
            }
            ++count_fields;
        }
        ++count;
    }
    free(line);
    fclose(database);
    data->content = content;
    return data;
}


void print_database(data_t *database)
{
    for (int icontent = 0; icontent < database->size; ++icontent)
    {
        print_vector(database->content[icontent].v, database->dim, database->content[icontent].label);
    }
}


void free_database(data_t *database)
{
    for (int icontent = 0; icontent < database->size; ++icontent)
    {
        free(database->content[icontent].v);
        database->content[icontent].v = NULL;
        free(database->content[icontent].label);
        database->content[icontent].label = NULL;
    }
    free(database->content);
    database->content = NULL;
    free(database);
    database = NULL;
}