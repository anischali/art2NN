#include "art2.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>




param_t * init_params(const char *filename)
{
    char *parsed[PARAM_END], *endptrs[PARAM_END] = {NULL};
    param_t *params = (param_t *) malloc(sizeof *params);
    
    FILE *params_file = fopen(filename, "r");

    size_t len = 0;
    char *line = NULL, *conf = NULL;
    ssize_t nread;
    int iconfig = 0;
    while ((nread = getline(&line, &len, params_file)) != -1)
    {
        if (iconfig < 3)
            line[strlen(line) - 1] = '\0';
        conf = line;
        while (*conf != ':')
            ++conf;
        parsed[iconfig++] = strdup((++conf));
    }
    free(line);
    line = NULL;
    params->db_learn = strdup(parsed[DB_LEARN]);
    params->db_test =  strdup(parsed[DB_TEST]);
    params->alpha = strtod(parsed[PARAM_ALPHA], &endptrs[PARAM_ALPHA]);
    params->sigma = strtod(parsed[PARAM_SIGMA], &endptrs[PARAM_SIGMA]);
    params->a = strtod(parsed[PARAM_A], &endptrs[PARAM_A]);
    params->b = strtod(parsed[PARAM_B], &endptrs[PARAM_B]);
    params->c = strtod(parsed[PARAM_C], &endptrs[PARAM_C]);
    params->d = strtod(parsed[PARAM_D], &endptrs[PARAM_D]);
    params->e = strtod(parsed[PARAM_E], &endptrs[PARAM_E]);
    params->vigilance = strtod(parsed[PARAM_VIGIL], &endptrs[PARAM_VIGIL]);
    params->nepoch = atoi(parsed[PARAM_NBEPOCH]);
    params->categories = atoi(parsed[PARAM_CAT]);
    params->dim = atoi(parsed[PARAM_DIM]);
    params->niter = atoi(parsed[PARAM_ITERATION]);
    params->db_size = atoi(parsed[PARAM_DBSIZE]);

    for (iconfig = 0; iconfig < PARAM_END; ++iconfig)
    {  
        if (parsed[iconfig] != NULL)
        {
            free(parsed[iconfig]);
            parsed[iconfig] = NULL;
        } 
    }
    fclose(params_file);
    return params;
}



void free_params(param_t *params)
{
    free(params->db_learn);
    params->db_learn = NULL;
    free(params->db_test);
    params->db_test = NULL;
    free(params);
    params = NULL;
}