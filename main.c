#include "art2.h"
#include "data.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

void freeman_skapura_test()
{
    param_t params;
    params.categories = 6;
    params.db_learn = "iris.data";
    params.db_size = 1;
    params.alpha = 0.4;
    params.dim = 5;
    params.a = 10.0f;
    params.b = 10.0f;
    params.c = 0.1f;
    params.sigma = 0.2;
    params.d = 0.9f;
    params.e = 0.0000000000000000000000001f;
    params.vigilance = 0.94f;
    params.nepoch = 1;
    params.niter = 1;
    params.buttom_up_init_value = 1.0f / ((1.0f - params.d) * sqrt((double)params.dim));

    double input[] = {0.2, 0.7, 0.1, 0.5, 0.4};
    char label[] = "Test vector";
    network_t *art2 = new_network(&params);
    //data_t *database = parse_database(params->db_learn, params);
    content_t content[2];
    content[0].v = input;
    content[0].label = label;
    quick_learn_art2(art2, content);
    print_matrix(art2->top_down, "Top down weights");
    print_matrix(art2->buttom_up, "Buttom up weights");
    print_clusters(art2, content);
    free_network(art2);
}

void test_other_data_slow()
{
    param_t *params = init_params("param.cfg");
    params->buttom_up_init_value = 1.0f / ((1.0f - params->d) * sqrt((double)params->dim));
    //params->sigma = 1.0 / sqrt((double) params->dim);
    network_t *art2 = new_network(params);
    data_t *database = parse_database(params->db_learn, params);
    slow_learn_art2(art2, database->content);
    print_clusters(art2, database->content);
    free_network(art2);
    free_database(database);
    free_params(params);
}

void test_other_data_quick()
{
    param_t *params = init_params("param.cfg");
    params->buttom_up_init_value = 1.0f / ((1.0f - params->d) * sqrt((double)params->dim));
    //params->sigma = 1.0 / sqrt((double) params->dim);
    network_t *art2 = new_network(params);
    data_t *database = parse_database(params->db_learn, params);
    quick_learn_art2(art2, database->content);
    print_clusters(art2, database->content);
    free_network(art2);
    free_database(database);
    free_params(params);
}


void usage_info()
{
    printf("Usage: ./art2 [option]\n");
    printf("options:\n\t--test : to test freeman & skapura processing example\n");
    printf("\t--quick : to test an other dataset in quick learning specified in param.cfg\n");
    printf("\t--slow : to test an other dataset in quick learning specified in param.cfg\n");
}

int main(int argc, char const *argv[])
{

    if (argc == 2)
    {
        if (strncmp(argv[1], "--test", 6) == 0)
        {
            freeman_skapura_test();
        }
        else if (strncmp(argv[1],"--quick", 7) == 0)
        {
            test_other_data_quick();
        }
        else if (strncmp(argv[1],"--slow",6) == 0)
        {
            test_other_data_slow();    
        }
    }
    else 
    {
        usage_info();
    }
    return 0;
}
