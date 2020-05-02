#ifndef DATA_H
#define DATA_H



typedef struct data_t data_t;
typedef struct content_t content_t;
struct param_t;
typedef struct param_t param_t;


struct data_t
{
    int size;
    int dim;
    content_t *content;
};



struct content_t
{
    double *v;
    double norm;
    char *label;   
};




void print_database(data_t *database);
void free_database(data_t *database);
data_t *parse_database(const char *db_filename, param_t *params);


#endif 