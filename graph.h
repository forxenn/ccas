#ifndef UTILS_GRAPH_H_
#define UTILS_GRAPH_H_

typedef struct {
    unsigned s;
    unsigned t;
} edge;

typedef struct {
    unsigned n; // number of nodes
    unsigned m; // number of edges
    edge *edges; 
    unsigned *adj;
    unsigned *sadj;

} graph;


typedef struct {
    int id;
    double r;
} rank;
graph *copy(graph *g);
graph *build_graph(char *fp, int k);
void reorder(graph *g, double *r);
graph* reserve(graph *g);
graph* shrink(graph *g, bool *is);
void core_order(graph *g, unsigned k);
void re_order(graph *g);
#endif