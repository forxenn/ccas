#ifndef UTILS_SCT_H_
#define UTILS_SCT_H_

struct SCT{
    unsigned long long size;
    std::vector<unsigned> id;
    std::vector<bool> type;
    unsigned long long *sadj;
    unsigned long long *adj;
    unsigned long long *srev;
    unsigned long long *rev;
    std::vector<unsigned long long> fa;
    std::vector<std::pair<unsigned long long, unsigned long long> > edges;
};

SCT *build_sct(graph *g, int k);
void update(graph *g, SCT *T, unsigned k, double *r, double *rr, double &rho, unsigned &sz, bool get_eps = false);
double combiner(unsigned n, unsigned m);
bool *check(graph *g, SCT *T, unsigned k, double *r, double bound);
SCT *build_mf(graph *g, int k);
void build_iter(graph *g, int k, double *r, bool flag);
void build_pava(graph *g, int k, double *r, double *rr);
#endif