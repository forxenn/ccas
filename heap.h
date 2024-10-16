#ifndef UTILS_HEAP_H_
#define UTILS_HEAP_H_


typedef struct {
    unsigned size;
    unsigned capacity;
    unsigned *id;
    unsigned *rev;
    double *val;      
} heap;

heap* init(unsigned n, double* r);
void heapify(heap* h, unsigned i);
std::pair<unsigned, double> extract_min(heap* h);
void decrease_key(heap* h, unsigned i, double new_val);
void destroy(heap* h);

#endif