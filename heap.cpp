#include<bits/stdc++.h>
#include"heap.h"

heap* init(unsigned n, double *r){
    heap* tmp = new heap();
    tmp->capacity = n;
    tmp->size = n;
    tmp->id = new unsigned[n + 1]();
    tmp->rev = new unsigned[n + 1]();
    tmp->val = new double[n + 1]();
    for(unsigned i = 1; i <= n; i++){
        tmp->id[i] = i - 1;
        tmp->rev[i - 1] = i;
        tmp->val[i] = r[i - 1];
    }
    for(unsigned i = n / 2; i >= 1; i--){
        heapify(tmp, i);
    }
    return tmp;
}

void heapify(heap* h, unsigned i) {
    unsigned smallest = i;
    unsigned l = 2 * i;
    unsigned r = 2 * i + 1;

    if (l <= h->size && h->val[l] < h->val[smallest])
        smallest = l;

    if (r <= h->size && h->val[r] < h->val[smallest])
        smallest = r;

    if (smallest != i) {
        // Swap values
        std::swap(h->id[i], h->id[smallest]);
        h->rev[h->id[i]] = i;
        h->rev[h->id[smallest]] = smallest;
        std::swap(h->val[i], h->val[smallest]);

        // Recursively heapify the affected sub-tree
        heapify(h, smallest);
    }
}

std::pair<unsigned, double> extract_min(heap* h) {
    if (h->size < 1) {
        std::cerr << "Heap underflow" << std::endl;
        return std::make_pair(-1,-1.0); // or throw an exception
    }
    std::pair<unsigned, double> res = std::make_pair(h->id[1], h->val[1]);
    h->id[1] = h->id[h->size];
    h->rev[h->id[1]] = 1;
    h->val[1] = h->val[h->size];
    h->size--;
    

    heapify(h, 1);

    return res;
}

void decrease_key(heap* h, unsigned i, double new_val) {
    unsigned idx = h->rev[i];
    h->val[idx] = new_val;
    while (idx > 1 && h->val[idx / 2] > h->val[idx]) {
        // Swap with parent
        std::swap(h->id[idx], h->id[idx / 2]);
        h->rev[h->id[idx]] = idx;
        h->rev[h->id[idx / 2]] = idx / 2;
        std::swap(h->val[idx], h->val[idx / 2]);

        idx /= 2;
    }
}

void destroy(heap* h) {
    delete[] h->id;
    delete[] h->rev;
    delete[] h->val;
    delete h;
}