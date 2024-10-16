#include <bits/stdc++.h>
#include "graph.h"
#include "sct.h"
#include <omp.h>

bool cmp_y(rank a, rank b){
    return a.r == b.r ? a.id > b.id : a.r > b.r;
}


std::pair<double, double> pava(graph *g, double *r,unsigned k, unsigned t, double rho){
    rank *yorder = new rank[g->n];
    unsigned i;
    for(i = 0; i < g->n; i++){
        yorder[i].id = i;   
        yorder[i].r = r[i];
    }
    
    std::sort(yorder, yorder + g->n, cmp_y);
    double r_max = 0;
    double ssum = 0;
    bool flag = 0;
    for(i = 0; i < g->n; i++){
        ssum += yorder[i].r;
        if(flag == 0 && 1.0 * combiner(i + 1, k) / (i + 1) > ssum / (i + 1) / t){
            flag = 1;
        }
        if(flag == 0){
            r_max = std::max(r_max, 1.0 * combiner(i + 1, k) / (i + 1));
        }
        if(flag == 1){
            r_max = std::max(r_max, ssum / (i + 1) / t);
        }
    }
    std::cout<<"r_max =  "<<r_max<<::std::endl;
    printf("Maximum density is %.10lf\n", rho);
    printf("Approxiamtion rate is %.10lf\n", rho / r_max);
    delete[] yorder;
    return std::make_pair(r_max, rho);
}

std::pair<double, double> pava_iter(graph *g, double *r, double *rr, unsigned k, unsigned t){
    double rho = 0;
    rank *yorder = new rank[g->n];
    for(unsigned  i = 0; i < g->n; i++){
        yorder[i].id = i;   
        yorder[i].r = r[i];
        rr[i] = 0;
    }
    build_pava(g, k, r, rr);
    std::sort(yorder, yorder + g->n, cmp_y);
    double r_max = 0;
    double ssum = 0;
    double rhosum = 0;
    bool flag = 0;
    for(unsigned i = 0; i < g->n; i++){
        rhosum += rr[yorder[i].id];
        rho = std::max(rho, rhosum / (i + 1));
        ssum += yorder[i].r;
        if(flag == 0 && 1.0 * combiner(i + 1, k) / (i + 1) > ssum / (i + 1) / t){
            flag = 1;
        }
        if(flag == 0){
            r_max = std::max(r_max, 1.0 * combiner(i + 1, k) / (i + 1));
        }
        if(flag == 1){
            r_max = std::max(r_max, ssum / (i + 1) / t);
        }
    }
    std::cout<<"r_max =  "<<r_max<<::std::endl;
    printf("Maximum density is %.10lf\n", rho);
    printf("Approxiamtion rate is %.10lf\n", r_max / (rho));
    delete[] yorder;
    return std::make_pair(r_max, rho);
}