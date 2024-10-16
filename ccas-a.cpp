#include <bits/stdc++.h>
#include "graph.h"
#include "sct.h"
#include "pava.h"
#include "heap.h"
#include <omp.h>

double ccnm[5000][5000];
void init_cc(){
    unsigned i, j;
    for(i = 0; i < 5000; i++){
        ccnm[i][0] = 1;
        for(j = 1; j <= i; j++){
            ccnm[i][j] = ccnm[i - 1][j - 1] + ccnm[i - 1][j];
        }
    }
}
bool ccmp(rank a, rank b){
    return a.r > b.r;
}
int main(int argc, char** argv){
    std::srand(time(NULL));
    unsigned k = atoi(argv[2]);
    unsigned max_iter = atoi(argv[3]);
    float err = atof(argv[4]);
    graph *g = build_graph(argv[1],3);
    auto start_time = std::chrono::high_resolution_clock::now();
    unsigned max_kcore = atoi(argv[5]);
    init_cc();
    unsigned tmp = 0;
    double res = ccnm[max_kcore][3] / max_kcore;
    while(ccnm[tmp][2] <= res) tmp++;
    unsigned pre = g->n;
    while(1){
        graph *g_reserved = reserve(g);
        unsigned *vis = new unsigned[g->n]();
        double *count = new double[g->n]();
        bool *is = new bool[g->n]();
        for(unsigned u = 0; u < g_reserved->n; u++){
            for(unsigned i = g_reserved->sadj[u]; i < g_reserved->sadj[u + 1]; i++) vis[g_reserved->adj[i]] = 1;
            for(unsigned i = g_reserved->sadj[u]; i < g_reserved->sadj[u + 1]; i++){
                unsigned v = g_reserved->adj[i];
                for (unsigned j = g_reserved->sadj[v]; j < g_reserved->sadj[v + 1]; j++) {
                    unsigned w = g_reserved->adj[j];
                    if(vis[w]) count[u]++,count[v]++,count[w]++;
                }
            }
            for(unsigned i = g_reserved->sadj[u]; i < g_reserved->sadj[u + 1]; i++) vis[g_reserved->adj[i]] = 0;
        } 
        for(unsigned i = 0; i < g->n; i++){
            is[i] = (count[i] > ccnm[tmp - 1][2]);
        }
        g = shrink(g, is);
        delete[] vis;
        delete[] count;
        delete[] is;
        std::cout<<g->n<<" "<<g->m<<std::endl;
        pre = g->n;
        break;
    }
    re_order(g);
    auto build_end = std::chrono::high_resolution_clock::now();
    auto build_duration = std::chrono::duration_cast<std::chrono::milliseconds>(build_end - start_time);
    for(unsigned kk = 3; kk <= k; kk++){
        std::cout<<"k = "<<kk<<" "<<"|V| = "<<g->n<<"|E| = "<<g->m<<std::endl;
        SCT *TT = build_sct(g, kk);
        double rho = 0;
        unsigned sz = 0;
        unsigned t = 0;
        double *r = new double[g->n]();
        double *rr = new double[g->n]();
        if(max_iter == 0){
            while(1){
                t += 1;
                update(g, TT, kk, r, rr, rho, sz, true);
                std::pair<double,double> tmp = pava(g, rr, kk, t, rho);
                if(tmp.second / tmp.first > 1 - err) break;
            }
        }
        else{
            unsigned t;
            for(t = 0; t < max_iter; t++){
                update(g, TT, k, r, rr, rho, sz, 1);
            }
            pava(g, rr, k, t, rho);
        }
        double res = ccnm[max_kcore][kk + 1] / max_kcore;
        while(ccnm[tmp][kk] <= res) tmp++;
        for(unsigned i = 0; i < g->n; i++) r[i] = 0;
        bool *is = check(g, TT, kk, r, ccnm[tmp - 1][kk - 1]);
        g = shrink(g, is);
        re_order(g);
        delete TT;
        delete[] r;
        delete[] rr;
        delete[] is;
        if(g->n == max_kcore){
            std::cout<<"stop by redution"<<std::endl;
            break;
        }
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "Total running time  " << duration.count() << " milliseconds." << std::endl;
}