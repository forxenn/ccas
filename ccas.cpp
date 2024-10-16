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
    std::cout << "Computing for k = " << k << "  Max_iter is: " << max_iter << std:: endl;

    graph *g = build_graph(argv[1],k);
    auto start_time = std::chrono::high_resolution_clock::now();
    unsigned max_kcore = atoi(argv[5]);
    init_cc();
    double res = ccnm[max_kcore][k] / max_kcore;
    unsigned tmp = 0;
    while(ccnm[tmp][k - 1] <= res) tmp++;
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
        if(g->n == pre) break;
        pre = g->n;
        delete[] vis;
        delete[] count;
        delete[] is;
        std::cout<<g->n<<" "<<g->m<<std::endl;
        
    }
    if(g->n == max_kcore){
        std::cout<<"stop by redution";
        std::cout<<"density : "<<res<<" "<<"size : "<<max_kcore<<std::endl;
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        std::cout << "Total running time  " << duration.count() << " milliseconds." << std::endl;
        return 0;
    }
    re_order(g);
    std::cout<<g->n<<" "<<g->m<<std::endl;
    auto tmp_time = std::chrono::high_resolution_clock::now();
    auto duration3 = std::chrono::duration_cast<std::chrono::milliseconds>(tmp_time - start_time);
    std::cout << "3-clique  " << duration3.count() << " milliseconds." << std::endl;
    SCT *TT = build_sct(g, k);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - tmp_time);
    std::cout << "Indexing time  " << duration.count() << " milliseconds." << std::endl;
    
    
    double *r = new double[g->n]();
    double *rr = new double[g->n]();
     
    double rho = 0;
    unsigned sz = 0;
    if (max_iter == 0) {
        double rho = 0;
        auto start_time_1 = std::chrono::high_resolution_clock::now();
        unsigned T = 1, t = 0;
        while(1){
            t++;
            update(g, TT, k, r, rr, rho, sz, 1);
            std::pair<double, double> res = pava(g, r, k, t, rho);
            auto end_time_1 = std::chrono::high_resolution_clock::now();
            auto duration_1 = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_1 - start_time);
            std::cout << "Total running time after iteration: " << t << "  " <<  duration_1.count() << " milliseconds." << std::endl;
            std::cout<<res.second<<" "<<res.first<<std::endl;
            if(res.second / res.first > 1 - err) break;
        }
    } else {
        unsigned t;
        for(t = 0; t < max_iter; t++){
            update(g, TT, k, r, rr, rho, sz, 1);
        }
        pava(g, rr, k, t, rho);
    }
    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "Total running time  " << duration.count() << " milliseconds." << std::endl;
}