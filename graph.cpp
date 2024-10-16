#include <bits/stdc++.h>
#include "graph.h"
#include <omp.h>

graph* new_graph() {
    graph *g = new graph();
    g->n = g->m = 0;
    g->edges = nullptr;
    g->adj = nullptr;
    g->sadj = nullptr;

    return g;
}
graph *copy(graph *g){
    graph *gtmp = new_graph();
    gtmp->n = g->n;
    gtmp->m = g->m;
    gtmp->edges = new edge[gtmp->m];
    gtmp->adj = new unsigned[gtmp->m];
    gtmp->sadj = new unsigned[gtmp->n + 1];
    for(unsigned i = 0; i < gtmp->m; i++){
        gtmp->edges[i] = g->edges[i];
        gtmp->adj[i] = g->adj[i];
    }
    for(unsigned i = 0; i < gtmp->n + 1; i++){
        gtmp->sadj[i] = g->sadj[i];
    }
    return gtmp;
}
graph* reserve(graph *g){
    graph *gg = new graph();
    gg->n = g->n;
    gg->m = g->m;
    gg->adj = new unsigned[gg->m]();
    gg->sadj = new unsigned[gg->n + 1]();
    for(unsigned i = 0; i < g->n; i++){
        for(unsigned j = g->sadj[i]; j < g->sadj[i + 1]; j++){
            gg->sadj[g->adj[j]]++;
        }
    }
    for(unsigned i = 1; i <= gg->n; i++){
        gg->sadj[i] += gg->sadj[i - 1];
    }
    for(unsigned i = 0; i < g->n; i++){
        for(unsigned j = g->sadj[i]; j < g->sadj[i + 1]; j++){
            gg->adj[--gg->sadj[g->adj[j]]] = i;
        }
    }
    return gg;
}
graph* shrink(graph *g, bool *is){
    graph *gg = new graph();
    unsigned *id = new unsigned[g->n]();
    for(unsigned i = 0; i < g->n; i++) if(is[i]) id[i] = gg->n++;
    gg->sadj = new unsigned[gg->n + 1]();

    for(unsigned i = 0; i < g->n; i++){
        if(!is[i]) continue;
        for(unsigned j = g->sadj[i]; j < g->sadj[i + 1]; j++){
            if(!is[g->adj[j]]) continue;
            gg->sadj[id[i]]++;
            gg->m++;
        }
    }

    
    for(unsigned i = 0; i < gg->n; i++){
        gg->sadj[i + 1] += gg->sadj[i];
    }
    gg->adj = new unsigned[gg->m];
    gg->edges = new edge[gg->m]();
    unsigned count = 0;
    for(unsigned i = 0; i < g->n; i++){
        if(!is[i]) continue;
        for(unsigned j = g->sadj[i]; j < g->sadj[i + 1]; j++){
            if(!is[g->adj[j]]) continue;
            gg->adj[--gg->sadj[id[i]]] = id[g->adj[j]]; 
            gg->edges[count].s = id[i];
            gg->edges[count].t = id[g->adj[j]];
            count++;
        }
    }
    delete[] id;
    return gg;
}
graph* read_graph(char *fp) {
    graph *g = new_graph();
    FILE *file = fopen(fp,"r");
    unsigned i;

    fscanf(file, "%u %u", &(g->n), &(g->m));
    g->edges = new edge[g->m]; 
    for(i = 0; i < g->m; i++) {
        fscanf(file, "%u %u", &(g->edges[i].s),&(g->edges[i].t));
        //assert(g->edges[i].s < g->n && g->edges[i].s >= 0);
        //assert(g->edges[i].t < g->n && g->edges[i].t >= 0);
    }

    fclose(file);
    return g;
}
bool e_cmp(edge a, edge b){
    return a.s == b.s ? a.t < b.t : a.s < b.s; 
}
void core_order(graph *g, unsigned k) {
    unsigned i, u, tmp;
    int head = 0,tail = -1;

    unsigned *deg = new unsigned[g->n]();
    unsigned *adj = new unsigned[(1LL * g->m) << 1];
    unsigned long long *sadj = new unsigned long long[g->n + 1]();
    unsigned *queue = new unsigned[g->n];
    unsigned *rev = new unsigned[g->n];
    unsigned *count = new unsigned[g->n]();
    char *remove = new char[g->n]();

    unsigned max_deg = 0;
    for(i = 0; i < g->m; i++) {
        deg[g->edges[i].s]++;
        deg[g->edges[i].t]++;
    }
    for(i = 0; i < g->n; i++) {
        sadj[i + 1] = sadj[i] + deg[i];
        max_deg = std::max(max_deg,deg[i]);
        deg[i] = 0;
    }
    for(i = 0; i < g->m; i++) {
        adj[sadj[g->edges[i].s] + deg[g->edges[i].s]++] = g->edges[i].t;
        adj[sadj[g->edges[i].t] + deg[g->edges[i].t]++] = g->edges[i].s;
    }
    //(k - 1) core decomposition
    k -= 2;
    for(int i = 0; i < g->n; i++){
        if(deg[i] <= k) {
            queue[++tail] = i;
        }
    }     
    unsigned long long ii;
    while(head <= tail){
        u = queue[head++];
        remove[u] = 1;
        for(ii = sadj[u]; ii < sadj[u + 1]; ii++) {
            if((--deg[adj[ii]]) == k) {
                queue[++tail] = adj[ii];
            }
        }
    }
    // building core ordered graph
    tmp = 0;
    max_deg = 0;
    for(i = 0; i < g->n; i++) {
        if(remove[i]) continue;
        deg[tmp] = deg[i];
        max_deg = std::max(max_deg,deg[i]);
        count[deg[tmp]]++;
        rev[i] = tmp++;
    }
    g->n = tmp;
    unsigned *rev_new = new unsigned[g->n];    
    for(int i = 1; i < g->n; i++){
        count[i] += count[i - 1];
    }
    max_deg = 0;
    for(int i = 0; i < g->n; i++){
        rev_new[i] = --count[deg[i]];
        deg[i] = 0;
    }
    tmp = 0;
    for(i = 0; i < g->m; i++){
        if(remove[g->edges[i].s] || remove[g->edges[i].t]) continue;
        tmp++;
    }
    edge *edges_new = new edge[tmp];
    tmp = 0;
    for(i = 0; i < g->m; i++){
        if(remove[g->edges[i].s] || remove[g->edges[i].t]) continue;
        g->edges[i].s = rev_new[rev[g->edges[i].s]];
        g->edges[i].t = rev_new[rev[g->edges[i].t]];
        if(g->edges[i].s > g->edges[i].t){
            g->edges[i].t ^= g->edges[i].s;
            g->edges[i].s ^= g->edges[i].t;
            g->edges[i].t ^= g->edges[i].s; 
        }
        edges_new[tmp++] = g->edges[i];
        deg[g->edges[i].s]++;
    }
    g->m = tmp;
    delete[] g->edges;
    g->edges = edges_new;
    std::sort(g->edges, g->edges + g->m, e_cmp);
    g->sadj = new unsigned[g->n + 1]();
    g->adj = new unsigned[g->m];
    max_deg = 0;
    unsigned sum_deg = 0;
    for(i = 0; i < g->n; i++) {
        g->sadj[i + 1] = g->sadj[i] + deg[i];
        max_deg = std::max(max_deg,deg[i]);
        sum_deg += deg[i];
        deg[i] = 0;
    }
    for(i = 0; i < g->m; i++){
        g->adj[g->sadj[g->edges[i].s] + deg[g->edges[i].s]++] = g->edges[i].t;
    }
    
    delete[] deg;
    delete[] adj;
    delete[] sadj;
    delete[] queue;
    delete[] rev;
    delete[] count;
    delete[] remove;
}
void re_order(graph *g){
    int* id = new int[g->n]();
    unsigned* deg = new unsigned[g->n]();
    unsigned* count = new unsigned[g->n]();
    unsigned* p = new unsigned[g->n]();
    unsigned** bucket = new unsigned*[g->n]();
    for(unsigned i = 0; i < g->m; i++){
        deg[g->edges[i].s]++;
        deg[g->edges[i].t]++;
    }
    for(unsigned i = 0; i < g->n; i++) count[deg[i]]++, id[i] = -1;
    for(unsigned i = g->n - 1; i >= 1; i--) count[i - 1] += count[i];
    for(unsigned i = 0; i < g->n; i++) bucket[i] = new unsigned[count[i]]();
    for(unsigned i = 0; i < g->n; i++) bucket[deg[i]][p[deg[i]]++] = i;
    unsigned cnt = 0;
    unsigned max_core = 0;
    std::vector<unsigned> *e = new std::vector<unsigned>[g->n]();
    for(unsigned u = 0; u < g->n; u++){
        for(unsigned k = g->sadj[u]; k < g->sadj[u + 1]; k++){
            unsigned v = g->adj[k];
            e[u].push_back(v), e[v].push_back(u);
        }
    }
    for(unsigned i = 0; i < g->n; i++){
        for(unsigned j = 0; j < p[i]; j++){
            unsigned u = bucket[i][j];
            if(id[u] != -1) continue;
            max_core = std::max(max_core, i);
            id[u] = cnt++;
            deg[u] = 0;
            for(auto v: e[u]){
                if(id[v] == -1){
                    deg[v]--;
                    bucket[deg[v]][p[deg[v]]++] = v;
                }
            }
        }
    }
    std::cout<<"core value: "<<max_core<<std::endl;
    for(unsigned i = 0; i < g->m; i++){
        g->edges[i].s = id[g->edges[i].s];
        g->edges[i].t = id[g->edges[i].t];
        if(g->edges[i].s > g->edges[i].t){
            g->edges[i].t ^= g->edges[i].s;
            g->edges[i].s ^= g->edges[i].t;
            g->edges[i].t ^= g->edges[i].s; 
        }
        deg[g->edges[i].s]++;
    }
    std::sort(g->edges, g->edges + g->m, e_cmp);
    for(unsigned i = 0; i < g->n; i++) {
        g->sadj[i + 1] = g->sadj[i] + deg[i];
        deg[i] = 0;
    }
    for(unsigned i = 0; i < g->m; i++){
        g->adj[g->sadj[g->edges[i].s] + deg[g->edges[i].s]++] = g->edges[i].t;
    }
    delete[] e;
    delete[] id;
    delete[] count;
    delete[] p;
    for(unsigned i = 0; i < g->n; i++){
        delete[] bucket[i];
    }
    delete[] bucket;
}
graph * build_graph(char *fp, int k) {
    printf("Reading graph from file %s\n", fp);
    graph *g = read_graph(fp);
    auto start_time = std::chrono::high_resolution_clock::now();
    printf("Graph size |V| = %u, |E| = %u\n", g->n, g->m);
    printf("Building graph structure\n");
    core_order(g, k);
    printf("Building finished\n");
    printf("Graph size after decompostion |V| = %u, |E| = %u\n", g->n, g->m);
    core_order(g, k);
    std::cout<<g->n<<" "<<g->m<<std::endl;
    return g;
    
    /*
    g_reserved = reserve(g);
    for(unsigned i = 0; i < g_reserved->n; i++) count[i] = 0;
        for(unsigned u = 0; u < g_reserved->n; u++){
        for(unsigned i = g_reserved->sadj[u]; i < g_reserved->sadj[u + 1]; i++) vis[g_reserved->adj[i]] = 1;
        for(unsigned i = g_reserved->sadj[u]; i < g_reserved->sadj[u + 1]; i++){
            unsigned v = g_reserved->adj[i];
            for(unsigned j = g_reserved->sadj[v]; j < g_reserved->sadj[v + 1]; j++) vis[g_reserved->adj[j]] += 2;
            for(unsigned j = g_reserved->sadj[v]; j < g_reserved->sadj[v + 1]; j++){
                unsigned w = g_reserved->adj[j];
                if(!vis[w]) continue;
                for(unsigned k = g_reserved->sadj[w]; k < g_reserved->sadj[w + 1]; k++){
                    unsigned x = g_reserved->adj[k];
                    if(vis[x] != 3) continue;
                    count[u]++,count[v]++,count[w]++,count[x]++;
                }

            }
            for(unsigned j = g_reserved->sadj[v]; j < g_reserved->sadj[v + 1]; j++) vis[g_reserved->adj[j]] -= 2;
        }
        for(unsigned i = g_reserved->sadj[u]; i < g_reserved->sadj[u + 1]; i++) vis[g_reserved->adj[i]] = 0;
    } 
    for(unsigned i = 0; i < g->n; i++){
        is[i] = (count[i] > ccnm[tmp - 1][3]);
    }
    g = shrink(g, is);
    */
    std::cout<<g->n<<" "<<g->m<<std::endl;
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "Running time for building graphs: " << duration.count() << " milliseconds." << std::endl;
    

    
}
bool cmp(rank x, rank y){
    return x.r == y.r ? x.id < y.id : x.r < y.r;
}
void reorder(graph *g, double *r) {
    int i;
    
    rank *rorder = new rank[g->n];
    unsigned *rev = new unsigned[g->n];
    unsigned *deg = new unsigned[g->n]();

    for(i = 0; i < g->n ; i++) {
        rorder[i].id = i;
        rorder[i].r = r[i];
    }

    unsigned max_deg = 0;
    std::sort(rorder, rorder + g->n, cmp);
    for(i = 0; i < g->n ; i++) {
        rev[rorder[i].id] = i;
    }
    for(i = 0; i < g->m; i++) {
        if(rev[g->edges[i].s] > rev[g->edges[i].t]) {
            g->edges[i].s ^= g->edges[i].t;
            g->edges[i].t ^= g->edges[i].s;
            g->edges[i].s ^= g->edges[i].t;
        }
        deg[g->edges[i].s]++;
    }   
    for(i = 0; i < g->n; i++) {
        g->sadj[i + 1] = g->sadj[i] + deg[i];
        max_deg = std::max(max_deg,deg[i]);
        deg[i] = 0;
    }
    for(i = 0; i < g->m; i++) {
        g->adj[g->sadj[g->edges[i].s] + deg[g->edges[i].s]++] = g->edges[i].t;
    }

    printf("max deg = %u\n",max_deg);

    delete[] rorder;
    delete[] rev;
    delete[] deg;
    

}