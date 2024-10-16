#include <bits/stdc++.h>
#include "graph.h"
#include "heap.h"
#include "sct.h"
#include <omp.h>
class CuckooHash
{
    using int32 = int32_t;
    const int unfilled = -1;

private:
	/* data */
	int32 capacity;
	int32 mask;
	int32 size;
	int32 buff_size = sizeof(int32);
	int32 *hashtable = nullptr;

	void rehash(int32 **_table) {
		int32 oldcapacity = capacity;
		mask = mask == 0 ? 1 : ((mask << 1) | 1);
		capacity = (mask + 1) * buff_size;
		int32 *newhash = new int32[capacity];
		memset((newhash), unfilled, sizeof(int32) * capacity);
		for (int32 i = 0; i < oldcapacity; ++i){
			if ((*_table)[i] != unfilled) insert((*_table)[i], &newhash);
		}
		std::swap((*_table), newhash);
		delete[] newhash;
	}
	void insert(const int32 &_u, int32 **_table) {
		
		int32 hs = hash1(_u);
		for (int32 i = 0; i < buff_size; ++i) {
			if ((*_table)[hs * buff_size + i] == unfilled){
				(*_table)[hs * buff_size + i] = _u;
				return;
			}
		}
		hs = hash2(_u);
		for (int32 i = 0; i < buff_size; ++i) {
			if ((*_table)[hs * buff_size + i] == unfilled){
				(*_table)[hs * buff_size + i] = _u;
				return;
			}
		}

		bool use_hash1 = true;
		int32 u = _u;
		for (int32 i = 0; i < mask; ++i) {
			int32 replaced;
			if (use_hash1) hs = hash1(u);
			else hs = hash2(u);
			int32 j = 0;
			for (; j < buff_size; ++j) {
				if ((*_table)[hs * buff_size + j] == unfilled) break;
			}
			if (buff_size == j) {
				replaced = std::move((*_table)[hs * buff_size]);
				j = 1;
				for (; j < buff_size; j++) {
					(*_table)[hs * buff_size + j - 1] =
						std::move((*_table)[hs * buff_size + j]);
				}
				(*_table)[hs * buff_size + j - 1] = u;
			}
			else {
				replaced = std::move((*_table)[hs * buff_size + j]);
				(*_table)[hs * buff_size + j] = u;
			}
			use_hash1 = hs == hash2(replaced);
			u = std::move(replaced);
			if (u == unfilled) return;
		}
		rehash(_table);
		insert(u, _table);
	}

	int32 hash1(const int32 x) { return x & mask;}
	int32 hash2(const int32 x) { return ~x & mask;}

public:
	CuckooHash(/* args */) {
		capacity = 0;
		hashtable = nullptr;
		mask = 0;
		size = 0;
	}
	~CuckooHash() {
		if (hashtable != nullptr) {
			delete[] hashtable;
			hashtable = nullptr;
		}
	}

	void reserve(int32 _size) {
		// mask = mask == 0 ? 1 : ((mask << 1) | 1);
		// if(_size < mask * buff_size) mask = 1;
		size = 0;
		while (_size >= mask * buff_size) mask = (mask << 1) | 1;
		int32 newCapacity = (mask + 1) * buff_size;

		if (capacity >= newCapacity) {
			memset(hashtable, unfilled, sizeof(int32) * capacity);
			return;
		}
//printf("heret2 %u %u\n", capacity, newCapacity);fflush(stdout);
		capacity = newCapacity;
		if (hashtable != nullptr) {
			delete [] hashtable;
			hashtable = nullptr;
		}

		hashtable = new int32[capacity];
		memset(hashtable, unfilled, sizeof(int32) * capacity);

	}

	void insert(const int32 &_u) {
		if (find(_u)) return;
		insert(_u, &hashtable);
		size++;
	}

	bool find(const int32 &_u) {

		int32 hs1 = hash1(_u);

		int32 hs2 = hash2(_u);

// assert(buff_size == 4 && sizeof (int32) == 4);
		__m128i cmp = _mm_set1_epi32(_u);
// if(buff_size*hs1 >= capacity) {
// 	printf("hs1 %d, cap %d\n", hs1, capacity);fflush(stdout);
// }
// assert(buff_size*hs1 < capacity);
		__m128i b1 = _mm_load_si128((__m128i*)&hashtable[buff_size * hs1]);
		__m128i b2 = _mm_load_si128((__m128i*)&hashtable[buff_size * hs2]);
		__m128i flag = _mm_or_si128(_mm_cmpeq_epi32(cmp, b1), _mm_cmpeq_epi32(cmp, b2));

		return _mm_movemask_epi8(flag) != 0;
	}
	int32 getcapacity() {return capacity;}
	int32 getmask() {return mask;}
	int32 *gethashtable() {return hashtable;}
};


class LinearSet {
private:
    uint32_t * vSet = nullptr;
    uint32_t * fIndex = nullptr;
    uint32_t sz;

public:
    LinearSet() {}
    LinearSet(uint32_t sz_) {
        resize(sz_);
    }
    void resize(uint32_t sz_) {
        sz = sz_;
        vSet = new uint32_t[sz];
        fIndex = new uint32_t[sz];
        for(uint32_t i = 0; i < sz; i++) {
            vSet[i] = fIndex[i] = i;
        }
    }

    ~LinearSet() {
        if(fIndex != nullptr) {
            delete [] fIndex; delete [] vSet; 
            fIndex = nullptr;
        } 
    }

    uint32_t * begin() {
        return vSet;
    }

    uint32_t operator [] (uint32_t i) {
        // if(i >= g->maxSize()) {
        //     printf("error index\n"); return -1;
        // }
        return vSet[i];
    }

    void changeTo(uint32_t u, uint32_t p) {
        uint32_t pU = fIndex[u];
        std::swap(fIndex[u], fIndex[vSet[p]]);
        std::swap(vSet[pU], vSet[p]);
    }
    uint32_t idx(uint32_t u) {
        return fIndex[u];
    }

    bool isin(uint32_t u, uint32_t st, uint32_t ed) {
        return st <= fIndex[u] && fIndex[u] < ed;
    }
    bool isIn(uint32_t u, uint32_t st, uint32_t ed) {
        return st <= fIndex[u] && fIndex[u] < ed;
    }

    void changeToByPos(uint32_t pU, uint32_t p) {
        std::swap(fIndex[vSet[pU]], fIndex[vSet[p]]);
        std::swap(vSet[pU], vSet[p]);
    }

    void copy(uint32_t * p, uint32_t r, uint32_t st = 0) {
        assert(st <= r);
        assert(r <= sz);

        if(r - st >= 4) {
            memcpy(p, vSet + st, sizeof(uint32_t) * (r - st));
        }
        else {
            for(uint32_t i = st; i < r; i++) {
                p[i - st] = vSet[i];
            }
        }
    }
};

double cnm[5000][5000];
void init_c(){
    unsigned i, j;
    for(i = 0; i < 5000; i++){
        cnm[i][0] = 1;
        for(j = 1; j <= i; j++){
            cnm[i][j] = cnm[i - 1][j - 1] + cnm[i - 1][j];
        }
    }
}
double combiner(unsigned n, unsigned m){
    return cnm[n][m];
}
unsigned Hsize = 0, Psize = 0;
LinearSet C;
void pivoter(SCT *T, unsigned sz, unsigned k,  unsigned *id, bool **edges, unsigned long long fa) {
    if(sz + Hsize + Psize < k) return;
    if(Hsize == k - 1){
        for(unsigned i = 0; i < sz; i++){
            T->id.push_back(C[i]);
            T->type.push_back(0);
            T->fa.push_back(i == 0 ? fa : T->size - 1);
            T->edges.push_back(std::make_pair(T->fa[T->size], T->size));
            T->size++;
        }
        return;
    }
    if(sz == 0) return;

    unsigned pivot = C[0], pivotDeg = 0;
    for(unsigned i = 0; i < sz; i++) {
        unsigned d = 0;
        for(unsigned j = 0; j < sz; j++) {
            if(edges[id[C[i]]][id[C[j]]]) d++;
        }
        if(d > pivotDeg) {
            pivot = C[i]; pivotDeg = d; 
        }
    }

    C.changeTo(pivot, --sz);

    for(unsigned i = 0, j = 0; i < sz; i++) {
        if(edges[id[pivot]][id[C[i]]]) C.changeToByPos(i, j++);
    }

    unsigned candSize = sz - pivotDeg;
    unsigned *cands = new unsigned[candSize];
    for(unsigned i = pivotDeg; i < sz; i++) cands[i-pivotDeg] = C[i];
    
    T->id.push_back(pivot);
    T->type.push_back(0);
    T->fa.push_back(fa);
    T->edges.push_back(std::make_pair(T->fa[T->size], T->size));
    Psize += 1;
    pivoter(T, pivotDeg, k, id, edges, T->size++);
    Psize -= 1;

    for(unsigned i = 0; i < candSize; i++) {
        unsigned u = cands[i];
        C.changeTo(u, --sz);
        unsigned newSz = 0;
        for(unsigned j = 0; j < sz; j++) {
            if(edges[id[u]][id[C[j]]]) C.changeToByPos(j, newSz++);
        }
        T->id.push_back(u);
        T->type.push_back(1);
        T->fa.push_back(fa);
        T->edges.push_back(std::make_pair(T->fa[T->size], T->size));
        Hsize += 1;
        pivoter(T, newSz, k, id, edges, T->size++);
        Hsize -= 1;
    }
    delete[] cands;
}
void pivoter_mf(SCT *T, unsigned sz, unsigned k,  unsigned *id, bool **edges, unsigned long long fa) {
    if(sz + Hsize + Psize < k) return;
    if(Hsize == k - 1){
        for(unsigned i = 0; i < sz; i++){
            T->size++;
        }
        return;
    }

    if(sz == 0) return;

    unsigned pivot = C[0], pivotDeg = 0;
    for(unsigned i = 0; i < sz; i++) {
        unsigned d = 0;
        for(unsigned j = 0; j < sz; j++) {
            if(edges[id[C[i]]][id[C[j]]]) d++;
        }
        if(d > pivotDeg) {
            pivot = C[i]; pivotDeg = d; 
        }
    }

    C.changeTo(pivot, --sz);

    for(unsigned i = 0, j = 0; i < sz; i++) {
        if(edges[id[pivot]][id[C[i]]]) C.changeToByPos(i, j++);
    }

    unsigned candSize = sz - pivotDeg;
    unsigned *cands = new unsigned[candSize];
    for(unsigned i = pivotDeg; i < sz; i++) cands[i-pivotDeg] = C[i];
    
    Psize += 1;
    pivoter_mf(T, pivotDeg, k, id, edges, T->size++);
    Psize -= 1;

    for(unsigned i = 0; i < candSize; i++) {
        unsigned u = cands[i];
        C.changeTo(u, --sz);
        unsigned newSz = 0;
        for(unsigned j = 0; j < sz; j++) {
            if(edges[id[u]][id[C[j]]]) C.changeToByPos(j, newSz++);
        }
        Hsize += 1;
        pivoter_mf(T, newSz, k, id, edges, T->size++);
        Hsize -= 1;
        std::cout<<T->size<<std::endl;
    }
    delete[] cands;
}
unsigned *Hnode, *Pnode;
void iter_kcca(unsigned k, double *r){
    std::pair<double, int> *hvertex = new std::pair<double, int>[Hsize];
    std::pair<double, int> *pvertex = new std::pair<double, int>[Psize];
    unsigned tmpp = 0, tmph = 0;
    if(Hsize != k) for(unsigned i = 0; i < Psize; i++) pvertex[tmpp++] = std::make_pair(r[Pnode[i]], Pnode[i]);
    for(unsigned i = 0; i < Hsize; i++) hvertex[tmph++] = std::make_pair(r[Hnode[i]], Hnode[i]);
    std::sort(hvertex, hvertex + tmph);
    std::sort(pvertex, pvertex + tmpp);
    unsigned left = tmpp;
    if(k - Hsize != 0){
        for(unsigned  i = 0; i < tmpp; i++){
            left--;
            if(r[pvertex[i].second] < r[hvertex[0].second] ){
                r[pvertex[i].second] += cnm[left][k - Hsize - 1];
            }
            else{
                r[hvertex[0].second] += cnm[left][k - Hsize - 1];
            }
        }
    }
    r[hvertex[0].second] += cnm[left][k - Hsize];
    delete[] hvertex;
    delete[] pvertex;
}
void iter_pava(unsigned k, double *r, double *rr){
    std::pair<double, int> *hvertex = new std::pair<double, int>[Hsize];
    std::pair<double, int> *pvertex = new std::pair<double, int>[Psize];
    unsigned tmpp = 0, tmph = 0;
    if(Hsize != k) for(unsigned i = 0; i < Psize; i++) pvertex[tmpp++] = std::make_pair(r[Pnode[i]], Pnode[i]);
    for(unsigned i = 0; i < Hsize; i++) hvertex[tmph++] = std::make_pair(r[Hnode[i]], Hnode[i]);
    std::sort(hvertex, hvertex + tmph);
    std::sort(pvertex, pvertex + tmpp);
    unsigned left = tmpp;
    if(k - Hsize != 0){
        for(unsigned  i = 0; i < tmpp; i++){
            left--;
            if(r[pvertex[i].second] < r[hvertex[0].second] ){
                rr[pvertex[i].second] += cnm[left][k - Hsize - 1];
            }
            else{
                rr[hvertex[0].second] += cnm[left][k - Hsize - 1];
            }
        }
    }
    rr[hvertex[0].second] += cnm[left][k - Hsize];
    delete[] hvertex;
    delete[] pvertex;
}
void iter_psctl(unsigned k, double *rr){
    std::pair<double, int> *hvertex = new std::pair<double, int>[Hsize];
    std::pair<double, int> *pvertex = new std::pair<double, int>[Psize];
    double *upp = new double[Psize];
    unsigned tmpp = 0, tmph = 0;
    if(Hsize != k) for(unsigned i = 0; i < Psize; i++) pvertex[tmpp] = std::make_pair(rr[Pnode[i]], Pnode[i]), upp[tmpp++] = cnm[Psize - 1][k - Hsize - 1];
    for(unsigned i = 0; i < Hsize; i++) hvertex[tmph++] = std::make_pair(rr[Hnode[i]], Hnode[i]);
    std::sort(hvertex, hvertex + tmph);
    std::sort(pvertex, pvertex + tmpp);
    double up = cnm[tmpp][k - tmph];
    unsigned idxp = 0;
    unsigned idxh = 0;
    unsigned remove = 0;
    while(up > 1e-6){
        double gap = 0;
        if(idxp == tmpp && idxh == tmph){
            gap = 1e100;
        }
        else if(idxp == tmpp){
            idxh += 1;
            gap = idxh == tmph ? 1e100 : rr[hvertex[idxh].second] - rr[hvertex[idxh - 1].second];
        }
        else if(idxh == tmph){
            idxp += 1;
            gap = idxp == tmpp ? 1e100 : rr[pvertex[idxp].second] - rr[pvertex[idxp - 1].second];
        }
        else{
            if(rr[pvertex[idxp].second] < rr[hvertex[idxh].second]){
                idxp += 1;
                gap = ((idxp == tmpp) ? rr[hvertex[idxh].second] : std::min(rr[pvertex[idxp].second], rr[hvertex[idxh].second])) - rr[pvertex[idxp - 1].second];
            }
            else{
                idxh += 1;
                gap = ((idxh == tmph) ? rr[pvertex[idxp].second] : std::min(rr[hvertex[idxh].second], rr[pvertex[idxp].second])) - rr[hvertex[idxh - 1].second];
            }
        }
        for(unsigned j = 0; j < idxp; j++){
            if(upp[j] < -0.5) continue; 
            double tmp = cnm[tmpp - remove - 1][k - Hsize - 1];
            if(upp[j] > tmp){
                upp[j] = tmp;
            }
        }
        for(unsigned j = 0; j < idxp; j++){
            if(upp[j] < -0.5) continue;
            gap = std::min(gap, upp[j]);
        }
        gap = std::min(up / (idxp + idxh - remove), gap);

        for(unsigned j = 0; j < idxp; j++){
            if(upp[j] < -0.5) continue; 
            if(upp[j] < 1e-2){
                upp[j] = -1;
                remove += 1;
            }
            if(upp[j] < -0.5) continue; 
            rr[pvertex[j].second] += gap, upp[j] -= gap;
            up -= gap;
            if(upp[j] < 1e-2){
                upp[j] = -1;
                remove += 1;
            }
        }
        for(unsigned j = 0; j < idxh; j++){ 
            rr[hvertex[j].second] += gap;
            up -= gap;
        }

    }
    delete[] hvertex;
    delete[] pvertex;
    delete[] upp;
}
void iter_mf(unsigned sz, unsigned k, unsigned *id, bool **edges, double *r, bool flag) {
    if(sz + Hsize + Psize < k) return;
    if(Hsize == k - 1){
        for(unsigned i = 0; i < sz; i++){
            Pnode[Psize++] = C[i];
        }
        Psize -= sz;
        flag ? iter_kcca(k, r) : iter_psctl(k, r);
        return;
    }

    if(sz == 0){
        flag ? iter_kcca(k, r) : iter_psctl(k, r);
        return;
    }

    unsigned pivot = C[0], pivotDeg = 0;
    for(unsigned i = 0; i < sz; i++) {
        unsigned d = 0;
        for(unsigned j = 0; j < sz; j++) {
            if(edges[id[C[i]]][id[C[j]]]) d++;
        }
        if(d > pivotDeg) {
            pivot = C[i]; pivotDeg = d; 
        }
    }

    C.changeTo(pivot, --sz);

    for(unsigned i = 0, j = 0; i < sz; i++) {
        if(edges[id[pivot]][id[C[i]]]) C.changeToByPos(i, j++);
    }

    unsigned candSize = sz - pivotDeg;
    unsigned *cands = new unsigned[candSize];
    for(unsigned i = pivotDeg; i < sz; i++) cands[i-pivotDeg] = C[i];
    
    Pnode[Psize++] = pivot;
    iter_mf(pivotDeg, k, id, edges, r, flag);
    Psize--;

    for(unsigned i = 0; i < candSize; i++) {
        unsigned u = cands[i];
        C.changeTo(u, --sz);
        unsigned newSz = 0;
        for(unsigned j = 0; j < sz; j++) {
            if(edges[id[u]][id[C[j]]]) C.changeToByPos(j, newSz++);
        }
        Hnode[Hsize++] = u;
        iter_mf(newSz, k, id, edges, r, flag);
        Hsize--;
    }
    delete[] cands;
}
void pava_mf(unsigned sz, unsigned k, unsigned *id, bool **edges, double *r, double *rr) {
    if(sz + Hsize + Psize < k) return;
    if(Hsize == k - 1){
        for(unsigned i = 0; i < sz; i++){
            Pnode[Psize++] = C[i];
        }
        Psize -= sz;
        iter_pava(k, r, rr);
        return;
    }

    if(sz == 0){
        iter_pava(k, r, rr);
        return;
    }

    unsigned pivot = C[0], pivotDeg = 0;
    for(unsigned i = 0; i < sz; i++) {
        unsigned d = 0;
        for(unsigned j = 0; j < sz; j++) {
            if(edges[id[C[i]]][id[C[j]]]) d++;
        }
        if(d > pivotDeg) {
            pivot = C[i]; pivotDeg = d; 
        }
    }

    C.changeTo(pivot, --sz);

    for(unsigned i = 0, j = 0; i < sz; i++) {
        if(edges[id[pivot]][id[C[i]]]) C.changeToByPos(i, j++);
    }

    unsigned candSize = sz - pivotDeg;
    unsigned *cands = new unsigned[candSize];
    for(unsigned i = pivotDeg; i < sz; i++) cands[i-pivotDeg] = C[i];
    
    Pnode[Psize++] = pivot;
    pava_mf(pivotDeg, k, id, edges, r, rr);
    Psize--;

    for(unsigned i = 0; i < candSize; i++) {
        unsigned u = cands[i];
        C.changeTo(u, --sz);
        unsigned newSz = 0;
        for(unsigned j = 0; j < sz; j++) {
            if(edges[id[u]][id[C[j]]]) C.changeToByPos(j, newSz++);
        }
        Hnode[Hsize++] = u;
        pava_mf(newSz, k, id, edges, r, rr);
        Hsize--;
    }
    delete[] cands;
}
unsigned long long npr = 1LL << 62;
SCT *build_sct(graph *g, int k){
    SCT *T = new SCT;
    T->size = 0;
    init_c();
    C.resize(g->n);
    std::vector<CuckooHash> mp;
    mp.resize(g->n);
    unsigned *deg = new unsigned[g->n]();
    for(unsigned i = 0; i < g->m; i++){
        deg[g->edges[i].s]++;
        deg[g->edges[i].t]++;
    }
    for(unsigned u = 0; u < g->n; u++) mp[u].reserve(deg[u]);
     for(unsigned i = 0; i < g->m; i++){
        mp[g->edges[i].s].insert(g->edges[i].t);
        mp[g->edges[i].t].insert(g->edges[i].s);
    }
    unsigned *id = new unsigned[g->n]();
    T->size = 0;
    for(unsigned u = 0; u < g->n; u++) {
        unsigned sz = 0;
        for(unsigned i = g->sadj[u]; i < g->sadj[u + 1]; i++) {
            unsigned v = g->adj[i];
            id[v] = sz;
            C.changeTo(v, sz++);
        }
        T->id.push_back(u);
        T->type.push_back(1);
        T->fa.push_back(npr);
        T->edges.push_back(std::make_pair(0, 0));
        Hsize = 1, Psize = 0;
        bool **edges = new bool*[sz];
        for(unsigned i = 0; i < sz; i++) edges[i] = new bool[sz];
        for(unsigned i = 0; i < sz; i++) for(unsigned j = 0; j < sz; j++) edges[i][j] = mp[C[i]].find(C[j]);
        pivoter(T, sz, k, id, edges, T->size++);
        for(unsigned i = 0; i < sz; i++) delete[] edges[i];
        delete[] edges;
    }
    std::cout<<T->size<<std::endl;
    T->sadj = new unsigned long long[T->size + 1]();
    T->adj = new unsigned long long[T->size]();
    T->srev = new unsigned long long[g->n + 1]();
    T->rev = new unsigned long long[T->size]();
    for(unsigned long long i = 0; i < T->size; i++){
        if(T->edges[i].first != 0 || T->edges[i].second != 0){
            T->sadj[T->edges[i].first + 1]++;
        }
        T->srev[T->id[i] + 1]++; 
    }
    for(unsigned long long i = 0; i < T->size; i++) T->sadj[i + 1] += T->sadj[i];
    for(unsigned i = 0; i < g->n; i++) T->srev[i + 1] += T->srev[i];
    for(unsigned i = g->n; i >= 1; i--) T->srev[i] = T->srev[i - 1];
    for(unsigned long long i = T->size; i >= 1; i--) T->sadj[i] = T->sadj[i - 1];
    for(unsigned long long i = 0; i < T->size; i++){
        if(T->edges[i].first != 0 || T->edges[i].second != 0){
            T->adj[T->sadj[T->edges[i].first + 1]++] = T->edges[i].second;
        }
        T->rev[T->srev[T->id[i] + 1]++] = i; 
    }

    return T;
}
SCT *build_mf(graph *g, int k){
    SCT *T = new SCT;
    T->size = 0;
    T->id.clear();
    T->type.clear();
    T->edges.clear();
    T->fa.clear();
    init_c();
    C.resize(g->n);
    std::vector<CuckooHash> mp;
    mp.resize(g->n);
    unsigned *deg = new unsigned[g->n]();
    for(unsigned i = 0; i < g->m; i++){
        deg[g->edges[i].s]++;
        deg[g->edges[i].t]++;
    }
    for(unsigned u = 0; u < g->n; u++) mp[u].reserve(g->sadj[u + 1] - g->sadj[u]);
     for(unsigned i = 0; i < g->m; i++){
        mp[g->edges[i].s].insert(g->edges[i].t);
        mp[g->edges[i].t].insert(g->edges[i].s);
    }
    unsigned *id = new unsigned[g->n];
    T->size = 0;
    for(unsigned u = 0; u < g->n; u++) {
        unsigned sz = 0;
        for(unsigned i = g->sadj[u]; i < g->sadj[u + 1]; i++) {
            unsigned v = g->adj[i];
            id[v] = sz;
            C.changeTo(v, sz++);
        }
        Hsize = 1;
        bool **edges = new bool*[sz];
        for(unsigned i = 0; i < sz; i++) edges[i] = new bool[sz];
        for(unsigned i = 0; i < sz; i++) for(unsigned j = 0; j < sz; j++) edges[i][j] = mp[C[i]].find(C[j]);
        pivoter_mf(T, sz, k, id, edges, T->size++);
        for(unsigned i = 0; i < sz; i++) delete[] edges[i];
        delete[] edges;
    }
    std::cout<<T->size<<std::endl;
    return T;
}
void build_iter(graph *g, int k, double *r, bool flag){
    init_c();
    Hnode = new unsigned[g->n];
    Pnode = new unsigned[g->n];
    C.resize(g->n);
    std::vector<CuckooHash> mp;
    mp.resize(g->n);
    unsigned *deg = new unsigned[g->n]();
    for(unsigned i = 0; i < g->m; i++){
        deg[g->edges[i].s]++;
        deg[g->edges[i].t]++;
    }
    for(unsigned u = 0; u < g->n; u++) mp[u].reserve(g->sadj[u + 1] - g->sadj[u]);
     for(unsigned i = 0; i < g->m; i++){
        mp[g->edges[i].s].insert(g->edges[i].t);
        mp[g->edges[i].t].insert(g->edges[i].s);
    }
    unsigned *id = new unsigned[g->n];
    for(unsigned u = 0; u < g->n; u++) {
        unsigned sz = 0;
        for(unsigned i = g->sadj[u]; i < g->sadj[u + 1]; i++) {
            unsigned v = g->adj[i];
            id[v] = sz;
            C.changeTo(v, sz++);
        }
        Psize = Hsize = 0;
        Hnode[Hsize++] = u;
        bool **edges = new bool*[sz];
        for(unsigned i = 0; i < sz; i++) edges[i] = new bool[sz];
        for(unsigned i = 0; i < sz; i++) for(unsigned j = 0; j < sz; j++) edges[i][j] = mp[C[i]].find(C[j]);
        iter_mf(sz, k, id, edges, r, flag);
        for(unsigned i = 0; i < sz; i++) delete[] edges[i];
        delete[] edges;
    }
}
void build_pava(graph *g, int k, double *r, double *rr){
    init_c();
    Hnode = new unsigned[g->n];
    Pnode = new unsigned[g->n];
    C.resize(g->n);
    std::vector<CuckooHash> mp;
    mp.resize(g->n);
    unsigned *deg = new unsigned[g->n]();
    for(unsigned i = 0; i < g->m; i++){
        deg[g->edges[i].s]++;
        deg[g->edges[i].t]++;
    }
    for(unsigned u = 0; u < g->n; u++) mp[u].reserve(g->sadj[u + 1] - g->sadj[u]);
     for(unsigned i = 0; i < g->m; i++){
        mp[g->edges[i].s].insert(g->edges[i].t);
        mp[g->edges[i].t].insert(g->edges[i].s);
    }
    unsigned *id = new unsigned[g->n];
    for(unsigned u = 0; u < g->n; u++) {
        unsigned sz = 0;
        for(unsigned i = g->sadj[u]; i < g->sadj[u + 1]; i++) {
            unsigned v = g->adj[i];
            id[v] = sz;
            C.changeTo(v, sz++);
        }
        Psize = Hsize = 0;
        Hnode[Hsize++] = u;
        bool **edges = new bool*[sz];
        for(unsigned i = 0; i < sz; i++) edges[i] = new bool[sz];
        for(unsigned i = 0; i < sz; i++) for(unsigned j = 0; j < sz; j++) edges[i][j] = mp[C[i]].find(C[j]);
        pava_mf(sz, k, id, edges, r, rr);
        for(unsigned i = 0; i < sz; i++) delete[] edges[i];
        delete[] edges;
    }
}
bool cmpT(std::pair<double, unsigned> a, std::pair<double, unsigned> b){
    return a.first == b.first ? a.second > b.second : a.first < b.first;
}
unsigned *pnode, *hnode, pnum, hnum;
void dfs_init(SCT *T, unsigned long long u, unsigned k, double *r, bool *del, long double &sum){
    del[u] = 1;
    if(T->type[u] == 0) pnode[pnum++] = T->id[u];
    if(T->type[u] == 1) hnode[hnum++] = T->id[u];
    if(T->sadj[u] == T->sadj[u + 1]){
        if(hnum != k){
            for(unsigned i = 0; i < pnum; i++) r[pnode[i]] += cnm[pnum - 1][k - hnum - 1];
        }
        for(unsigned i = 0; i < hnum; i++) r[hnode[i]] += cnm[pnum][k - hnum];
        sum += cnm[pnum][k - hnum];
    }
    else{
        for(unsigned long long i = T->sadj[u]; i < T->sadj[u + 1]; i++){
            dfs_init(T, T->adj[i], k, r, del, sum);
        }
    }
    if(T->type[u] == 0) pnum--;
    if(T->type[u] == 1) hnum--;
}
void dfs_clear(SCT *T, unsigned long long u, bool *del){
    if(del[u] == 1 && T->type[u] == 1) return;
    del[u] = 1;
    for(unsigned long long i = T->sadj[u]; i < T->sadj[u + 1]; i++){
        dfs_clear(T, T->adj[i], del);
    }
}
bool dfs_up(SCT *T, unsigned long long u, unsigned k, double *r, bool *del,long double &sum, unsigned *id, unsigned &top, bool *moved, unsigned long long escape){
    if(del[u]){
        if(T->type[u] == 1) return false; 
    }
    else if(u != escape){
        if(T->type[u] == 0) pnode[pnum++] = T->id[u];
        if(T->type[u] == 1) hnode[hnum++] = T->id[u];
    }
    bool flag = 0;
    if(T->sadj[u] == T->sadj[u + 1]){
        if(hnum < k){
            for(unsigned i = 0; i < pnum; i++){
                r[pnode[i]] -= cnm[pnum - 1][k - hnum - 1];
                if(!moved[pnode[i]]){
                    moved[pnode[i]] = 1;
                    id[top++] = pnode[i];
                }
            }
        }
        if(hnum <= k){
            for(unsigned i = 1; i < hnum; i++){
                r[hnode[i]] -= cnm[pnum][k - hnum];
                if(!moved[hnode[i]]){
                    moved[hnode[i]] = 1;
                    id[top++] = hnode[i];
                }
            }
            sum -= cnm[pnum][k - hnum];
            flag = 1;
        }
    }
    else{
        for(unsigned long long i = T->sadj[u]; i < T->sadj[u + 1]; i++){
            flag |= dfs_up(T, T->adj[i], k, r, del, sum, id, top, moved, escape);
        }
    }
    if(!del[u] && u != escape){
        if(T->type[u] == 0) pnum--;
        if(T->type[u] == 1) hnum--;
    }
    if(flag == 0){
        del[u] = 1;
        return false;
    }
    return flag;
    
}
void up(SCT *T, unsigned long long u, unsigned k, double *r, bool *del, long double &sum, unsigned *id, unsigned &top, bool *moved){
    if(del[u]) return;
    pnum = 0;
    hnum = 1;
    unsigned long long v = u;
    while(T->fa[v] != npr){
        v = T->fa[v];
        if(del[v]) continue;
        if(T->type[v] == 0) pnode[pnum++] = T->id[v]; 
        if(T->type[v] == 1) hnode[hnum++] = T->id[v]; 
    }
    dfs_up(T, u, k, r, del, sum, id, top, moved, u);
    if(T->type[u] == 1) dfs_clear(T, u, del);
    del[u] = 1;
}
void dfs_r(SCT *T, unsigned long long u, unsigned k, double *rr){
    if(T->type[u] == 0) pnode[pnum++] = T->id[u];
    if(T->type[u] == 1) hnode[hnum++] = T->id[u];
    if(T->sadj[u] == T->sadj[u + 1]){
        std::pair<double, int> *hvertex = new std::pair<double, int>[hnum];
        std::pair<double, int> *pvertex = new std::pair<double, int>[pnum];
        double *upp = new double[pnum];
        unsigned tmpp = 0, tmph = 0;
        if(hnum != k) for(unsigned i = 0; i < pnum; i++) pvertex[tmpp] = std::make_pair(rr[pnode[i]], pnode[i]), upp[tmpp++] = cnm[pnum - 1][k - hnum - 1];
        for(unsigned i = 0; i < hnum; i++) hvertex[tmph++] = std::make_pair(rr[hnode[i]], hnode[i]);
        std::sort(hvertex, hvertex + tmph);
        std::sort(pvertex, pvertex + tmpp);
        double up = cnm[tmpp][k - tmph];
        unsigned idxp = 0;
        unsigned idxh = 0;
        unsigned remove = 0;
        while(up > 1e-6){
            double gap = 0;
            if(idxp == tmpp && idxh == tmph){
                gap = 1e100;
            }
            else if(idxp == tmpp){
                idxh += 1;
                gap = idxh == tmph ? 1e100 : rr[hvertex[idxh].second] - rr[hvertex[idxh - 1].second];
            }
            else if(idxh == tmph){
                idxp += 1;
                gap = idxp == tmpp ? 1e100 : rr[pvertex[idxp].second] - rr[pvertex[idxp - 1].second];
            }
            else{
                if(rr[pvertex[idxp].second] < rr[hvertex[idxh].second]){
                    idxp += 1;
                    gap = ((idxp == tmpp) ? rr[hvertex[idxh].second] : std::min(rr[pvertex[idxp].second], rr[hvertex[idxh].second])) - rr[pvertex[idxp - 1].second];
                }
                else{
                    idxh += 1;
                    gap = ((idxh == tmph) ? rr[pvertex[idxp].second] : std::min(rr[hvertex[idxh].second], rr[pvertex[idxp].second])) - rr[hvertex[idxh - 1].second];
                }
            }
            for(unsigned j = 0; j < idxp; j++){
                if(upp[j] < -0.5) continue; 
                double tmp = cnm[tmpp - remove - 1][k - hnum - 1];
                if(upp[j] > tmp){
                    upp[j] = tmp;
                }
            }
            for(unsigned j = 0; j < idxp; j++){
                if(upp[j] < -0.5) continue;
                gap = std::min(gap, upp[j]);
            }
            gap = std::min(up / (idxp + idxh - remove), gap);

            for(unsigned j = 0; j < idxp; j++){
                if(upp[j] < -0.5) continue; 
                if(upp[j] < 1e-2){
                    upp[j] = -1;
                    remove += 1;
                }
                if(upp[j] < -0.5) continue; 
                rr[pvertex[j].second] += gap, upp[j] -= gap;
                up -= gap;
                if(upp[j] < 1e-2){
                    upp[j] = -1;
                    remove += 1;
                }
            }
            for(unsigned j = 0; j < idxh; j++){ 
                rr[hvertex[j].second] += gap;
                up -= gap;
            }

        }
        delete[] hvertex;
        delete[] pvertex;
        delete[] upp;
    }
    else{
        for(unsigned long long i = T->sadj[u]; i < T->sadj[u + 1]; i++){
            dfs_r(T, T->adj[i], k, rr);
        }
    }
    if(T->type[u] == 0) pnum--;
    if(T->type[u] == 1) hnum--;
}
void dfs_rr(SCT *T, unsigned long long u, unsigned k, double *rr, bool *del, unsigned long long escape){
    if(del[u] == 0 && u != escape){
        if(T->type[u] == 0) pnode[pnum++] = T->id[u];
        if(T->type[u] == 1) hnode[hnum++] = T->id[u];
    }
    if(T->sadj[u] == T->sadj[u + 1]){
        if(hnum > k) return;
        if(hnum <= k){
            std::pair<double, int> *hvertex = new std::pair<double, int>[hnum];
            std::pair<double, int> *pvertex = new std::pair<double, int>[pnum];
            double *upp = new double[pnum];
            unsigned tmpp = 0, tmph = 0;
            if(hnum != k) for(unsigned i = 0; i < pnum; i++) pvertex[tmpp] = std::make_pair(rr[pnode[i]], pnode[i]), upp[tmpp++] = cnm[pnum - 1][k - hnum - 1];
            for(unsigned i = 0; i < hnum; i++) hvertex[tmph++] = std::make_pair(rr[hnode[i]], hnode[i]);
            std::sort(hvertex, hvertex + tmph);
            std::sort(pvertex, pvertex + tmpp);
            double up = cnm[tmpp][k - tmph];
            unsigned idxp = 0;
            unsigned idxh = 0;
            unsigned remove = 0;
            while(up > 1e-6){
                double gap = 0;
                if(idxp == tmpp && idxh == tmph){
                    gap = 1e100;
                }
                else if(idxp == tmpp){
                    idxh += 1;
                    gap = idxh == tmph ? 1e100 : rr[hvertex[idxh].second] - rr[hvertex[idxh - 1].second];
                }
                else if(idxh == tmph){
                    idxp += 1;
                    gap = idxp == tmpp ? 1e100 : rr[pvertex[idxp].second] - rr[pvertex[idxp - 1].second];
                }
                else{
                    if(rr[pvertex[idxp].second] < rr[hvertex[idxh].second]){
                        idxp += 1;
                        gap = ((idxp == tmpp) ? rr[hvertex[idxh].second] : std::min(rr[pvertex[idxp].second], rr[hvertex[idxh].second])) - rr[pvertex[idxp - 1].second];
                    }
                    else{
                        idxh += 1;
                        gap = ((idxh == tmph) ? rr[pvertex[idxp].second] : std::min(rr[hvertex[idxh].second], rr[pvertex[idxp].second])) - rr[hvertex[idxh - 1].second];
                    }
                }
                for(unsigned j = 0; j < idxp; j++){
                    if(upp[j] < -0.5) continue; 
                    double tmp = cnm[tmpp - remove - 1][k - hnum - 1];
                    if(upp[j] > tmp){
                        upp[j] = tmp;
                    }
                }
                for(unsigned j = 0; j < idxp; j++){
                    if(upp[j] < -0.5) continue;
                    gap = std::min(gap, upp[j]);
                }
                gap = std::min(up / (idxp + idxh - remove), gap);

                for(unsigned j = 0; j < idxp; j++){
                    if(upp[j] < -0.5) continue; 
                    if(upp[j] < 1e-2){
                        upp[j] = -1;
                        remove += 1;
                    }
                    if(upp[j] < -0.5) continue; 
                    rr[pvertex[j].second] += gap, upp[j] -= gap;
                    up -= gap;
                    if(upp[j] < 1e-2){
                        upp[j] = -1;
                        remove += 1;
                    }
                }
                for(unsigned j = 0; j < idxh; j++){ 
                    rr[hvertex[j].second] += gap;
                    up -= gap;
                }

            }
            delete[] hvertex;
            delete[] pvertex;
            delete[] upp;
        } 
    }
    else{
        for(unsigned long long i = T->sadj[u]; i < T->sadj[u + 1]; i++){
            dfs_rr(T, T->adj[i], k, rr, del, escape);
        }
    }
    if(del[u] == 0 && u != escape){
        if(T->type[u] == 0) pnum--;
        if(T->type[u] == 1) hnum--;
    }
}
void get_rr(SCT *T, unsigned long long u, unsigned k, double *r, bool *del){
    pnum = 0;
    hnum = 0;
    del[u] = 0;
    hnode[hnum++] = T->id[u];
    unsigned long long v = u;
    while(T->fa[v] != npr){
        v = T->fa[v];
        if(del[v]) continue;
        if(T->type[v] == 0) pnode[pnum++] = T->id[v]; 
        if(T->type[v] == 1) hnode[hnum++] = T->id[v]; 
    }
    dfs_rr(T, u, k, r, del, u);
}
void update(graph *g, SCT *T, unsigned k, double *r, double *rr, double &rho, unsigned &sz, bool get_eps){
    long double sum = 0;
    bool *del = new bool[T->size]();
    bool *moved = new bool[g->n]();
    unsigned *id = new unsigned[g->n]();
    unsigned top = 0;
    hnode = new unsigned[g->n]();
    pnode = new unsigned[g->n]();
    pnum = hnum = 0;
    for(unsigned long long i = 0; i < T->size; i++){
        if(!del[i]) dfs_init(T, i, k, r, del, sum);
    }
    for(unsigned long long i = 0; i < T->size; i++) del[i] = 0;
    heap *H = init(g->n, r);
    if(sum / g->n > rho){
        rho = sum / g->n;
        sz = g->n;
    }
    while(H->size){
        std::pair<unsigned, double> tmp = extract_min(H);
        for(unsigned long long s = T->srev[tmp.first]; s < T->srev[tmp.first + 1]; s++){
            unsigned long long i = T->rev[s];
            up(T, i, k, r, del, sum, id, top, moved);
        }
        for(unsigned i = 0; i < top; i++){
            moved[id[i]] = 0;
            decrease_key(H, id[i], r[id[i]]);
        }
        top = 0;
        if(H->size != 0){
            if(sum / H->size > rho){
                rho = sum / H->size;
                sz = H->size;
            }
        }
    }
    double tmp_sum  = 0;
    for(unsigned i = 0; i < g->n; i++) tmp_sum += r[i];
    destroy(H);
    if(get_eps){
        for(long long i = T->size - 1; i >= 0; i--){
            if(T->fa[i] == npr){
                hnum = pnum = 0;
                dfs_r(T, i, k, rr);
            }
        }
        /*
        for(unsigned long long i = 0; i < pos.size(); i++){
            assert(del[pos[i]] == 1);
            get_rr(T, pos[i], k, rr, del);
        }
        */
    }
    delete[] pnode;
    delete[] hnode;
    delete[] del;
    delete[] id;
    delete[] moved;
}

bool *check(graph *g, SCT *T, unsigned k, double *r, double bound){
    long double sum = 0;
    bool *del = new bool[T->size]();
    bool *moved = new bool[g->n]();
    unsigned *id = new unsigned[g->n]();
    unsigned top = 0;
    pnum = hnum = 0;
    hnode = new unsigned[g->n]();
    pnode = new unsigned[g->n]();
    for(unsigned long long i = 0; i < T->size; i++){
        if(!del[i]) dfs_init(T, i, k, r, del, sum);
    }
    for(unsigned long long i = 0; i < T->size; i++) del[i] = 0;
    heap *H = init(g->n, r);
    bool *is = new bool[g->n]();
    for(unsigned i = 0; i < g->n; i++) is[i] = 1;
    while(H->size){
        std::pair<unsigned, double> tmp = extract_min(H);
        if(tmp.second < bound) is[tmp.first] = 0;
        else break;
        for(unsigned long long s = T->srev[tmp.first]; s < T->srev[tmp.first + 1]; s++){
            unsigned long long i = T->rev[s];
            up(T, i, k, r, del, sum, id, top, moved);
        }
        for(unsigned i = 0; i < top; i++){
            moved[id[i]] = 0;
            decrease_key(H, id[i], r[id[i]]);
        }
        top = 0;
    }
    destroy(H);
    delete[] pnode;
    delete[] hnode;
    delete[] moved;
    delete[] id;
    delete[] del;
    return is;
}

