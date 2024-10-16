#ifndef UTILS_PAVA_H_
#define UTILS_PAVA_H_

std::pair<double, double> pava(graph *g, double *r,unsigned k, unsigned t, double rho);
bool cmp_y(rank a, rank b);
std::pair<double, double> pava_iter(graph *g, double *r, double *rr, unsigned k, unsigned t);
#endif