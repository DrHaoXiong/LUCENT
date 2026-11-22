#include <math.h>
#include <stdint.h>
#include <stddef.h>

typedef double (*unary)(double, int);
typedef double (*binary)(double, double, int);
static unary phi;
static binary owen;

void bind_special(unary normal_cdf, binary owens_t) {
    phi = normal_cdf;
    owen = owens_t;
}

static double corner(double h, double k, double ph, double pk, double r, double sd) {
    if (h == -INFINITY || k == -INFINITY) return 0.;
    if (h == INFINITY) return pk;
    if (k == INFINITY) return ph;
    if (h == 0. && k == 0.) return 0.25 + asin(r) / (2.*acos(-1.));
    if (h == 0.) return .5*pk - owen(k, -r/sd, 1);
    if (k == 0.) return .5*ph - owen(h, -r/sd, 1);
    return .5*(ph+pk) - owen(h,(k/h-r)/sd,1) - owen(k,(h/k-r)/sd,1)
           - ((signbit(h) != signbit(k)) ? .5 : 0.);
}
