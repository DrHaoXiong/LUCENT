"""Observed-only R2 discrete Gaussian-copula pair-likelihood / raw-T MC.

No generator setting, detection probability, histogram or true covariance input.
The fitted marginal and spatial shape are plug-in estimates, not an exact test.
"""
from functools import lru_cache
import hashlib

import numpy as np
from scipy.integrate import quad
from scipy.spatial.distance import cdist
from scipy.special import ndtr, ndtri

GRID_FRACTIONS = (0., .25, .5, .75, .9, .97, .99, .999)
GRID_RANGES = (.5, 1., 2., 4., 8.)

@lru_cache(maxsize=8)
def quadrature(order):
    node, weight = np.polynomial.legendre.leggauss(order)
    return (node+1)/2, weight/2

def normal_interval(a, b):
    """Stable Phi(b)-Phi(a) without an arbitrary probability floor."""
    a, b = np.broadcast_arrays(a, b)
    return np.where(a > 0, ndtr(-a)-ndtr(-b), ndtr(b)-ndtr(a))

def rectangle(pa_lo, pa_hi, pb_lo, pb_hi, rho, order=128):
    """Bivariate normal rectangle via deterministic 1-D probability quadrature.

    Integrate over the rarer marginal interval; infinities handled by ndtri.
    Independence is exact. Never clip a rectangle probability to a fitted floor.
    """
    al, ah, bl, bh, r = np.broadcast_arrays(pa_lo, pa_hi, pb_lo, pb_hi, rho)
    shape = r.shape
    al, ah, bl, bh, r = [x.ravel() for x in (al, ah, bl, bh, r)]
    switch = ah-al > bh-bl
    al, ah, bl, bh = (np.where(switch, bl, al), np.where(switch, bh, ah),
                      np.where(switch, al, bl), np.where(switch, ah, bh))
    node, weight = quadrature(order)
    latent = ndtri(al[:, None]+(ah-al)[:, None]*node[None, :])
    sigma = np.sqrt(1-r*r)
    lo = (ndtri(bl)[:, None]-r[:, None]*latent)/sigma[:, None]
    hi = (ndtri(bh)[:, None]-r[:, None]*latent)/sigma[:, None]
    prob = (ah-al)*(normal_interval(lo, hi)@weight)
    prob[r == 0] = (ah-al)[r == 0]*(bh-bl)[r == 0]
    if not np.isfinite(prob).all() or (prob <= 0).any():
        raise ValueError('nonpositive/nonfinite rectangle integral; no probability floor')
    return prob.reshape(shape)

def rectangle_reference(al, ah, bl, bh, rho):
    """Independent adaptive scipy quadrature for numerical tests, not fitting."""
    if ah-al > bh-bl:
        al, ah, bl, bh = bl, bh, al, ah
    if rho == 0:
        return (ah-al)*(bh-bl)
    lower, upper = ndtri(bl), ndtri(bh)
    sd = np.sqrt(1-rho*rho)
    def fun(u):
        x = ndtri(u)
        return float(normal_interval((lower-rho*x)/sd, (upper-rho*x)/sd))
    value, error = quad(fun, al, ah, epsabs=1e-12, epsrel=2e-8, limit=250)
    return value

class CountGeometry:
    def __init__(self, coordinates, gaussian_model):
        self.coordinates = np.asarray(coordinates, float)
        self.model = gaussian_model
        self.n = len(coordinates)
        distance = cdist(self.coordinates, self.coordinates)
        i, j = np.triu_indices(self.n, 1)
        order = np.argsort(distance[i, j], kind='stable')
        if len(i) <= 2048:
            chosen = np.arange(len(i))
        else:
            near = order[:1024]
            remaining = order[1024:]
            fingerprint = hashlib.sha256(self.coordinates.astype('<f8').tobytes()).digest()
            gen = np.random.Generator(np.random.PCG64DXSM(int.from_bytes(fingerprint[:16], 'little')))
            chosen = np.r_[near, gen.choice(remaining, size=1024, replace=False)]
        self.pairs = np.column_stack([i[chosen], j[chosen]])
        self.distance = distance[i[chosen], j[chosen]]
        self.rho = np.array([np.zeros(len(chosen)) if s.fraction == 0 else
                             s.fraction*np.exp(-self.distance/(gaussian_model.scale*s.range_multiplier))
                             for s in gaussian_model.shapes])
        assert len({tuple(x) for x in self.pairs}) == len(chosen)

    def fit(self, response):
        y = np.asarray(response, float)
        if y.shape != (self.n,) or not np.isfinite(y).all():
            raise ValueError('invalid complete response pool')
        values, labels, counts = np.unique(y, return_inverse=True, return_counts=True)
        if len(values) < 2:
            raise ValueError('RESPONSE_CONSTANT')
        cumulative = np.r_[0., np.cumsum(counts)/self.n]
        left, right = labels[self.pairs[:, 0]], labels[self.pairs[:, 1]]
        al, ah = cumulative[left], cumulative[left+1]
        bl, bh = cumulative[right], cumulative[right+1]
        from .rectangles import fast_rectangle
        likelihood = np.log(fast_rectangle(al[None, :], ah[None, :], bl[None, :], bh[None, :], self.rho)).sum(axis=1)
        chosen = int(np.argmax(likelihood))
        # Only numerical convergence, not alternative fitted methods.
        best = np.argsort(-likelihood, kind='stable')[:2]
        refined = np.log(rectangle(al[None, :], ah[None, :], bl[None, :], bh[None, :], self.rho[best], 256)).sum(axis=1)
        drift = np.abs(refined-likelihood[best])
        if drift.max() > .02 or int(best[np.argmax(refined)]) != chosen:
            likelihood = np.log(rectangle(al[None, :], ah[None, :], bl[None, :], bh[None, :], self.rho, 256)).sum(axis=1)
            chosen = int(np.argmax(likelihood))
            order_used = 256
        else:
            order_used = 0  # analytic Owen-T; adaptive quadrature only for tiny rectangles
        shape = self.model.shapes[chosen]
        return {'shape': chosen, 'fraction': shape.fraction, 'range_multiplier': shape.range_multiplier,
                'values': values, 'cumulative': cumulative, 'counts': counts,
                'log_likelihood': float(likelihood[chosen]), 'pair_count': len(self.pairs),
                'quadrature_order': order_used, 'quadrature_loglik_drift': float(drift.max())}


