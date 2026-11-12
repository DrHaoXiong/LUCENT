"""Algebraically identical bivariate rectangle implementation performance probe."""
import json
import time
import numpy as np
from scipy.special import ndtr, ndtri, owens_t
from .copula import rectangle, rectangle_reference

def cdf(h, k, r):
    h, k, r = np.broadcast_arrays(h, k, r)
    out = np.empty(h.shape)
    low = np.isneginf(h) | np.isneginf(k)
    out[low] = 0.
    hi_h = np.isposinf(h) & ~low
    out[hi_h] = ndtr(k[hi_h])
    hi_k = np.isposinf(k) & ~low & ~hi_h
    out[hi_k] = ndtr(h[hi_k])
    finite = ~(low | hi_h | hi_k)
    hh, kk, rr = h[finite], k[finite], r[finite]
    z = np.empty(hh.shape)
    sd = np.sqrt(1-rr*rr)
    both_zero = (hh == 0) & (kk == 0)
    z[both_zero] = .25+np.arcsin(rr[both_zero])/(2*np.pi)
    hz = (hh == 0) & ~both_zero
    kz = (kk == 0) & ~both_zero
    z[hz] = .5*ndtr(kk[hz])-owens_t(kk[hz], -rr[hz]/sd[hz])
    z[kz] = .5*ndtr(hh[kz])-owens_t(hh[kz], -rr[kz]/sd[kz])
    other = ~(both_zero | hz | kz)
    a, b, rho, sig = hh[other], kk[other], rr[other], sd[other]
    z[other] = .5*(ndtr(a)+ndtr(b))-owens_t(a, (b/a-rho)/sig)-owens_t(b, (a/b-rho)/sig)-(1-np.sign(a)*np.sign(b))/4
    out[finite] = z
    return out

def fast_rectangle(al, ah, bl, bh, r):
    al, ah, bl, bh, r = np.broadcast_arrays(al, ah, bl, bh, r)
    lower_a, upper_a, lower_b, upper_b = ndtri(al), ndtri(ah), ndtri(bl), ndtri(bh)
    p = cdf(upper_a, upper_b, r)-cdf(lower_a, upper_b, r)-cdf(upper_a, lower_b, r)+cdf(lower_a, lower_b, r)
    fallback = (p < 1e-10) | ~np.isfinite(p)
    if fallback.any():
        p[fallback] = rectangle(al[fallback], ah[fallback], bl[fallback], bh[fallback], r[fallback], 256)
    return p

