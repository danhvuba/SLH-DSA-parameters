from sage.all import *
from functools import lru_cache

# ------------------------------------------------
# Cached probability functions
# ------------------------------------------------

# Pr[exactly r sigs hit the leaf targeted by this forgery attempt] 
# = (qs,r) * 1/2^(hr) * (1-1/2^h)^(qs-r) 
@lru_cache(None)
def qhitprob(leaves, qs, r):
    p = 1 / F(leaves)
    return binomial(qs, r) * p^r * (1 - p)^(qs - r)


# Pr[FORS forgery given that exactly r sigs hit the leaf] 
# = (1-(1-1/(2^a))^r)^k = (1-(1-1/t)^r)^k
@lru_cache(None)
def forgeryprob(a, r, k):
    base = 1 - (1 - 1 / F(2^a))^r
    return base^k


# Compute WOTS+ length
#
# len = len1 + len2
#
# len1 = ceil(m / log2(w))
# len2 = floor(log2(len1*(w-1)) / log2(w)) + 1
@lru_cache(None)
def wotslen(m, w):
    len1 = ceil(m / log(w, 2))
    len2 = floor(log(len1 * (w - 1), 2) / log(w, 2)) + 1
    return len1 + len2


# ------------------------------------------------
# Print header
# ------------------------------------------------

print(f"    n     h     d     h'    a     k     m    sigsize         p         Y"
      "           NIST          Loose      This work           "
      "NIST          Loose      This work")


# ------------------------------------------------
# Parameter sets
# ------------------------------------------------

maxsigs = 2^64                 # Maximum number of signatures allowed per seed
w = 16                         # Winternitz parameter 
for n,h,a,k,d in [
    (16,63,12,14,7),
    (16,66,6,33,22),
    (24,63,14,17,7),
    (24,66,8,33,22),
    (32,64,14,22,8),
    (32,68,9,35,17)
]:
    # Security parameter tsec = 8*n
    tsec = 8*n
    
    # Number of leaves in the hypertree
    leaves = 2^h
    
    # Number of leaves in each FORS tree
    t = 2^a

    # Output of hash functions
    hashbytes = n
    
    # Real interval field with precision tsec + 100 bits
    # used to compute probabilities with high precision
    F = RealIntervalField(tsec + 100)
    
    # Stop condition for sigma summation
    # when the probability contribution becomes extremely small
    sigma_stop = F(2)**(-20 * tsec)

    # ------------------------------------------------
    # Compute sigma
    # σ = Σ_{r ≥ 1} Pr[r hits] · Pr[FORS forgery | r hits]
    #
    # Pr[r hits] =
    #   C(qs, r) · (1/2^h)^r · (1 − 1/2^h)^(qs − r)
    #
    # Pr[FORS forgery | r hits] =
    #   (1 − (1 − 1/2^a)^r)^k
    # ------------------------------------------------
    sigma = F(0)
    r = 1

    while True:

        p = min(1, forgeryprob(a,r,k))
        q = qhitprob(leaves,maxsigs,r) * p

        sigma += q
        r += 1

        if q < sigma_stop:
            break

    # WOTS+ chain length (constant because n and w are fixed)
    wots_len = wotslen(8 * n, w)
    
    # Signature size in bytes
    sigsize = (h + k * (a + 1) + d * wots_len + 1) * n

    # Message digest size
    m = (
        floor((k*a + 7)/8)
        + floor((h - h/d + 7)/8)
        + floor((h/d + 7)/8)
    )

    p = (2**h) * k * (2**a)
    p_1 = maxsigs * k
    c = p/p_1
    y = p_1 /(1- p_1/(p-c+1))**c
    # ------------------------------------------------
    # Print parameters
    # ------------------------------------------------

    print(f"{int(n):5d} {h:5d} {d:5d} {int(h/d):5d} {a:5d} {k:5d} {m:5d} ", end="")
    print(f"{int(sigsize):10d}     {float(log(p,2)):.2f}     {float(log(y,2)):.2f}",end="")

    # ------------------------------------------------
    # Classical security 
    # ------------------------------------------------

    # NIST
    bitSec = -log(1/(2**tsec) + sigma, 2)
    print(f"{float(bitSec):15.3f}", end="")

    # Loose bound
    bitSec = -log((((2**h)*k*t)*k + w + 7)/(2**tsec) + sigma, 2)
    print(f"{float(bitSec):15.3f}", end="")

    # This work
    bitSec = -log((maxsigs*k + w + 7)/(2**tsec) + sigma, 2)
    print(f"{float(bitSec):15.3f}", end="")


    # ------------------------------------------------
    # Quantum security 
    # ------------------------------------------------

    # NIST
    bitSec = -log(1/(2**tsec) + sigma, 2)
    print(f"{float(bitSec/2):15.3f}", end="")

    # Loose bound
    A = sigma + (((2**h)*k*t)*k + 8) * 32 / (2**tsec)

    B = (
        (w + 2) * 12 / (2**(tsec/2))
        + (((2**h)*k*t)*k + 8) * 32 / (2**tsec)
    )

    qh = (-B + sqrt(B**2 + 4*A)) / (2*A)
    print(f"{float(log(qh,2)):15.3f}", end="")


    # This work
    A = sigma + (maxsigs*k + 8) * 32 / (2**tsec)

    B = (
        (w + 2) * 12 / (2**(tsec/2))
        + (maxsigs*k + 8) * 32 / (2**tsec)
    )

    qh = (-B + sqrt(B**2 + 4*A)) / (2*A)
    print(f"{float(log(qh,2)):15.3f}")
