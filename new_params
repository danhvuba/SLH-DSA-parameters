from sage.all import *
from functools import lru_cache

# ------------------------------------------------
# Fixed parameters
# ------------------------------------------------

tsec = 328                     # Security parameter tsec = 8*n
maxsigs = 2^64                 # Maximum number of signatures allowed per seed
maxsigbytes = 200000           # Maximum signature size we allow to print (bytes)

w = 16                         # Winternitz parameter (fixed)

# Real interval field with precision tsec + 100 bits
# used to compute probabilities with high precision
F = RealIntervalField(tsec + 100)

# if sigma exceeds this value we discard the parameter set
sigmalimit = F((2 * 64) * 2^(-tsec))

# Stop condition for sigma summation
# when the probability contribution becomes extremely small
sigma_stop = F(2)^(-20 * tsec)

# log2(max signatures)
s = int(log(maxsigs, 2))

# n = security bytes
n = tsec // 8


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


# WOTS+ chain length (constant because n and w are fixed)
wots_len = wotslen(8 * n, w)


# ------------------------------------------------
# Print table header
# ------------------------------------------------

print(f"                         params                     byte            "
      "  Keccak calls (F,H,PRF,T_len)                      classical security                          quantum security")

print(f"    n     h     d     h'   a     k     w     m   "
      "sigsize        KeyGen        Signing        Verify      "
      "NIST          Loose          This work     "
      "NIST          Loose          This work")


# =================================================
# Parameter search
# =================================================

# Iterate possible hypertree heights around log2(maxsigs)
for h in range(s - 10, s + 10):

    # Number of leaves in the hypertree
    leaves = 2^h

    # Height of each FORS Merkle tree
    for a in range(1, 30):

        # Number of leaves in each FORS tree
        t = 2^a

        # Number of FORS trees
        for k in range(1, 64):

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

                # Probability that r signatures allow forgery
                p = min(1, forgeryprob(a, r, k))

                # Combined probability of r hits and successful forgery
                q = qhitprob(leaves, maxsigs, r) * p

                sigma += q
                r += 1

                # Stop when contribution becomes negligible
                if q < sigma_stop:
                    break

            # If sigma exceeds the security threshold
            # discard this parameter set
            if sigma > sigmalimit * k:
                continue


            # ------------------------------------------------
            # Hypertree layers
            # ------------------------------------------------

            for d in range(2, h):

                # d must divide h
                if h % d != 0:
                    continue

                if h > 64 + (h / d):
                    continue

                # Signature size in bytes
                sigsize = (h + k * (a + 1) + d * wots_len + 1) * n

                # Discard if signature becomes too large
                if sigsize >= maxsigbytes:
                    continue

                # ------------------------------------------------
                # Message digest size
                # ------------------------------------------------

                # m = number of bytes required for message digest (output of H_msg)
                m = (
                    floor((k * a + 7) / 8)
                    + floor((h - h/d + 7) / 8)
                    + floor((h/d + 7) / 8)
                )

                # SHA3-512 has a 64-byte output.
                if m >= 64:
                    continue

                # ------------------------------------------------
                # Keccak permutation cost estimation
                # ------------------------------------------------

                # PRF and F cost
                keccakPRF = keccakF = (
                    ceil(((2 * n + 32) * 8 + 4) / (1600 - 1024)) + 1
                )

                # Hash H cost
                keccakH = (
                    ceil(((3 * n + 32) * 8 + 4) / (1600 - 1024)) + 1
                )

                # T_len cost
                keccakT = (
                    ceil(((n + 32 + wots_len * n) * 8 + 4) / (1600 - 1024)) + 1
                )

                # ------------------------------------------------
                # Key generation cost
                # ------------------------------------------------

                keyGen = (
                    (2^(h/d) * w * wots_len) * keccakF
                    + (2^(h/d) - 1) * keccakH
                    + (2^(h/d) * wots_len) * keccakPRF
                    + (2^(h/d)) * keccakT
                )

                # ------------------------------------------------
                # Signing cost
                # ------------------------------------------------

                signSpeed = (
                    (k * t + d * 2^(h/d) * w * wots_len) * keccakF
                    + (k * (t - 1) + d * (2^(h/d) - 1)) * keccakH
                    + (k * t + d * 2^(h/d) * wots_len) * keccakPRF
                    + (d * 2^(h/d)) * keccakT
                )

                # ------------------------------------------------
                # Verification cost
                # ------------------------------------------------

                verifySpeed = (
                    (k + d * w * wots_len) * keccakF
                    + (k * a + h) * keccakH
                    + d * keccakT
                )

                # ------------------------------------------------
                # Print parameter set
                # ------------------------------------------------

                print(
                    f"{n:5d} {h:5d} {d:5d} {int(h/d):5d}"
                    f"{a:5d} {k:5d} {w:5d} {m:5d}"
                    f"{int(sigsize):10d}",
                    end=""
                )

                print(
                    f"{int(keyGen):15d}{int(signSpeed):15d}{int(verifySpeed):15d}",
                    end=""
                )


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
