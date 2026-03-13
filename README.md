# SLH-DSA Parameter Evaluation (Sage)

This repository contains two Sage scripts used in the paper **Another Look at the Security Bounds of the SPHINCS+ Framework** to evaluate the parameter sets of **SLH-DSA**, the stateless hash-based digital signature scheme standardized by NIST.

## Overview

The provided Sage code implements the security evaluation and parameter analysis for SLH-DSA parameter sets.

These scripts reproduce the evaluation methodology used in **Another Look at the Security Bounds of the SPHINCS+ Framework** to assess the security bits and signature sizes of the scheme.

The implementation follows the evaluation approach described in:

- *Another Look at the Security Bounds of the SPHINCS+ Framework* (evaluation methodology used in the paper)
- *The SPHINCS+ Signature Framework*
- *SPHINCS+: Submission to the NIST Post-Quantum Cryptography Project*, v3.1

## Files

- **NIST_params.sage**  
  Sage script for evaluating SLH-DSA parameter sets and computing the corresponding security bits.

- **new_params.sage**  
  Sage script used to generate and test candidate parameter configurations.

## Requirements

- **SageMath**

## Usage

Run the scripts using Sage:

```bash
sage script_name.sage
```

The scripts will output tables containing the evaluated parameter sets and their estimated security bits.

## References

1. *Another Look at the Security Bounds of the SPHINCS+ Framework*
2. *The SPHINCS+ Signature Framework*
3. *SPHINCS+: Submission to the NIST Post-Quantum Cryptography Project*, Version 3.1
