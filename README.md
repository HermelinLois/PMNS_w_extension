# PMNS with extension Rings

This project explores an extension of the Polynomial Modular Number System (PMNS) to extension rings with prime $p$ and extension degree $k$.

## Overview
The goal is to adapt PMNS arithmetic to work in extension rings while preserving efficient modular reduction.

## Methods
To achieve this, we implement reduction techniques based on:
- Montgomery reduction
- Babai rounding

Additionally, we extend the polynomial form of \(E\), which is used as the reduction polynomial in PMNS:
- $E = X^n -\lambda$

## Dependencies
- SageMath (for the Python implementation)