# PMNS with extension Fields

This project explores an extension of the Polynomial Modular Number System (PMNS) to extension fields with prime $p$ and extension degree $k$.

## Overview
The goal is to adapt PMNS arithmetic to work in extension fields while preserving efficient modular reduction.

## Methods
To achieve this, we implement reduction techniques based on:
- Montgomery reduction
- Babai rounding

Additionally, we extend the polynomial form of \(E\), which is used as the reduction polynomial in PMNS:
- $E = X^n -\lambda$
- $E = X^n -\alpha X^k -\beta$
And we try an optimised research of parameters using specific irreducible modulus to construct extension field.

## Dependencies
- SageMath (for the Python implementation)
- Jinja2 (for c code writing) => you can overwrite this dependance by using directly c code stock in template