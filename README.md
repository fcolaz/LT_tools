# LT_Tools

This repository provides MATLAB scripts for the finite-length analysis of LT codes under **peeling decoding** [1] and **inactivation decoding** [2].  

## Features

- **Peeling decoding analysis:** Computes the probability of decoding failure, stopping probability, expected cloud size, and ripple size.  
- **Inactivation decoding analysis:** Estimates the expected number of inactivations, which serves as a measure of decoding complexity.  

## Usage

### Peeling Decoding Analysis

The script **`run_finite_length_analysis_peeling.m`** runs the finite-length analysis introduced in [1] for an LT code under peeling decoding.  

- The core function **`peeling_decoding_analysis.m`** implements the analysis using dynamic programming.  
- The script returns:  
  - Probability of decoding failure.  
  - Probability that decoding stops with `u` unrecovered input symbols.  
  - Expected cloud and ripple sizes as functions of `u`.  

### Inactivation Decoding Analysis  

The script **`run_finite_length_analysis_inactivation.m`** runs the finite-length analysis introduced in [2] for an LT code under inactivation decoding.  

- The core function **`get_n_inact.m`** implements the analysis via dynamic programming.  
- The script returns:  
  - The expected number of inactivations, which reflects decoding complexity (see [2] for details).  

## Citation  

If you use this implementation in your research, please cite it using the DOI:  

[![DOI](https://zenodo.org/badge/DOI/YOUR_DOI_HERE.svg)](https://doi.org/10.5281/zenodo.14889584)  

```bibtex
@misc{your_repo_name,
  author    = {Lázaro, Francisco},
  title     = {LT_Tools: Finite-Length Analysis of LT Codes},
  year      = {2025},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.14889584}
}
```

## References  

- **[1]** Richard Karp, Michael Luby, and Amin Shokrollahi. *"Finite length analysis of LT codes."* Proceedings of the International Symposium on Information Theory (ISIT), 2004.  

- **[2]** Francisco Lázaro, Gianluigi Liva, and Gerhard Bauch. *"Inactivation decoding of LT and Raptor codes: Analysis and code design."* IEEE Transactions on Communications, vol. 65, no. 10, pp. 4114-4127, 2017.  

