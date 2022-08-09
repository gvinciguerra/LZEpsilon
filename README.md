# $\text{LZ}_\varepsilon$

This repository implements a rank/select dictionary, described [here](https://doi.org/10.4230/LIPIcs.ISAAC.2021.64), based on a combination of Lempel-Ziv and [LA-vector](https://github.com/gvinciguerra/la_vector) compression.

## Usage

This is a header-only library. To compile the [example](example.cpp), use the following commands:

```sh
git clone https://github.com/gvinciguerra/LZEpsilon.git
cd LZEpsilon
cmake . -DCMAKE_BUILD_TYPE=Release
make -j8
```

## License

This project is licensed under the terms of the Apache License 2.0.

If you use this code for your research, please cite:

> Paolo Ferragina, Giovanni Manzini, and Giorgio Vinciguerra. Repetition- and linearity-aware rank/select dictionaries.  In: Proceedings of the 32nd International Symposium on Algorithms and Computation (ISAAC), 2021.

```bibtex
@inproceedings{Ferragina:2021isaac,
  author = {Ferragina, Paolo and Manzini, Giovanni and Vinciguerra, Giorgio},
  booktitle = {Proceedings of the 32nd International Symposium on Algorithms and Computation (ISAAC)},
  doi = {10.4230/LIPIcs.ISAAC.2021.64},
  title = {Repetition- and linearity-aware rank/select dictionaries},
  year = {2021}}
```
