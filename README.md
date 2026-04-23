# ToricVerticalSystems.jl

[![CI](https://github.com/oskarhenriksson/ToricVerticalSystems.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/oskarhenriksson/ToricVerticalSystems.jl/actions/workflows/ci.yml)

This repository contains a proof-of-concept implementation of the methods described in the  paper [Toric invariance of vertically parametrized systems](https://arxiv.org/abs/2411.15134) by Elisenda Feliu and Oskar Henriksson. It also contains data from the database experiments described in the paper. 

## File descriptions
The repository contains the following files:
* A directory `src` that contains Julia functions for testing whether a network satisfies the various notions of toricity treated in the paper, as well as properties that can easily be checked in the presence of toricity such as (local) ACR and multistationarity. (These functions also come with tests in a directory `test`.)
* A directory `results` that contains list of networks in [ODEbase](https://www.odebase.org/) (as of November 2, 2023) with various properties. All the results are summarized in a text file `report.txt` and the table `results.csv` (see the file `legend.md` for explanation of the column headings).

## Examples
See the notebook `IDH_example.ipynb` for an illustration of how our functions can be applied to the IDHKP-IDH network (see Example 1.1 in the paper and Figure 8.1(a)).

## Dependencies
The code is based on `Oscar` and `DataStructures`, `Graphs` and `HomotopyContinuation`. 