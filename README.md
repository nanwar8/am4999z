# AM4999Z Honours Project

The project uses stochastic simulation methods to study a two-box gene regulatory model with protein production, degradation, diffusion between compartments, and mRNA dynamics. 
Our focus is on integrating stochasticity and time delays into one central model using the Gillepsie's algorithm.
The model is based off of the one presented in the following paper: https://doi.org/10.1039/C1MB05086A

## Project Overview
The main goal of this project is to simulate the time evolution of three variables:

- `p1`: protein population in box 1
- `p2`: protein population in box 2
- `m`: mRNA population

The model is simulated using the Gillespie algorithm, which generates exact stochastic trajectories based on reaction rates. 
I also implemented methods for sampling the system on a uniform time grid so that trajectories can be compared more easily and used for plotting and histogram analysis.
Our main changes from the paper includes p1, which serves as the inactivated version of the protein, and must be converted to p2, which is the activated protein interacting with m. 

## Features

The code currently includes:
- exponential waiting-time sampling for stochastic events
- event selection based on normalized reaction probabilities
- reaction-rate calculations for the two-box system
- state updates for protein and mRNA copy numbers
- a standard Gillespie simulation with non-uniform event times
- interpolation/sampling onto a uniform time grid
- a direct Gillespie-on-grid implementation
- time-series plots and uniform time histograms for `p1`, `p2`, and `m`
- parameter sweeps over diffusion values `j`

## Model components

The simulation includes the following types of events:
1. production of `p1`
2. degradation of `p1`
3. production of `p2`
4. basal production of `m`
5. `p2`-dependent production of `m`
6. degradation of `m`
7. degradation of `p2`
8. diffusion of `p1` into `p2`

## File Contents
This repository includes code (plotfit_p1_p2_m_timestepping_gillepsie.py) used to:
- run stochastic simulations of the model
- compare non-uniform and uniform-time implementations
- generate time-series plots
- generate histograms of molecule-number distributions
- explore how changing diffusion parameters affects system behavior

## Methods
The project is based on stochastic simulation using the Gillespie algorithm. Two implementations are included:
### 1. Standard Gillespie simulation
This version computes reaction times directly from the total rate and stores system states at event times.
### 2. Uniform-grid sampling
This version converts the stochastic trajectory into values on a fixed time grid, making it easier to compare trajectories and construct histograms.

## Parameters

The script currently uses fixed parameter values for:

- production rates
- degradation rates
- feedback terms
- diffusion strength `j`
- simulation end time
- uniform sampling interval

These can be modified directly in the script to explore different regimes of the model, and are defined in the script as well.
We work off of the case 3 parameters presented in the paper. 'j' is varied as we are interested in observing what happens when we effect the rate at which p1 becomes p2.

## Outputs

The code produces:

- time-series plots for `p1`, `p2`, and `m`
- histograms of simulated values after an initial transient period
- comparisons across different values of the diffusion parameter `j`

## Notes

This repository reflects work in progress for my honours project. The code is being actively organized and updated as the project develops.
