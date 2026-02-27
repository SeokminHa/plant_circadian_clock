## Phase-Aligned Synthesis and Stability Modulation Amplify Circadian Protein Oscillations in Plants

This repository contains the MATLAB code used to simulate and analyze the ODE-based models described in the Supplementary Information of our paper: "Phase-aligned synthesis and stability modulation amplify circadian protein oscillations in plants."


## [Folder: Figure1_to_3] Mathematical Model and Parameter Fitting

Objective: To determine whether previously identified Protein-Protein Interactions (PPIs) are sufficient to recapitulate TOC1 and ZTL protein dynamics. If they are not, we investigate whether adding new PPIs can help fit the time-course data, and predict how novel interactions lead to reciprocal stabilization. (See Supplementary Methods for detailed methodology)

# Repository Structure
- multi_degradation_ODE_v6_5.m: Defines the ODE system describing rhythmic synthesis and degradation dynamics.
- solve_ODE.m: Numerically solves the ODE system using MATLAB’s built-in solvers.
- sa.m: Simulated annealing (SA) algorithm for parameter estimation.
- adjust_range.m: Utility function for adjusting parameter search ranges (see Supplementary Tables 4 and 5 for ranges).
- save_to_csv.m: Utility function for saving outputs to a CSV file.
- time_courses_1110.m: Plots the time courses of proteins and the stability of TOC1.

# key variable: bd_type
In the code, the core variable is bd_type. It is a 4-element binary array (e.g., [0,0,0,0]) where each component indicates the presence (1) or absence (0) of a specific protein interaction:
- First element: GI-PRR3 interaction
- Second element: ZTL_light-PRR3 interaction
- Third element: ZTL_dark-PRR3 interaction
- Fourth element: TOC1-GI interaction

# Execution Process
# A. Testing the baseline model (known PPIs only)
Run sa([0,0,0,0], ran_num) using a fixed random seed (e.g., sa([0,0,0,0], 13)).

- This procedure tests whether the model based strictly on previously known PPIs can simulate the experimentally measured time courses of TOC1, ZTL, PRR3, and GI proteins.
- Note: You can use various random numbers to run multiple instances and speed up the fitting process.
- Result: After running 1,000 different parameter fits using the SA algorithm (Supplementary Table 1), none were able to recapitulate the experimental time courses for TOC1 and ZTL. This confirms that previously identified PPIs are insufficient.

# B. Testing the extended model (novel PPIs added)
Run sa([1,1,1,1], ran_num) with a fixed random seed (e.g., sa([1,1,1,1], 13)).

- This introduces four new possible interactions into the model.
- Result: Using this code, we found 2,000 different parameter sets (Supplementary Table 2) capable of simulating the experimental time courses (Fig. 1h, right).
- When classifying these 2,000 successful models, ZTL-PRR3 and GI-PRR3 interactions—leading specifically to reciprocal stabilization—emerged as a consistent feature across the five dominant scenarios. Parameters corresponding to this scenario are saved in para_1110_full.csv.

# C. Plotting the results
Run time_courses_1110.m. This generates time-course plots for individual proteins and their respective complexes (Extended Data Fig. 3).
- Figure 11: TOC1 time course and specific complexes (e.g., TOC1-PRR3).
- Figure 12: ZTL time course and specific complexes (e.g., ZTL_light-PRR3).
- Figure 13: GI time course and specific complexes (e.g., ZTL_light-GI).
- Figure 14: PRR3 time course and specific complexes (e.g., GI-PRR3).

This script also calculates TOC1 stability and its amplitude. The WT TOC1 stability amplitude data is exported to TOC1 stability_amplitude.csv.

## [Folder Figure4_and_4S] Rhythmic Changes in Protein Stability

Objective: Rhythmic changes in protein stability over time play a pivotal role in generating strong circadian rhythms of proteins in other organisms, including mice and Drosophila. For a given protein $p(t)$ driven by mRNA $m(t)$, we assume the following dynamics:

$$\frac{dp}{dt}(t)=Am(t)-K\left(1+B\cos\left(\frac{2\pi}{24}(t-\varphi)\right)\right)p(t)$$

# Repository Structure
- `Fig4_rhythmic_stability.m` : Code for parameter estimation.
- `ODE_rhythmic_degra.m` : Numerically solves the ODE system using MATLAB’s built-in solvers based on the equation above.

# Execution process:

# A. Data Preparation
For the desired protein, save the normalized mRNA and protein time-course data as mRNA_time_course.csv and protein_time_course.csv, respectively.
- Row 1: Time in hours (a number between 0 and 24).
- Row 2: Relative abundance for each time point.
(See the provided CSV files for formatting examples).

# B. Parameter Estimation
Run 'Fig4_rhythmic_stability(ran_num).m' with a fixed random seed (e.g., Fig4_rhythmic_stability(1)). This will save the randomly selected parameters ($A$, $K$, $B$, $\varphi$) and calculate the Mean Square Error (MSE) between the simulated and input protein time courses.
- Note: Use various random numbers to run multiple instances in parallel and speed up the search process.

C. Parameter Selection
Select the top 100 parameter sets that induced the lowest MSE for downstream analysis.
