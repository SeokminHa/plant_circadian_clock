## Phase-Aligned Synthesis and Stability Modulation Amplify Circadian Protein Oscillations in Plants

This repository contains the MATLAB code used to simulate and analyze the ODE-based models described in the Supplementary Information of our paper: "Phase-aligned synthesis and stability modulation amplify circadian protein oscillations in plants."


## [Folder: Figure1] Mathematical Model and Parameter Fitting

Objective: To determine whether previously identified Protein-Protein Interactions (PPIs) are sufficient to recapitulate TOC1 and ZTL protein dynamics. If they are not, we investigate whether adding new PPIs can help fit the time-course data, and predict how novel interactions lead to reciprocal stabilization (See Supplementary Methods for detail).

# Repository Structure
- `multi_degradation_ODE_v6_5.m`: Defines the ODE system describing rhythmic synthesis and degradation dynamics.
- `solve_ODE.m`: Numerically solves the ODE system using MATLAB’s built-in solvers.
- `sa.m`: Simulated annealing (SA) algorithm for parameter estimation.
- `adjust_range.m`: Utility function for adjusting parameter search ranges (see Supplementary Tables 4 and 5 for ranges).
- `save_to_csv.m`: Utility function for saving outputs to a CSV file.
- `classify_1111_perturb_add.m` : Performs post hoc PPI-removal perturbation analysis and classifies the 2,000 successful models into essential PPI scenarios.
- `time_courses_1110.m`: Plots the time courses of proteins and the stability of TOC1.

# key variable: bd_type
In the code, the core variable is bd_type. It is a 4-element binary array (e.g., [0,0,0,0]) where each component indicates the presence (1) or absence (0) of a specific protein interaction:
- First element: GI-PRR3 interaction
- Second element: ZTL_light-PRR3 interaction
- Third element: ZTL_dark-PRR3 interaction
- Fourth element: TOC1-GI interaction

# Execution Process
# A. Testing the baseline model (known PPIs only)
Run `sa([0,0,0,0], ran_num)` using a fixed random seed (e.g., sa([0,0,0,0], 13)).

- This procedure tests whether the model based strictly on previously known PPIs can simulate the experimentally measured time courses of TOC1, ZTL, PRR3, and GI proteins.
- Note: You can use various random numbers to run multiple instances and speed up the fitting process.
- Result: After running 1,000 different parameter fits using the SA algorithm (Supplementary Table 1), none were able to recapitulate the experimental time courses for TOC1 and ZTL. This confirms that previously identified PPIs are insufficient.

# B. Testing the extended model (novel PPIs added)
Run `sa([1,1,1,1], ran_num)` with a fixed random seed (e.g., sa([1,1,1,1], 13)).

- This introduces four new possible interactions into the model.
- Result: Using this code, we found 2,000 different parameter sets (Supplementary Table 2) capable of simulating the experimental time courses (Figure 1H, right).

# C. Filtering successful models to identify essential PPIs

Run `classify_1111_perturb_add.m`. For each of the 2,000 fitted models, this selects the PPIs that are essential for maintaining the fitted protein time-course dynamics — an interaction is retained as essential if deleting it substantially increases the fitting error.

- Output: `classify_result_1111_essential.csv`, a 2,000 × 8 matrix. Each row is one model; each column pair is one interaction (two elements because A→B and B→A effects are evaluated separately).

| Columns | Interaction |
|---|---|
| 1–2 | TOC1 – GI |
| 3–4 | ZTL_light – PRR3 |
| 5–6 | ZTL_dark – PRR3 |
| 7–8 | GI – PRR3 |

Each element is a number from 1 to 4 (see Figure S2A):
1 = no interaction, 2 = stabilize, 3 = destabilize, 4 = no stability change

 The distribution of these scenarios was highly skewed (see `classify_result_sorted_1111_essential.csv`), with the three most frequent scenarios, R1, R2, and R3, dominating the results. ZTL-PRR3 and GI-PRR3 interactions—leading specifically to reciprocal stabilization—emerged as a consistent feature. Parameters corresponding to this scenario are saved in `para_1110_full.csv`.


# D. Plotting the results
Run `time_courses_1110.m`. This generates time-course plots for individual proteins and their respective complexes (Figure S4).
- Figure 11: TOC1 time course and specific complexes (e.g., TOC1-PRR3).
- Figure 12: ZTL time course and specific complexes (e.g., ZTL_light-PRR3).
- Figure 13: GI time course and specific complexes (e.g., ZTL_light-GI).
- Figure 14: PRR3 time course and specific complexes (e.g., GI-PRR3).



## [Folder: Figure3A, 3B] Stability Calculation

Calculates TOC1 degradation rate and stability over time, and the amplitude of TOC1 stability, for WT and mutant conditions.

# Repository Structure

- `multi_degradation_ODE_v6_5.m`: Defines the ODE system describing rhythmic synthesis and degradation dynamics.
- `degradation_rate_WT.m`: Calculates the degradation rate of TOC1 over time.
- `stability_TOC1.m`: Calculates the stability of TOC1 over time.

# A. Calculate the degradation rate of TOC1 over time (WT and mutants)

Run `degradation_rate_WT.m`. To simulate the mutants, edit line 64:

| Condition | Line 64 | Output |
|---|---|---|
| WT | (unchanged) | `TOC1_degra_wt.csv` |
| *prr3-1* | `tp=0;` | `TOC1_degra_p.csv` |
| *ztl103* | `tz=0;` | `TOC1_degra_z.csv` |
| *prr3-1ztl103* | `tz=0; tp=0;` | `TOC1_degra_pz.csv` |

- Output format: the first row is time; the following 257 rows are the TOC1 degradation rate over time, one row per model.

# B. Calculate TOC1 stability (WT and mutants)

Run stability_TOC1.m. This computes TOC1 stability from the degradation rates above and extracts its amplitude.

Output:

- `TOC1_raw_stab_wt_summary.csv` — TOC1 stability in WT. Row 1: time. Row 2: mean stability across the 257 models. Row 3: standard deviation. Equivalent files are produced for the three mutant conditions.
- `TOC1_stab_amp_raw_all.csv` — amplitude of TOC1 stability for each of the 257 models (one column per model).

| Row | Condition |
|---|---|
| 1 | WT |
| 2 | *prr3-1* |
| 3 | *ztl103* |
| 4 | *prr3-1ztl103* |


## [Folder Figure4_and_S5] Rhythmic Changes in Protein Stability

Objective: Rhythmic changes in protein stability over time play a pivotal role in generating strong circadian rhythms of proteins in other organisms, including mice and Drosophila. For a given protein $p(t)$ driven by mRNA $m(t)$, we assume the following dynamics:

$$\frac{dp}{dt}(t)=Am(t)-K\left(1+B\cos\left(\frac{2\pi}{24}(t-\varphi)\right)\right)p(t)$$

# Repository Structure
- `Fig4_rhythmic_stability.m`: Code for parameter estimation.
- `ODE_rhythmic_degra.m`: Numerically solves the ODE system using MATLAB’s built-in solvers based on the equation above.

# Execution process:

# A. Data Preparation
For the desired protein, save the normalized mRNA and protein time-course data as mRNA_time_course.csv and protein_time_course.csv, respectively.
- Row 1: Time in hours (a number between 0 and 24).
- Row 2: Relative abundance for each time point.
(See the provided CSV files for formatting examples).

# B. Parameter Estimation
Run `Fig4_rhythmic_stability(ran_num).m` with a fixed random seed (e.g., Fig4_rhythmic_stability(1)). This will save the randomly selected parameters ($A$, $K$, $B$, $\varphi$) and calculate the Mean Square Error (MSE) between the simulated and input protein time courses.
- Note: Use various random numbers to run multiple instances in parallel and speed up the search process.

C. Parameter Selection
Select the top 100 parameter sets that induced the lowest MSE for downstream analysis.
