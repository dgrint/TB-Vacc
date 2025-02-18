# TB-Vacc

This repository contains code to run a simulation study of the impact of imperfect sensitivity and specificity on the interpretation of a non-inferiority outcome in a trial of a new TB vaccine.

The file 'R scripts/Impact of specificity simulation.R' runs the simulation. For ease of readability, the code is duplicated under section headings for each combination of input parameters.

Simulations are run for three 2-year TB infection risk scenarios:
* 2%
* 5%
* 8%

For each scenario, simulations are run for the following combinations of sensitivity and specificity of the tool used to measure TB infection:
* 100%; 100;
* 100%; 95%;
* 95%; 100%;
* 95%; 95%;
* 95%; 98%;
* 64%; 85%;
 