# Joint_Modeling_Mediation
# The three folders provide codes for the simulation analysis (four settings) and data analysis (CD4 and PI)

# Within each simulation setting, it contains (1)data generation code; (2) joint modeling fitting code;
# (3) NDE, NIE and CI computation code for each time; (4) parallel code for multiple time points; (5) summary code
#For simulatino data, to repeat the result, the following step is needed (1) run data generation code with R and save the simulated datasets; (2) run the model #fitting code with SAS to obtain fitting coefficient and variance-covariance matrix for each datasets; (3) using parallel code to run NDE, NIE computation code to #obtain point estimate and CI; (4) using summary result to generate simulation result table for each setting

# For CD4 example, it contains (1) Model fitting for model I-II and separate model; (2) NDE, NIE and CI computation code for each time for model I-II;
#(3) NDE, NIE and CI computation code for each time for model I-II; (4) parallel codes for multiple time points; (5) code to generate figure

# For PI example, it contains (1) Model fitting for model I-II; (2) Model fitting for separate model (3) NDE, NIE and CI computation code for each time for model I-II;
#(4) NDE, NIE and CI computation code for each time for model I-II; (5) parallel codes for multiple time points; (5) code to generate figure


