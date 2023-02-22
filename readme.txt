I. R Based Methods

  To run the R based methods, make sure you have R and all the required dependency installed.
  You can either run them on command line or on your console.

  For example, to run our LE method on the simulation SBM dataset, use the following code on bash:
  # Set directory 
    cd Simulation_Experiments
  # Run the Rscript on command line
    Rscript sim-LE.R SBM

  To run our LE method on the real enron dataset, use the following code on bash:
  # Set directory 
    cd Real_Data_Examples
  # Run the Rscript on command line
    Rscript Real-LE.R enron

  It is possible to set your directory and edit the sequences of rho, degrees
  and the number of iterations in the sim-${method}.R and Real-${method}.R files.

II. Python Based Method (NS)

  To run the Python based method NS, make sure you have python3 and all the required dependency installed.
  You can either run them on command line or on your console.

  For example, to run our LE method on the simulation SBM dataset, use the following code on bash:
  # Set directory 
    cd Simulation_Experiments
  # Run the Rscript on command line
    python3 sim-NS.py SBM

  To run our LE method on the real enron dataset, use the following code on bash:
  # Set directory 
    cd Real_Data_Examples
  # Run the Rscript on command line
    python3 Real-NS.python enron

  p.s. For the real data version, there is no argument for degrees because it cannot be specified for real networks.

  It is possible to edit the sequences of rho, degrees and the number of iterations in the sim-NS.py and Real-NS.py files.
  
 Referecnce:
 - Tianxi Li, Yun-Jhong Wu, Elizaveta Levina & Ji Zhu (2023) Link Prediction for Egocentrically Sampled Networks, Journal of Computational and Graphical Statistics, 
   DOI: 10.1080/10618600.2022.2163648(https://www.tandfonline.com/doi/suppl/10.1080/10618600.2022.2163648)
