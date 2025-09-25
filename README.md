Data and Scripts for 'Ionic Liquid–Electrode Interface: From One Law to Fit Them All to One Model to Predict Them All'

----------Description-----------

This project contains the complete dataset and all supporting scripts for the article, 'Ionic Liquid–Electrode Interface: From One Law to Fit Them All to One Model to Predict Them All.' Included are the empirical data, scripts to reproduce all figures in the publication, and a toolkit for data analysis and energy density prediction in ionic liquid systems.

----------Project Structure----------

    • Analysis.ipynb: A Jupyter Notebook to automatically read and analyze experimental or simulated .csv data. They can conduct either capacitance-potential C(U) or surface charge density-potential Q(U) fitting on experimental or simulated .csv data.
    
    • Prediction.py: A Python script to predict energy density by running a parameter sweep.
    
    • requirements.txt: A list of required Python packages to run these scripts.
    
    • exampledatastructure.csv: Example data structure that can be analyzed by Analysis.ipynb
    
    • “Data” folder: .csv files containing  empirical data used in the article.
    
    • “Figure” folder: scripts for reproducing the figures in the article.
    
    
----------Installation of the required packages----------

pip install -r requirements.txt

----------Usage & Workflow of Analysis.ipynb and Prediction.py----------

    1. Prepare Data: Place your .csv data files in the root directory of this project. Ensure they follow the format described below.
    
    2. Analyze Experimental Data:
    
        ◦ Open and run the Analysis.ipynb notebook.
        
        ◦ The script will automatically find all .csv files, perform the fitting, print the resulting model parameters to the screen, and save a plot as fitting.png.
        
    3. Predict Performance:
    
        ◦ Run the Prediction.py script from your terminal.
        
        ◦ This script runs a parameter sweep to find the optimal theoretical parameters for energy density (with quantum effect taken into account).
        
        ◦ It will save its output to bestcomb.txt and generate a plot named predictionfinal.png showing dependence of surface charge density, capacitance, and energy density on applied voltage.
        
        
----------Data Format----------

The input .csv files should be structured as follows:

    • Row 1: Metadata keys (e.g., uM_pos, sM_pos). sM (surface charge density at PMC) is mandatory. uM (PMC) can be left blank to be treated as a fitted parameter.
    
    • Row 2: Metadata values.
    
    • Data Table: The numerical data should follow the metadata, with an optional header row.
    
    • Default Column Order (if no header): U, dU, S, C, dC.
    
    • Units:
    
        ◦ Potential (U): V
        
        ◦ Capacitance (C): μF/cm²
        
        ◦ Surface Charge Density (S): μC/cm²
        
        
----------License----------

    • The code in this repository is licensed under the MIT License.
    
    • The data is licensed under the Creative Commons Attribution 4.0 International (CC-BY 4.0) license.
    
    
----------Authors & Citation----------

    • Ba Long Nguyen (ORCID: 0009-0003-7682-9851): balongn99@gmail.com
    
    • Vladislav Ivanistsev (ORCID: 0000-0003-4517-0540): vladislav.ivanistsev@gmail.com
    
    
How to Cite

If you use this code or data in your research, please cite both the associated article and the project itself.

    1. Article DOI: https://doi.org/10.1016/j.elecom.2025.108049

    2. Scripts and data Zenodo DOI: to be updated.


----------Acknowledgments----------

This work was supported by the Estonian Ministry of Education and Research (TK210) and the Estonian Research Council (grant STP52). Results were obtained using the High Performance Computing Center of the University of Tartu. We are grateful to N. Nishi, S. Katakura, and R. Costa for sharing their empirical data.

