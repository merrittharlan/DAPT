# DAPT (Discharge Estimation via Arrays of Pressure Transducers)
Code to estimate river discharge from pressure transducers

Included in this git repository is results (and scripts to process the data and produce results) from the study 'Discharge Estimation from Dense Arrays of Pressure Transducers.'

Data: Included in the "input" folder are pressure transducer and river centerline data. In the "TOWNS" folder, there are data from the Tanana, Olentangy, Willamette, North Saskatchewan, and Sacramento rivers. In the "PAD"" folder, there are data from river reaches in the Peace-Athabasca Delta, including discharge from ADCP transects and GNSS data. The "Raw" folder contains data from the entire 2018-2019 field campaign, including several pressure transducers and discharge measurements not explicitly included in this study. 

Code: To reproduce results from the paper, and estimate discharge using the Bayesian-AMHG-Manning algorithm (BAM), follow the scripts in this order:

1) Format the water surface elevation data in both the gauged rivers (TOWNS), and ungauged rivers (PAD), with TOWNS_WSE.Rmd, and PAD_WSE.Rmd. For the PAD, format the ADCP discharge data with PAD_Discharge.Rmd.

2) Estimate width from the digital elevation model (downloaded either from the National Elevation Dataset, or Canadian Digital Elevation Model), and create BAM data files in order to run BAM with TOWNS_BAMDATA.Rmd, and PAD_BAMDATA.Rmd

3) Run sensitivity analyses, and compare with a temporary gauge using scripts DAPT_BAM_runs.R with functions in DAPT_Sensitivity.R and DAPT_Plotting.R for the TOWNS dataset, and DAPT_BAM_runs.R for the PAD dataset. Instead of generating the data in steps 1 and 2, the R BAMdata files needed to run these scripts are under the folder "output/BAMdata". All sensitivity results are in the "Sensitivity" folder. 

Please refer to the citation below, or send me an email at mharlan@umass.edu if you have any questions!


Harlan, M. E., Gleason, C. J., Altenau, E. H., Butman, D., Carter, T., Chu, V. W., et al. (2021). Discharge Estimation from Dense Arrays of Pressure Transducers. Water Resources Research, 57, e2020WR028714. https://doi.org/10.1029/2020WR028714