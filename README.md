# Replication codes
Codes to replicate the results in "Oil shocks and stock market: Revisiting the dynamics" https://doi.org/10.1016/j.eneco.2021.105111
data.csv 
----------------------
Contains the data used to estimate the models
 - Period in months
 - op: Monthly returns of crude oil prices
 - xoi: Monthly returns of NYSE Arca Oil Index 
 - EPI: Economic Policy Uncertainty Index ( Scaled the series by dividing it by 100)
 - rv : Realised volatility constructed based on daily BSE Sensex
 - ret : Monthly returns of BSE Sensex

datalab.csv
-----------------------
Period in months

TVPcodes.r
-----------------------
R codes to estimate TVP VAR and to save the IRF and Accumulated IRF for generating plots as csv files
The paper estimates two models. The codes given in the files estimates firmst model1 -{op,xoi,epi,ret}.
The model 2 -{op,xoi,epi,rv}can be estimated easily by modifying the codes

TVP_funs.r
-----------------------
Codes for functions used in TVPcodes.r

3dgraphs.r
----------------------------

Codes to generate 3D plots given in Appendix
cirf_plots.r
-------------------------
Codes to generate the accumulated IRF plots given in the main text
