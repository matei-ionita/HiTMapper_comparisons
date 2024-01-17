# HiTMapper Comparisons
Comparison of HiTMapper with manual gating and other algorithms as well as several analysis perspectives

| Dataset  |  Cytometry Modality  |  Description  |  Reference  |  
|---------------|---------------|-----------|-----------|
|  Levine Benchmark Data  |  Mass Cytometry (CyTOF)  |  The data set is a 32-dimensional mass cytometry (CyTOF) data set, consisting of protein expression levels for n = 265,627 cells, p = 32 protein markers (dimensions), and k = 14 manually gated cell populations (clusters), from h = 2 individuals. Cluster labels are available for 39% (104,184) of the cells.  | Levine et al. (2015)  | 
|  COVID Vaccine (Day ?)  |  Mass Cytometry (CyTOF)   |  Description  |  Mathew, Giles, Baxter, Oldridge,  Greenplate, Wu, Alanio, et al. (2020)  |  
|  Acute COVID (Day 0 and Day 7)  | Mass Cytometry (CyTOF)  |  Description  |  Mathew, Giles, Baxter, Oldridge,  Greenplate, Wu, Alanio, et al. (2020)  |  
|  Acute COVID (Day ?)  |  Flow Cytometry   |  Small subset of the acute COVID flow cytometry data  |  Mathew, Giles, Baxter, Oldridge,  Greenplate, Wu, Alanio, et al. (2020) | 


## Comparison of HiTMapper, FlowSOM, Rphenograph, and FastPG using 32-dimensional CyTOF data from Levine et al


| Method  |  Description  |  Reference  |  GitHub |
|---------------|---------------|-----------|-----------|
|  HiTMapper  |  Description  |  Version  |  [HiTMapper Code](https://github.com/matei-ionita/HiTMapper/tree/master)  | 
|  FlowSOM  |  Description  |  Version  |  [FlowSOM Code](https://github.com/SofieVG/FlowSOM)  |  
|  Rphenograph  |  Description  |  Version  |  [Rphenograph Code](https://github.com/JinmiaoChenLab/Rphenograph)  | 
|  FastPG  |  Description  |  Version  |  [FastPG Code](https://github.com/sararselitsky/FastPG)  |

- Additionally, there's a command-line implementation of FastPG by MCMICRO, https://github.com/labsyspharm/mcmicro-fastPG 

## Mass Cytometry (CyTOF) data analysis using HiTMapper
* Led by Matei Ionita
* Acute COVID (Day ?)
* COVID Vaccine (Day ?)

## Flow Cytometry data analysis using HiTMapper
* Led by Van Truong
* Acute Covide (Day 0 and Day 7)


