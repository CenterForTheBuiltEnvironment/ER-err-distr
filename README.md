# Evaluating Carbon Emission Factors for Accurate Accounting in Commercial Buildings #
## Project overview ##
Accurately quantifying operational carbon emissions from buildings is essential to verify that decarbonization targets are met. Yet, many current practices rely on overly simplified yearly average emission factors that overlook the temporal variability of the electrical grid generation. This issue becomes increasingly critical with the rising penetration of renewable energy and the adoption of demand-side strategies in buildings, such as load shifting. We evaluated how the temporal resolution (annual, season, time-of-day, season-hour, month-hour) of grid emission factors impacts the accuracy of carbon emissions accounting in commercial buildings across 18 grid regions in the United States. We applied emission factors ranging from annual averages to hourly intervals, sourced from a publicly available dataset, to hourly electricity consumption profiles from over 600 real commercial buildings. By considering hourly emission factors as the ground truth reference, we then quantified the resulting errors and uncertainties due to the reduced temporal resolution of the emission factors. Additionally, we assessed how buildings’ on-site solar PV generation affects avoided emissions accounting from grid exports when electricity generation exceeds demand. We found that annual and seasonal averages are not sufficiently accurate and they should not be used. Season-hour and month-hour emission factors consistently deliver sufficient accuracy across diverse U.S. grid regions, with median errors typically less than 5%. This is also the case when quantifying the avoided emissions when utility export is available, yet the results are more variable and dependent on the grid generation mix. We also found that as on-site solar generation capacity increases, the median normalized fractional error of the building operational carbon emissions can be between 5% to 30% when assessed using coarse emission factors such as annual or seasonal emission factors. These findings highlight the need for future standards and guidelines to require at least the use of season-hour or month-hour emission factors’ resolution to achieve acceptable emissions accounting accuracy.
## Highlights ##
*	Annual average emission factors are inaccurate and should not be used for assessing building operational or avoided carbon emissions. 
*	Season-hour and month-hour average emission factors provide the most reliable balance between accuracy and practicality for carbon emissions accounting.
*	The importance of temporal resolution increases significantly in grid regions with higher solar PV generation. 
*	Coarse-resolution emission factors risk underestimating the true impact of retrofits such as efficiency upgrades or thermal storage strategies.

## Graphical summary ##
<img width="2000" height="1125" alt="image" src="https://github.com/user-attachments/assets/94c4d04f-2006-4d92-bf58-67e02353cce4" />

## Emission factors
<img width="2000" height="1125" alt="image" src="https://github.com/user-attachments/assets/ae2d40a4-f912-4c59-9c3a-e26389de9f38" />

## Structure ##
This repository includes: 
* [code](code/): R code for data analysis
* [readfiles](readfiles/): all input files and cleaned data
* [docs](docs/): manuscript
* [figs](figs/): generated figures for visualization

## Policy implication
There is a growing interest and effort in decarbonizing the building industry and the power sector globally. Apart from technical development and financial incentives for carbon-free energy deployment, we believe that a more comprehensive assessment method is required to ensure the accuracy of carbon emissions accounting associated with building operations, especially on the temporal resolution of the emission factors. Therefore, we summarized the key policy implications of this study as:
* Differentiation between NZCBs and NZEBs. The operational carbon emissions associated with a building should not be considered as a simple multiplication of the annual energy consumption and a static emission factor, which falsely recognizes a net-zero energy building as a net-zero carbon building in operation (if excluding embodied carbon emissions assessment). The environmental impact of each kWh of energy needs to be assessed with improved temporal accuracy.
* Season-hour and month-hour average emission factors should be used for emissions accounting purposes instead of annual or seasonal average emission factors, since they are sufficiently accurate to reflect the seasonal and diurnal variations of the energy generation from the grid. Exceptionally, those emission factors may not be suitable for real-time applications due to weather synchronization.

