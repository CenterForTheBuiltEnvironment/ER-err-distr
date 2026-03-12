# Evaluating Carbon Emission Factors for Accurate Accounting in Commercial Buildings #
## Project overview ##
Accurately quantifying operational carbon emissions from buildings is essential to verify that decarbonization targets are met. Yet many current practices rely on overly simplified yearly average emission factors that overlook the temporal variability of the electrical grid generation. This issue becomes increasingly critical with the rising penetration of renewable energy and the adoption of demand-side strategies in buildings, such as load shifting. Existing literature on net-zero carbon buildings rarely examines how the temporal resolution of emission factors affects accounting accuracy. In this study, we evaluated a range of grid emission factors (annual, seasonal, time-of-day, season-hour, and month-hour) and quantified their impact on carbon emissions accounting accuracy for commercial buildings across 18 U.S. grid regions. These emission factors, derived from publicly available datasets, are applied to measured hourly electricity consumption profiles from over 600 real commercial buildings. By considering hourly emission factors as the benchmark reference, we then quantified the resulting errors and uncertainties due to the reduced temporal resolution of the emission factors. Additionally, we assessed how buildings' on-site solar PV generation affects avoided emissions accounting from grid exports when electricity generation exceeds demand. As a result, we found that annual and seasonal averages are unreliable and should not be used for net-zero target assessments. Instead, we recommend incentivizing the adoption of season-hour and month-hour emission factors, which consistently deliver sufficient accuracy across diverse U.S. grid regions, with median errors typically less than 10%. This is also the case when quantifying the avoided emissions when utility export is available, yet the results are more variable and dependent on the grid generation mix. In particular, when coarse emission factors such as annual or seasonal averages are used, increasing onsite solar capacity can raise the median normalized fractional error to approximately 15% for operational emissions and to more than 100% for avoided emissions. These findings highlight the need for future standards and guidelines to consider at least the use of season-hour or month-hour emission factors' resolution to achieve acceptable emissions accounting accuracy.
## Highlights ##
*	Annual average emission factors are inaccurate and should not be used for assessing building operational or avoided carbon emissions. 
*	Season-hour and month-hour average emission factors deliver sufficient accuracy with median errors typically below 10%, providing the best practical balance for accounting.
*	The importance of temporal resolution increases significantly in grid regions with higher solar PV generation. 
*	Coarse-resolution emission factors risk underestimating the true impact of retrofits such as efficiency upgrades or thermal storage strategies.

## Graphical summary ##
<img width="2000" height="1125" alt="image" src="<img width="2000" height="1125" alt="image" src="https://github.com/user-attachments/assets/be4cc11c-8cd0-40eb-95d1-0818aa895449" />

## Emission factors
<img width="2000" height="1125" alt="image" src="https://github.com/user-attachments/assets/ae2d40a4-f912-4c59-9c3a-e26389de9f38" />

## Structure ##
This repository includes: 
* [code](code/): R code for data analysis
* [readfiles](readfiles/): all input files and cleaned data
* [docs](docs/): manuscript and supplementary material
* [figs](figs/): generated figures for visualization

## Policy implication
There is a growing interest and effort in decarbonizing the building industry and the power sector globally. Apart from technical development and financial incentives for carbon-free energy deployment, we believe that a more comprehensive assessment method is required to ensure the accuracy of carbon emissions accounting associated with building operations, especially on the temporal resolution of the emission factors. Therefore, we summarized the key policy implications of this study as:
* We recommend standards such as ASHRAE 228, EU EPBD and IEA explicitly incentivize the use of season-hour or month-hour average emission factors for operational carbon accounting and disincentivize the use of annual averages (e.g., add a penalty to discourage this practice). Exceptionally, those emission factors may not be suitable for real-time applications due to weather synchronization. In these cases, actual emissions data (e.g., from [WattTime](https://docs.watttime.org/)) are required.
* We further recommend that those standards explicitly decouple net-zero carbon buildings (NZCBs) from net-zero energy buildings (NZEBs). This decoupling action should include: (1) clear separation in reporting between lifecycle impacts associated with operational energy use (module B6) and those associated with exported utilities (module D2); and similarly (2) incentivize the use of season-hour or month-hour emission factors to assess environmental impact of each kWh of exported utility.

