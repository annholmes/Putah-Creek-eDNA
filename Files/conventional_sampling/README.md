# 🐟 Conventional Fish Sampling Data (2020–2021)

This folder contains data from conventional sampling of fish assemblages in Putah Creek (2020-2021) that was compared to eDNA results in the paper *Spatiotemporal stability of fish communities in a regulated stream: Insights from environmental DNA*.


## 🗂️ Files

| Filename                         | Method              | Year | Location        | Season | Notes |
|----------------------------------|---------------------|------|-----------------|--------|-------|
| `PutahCk_electrofishing_comp_2020.csv`        | Electrofishing      | 2020 | Multiple sites  | Fall   | Partial overlap with eDNA sites |
| `PutahCk_electrofishing_comp_2021.csv`        | Electrofishing      | 2021 | Multiple sites  | Fall   | Partial overlap with eDNA sites |
| `PutahCk_fyke_trap_comp_2020.csv`             | Fyke Trap           | 2020 | Russell Ranch   | Spring | Single trap |
| `PutahCk_rotary_screw_trap_comp_2021.csv`     | Rotary Screw Trap   | 2021 | Step Weir 2     | Spring | Single trap |

### Data Fields

- **collecting_event**: Unique identifier for each sampling instance, constructed from:  
  • Site code (see [Table 1](../../Tables/Table1_Putah_Creek_sites.md)) and Jacinto et al. 2022, 2023)  
  • Season: `SP` (Spring) or `FA` (Fall)  
  • Year: `20` or `21`  
  • For rotary screw trap data, the month is included (`APR` or `MAY`) due to multiple sampling events  
    (e.g., `SW2DSP21APR`)

- **method**: Sampling method and measurement type. One of:  
  • `eDNA (reads)` – sequence reads from environmental DNA  
  • `Fyke (individuals)` – individual fish count from fyke trap  
  • `Fyke (biomass)` – biomass (g) from fyke trap  
  • `Rotary Screw Trap (individuals)` – individual fish count  
  • `Rotary Screw Trap (biomass)` – biomass (g)  
  • `Electrofishing (individuals)` – individual fish count; biomass was not measured

- **species**: Binomial name of the taxon

- **common_name**: Common name of the taxon

- **common_name_label**: Field used for labeling selected species in plots

- **names**: Formatted species label for figures and tables

- **names_by_family**: Label including family grouping, used for organizing species in legends

- **frequency**: Count of individuals (trap/electrofishing) or sequence reads (eDNA); biomass in grams (when available) 

- **relative_abundance**: Proportion within the sample (based on the `frequency` field)


## 📊 Data Collection and Funding

• Fyke trap and rotary screw trap data were collected and provided by Mackenzie Miner and Andrew Rypel.  
• Electrofishing data were collected and provided by Tim Salamunovich and TRPA Fish Biologists.  
• Multi-method sampling at site FRPG (fall 2021) is included in the electrofishing 2021 data file; these data were collected by the UC Davis Fish Class (WFC 120), and curated and provided by Teejay O'Rear and John Durand. The sampling methods were electrofishing, traps, gill nets, and beach seine.  
• Funding for annual conventional surveys was provided by the Solano County Water Agency.

## 📄 Methods and Related Work

### Electrofishing methods

**Jacinto, E., Fangue, N. A., Cocherell, D. E., Kiernan, J. D., Moyle, P. B., & Rypel, A. L.** (2022).  
*Increasing Stability of a Native Freshwater Fish Assemblage Following Flow Rehabilitation* [Dataset]. **Zenodo**. https://doi.org/10.5281/zenodo.7822308  
> This Zenodo archive includes electrofishing data from 1993–2017 used for analysis by Jacinto et al. (2023). The electrofishing data from 2020–2021 in this repository are part of the same long-term survey.

**Jacinto, E., Fangue, N. A., Cocherell, D. E., Kiernan, J. D., Moyle, P. B., & Rypel, A. L.** (2023).  
*Increasing stability of a native freshwater fish assemblage following flow rehabilitation.* *Ecological Applications, 33*(5), e2868.
https://doi.org/10.1002/eap.2868

### Trap methods

**Miner, M. C., Jacinto, E., Holmes, A., Cocherell, D. E., Baird, S. E., Schreier, A. M., Fangue, N. A., & Rypel, A. L.** (2020).  
*Origin and abundance of Chinook salmon in Putah Creek.* Annual report to the Solano County Water Agency, September 2020.

**Miner, M. C.** (2022).  
*Migratory phenology and spatial distributions of a recovering Chinook salmon run in a flow-regulated creek: Considerations for management.* MS Thesis, University of California, Davis.
https://escholarship.org/uc/item/93z321vk

**Rypel, A. L., Miner, M. C., Hitt, L., Holmes, A., Baird, S. E., Cocherell, D. E., Schreier, A. M., & Fangue, N. A.** (2022).  
*Origin and abundance of Chinook salmon in Putah Creek.* Annual report to the Solano County Water Agency 2021/2022.

## 📌 Notes

• Electrofishing data for sites not paired with eDNA are not provided here.  
• Some eDNA sites were not sampled using conventional methods.
