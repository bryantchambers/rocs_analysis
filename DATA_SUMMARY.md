# Data Summary for Taxon Importance Analysis

## Prokaryote VST Data
- Dimensions: 214 x 1797 (Samples x Taxa)
- Range of values: -9.59834453381046 to 11.3084352487707 
- Any NAs: FALSE 

## WGCNA Module Assignments
- Total taxa assigned: 1797 
- Module Distribution:



|module    |    N|
|:---------|----:|
|grey      | 1196|
|turquoise |  237|
|blue      |  138|
|brown     |   85|
|yellow    |   76|
|green     |   65|

## HMM States
- Total samples with state labels: 214 
- State Label Distribution:



|label |  N|
|:-----|--:|
|G-A   | 70|
|IG-B  | 53|
|IG-C  | 35|
|IG-E  | 56|

## Sample Alignment
- Samples in both VST and HMM: 214 
- Samples missing from VST: 0 
- Samples missing from HMM: 0 

## FuzzyForest Requirements
- Features with zero variance:
  - None

## Taxon ID Consistency
- Intersection VST & Modules: 1797 
- Intersection VST & Metadata: 1797 
- Intersection VST & Raw (subspecies): 1797 

## Fuzzy Forest Results
- Model Accuracy: 0.9032258 
- Model Kappa: 0.8673797 

- Top 5 Fuzzy Forest Drivers:



|taxon              | variable_importance|module    |
|:------------------|-------------------:|:---------|
|S__GCA_905182465.1 |           0.0457578|green     |
|S__3300020317_19   |           0.0383628|green     |
|S__GCA_905478185.1 |           0.0362573|turquoise |
|S__GCA_012964335.1 |           0.0351686|turquoise |
|S__3300027779_22   |           0.0330220|turquoise |

## Network Statistics Results
- Total taxa analyzed: 1797 
- Potential keystones identified: 50 
- Importance Type Distribution:



|importance_type |    N|
|:---------------|----:|
|Peripheral      | 1280|
|Connector       |  436|
|Keystone        |   50|
|Hub             |   27|
|Hidden Gem      |    4|

- Top 5 Keystone Taxa:



|taxon                     |module    |functional_group             |species                      |
|:-------------------------|:---------|:----------------------------|:----------------------------|
|S__GCA_018672475.1        |turquoise |Core_heterotrophy            |s__GCA-2724215 sp018672475   |
|S__GCA_012960835.1        |turquoise |Particle_heterotroph         |s__DUAL01 sp012960835        |
|S__TARA_ARC_108_MAG_00128 |turquoise |Particle_heterotroph         |s__GCA-002697625 sp905478555 |
|S__GCA_905479195.1        |yellow    |Sulfate_reduction_diagenetic |s__GCA-2698665 sp905479195   |
|S__GCA_014237975.1        |turquoise |Particle_heterotroph         |s__Mariniblastus sp014237975 |

