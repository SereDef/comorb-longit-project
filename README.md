# <ins>Dash application</ins>

**Longitudinal modelling of the co-development of depression and cardio-metabolic risk from childhood to young adulthood**

The study is based on data from the ALSPAC study.

This repo includes: 
  
  0. data preparation: i.e., extraction and transformation of all relevant variables.
  1. generalized cross-lagged panel models using lavaan
  2. (cross-lagged panel) network analysis using qgraph and psychonetrics (tbd)
  3. descriptive figures
  4. interctive dashboard to navigate results

## Dependencies

The python application:

| _Dependencies_            | version |
| ------------------------- | ------- |
| pandas                    | 2.0.3   |
| numpy                     | 1.25.2  |
| dash                      | 2.12.1  |
| dash_bootstrap_components | 1.4.2   |
| dash_cytoscape            | 0.3.0   |
| plotly                    | 5.16.1  |
| pyreadr                   | 0.4.9   |
