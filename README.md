# <ins>Dash application</ins>

**Longitudinal modelling of the co-development of depression and cardio-metabolic risk from childhood to young adulthood**

The study is based on data from the ALSPAC study.

This repo includes: 
  
  0. data preparation: i.e., extraction and transformation of all relevant variables.
  1. generalized cross-lagged panel models using lavaan
  2. (cross-lagged panel) network analysis using qgraph and psychonetrics (tbd)
  3. descriptive figures
  4. interctive dashboard to navigate results. The application is live at https://longit-comorbidity.onrender.com

## Dependencies

R scripts:

| _Dependencies_            | version | script |
| ------------------------- | ------- | ------ |
| foreign                   | 0.8-84  | 0      |
| lavaan                    | 0.6-16  | 1      |
| tidySEM                   | 0.2.4   | 1      |
| foreach                   | 1.5.2   | 1      |
| psychonetrics             | 0.11    | 2      |
| qgraph                    | 1.9.5   | 2      |
| dplyr                     | 1.1.2   | 2      |
| tidyr                     | 1.3.0   | 2      |
| psych                     | 2.3.6   | 2      |
| gdata                     | 2.19.0  | 2      |

Python application (see also requirements.txt):

| _Dependencies_            | version |
| ------------------------- | ------- |
| pandas                    | 2.1.1   |
| numpy                     | 1.25.2  |
| dash                      | 2.12.1  |
| dash_bootstrap_components | 1.4.2   |
| dash_cytoscape            | 0.3.0   |
| plotly                    | 5.16.1  |
| pyreadr                   | 0.4.9   |

## Running the application locally
It may be faster play around with the application on your local machine. To run the application locally (assuming a UNIX system with Python 3 installed):

```
git clone https://github.com/SereDef/comorb-longit-project
cd comorb-longit-project

pip install -r requirements.txt
python3 src/app.py
```
