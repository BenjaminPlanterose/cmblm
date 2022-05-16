# Combinatorial linear models: an imputation-free method for linear regression prediction on datasets with missing values 

#### Benjamin Planterose Jiménez<sup>1</sup>, Manfred Kayser<sup>1</sup>, Athina Vidaki<sup>1</sup>, Amke Caliebe<sup>2, 3</sup>

<sup>1</sup> Department of Genetic Identification, Erasmus MC University Medical Center Rotterdam, Rotterdam, the Netherlands.

<sup>2</sup> Institute of Medical Informatics and Statistics, Kiel University, Kiel, Germany.

<sup>3</sup> University Medical Centre Schleswig-Holstein, Kiel, Germany.


## Tested on:

    Operating system: Ubuntu 18.04.5 LTS
    R: R version 4.1.2 (2021-11-01) -- "Bird Hippie"

## Structure

All R-scripts employed are included and segmented as per sections in B. Planterose *et al* (**2022**).

  1) **Theory**: Implementations of the derived theoretical results
  
  2) **Simulations**: R-scripts employed, raw simulated data and summary statistics for all three simulation rounds. This includes 21 simulations 
for round 1, 6 simulations for round 2 and 6 simulations for round 3. Parameter space is described on *Table 1* in B. Planterose *et al* (**2022**).
  
  3) **Application**: R-script employed, processed data for EWAS datahub (Horvath´s CpG methylation + phenotypes: age and tissue), results from the missing value injection 
titration assays and the 4 trained models (2 cmb-lm, 2 cmb-ridge; see publication for details).
  

For any questions/concerns arising from this Github project, either [report an issue](https://github.com/BenjaminPlanterose/cmblm/issues) or simply contact b.planterosejimenez at erasmusmc.nl

### References and Supporting Information
B. Planterose *et al* (**2022**). Combinatorial linear models: an imputation-free method for linear regression prediction on datasets with missing values. Submitted to *Biostatistics*





