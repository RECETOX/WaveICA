# WaveICA
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/r-recetox-waveica/README.html)

Removal of batch effects for large-scale untargeted metabolomics data based on wavelet analysis.
The WaveICA R package provides a new algorithm to removing batch effects for metabolomics data. The details of this
method are in the papers "[WaveICA: a novel algorithm to remove batch effects for large-scale untargeted metabolomics
data based on wavelet analysis](https://doi.org/10.1016/j.aca.2019.02.010)" and "[WaveICA 2.0: a novel batch effect
removal method for untargeted metabolomics data without using batch information](https://doi.org/10.1007/s11306-021-01839-7)."

This fork is a slightly modified WaveICA by [RECETOX](https://github.com/RECETOX) that includes both original WaveICA
and WaveICA 2.0, and contains fewer dependencies.

You can install WaveICA Package Using the following commands:
    
    devtools::install_github("RECETOX/WaveICA",host="https://api.github.com")

To run WaveICA, you can use the following command:
    
    library(recetox.waveica)
    features <- waveica(features_table, wf, batch, factorization, group, K, t, t2, alpha)

Alternatively, if your data contains no batch information or the samples have been measured in a single batch you can run:

    library(recetox.waveica)
    features <- waveica_nonbatchwise(features_table, wf, injection_order, alpha, cutoff, K)

For further information on how to use the tool please refer to our docs in this repository.

When using **waveica**, please cite as: Deng K, Zhang F, Tan Q, Huang Y, Song W, Rong Z, Zhu ZJ, Li K, Li Z. 
WaveICA: A novel algorithm to remove batch effects for large-scale untargeted metabolomics data based on wavelet
analysis. Anal Chim Acta. 2019 Jul 11;1061:60-69. doi: 10.1016/j.aca.2019.02.010. Epub 2019 Feb 19. PMID: 30926040.

When using **waveica_nonbatchwise**, please cite as: Deng K, Zhao F, Rong Z, Cao L, Zhang L, Li K, Hou Y, Zhu ZJ.
WaveICA 2.0: a novel batch effect removal method for untargeted metabolomics data without using batch information.
Metabolomics. 2021 Sep 20;17(10):87. doi: 10.1007/s11306-021-01839-7. PMID: 34542717.

If you have any questions regarding this fork please contact RECETOX development team at [GalaxyToolsDevelopmentandDeployment@space.muni.cz](mailto:GalaxyToolsDevelopmentandDeployment@space.muni.cz)



