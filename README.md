# gCNV validation

[IntensityRankSumAnnotator](http://software.broadinstitute.org/software/genomestrip/org_broadinstitute_sv_annotation_IntensityRankSumAnnotator.html) tool is used to perform in-silico validation of Copy Number Variants (CNVs) in the UKBB dataset using SNP array intensity data.

## Usage

```
git clone
zip wdl/dependencies.zip wdl/*
cromshell submit wdl/ukbbValidation.wdl array-validation.json config.json wdl/dependencies.zip
```

### Contact and Credits

Copyright (c) 2021 Talkowski Lab and The Broad Institute of M.I.T. and Harvard  
Contact: [Alba Sanchis-Juan](mailto:asanchis-juan@mgh.harvard.edu)

SV aggregation team: Ryan Collins, Jack Fu, Isaac Wong, Alba Sanchis-Juan and Harrison Brand on behalf of the Talkowski Laboratory
