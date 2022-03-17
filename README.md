# gCNV validation

[IntensityRankSumAnnotator](http://software.broadinstitute.org/software/genomestrip/org_broadinstitute_sv_annotation_IntensityRankSumAnnotator.html) tool is used to perform in-silico validation of Copy Number Variants (CNVs) in the UKBB dataset using SNP array intensity data.

## Table of Contents
* [Requirements](#requirements)
* [Usage](#usage)
* [Contact and credits](#contact)

## <a name="requirements">Requirements</a>
### Deployment and execution:
* A [Google Cloud](https://cloud.google.com/) account.
* A workflow execution system supporting the [Workflow Description Language](https://openwdl.org/) (WDL), such as:
  * [Cromwell](https://github.com/broadinstitute/cromwell) (v36 or higher). A dedicated server is highly recommended.
* [cromshell](https://github.com/broadinstitute/cromshell) for interacting with a dedicated Cromwell server.

### Data:
* Directory path with bed files containing UKBB gCNV output (per chromosome).
  * VCF header template
* List of UKBB SNP array files in VCF format: 
* List of samples on which to run GenomeStrip IRS

## <a name="usage">Usage</a>
The main scripts to run this analysis are:
* `ukbbValidation.wdl`: this workflow reformats SNP array and gCNV data from the UKBB and calls GenomeStrip IRS for <i>in-silico</i> CNV validation.   
* `genomeStripIRS.wdl`: runs GenomeStrip IRS and can be executed on its own. 

### Execution

```
> git clone https://github.com/talkowski-lab/cnv-validation.git
> cd cnv-validation/wdl
> zip dependencies.zip *

> cromshell submit ukbbValidation.wdl /path/to/array-validation.json /path/to/config.json dependencies.zip
```

## <a name="contact">Contact and credits</a>

Copyright (c) 2022 Talkowski Lab and The Broad Institute of M.I.T. and Harvard  
Contact: [Alba Sanchis-Juan](mailto:asanchis-juan@mgh.harvard.edu)

SV aggregation team: Ryan Collins, Jack Fu, Isaac Wong, Alba Sanchis-Juan and Harrison Brand on behalf of the Talkowski Laboratory
