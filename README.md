# Introduction

This package is used to analyze CRISPR editing data from chip synthesis. It is based on the following assumptions:
  - The synthesized plasmids can be very different from the designed plasmids. We classify the synthesized plasmids as functional and non-functional. Only functional plasmids can induce CRISPR editing.
  - The top common sythesized functional plasmids instead of the designed plasmids should be used as editing reference.
  - If a treat read has more than one reference functional plasmids (based on barcode), we distribute the read count as follows.
    - Normalize the reference count to get the priori distribution for read count.
    - Normalize the alignment score by a temperature and use softmax to calculate the conditional probability.
    - Compose the priori distribution and conditional probability to get the posteriori distribution of read count across all references.
  - Both functional and non-functional plasmids are transferred to treat samples.
  - With the functional plasmids as reference, the called mutants in treat samples comes either from edited functional plasmids or non-functional plasmids.
  - The abundance of non-functional plasmids is similar in treat and control samples. Therefore, one may substract the mutant frequency in control from that in treat. The remained mutants are expected to be edited reference functional plasmids.

The previous analysis piplines either use designed plasmids as reference or does not substract non-functional plasmids from total mutants. We compose both methods. We also use score as energy and apply an energy based method to distribute the treat read count to multiple references. We use the bioconda package rearr backend by an efficient and accurate chimeric alignment engine to call mutants from treat reads. Rearr allow us to discriminate templated insertions. 


# Install

```shell
$ pip install naapam
```

# TODO

- [ ] 文档
  - 用法说明
  - 例子
  - mermaid
- [ ] 发布包

# Dependencies

- bowtie2
- gawk
