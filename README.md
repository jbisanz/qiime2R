# qiime2R
A package for importing qiime artifacts into an R session.

### Background
The [qiime artifact](https://docs.qiime2.org/2018.4/concepts/#data-files-qiime-2-artifacts) is a method for storing the input and outputs for [QIIME2](https://qiime2.org/) along with associated metadata and provenance information about how the object was formed. This method of storing objects has a number of obvious advantages; however, on the surface it does not lend itself to easy import to R for the R-minded data scientist. In reality, the .qza file is a compressed directory with an intuitive structure. This package is trying to simplify the process through a simple `read_qza` function. The artifact is unpacked in to /tmp (or another directory if specified using `tmp="/yourdirhere`) and the raw data and associated metadata are read into a named list (see below). The object is then removed from the tmp dir (unless user specifies `rm=F`). Data are typically returned as either a matrix, data.frame, phylo object (trees), or DNAStringSets (nucleic acid sequences).

### Dependencies
* [biomformat](https://www.bioconductor.org/packages/release/bioc/html/biomformat.html)
* [yaml](https://cran.r-project.org/web/packages/yaml/index.html)
* [BioStrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html)

### Supported semantic types

* FeatureTable[Balance]
* FeatureTable[Composition]
* FeatureTable[Frequency]
* FeatureTable[PresenceAbsence]
* FeatureTable[RelativeFrequency]
* Phylogeny[Rooted]
* Phylogeny[Unrooted]
* DistanceMatrix
* DeblurStats
* QualityFilterStats
* FeatureData[Taxonomy]
* PCoAResults
* FeatureData[AlignedSequence]
* FeatureData[Sequence]

### To do
Improve support for more/new semanatic types and include a function for visualizing and/or summarizing providence in an intelligble format.

### Example usage

This example is using data derived from the [moving pictures tutorial](https://docs.qiime2.org/2018.4/tutorials/moving-pictures/).

```
otus<-read_qza("~/QIIME2/mvpics/table.qza")
```

To access the actual data stored within the object, use the .$data

```
otus$data[1:5,1:5] #show first 5 samples and first 5 taxa
                                 L1S105 L1S140 L1S208 L1S257 L1S281
4b5eeb300368260019c1fbc7a3c718fc   2222      0      0      0      0
fe30ff0f71a38a39cf1717ec2be3a2fc      5      0      0      0      0
d29fe3c70564fc0f69f2c03e0d1e5561      0      0      0      0      0
868528ca947bc57b69ffdf83e6b73bae      0   2276   2156   1205   1772
154709e160e8cada6bfb21115acc80f5    812   1176    713    407    242
```

We can also get information about the object:
```
names(otus)
[1] "uuid"       "type"       "format"     "contents"   "version"   
[6] "data"       "provenance"
```

```
otus$uuid
[1] "6a560288-898e-4c1d-92ac-dd8d7822dcc9"
otus$type
[1] "FeatureTable[Frequency]"
```

We can also get a complete list of the files within the artifact and their sizes:
```
files	size
6a560288-898e-4c1d-92ac-dd8d7822dcc9/metadata.yaml	96
6a560288-898e-4c1d-92ac-dd8d7822dcc9/VERSION	39
6a560288-898e-4c1d-92ac-dd8d7822dcc9/data/feature-table.biom	86181
6a560288-898e-4c1d-92ac-dd8d7822dcc9/provenance/citations.bib	2387
6a560288-898e-4c1d-92ac-dd8d7822dcc9/provenance/metadata.yaml	96
6a560288-898e-4c1d-92ac-dd8d7822dcc9/provenance/VERSION	39
6a560288-898e-4c1d-92ac-dd8d7822dcc9/provenance/action/action.yaml	4969
6a560288-898e-4c1d-92ac-dd8d7822dcc9/provenance/artifacts/640f363f-519b-4152-9664-07f9a33e9f25/citations.bib	856
6a560288-898e-4c1d-92ac-dd8d7822dcc9/provenance/artifacts/640f363f-519b-4152-9664-07f9a33e9f25/metadata.yaml	130
6a560288-898e-4c1d-92ac-dd8d7822dcc9/provenance/artifacts/640f363f-519b-4152-9664-07f9a33e9f25/VERSION	39
6a560288-898e-4c1d-92ac-dd8d7822dcc9/provenance/artifacts/640f363f-519b-4152-9664-07f9a33e9f25/action/action.yaml	4650
6a560288-898e-4c1d-92ac-dd8d7822dcc9/provenance/artifacts/640f363f-519b-4152-9664-07f9a33e9f25/action/barcodes.tsv	756
6a560288-898e-4c1d-92ac-dd8d7822dcc9/provenance/artifacts/bbb00eed-8dec-475e-a07b-4a527a2f1863/citations.bib	856
6a560288-898e-4c1d-92ac-dd8d7822dcc9/provenance/artifacts/bbb00eed-8dec-475e-a07b-4a527a2f1863/metadata.yaml	98
6a560288-898e-4c1d-92ac-dd8d7822dcc9/provenance/artifacts/bbb00eed-8dec-475e-a07b-4a527a2f1863/VERSION	39
6a560288-898e-4c1d-92ac-dd8d7822dcc9/provenance/artifacts/bbb00eed-8dec-475e-a07b-4a527a2f1863/action/action.yaml	4220
```
The provenance is supplied as a nested set of lists but it can be difficult to interpret at this point: `otus$provenance`.

Here is how it could be imported to form a phyloseq object
```
library(phyloseq)
tree<-read_qza("~/QIIME2/mvpics/rooted-tree.qza")

taxonomy<-read_qza("~/QIIME2/mvpics/taxonomy.qza")
  tax_table<-do.call(rbind, strsplit(as.character(taxonomy$data$Taxon), "; "))
  colnames(tax_table)<-c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  rownames(tax_table)<-taxonomy$data$Feature.ID
  
metadata<-read.table("~/QIIME2/mvpics/sample-metadata.tsv", sep='\t', header=T, row.names=1, comment="")
metadata<-metadata[-1,]#remove the second line that specifies the data type

physeq<-phyloseq(otu_table(otus$data, taxa_are_rows = T), phy_tree(tree$data), tax_table(tax_table), sample_data(metadata))

physeq
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 759 taxa and 34 samples ]
sample_data() Sample Data:       [ 34 samples by 10 sample variables ]
tax_table()   Taxonomy Table:    [ 759 taxa by 7 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 759 tips and 757 internal nodes ]
```

You could then do whatever analysis you wanted (in this case we are ignoring normalization/scaling considerations):
```
plot_ordination(physeq, ordinate(physeq, "DCA"), type="samples", color="BodySite")
```
![physeq](https://github.com/jbisanz/qiime2R/raw/master/images/physeq.png)

You could also try [MicrobeR](https://github.com/jbisanz/MicrobeR)
```
devtools::install_github("jbisanz/MicrobeR")
Microbiome.Barplot(Summarize.Taxa(otus$data, as.data.frame(tax_table))$Family, metadata, CATEGORY="BodySite")
```
![microber](https://github.com/jbisanz/qiime2R/raw/master/images/microber.png)
