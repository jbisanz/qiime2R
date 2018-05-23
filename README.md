# qiime2R
A package for importing qiime artifacts into an R session.

### Background
The [qiime artifact](https://docs.qiime2.org/2018.4/concepts/#data-files-qiime-2-artifacts) is a method for storing the input and outputs for [QIIME2](https://qiime2.org/) along with associated metadata and provenance information about how the object was formed. This method of storing objects has a number of obvious advantages; however, on the surface it does not lend itself to easy import to R for the R-minded data scientist. In reality, the .qza file is a compressed directory with an intuitive structure.

While it is possible to export data on a one-by-one basis from the qiime artifacts using qiime's built in suite of export features this is problematic and runs antithetical to the purpose of the artifact format for the following reasons:

* Export of data from the artifact using QIIME2 requires an installation which may not be available on the user's computer and may not be trivial to install for a novice user
* Export of the data will loose the associated provenance information. Now the origin of the data can't be traced and the parameters that led to its generation has been lost.
* Export of the data on a one-by-one basis is tedious and creates multiple copies of intermediate files
* R has many options for advanced data analysis and/or visualization that may not natively supported in QIIME or python environments

This package is trying to simplify the process of getting the artifact into R without discarding any of the associated data through a simple `read_qza` function. The artifact is unpacked in to /tmp (or another directory if specified using `tmp="/yourdirhere`) and the raw data and associated metadata are read into a named list (see below). The object is then removed from the tmp dir (unless user specifies `rm=F`). Data are typically returned as either a matrix, data.frame, phylo object (trees), or DNAStringSets (nucleic acid sequences).

### Dependencies
* [biomformat](https://www.bioconductor.org/packages/release/bioc/html/biomformat.html)
* [yaml](https://cran.r-project.org/web/packages/yaml/index.html)
* [BioStrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html)
* [Phyloseq](https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html)
* [Hmisc](https://www.google.com/search?q=Hmisc&rlz=1C5CHFA_enIN774IN774&oq=Hmisc&aqs=chrome.0.69i59j69i60l3j0l2.628j1j9&sourceid=chrome&ie=UTF-8)

### To do
Include a function for visualizing and/or summarizing providence in an intelligble format.

### Known issues
The artifact metadata import may fail on windows PCs due to an issue parsing `.yaml` manifests.

Please report any issues using the issues tab above. Community contributions to coding would be appreciated!


### Example usage

This example is using data derived from the [moving pictures tutorial](https://docs.qiime2.org/2018.4/tutorials/moving-pictures/).
```
?read_qza
Usage:

     read_qza(file, tmp, rm)
     
Arguments:

    file: path to the input file, ex:
          file="~/data/moving_pictures/table.qza"

     tmp: a temporary directory that the object will be decompressed to
          (default="/tmp")

      rm: should the decompressed object be removed at completion of
          function (T/F default=T)

Value:

     a named list of the following objects:

        • artifact$data - the raw data ex OTU table as matrix or tree
          in phylo format

        • artifact$uuid - the unique identifer of the artifact

        • artifact$type - the semantic type of the object (ex
          FeatureData[Sequence])

        • artifact$format - the format of the qiime artifact

        • artifact$provenance - information tracking how the object was
          created

        • artifact$contents - a table of all the files contained within
          the artifact and their file size

        • artifact$version - the reported version for the artifact, a
          warning error may be thrown if a new version is seen


```

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
The provenance is supplied as a nested set of lists but it can be difficult to interpret at this point: `otus$provenance`. A function for printing the provenance is `print_provenance()`

```
print_provenance(otus)
artifact$provenance = list 3 (103880 bytes)
.  6a560288-898e-4c1d-92ac-dd8d7822dcc9/provenance/action/action.yaml = list 4
. .  execution = list 2
. . .  uuid = character 1= 70ea1d17-ffcf-4113-8ce 
. . .  runtime = list 3
. . . .  start = character 1= 2018-04-26T09:43:13.48 
. . . .  end = character 1= 2018-04-26T09:47:32.33 
. . . .  duration = character 1= 4 minutes, 18 seconds, 
...
```

If you want to build a phyloseq object for further analysis you could use the function `qza_to_phyloseq()`
```
physeq<-qza_to_phyloseq("table.qza","rooted-tree.qza","taxonomy.qza", "sample-metadata.tsv")
physeq
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 759 taxa and 34 samples ]
## sample_data() Sample Data:       [ 34 samples by 10 sample variables ]
## tax_table()   Taxonomy Table:    [ 759 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 759 tips and 757 internal nodes ]
```

OR you could manually build it:

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
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 759 taxa and 34 samples ]
## sample_data() Sample Data:       [ 34 samples by 10 sample variables ]
## tax_table()   Taxonomy Table:    [ 759 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 759 tips and 757 internal nodes ]
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
