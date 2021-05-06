# Tutorial: Integrating QIIME2 and R for data visualization and analysis using qiime2R (*v0.99.6*)

## Background
The [qiime artifact](https://docs.qiime2.org/2018.4/concepts/#data-files-qiime-2-artifacts) is a method for storing the input and outputs for [QIIME2](https://qiime2.org/) along with associated metadata and provenance information about how the object was formed. This method of storing objects has a number of obvious advantages; however, on the surface it does not lend itself to easy import to R for the R-minded data scientist. In reality, the .qza file is a compressed directory with an intuitive structure.

While it is possible to export data on a one-by-one basis from the qiime artifacts using qiime's built in suite of export features this is problematic and runs antithetical to the purpose of the artifact format for the following reasons:

* Export of data from the artifact using QIIME2 requires an installation which may not be available on the user's computer and may not be trivial to install for a novice user
* Export of the data will loose the associated provenance information. Now the origin of the data can't be traced and the parameters that led to its generation have been lost.
* Export of the data on a one-by-one basis is tedious and creates multiple copies of intermediate files
* R has many options for advanced data analysis and/or visualization that may not natively supported in QIIME or python environments

This package is trying to simplify the process of getting the artifact into R without discarding any of the associated data through a simple `read_qza` function. The artifact is unpacked in to a temporary directory and the raw data and associated metadata are read into a named list (see below). Data are typically returned as either a data.frame, phylo object (trees), or DNAStringSets (nucleic acid sequences).

## Functions

* `read_qza()` - Function for reading artifacts (.qza). 
* `qza_to_phyloseq()` - Imports multiple artifacts to produce a phyloseq object.
* `read_q2metadata()` - Reads qiime2 metadata file (containing q2-types definition line)
* `write_q2manifest()` - Writes a read manifest file to import data into qiime2
* `theme_q2r()` - A ggplot2 theme for for clean figures.
* `print_provenance()` - A function to display provenance information.
* `is_q2metadata()` - A function to check if a file is a qiime2 metadata file.
* `parse_taxonomy()` - A function to parse taxonomy strings and return a table where each column is a taxonomic class.
* `parse_ordination()` - A function to parse the internal ordination format.
* `read_q2biom()` - A function for reading QIIME2 biom files in format v2.1
* `make_clr()` - Transform feature table using centered log2 ratio.
* `make_proportion()` - Transform feature table to proportion (sum to 1).
* `make_percent()` - Transform feature to percent (sum to 100).
* `interactive_table()` - Create an interactive table in Rstudio viewer or rmarkdown html.
* `summarize_taxa()`- Create a list of tables with abundances sumed to each taxonomic level.
* `taxa_barplot()` - Create a stacked barplot using ggplot2.
* `taxa_heatmap()` - Create a heatmap of taxonomic abundances using gplot2.
* `corner()` - Show top corner of a large table-like obejct.
* `min_nonzero()` - Find the smallest non-zero, non-NA in a numeric vector.
* `mean_sd()` - Return mean and standard deviation for plotting.
* `subsample_table()` - Subsample a table with or without replacement.
* `filter_features()` - Remove low abundance features by number of counts and number of samples they appear in.

## Installing qiime2R

qiime2R is currently available via github which can easily be installed in R via the following command:
```
if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R")
```

***

## Reading Artifacts (.qza)

This example is using data derived from the [moving pictures tutorial](https://docs.qiime2.org/2020.2/tutorials/moving-pictures/). The main function we will need is `read_qza()`:

```
?read_qza
read qiime2 artifacts (.qza)

Description

  extracts embedded data and object metadata into an R session

Usage

  read_qza(file, tmp, rm)

Arguments

  file - path to the input file, ex: file="~/data/moving_pictures/table.qza"
  tmp - a temporary directory that the object will be decompressed to (default="tempdir()")
  rm - should the decompressed object be removed at completion of function (T/F default=TRUE)

Value
  
  a named list of the following objects:
    artifact$data - the raw data ex OTU table as matrix or tree in phylo format
    artifact$uuid - the unique identifer of the artifact
    artifact$type - the semantic type of the object (ex FeatureData[Sequence])
    artifact$format - the format of the qiime artifact
    artifact$provenance - information tracking how the object was created
    artifact$contents - a table of all the files contained within the artifact and their file size
    artifact$version - the reported version for the artifact, a warning error may be thrown if a new version is seen
```

We will start be reading in a table of sequence variants (SVs):
```
SVs<-read_qza("table.qza")
```

When the artifact is imported, there are a number of pieces of information included. To see them we can use the names command:
```
names(SVs)
[1] "uuid"       "type"       "format"     "contents"   "version"   
[6] "data"       "provenance"
```

To access the actual data stored within the object, access the data as below:

```
SVs$data[1:5,1:5] #show first 5 samples and first 5 taxa
#                                 L1S105 L1S140 L1S208 L1S257 L1S281
#4b5eeb300368260019c1fbc7a3c718fc   2183      0      0      0      0
#fe30ff0f71a38a39cf1717ec2be3a2fc      5      0      0      0      0
#d29fe3c70564fc0f69f2c03e0d1e5561      0      0      0      0      0
#868528ca947bc57b69ffdf83e6b73bae      0   2249   2117   1191   1737
#154709e160e8cada6bfb21115acc80f5    802   1174    694    406    242
```

In the above file each row denotes a sequence variant where in the actual text is the hash of the complete sequence. See [here](https://www.md5hashgenerator.com/) for an example tool for generating hashes.

We can also look at the unique identifier for this object:
```
SVs$uuid
[1] "706b6bce-8f19-4ae9-b8f5-21b14a814a1b"
```

We can see the type of artifact:
```
SVs$type
[1] "FeatureTable[Frequency]"
```

We can also get a complete list of the files within the artifact and their sizes:
```
SVs$contents
#                                                                                                           files.Name files.Length          files.Date size
#1                                                                  706b6bce-8f19-4ae9-b8f5-21b14a814a1b/metadata.yaml           96 2020-02-28 10:21:00  224
#2                                                                  706b6bce-8f19-4ae9-b8f5-21b14a814a1b/checksums.md5         1341 2020-02-28 10:21:00  224
#3                                                                        706b6bce-8f19-4ae9-b8f5-21b14a814a1b/VERSION           39 2020-02-28 10:21:00  224
#4                                                       706b6bce-8f19-4ae9-b8f5-21b14a814a1b/provenance/metadata.yaml           96 2020-02-28 10:21:00  224
#5                                                       706b6bce-8f19-4ae9-b8f5-21b14a814a1b/provenance/citations.bib         4305 2020-02-28 10:21:00  224
#6                                                             706b6bce-8f19-4ae9-b8f5-21b14a814a1b/provenance/VERSION           39 2020-02-28 10:21:00  224
#7        706b6bce-8f19-4ae9-b8f5-21b14a814a1b/provenance/artifacts/4de0fc23-6462-43d3-8497-f55fc49f5db6/metadata.yaml          130 2020-02-28 10:20:00  224
#8        706b6bce-8f19-4ae9-b8f5-21b14a814a1b/provenance/artifacts/4de0fc23-6462-43d3-8497-f55fc49f5db6/citations.bib         3488 2020-02-28 10:20:00  224
#9              706b6bce-8f19-4ae9-b8f5-21b14a814a1b/provenance/artifacts/4de0fc23-6462-43d3-8497-f55fc49f5db6/VERSION           39 2020-02-28 10:20:00  224
#10  706b6bce-8f19-4ae9-b8f5-21b14a814a1b/provenance/artifacts/4de0fc23-6462-43d3-8497-f55fc49f5db6/action/action.yaml         5473 2020-02-28 10:20:00  224
#11 706b6bce-8f19-4ae9-b8f5-21b14a814a1b/provenance/artifacts/4de0fc23-6462-43d3-8497-f55fc49f5db6/action/barcodes.tsv          757 2020-02-28 10:20:00  224
#12       706b6bce-8f19-4ae9-b8f5-21b14a814a1b/provenance/artifacts/f5d67104-9506-4373-96e2-97df9199a719/metadata.yaml           98 2020-02-28 10:20:00  224
#13       706b6bce-8f19-4ae9-b8f5-21b14a814a1b/provenance/artifacts/f5d67104-9506-4373-96e2-97df9199a719/citations.bib         2774 2020-02-28 10:20:00  224
#14             706b6bce-8f19-4ae9-b8f5-21b14a814a1b/provenance/artifacts/f5d67104-9506-4373-96e2-97df9199a719/VERSION           39 2020-02-28 10:20:00  224
#15  706b6bce-8f19-4ae9-b8f5-21b14a814a1b/provenance/artifacts/f5d67104-9506-4373-96e2-97df9199a719/action/action.yaml         4893 2020-02-28 10:20:00  224
#16                                                 706b6bce-8f19-4ae9-b8f5-21b14a814a1b/provenance/action/action.yaml         5617 2020-02-28 10:21:00  224
#17                                                       706b6bce-8f19-4ae9-b8f5-21b14a814a1b/data/feature-table.biom       114900 2020-02-28 10:21:00  224
```

We can also print the providence; however, it is probably easier to use the [q2view tool](https://view.qiime2.org/) for a graphical aid in its interpretation.
```
print_provenance(SVs)
artifact$provenance = list 3 (118072 bytes)
.  706b6bce-8f19-4ae9-b8f5-21b14a814a1b/provenance/artifacts/4de0fc23-6462-43d3-8497-f55fc49f5db6/action/action.yaml = list 4
. .  execution = list 2
. . .  uuid = character 1= 8614aa95-3638-4e49-88f 
. . .  runtime = list 3
. . . .  start = character 1= 2020-02-28T10:19:02.86 
. . . .  end = character 1= 2020-02-28T10:19:24.46 
. . . .  duration = character 1= 21 seconds, and 605755 
...
```

***

## Reading Metadata

If you are using a qiime2 metadata file (outlined [here](https://docs.qiime2.org/2020.2/tutorials/metadata/)), you can use the supplied function `read_q2metadata()` as below. If using a standard tsv or csv file, you can use `read.table()`, `readr::read_tsv()`, or `readr::csv()`.

```
metadata<-read_q2metadata("sample-metadata.tsv")
head(metadata) # show top lines of metadata
#  SampleID barcode-sequence body-site year month day   subject reported-antibiotic-usage days-since-experiment-start
#2     L1S8     AGCTGACTAGTC       gut 2008    10  28 subject-1                       Yes                           0
#3    L1S57     ACACACTATGGC       gut 2009     1  20 subject-1                        No                          84
#4    L1S76     ACTACGTGTGGT       gut 2009     2  17 subject-1                        No                         112
#5   L1S105     AGTGCGATGCGT       gut 2009     3  17 subject-1                        No                         140
#6   L2S155     ACGATGCGACCA left palm 2009     1  20 subject-1                        No                          84
#7   L2S175     AGCTATCCACGA left palm 2009     2  17 subject-1                        No                         112
```

***

## Reading Taxonomy

Note, when taxonomy is imported, a single string is returned along with a confidence score. For many analysis we will want to break up this string and for that purpose the `parse_taxonomy()` function is provided:
```
taxonomy<-read_qza("taxonomy.qza")
head(taxonomy$data)
#                        Feature.ID                                                                                                                            Taxon Confidence
#1 4b5eeb300368260019c1fbc7a3c718fc                          k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Bacteroides; s__  0.9972511
#2 fe30ff0f71a38a39cf1717ec2be3a2fc                           k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Neisseriales; f__Neisseriaceae; g__Neisseria  0.9799427
#3 d29fe3c70564fc0f69f2c03e0d1e5561                                k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Streptococcaceae; g__Streptococcus  1.0000000
#4 868528ca947bc57b69ffdf83e6b73bae                          k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Bacteroides; s__  0.9955859
#5 154709e160e8cada6bfb21115acc80f5                               k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Bacteroides  1.0000000
#6 1d2e5f3444ca750c85302ceee2473331 k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pasteurellales; f__Pasteurellaceae; g__Haemophilus; s__parainfluenzae  0.9455365
taxonomy<-parse_taxonomy(taxonomy$data)
head(taxonomy)
#                                  Kingdom         Phylum               Class           Order           Family         Genus        Species
#4b5eeb300368260019c1fbc7a3c718fc Bacteria  Bacteroidetes         Bacteroidia   Bacteroidales   Bacteroidaceae   Bacteroides           <NA>
#fe30ff0f71a38a39cf1717ec2be3a2fc Bacteria Proteobacteria  Betaproteobacteria    Neisseriales    Neisseriaceae     Neisseria           <NA>
#d29fe3c70564fc0f69f2c03e0d1e5561 Bacteria     Firmicutes             Bacilli Lactobacillales Streptococcaceae Streptococcus           <NA>
#868528ca947bc57b69ffdf83e6b73bae Bacteria  Bacteroidetes         Bacteroidia   Bacteroidales   Bacteroidaceae   Bacteroides           <NA>
#154709e160e8cada6bfb21115acc80f5 Bacteria  Bacteroidetes         Bacteroidia   Bacteroidales   Bacteroidaceae   Bacteroides           <NA>
#1d2e5f3444ca750c85302ceee2473331 Bacteria Proteobacteria Gammaproteobacteria  Pasteurellales  Pasteurellaceae   Haemophilus parainfluenzae
```
***

## Creating a Phyloseq Object

A wrapper function called `qza_to_phyloseq()` is provided which links multiple `read_qza()` calls together to create a phyloseq object for subsequent analysis as per the [phyloseq tutorials](https://joey711.github.io/phyloseq/). An example usage is shown below:

```
physeq<-qza_to_phyloseq(
    features="inst/artifacts/2020.2_moving-pictures/table.qza",
    tree="inst/artifacts/2020.2_moving-pictures/rooted-tree.qza",
    taxonomy="inst/artifacts/2020.2_moving-pictures/taxonomy.qza",
    metadata = "inst/artifacts/2020.2_moving-pictures/sample-metadata.tsv"
    )
physeq
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 759 taxa and 34 samples ]
## sample_data() Sample Data:       [ 34 samples by 10 sample variables ]
## tax_table()   Taxonomy Table:    [ 759 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 759 tips and 757 internal nodes ]
```

***

## Example Visualizations

In this next section, methods for publication-ready figures will be demonstrated through extensive use of the [tidyverse](https://www.tidyverse.org/), and particularly ggplot2 and its many extensions. If starting with R, learning these packages early on will have a bit of a learning curve, but pay off massively in the long term. Each section contains a self-contained code snippets capable of replicating the figures and links to the data used. A lot can be learned from modifying parameters of the graphing commands to customize and play around with the visualizations. Note, the pipe operator (%>%) is being used frequently here and passes the output of one line to the next. Learn more about piping in R [here](https://www.datacamp.com/community/tutorials/pipe-r-tutorial).

All the data we will use is from the "Moving Pictures" tutorial and/or related content which can be found [here](https://docs.qiime2.org/2020.2/tutorials/moving-pictures/).

### Alpha Diversity Over Time

Data:
* [Metadata](https://data.qiime2.org/2020.2/tutorials/moving-pictures/sample_metadata.tsv)
* [Shannon Diversity](https://docs.qiime2.org/2020.2/data/tutorials/moving-pictures/core-metrics-results/shannon_vector.qza)

```
library(tidyverse)
library(qiime2R)

metadata<-read_q2metadata("sample-metadata.tsv")
shannon<-read_qza("shannon_vector.qza")

shannon<-shannon$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged
```
This dataset provides a good example of why you should always look at your data. There are an extra 3 samples in the metadata which do not have an assigned Shannon diversity value as seen below:
```
gplots::venn(list(metadata=metadata$SampleID, shannon=shannon$SampleID))
```
![sol](https://github.com/jbisanz/qiime2R/raw/master/images/Shannon_overlap.png)

If this was your dataset, you should look into why this was and try to see why these samples have been lost. We will continue after sorting them out:

```
metadata<-
metadata %>% 
  left_join(shannon)
head(metadata)
#  SampleID barcode-sequence body-site year month day   subject reported-antibiotic-usage days-since-experiment-start  shannon
#1     L1S8     AGCTGACTAGTC       gut 2008    10  28 subject-1                       Yes                           0 3.188316
#2    L1S57     ACACACTATGGC       gut 2009     1  20 subject-1                        No                          84 3.985702
#3    L1S76     ACTACGTGTGGT       gut 2009     2  17 subject-1                        No                         112 3.461625
#4   L1S105     AGTGCGATGCGT       gut 2009     3  17 subject-1                        No                         140 3.972339
#5   L2S155     ACGATGCGACCA left palm 2009     1  20 subject-1                        No                          84 5.064577
#6   L2S175     AGCTATCCACGA left palm 2009     2  17 subject-1                        No                         112 4.966707
```
Note how we have now added a column for the shannon diversity to our table which will make everything easy to plot. Given this dataset is longitudinal we can plot it as such. Note: the code I am using here would handle replicate data and plot error bars; however, this case does not have replicates. Also, a quick note is that R does generally not like `-`s in column names. To use these names we can surround them with ticks.

```
metadata %>%
  filter(!is.na(shannon)) %>%
  ggplot(aes(x=`days-since-experiment-start`, y=shannon, color=`body-site`)) +
  stat_summary(geom="errorbar", fun.data=mean_se, width=0) +
  stat_summary(geom="line", fun.data=mean_se) +
  stat_summary(geom="point", fun.data=mean_se) +
  xlab("Days") +
  ylab("Shannon Diversity") +
  theme_q2r() + # try other themes like theme_bw() or theme_classic()
  scale_color_viridis_d(name="Body Site") # use different color scale which is color blind friendly
  ggsave("Shannon_by_time.pdf", height=3, width=4, device="pdf") # save a PDF 3 inches by 4 inches
```
![sbt](https://github.com/jbisanz/qiime2R/raw/master/images/Shannon_by_time.png)

Imagine we wanted to instead ask if there was an effect of antibiotics on diversity. We can change the type of plot:

```
metadata %>%
  filter(!is.na(shannon)) %>%
  ggplot(aes(x=`reported-antibiotic-usage`, y=shannon, fill=`reported-antibiotic-usage`)) +
  stat_summary(geom="bar", fun.data=mean_se, color="black") + #here black is the outline for the bars
  geom_jitter(shape=21, width=0.2, height=0) +
  coord_cartesian(ylim=c(2,7)) + # adjust y-axis
  facet_grid(~`body-site`) + # create a panel for each body site
  xlab("Antibiotic Usage") +
  ylab("Shannon Diversity") +
  theme_q2r() +
  scale_fill_manual(values=c("cornflowerblue","indianred")) + #specify custom colors
  theme(legend.position="none") #remove the legend as it isn't needed
  ggsave("../../../images/Shannon_by_abx.pdf", height=3, width=4, device="pdf") # save a PDF 3 inches by 4 inches
```
![sba](https://github.com/jbisanz/qiime2R/raw/master/images/Shannon_by_abx.png)

Maybe instead the question is if the two individuals have different diversity?

```
metadata %>%
  filter(!is.na(shannon)) %>%
  ggplot(aes(x=subject, y=shannon, fill=`subject`)) +
  stat_summary(geom="bar", fun.data=mean_se, color="black") + #here black is the outline for the bars
  geom_jitter(shape=21, width=0.2, height=0) +
  coord_cartesian(ylim=c(2,7)) + # adjust y-axis
  facet_grid(~`body-site`) + # create a panel for each body site
  xlab("Antibiotic Usage") +
  ylab("Shannon Diversity") +
  theme_q2r() +
  scale_fill_manual(values=c("cornflowerblue","indianred")) + #specify custom colors
  theme(legend.position="none") #remove the legend as it isn't needed
  ggsave("Shannon_by_person.pdf", height=3, width=4, device="pdf") # save a PDF 3 inches by 4 inches
```
![sbp](https://github.com/jbisanz/qiime2R/raw/master/images/Shannon_by_person.png)

Note: It is important to remember that these a graphical tools, but to make a claim, you should also generate some form of statistical support taking the nature of the data into consideration (ie that it uses repeated sampling over time).

### Plotting PCoA

Here I am going to demo a PCoA with extra metadata mapped to it. I indicate the body site, antibiotic usage, and the shannon diversity by mapping to different aesthetic parameters.

Data:
* [Metadata](https://data.qiime2.org/2020.2/tutorials/moving-pictures/sample_metadata.tsv)
* [Unweighted UniFrac PCoA](https://docs.qiime2.org/2020.2/data/tutorials/moving-pictures/core-metrics-results/unweighted_unifrac_pcoa_results.qza)
* [Shannon Diversity](https://docs.qiime2.org/2020.2/data/tutorials/moving-pictures/core-metrics-results/shannon_vector.qza)

```{r}
library(tidyverse)
library(qiime2R)

metadata<-read_q2metadata("sample-metadata.tsv")
uwunifrac<-read_qza("unweighted_unifrac_pcoa_results.qza")
shannon<-read_qza("shannon_vector.qza")$data %>% rownames_to_column("SampleID") 

uwunifrac$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>%
  left_join(shannon) %>%
  ggplot(aes(x=PC1, y=PC2, color=`body-site`, shape=`reported-antibiotic-usage`, size=shannon)) +
  geom_point(alpha=0.5) + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  scale_shape_manual(values=c(16,1), name="Antibiotic Usage") + #see http://www.sthda.com/sthda/RDoc/figure/graphs/r-plot-pch-symbols-points-in-r.png for numeric shape codes
  scale_size_continuous(name="Shannon Diversity") +
  scale_color_discrete(name="Body Site")
  ggsave("PCoA.pdf", height=4, width=5, device="pdf") # save a PDF 3 inches by 4 inches
```
![pcoa](https://github.com/jbisanz/qiime2R/raw/master/images/PCoA.png)

### Plotting a Heatmap

Data:
* [Metadata](https://data.qiime2.org/2020.2/tutorials/moving-pictures/sample_metadata.tsv)
* [Taxonomy](https://docs.qiime2.org/2020.2/data/tutorials/moving-pictures/taxonomy.qza)
* [Feature Table](https://docs.qiime2.org/2020.2/data/tutorials/moving-pictures/table.qza)

Heat maps provide a very good overview of the composition of samples but tend to not perform well for complex samples or visualizing small differences. qiime2R includes a function for making heatmaps called `taxa_heatmap()`.

```
library(tidyverse)
library(qiime2R)

metadata<-read_q2metadata("sample-metadata.tsv")
SVs<-read_qza("table.qza")$data
taxonomy<-read_qza("taxonomy.qza")$data %>% parse_taxonomy()

taxasums<-summarize_taxa(SVs, taxonomy)$Genus

taxa_heatmap(taxasums, metadata, "body-site")

ggsave("heatmap.pdf", height=4, width=8, device="pdf") # save a PDF 4 inches by 8 inches
```
![heatmap](https://github.com/jbisanz/qiime2R/raw/master/images/heatmap.png)

### Making a taxonomic barplot

Data:
* [Metadata](https://data.qiime2.org/2020.2/tutorials/moving-pictures/sample_metadata.tsv)
* [Taxonomy](https://docs.qiime2.org/2020.2/data/tutorials/moving-pictures/taxonomy.qza)
* [Feature Table](https://docs.qiime2.org/2020.2/data/tutorials/moving-pictures/table.qza)

Heat maps provide a very good overview of the composition of samples but tend to not perform well for complex samples or visualizing small differences. qiime2R includes a function for making heatmaps called `taxa_heatmap()`.

```
library(tidyverse)
library(qiime2R)

metadata<-read_q2metadata("sample-metadata.tsv")
SVs<-read_qza("table.qza")$data
taxonomy<-read_qza("taxonomy.qza")$data %>% parse_taxonomy()

taxasums<-summarize_taxa(SVs, taxonomy)$Genus

taxa_barplot(taxasums, metadata, "body-site")

ggsave("barplot.pdf", height=4, width=8, device="pdf") # save a PDF 4 inches by 8 inches
```
![barplot](https://github.com/jbisanz/qiime2R/raw/master/images/barplot.png)


### Differential Abundance Analysis

There are many ways to find features that are deferentially abundant across conditions. In this case we are going to use the output of q2-aldex as described [here](https://library.qiime2.org/plugins/q2-aldex2/24/).

```
library(tidyverse)
library(qiime2R)
library(ggrepel) # for offset labels
library(ggtree) # for visualizing phylogenetic trees
library(ape) # for manipulating phylogenetic trees
metadata<-read_q2metadata("sample-metadata.tsv")
SVs<-read_qza("table.qza")$data
results<-read_qza("differentials.qza")$data
taxonomy<-read_qza("taxonomy.qza")$data
tree<-read_qza("rooted-tree.qza")$data
```

#### Volcano Plot
```
results %>%
  left_join(taxonomy) %>%
  mutate(Significant=if_else(we.eBH<0.1,TRUE, FALSE)) %>%
  mutate(Taxon=as.character(Taxon)) %>%
  mutate(TaxonToPrint=if_else(we.eBH<0.1, Taxon, "")) %>% #only provide a label to signifcant results
  ggplot(aes(x=diff.btw, y=-log10(we.ep), color=Significant, label=TaxonToPrint)) +
  geom_text_repel(size=1, nudge_y=0.05) +
  geom_point(alpha=0.6, shape=16) +
  theme_q2r() +
  xlab("log2(fold change)") +
  ylab("-log10(P-value)") +
  theme(legend.position="none") +
  scale_color_manual(values=c("black","red"))
  ggsave("volcano.pdf", height=3, width=3, device="pdf")
```
![heatmap](https://github.com/jbisanz/qiime2R/raw/master/images/volcano.png)

### Plotting Per-Feature Abundances
If we only wanted to plot one significantly different feature, say 4b5eeb300368260019c1fbc7a3c718fc, we can plot only its abundances. The abundances are not exported by ALDEx2 so instead we will recalculate the CLR-abundances using an approximation of ALDEx; however, in R this could be easily calculated using the package itself.

```
clr<-apply(log2(SVs+0.5), 2, function(x) x-mean(x))
clr %>%
  as.data.frame() %>%
  rownames_to_column("Feature.ID") %>%
  gather(-Feature.ID, key=SampleID, value=CLR) %>%
  filter(Feature.ID=="4b5eeb300368260019c1fbc7a3c718fc") %>%
  left_join(metadata) %>%
  filter(`body-site`=="gut") %>%
  ggplot(aes(x=subject, y=CLR, fill=subject)) +
  stat_summary(geom="bar", color="black") +
  geom_jitter(width=0.2, height=0, shape=21) +
  theme_q2r() +
  theme(legend.position="none")
  ggsave("aldexbar.pdf", height=2, width=1.5, device="pdf") 
```
![aldexbar](https://github.com/jbisanz/qiime2R/raw/master/images/aldexbar.png)

### Plotting a Phylogenetic Tree

Finally, we can visualize the results on a phylogenetic tree with the help of [ggtree](https://guangchuangyu.github.io/software/ggtree/).

```
results<-results %>% mutate(Significant=if_else(we.eBH<0.1,"*", ""))

tree<-drop.tip(tree, tree$tip.label[!tree$tip.label %in% results$Feature.ID]) # remove all the features from the tree we do not have data for
ggtree(tree, layout="circular") %<+% results +
  geom_tippoint(aes(fill=diff.btw), shape=21, color="grey50")  +
  geom_tiplab2(aes(label=Significant), size=10) +
  scale_fill_gradient2(low="darkblue",high="darkred", midpoint = 0, mid="white", name="log2(fold-change") +
  theme(legend.position="right")
ggsave("tree.pdf", height=10, width=10, device="pdf", useDingbats=F)
```
![tree](https://github.com/jbisanz/qiime2R/raw/master/images/tree.png)

***

Note: please report problems or suggestions on the [qiime2 forum](https://forum.qiime2.org/t/tutorial-integrating-qiime2-and-r-for-data-visualization-and-analysis-using-qiime2r/4121) or the [github issue tracker](https://github.com/jbisanz/qiime2R/issues).

