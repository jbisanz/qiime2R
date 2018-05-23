## --------------------------------------------------------------------------
library(qiime2R)
library(phyloseq)
library(tidyverse)

## --------------------------------------------------------------------------
download.file("https://docs.qiime2.org/2018.4/data/tutorials/moving-pictures/table.qza", "table.qza")
SequenceVariants<-read_qza("table.qza")

## --------------------------------------------------------------------------
names(SequenceVariants)

## --------------------------------------------------------------------------
print_provenance(SequenceVariants)

## --------------------------------------------------------------------------
download.file("https://data.qiime2.org/2018.4/tutorials/moving-pictures/sample_metadata.tsv", "sample_metadata.tsv")
download.file("https://docs.qiime2.org/2018.4/data/tutorials/moving-pictures/taxonomy.qza", "taxonomy.qza")
download.file("https://docs.qiime2.org/2018.4/data/tutorials/moving-pictures/rooted-tree.qza", "rooted-tree.qza")

phyobj<-qza_to_phyloseq(features="table.qza", taxonomy = "taxonomy.qza", tree = "rooted-tree.qza", metadata="sample_metadata.tsv")

phyobj

## --------------------------------------------------------------------------
download.file("https://docs.qiime2.org/2018.4/data/tutorials/moving-pictures/core-metrics-results/shannon_vector.qza", "shannon_vector.qza")

metadata<-
  read_tsv("sample_metadata.tsv", comment="#q2:types") #to exclude the column denoting the variable class

read_qza("shannon_vector.qza")$data %>%
  as.data.frame() %>%
  rownames_to_column("#SampleID") %>% # to allow a smooth joining with the metadata
  left_join(metadata) %>%
  ggplot(
    aes(
      x=DaysSinceExperimentStart, 
      y=shannon, 
      group=BodySite, 
      color=BodySite, 
      shape=ReportedAntibioticUsage)
        ) +
  geom_line() +
  geom_point() +
  facet_wrap(~Subject) + #make a separate plot for each subject
  theme_bw()

## --------------------------------------------------------------------------
download.file("https://docs.qiime2.org/2018.4/data/tutorials/moving-pictures/core-metrics-results/unweighted_unifrac_pcoa_results.qza", "unweighted_unifrac_pcoa_results.qza")

pcoa<-read_qza("unweighted_unifrac_pcoa_results.qza")$data
pcoa$Vectors %>%
  rename("#SampleID"=SampleID) %>% # to match the metadata
  left_join(metadata) %>%
  ggplot(
    aes(
      x=PC1, 
      y=PC2, 
      color=BodySite,
      shape=Subject
        )
    ) +
  geom_point() +
  theme_bw() +
  xlab(paste0("PC1: ", pcoa$ProportionExplained["PC1"])) + #add variance explained to axis
  ylab(paste0("PC2: ", pcoa$ProportionExplained["PC2"]))


