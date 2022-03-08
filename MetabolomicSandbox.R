library(phyloseq)


btab<-read.table(file.path("../../Data","Metabolomic Data", "METABOLOMICS-SNETS-V2-9965902b-download_cluster_buckettable-main.tsv"),sep="\t", comment.char = "", check.names=FALSE, header=TRUE, row.names=1)

#prepare temporary sample metadata (until ER provides)
tempsamplemetadata<-strsplit( names(btab), "_")
for (i in 1:length(tempsamplemetadata)) {
    tempsamplemetadata[[i]][3]<-substr(tempsamplemetadata[[i]][3], start=1, stop=2)
}
tempsamplemetadata[ sapply(tempsamplemetadata, length)<5][[1]]<- c("Blank", "Blank", "RA", "01", "2404")
tempsamplemetadata<-t(as.data.frame(sapply(tempsamplemetadata, unlist)))
colnames(tempsamplemetadata)<-c("PI", "SampleNumber", "SampleCode", "01", "MiscNumber")
tempsamplemetadata_df<-sample_data(data.frame(tempsamplemetadata))

# sort out names of btab
names(btab)<-rownames(tempsamplemetadata_df)

#names(btab)<-sapply( strsplit( names(btab), "_"),  '[', 3)


alternative_names<-sapply( # names of form ###_XX#
                    strsplit(
                            sapply(
                                    strsplit( names(btab), "_01"), 
                                    '[', 1),
                            "ER_"), 
                    '[', 2)

btab<-otu_table(btab, taxa_are_rows=TRUE)
btab<-phyloseq(btab, tempsamplemetadata_df)


#In order to conduct further analyses all samples are transformed to have an equal number of reads. 
#This is as an alternative to rarefying in order to avoid data loss.
btab = transform_sample_counts(btab, function(x) 10000 * x/sum(x))


btab_ord <- ordinate(btab, "NMDS", "bray")
plot_ordination(btab, btab_ord, type="sample", color="SampleNumber")

