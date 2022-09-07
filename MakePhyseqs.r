# 1. Packages and path to project directory ######
    set.seed(0.1)
    
    library(tidyverse)
    library(phyloseq)
    library(ggplot2)
    library(treemapify)
    library(viridis)

    library(breakaway)
    
    #params
        path<-"/home/j/Dropbox/PNG_Paper" #provide the path to the project folder
        
    #1.1 load functions
    #funcs
    function_files<-list.files(file.path(path, "Functions"))
    sapply(file.path(path, "Functions",function_files),source)
    
    function_files<-list.files(file.path("/home/j/Dropbox/BioinformaticPipeline_Env/BioinformaticPipeline/SupportFunctions"))
    sapply(file.path("/home/j/Dropbox/BioinformaticPipeline_Env/BioinformaticPipeline/SupportFunctions",function_files),source)

# 2. Data prep
    # 2a. Sequencing Data
        # format as phyloseq
            #16s
                metadata_16s<-read.table(file=file.path(path, "Data", "16s_metadata.csv"), sep=",", header=TRUE, row.names=1) %>% 
                                as_tibble %>%
                                mutate_if(is.character, as.factor) %>%
                                mutate  (pH=factor(pH, order = TRUE,  levels = c("Control pH", "Medium pH", "Low pH")) )  %>%
                                sample_data
                sample_names(metadata_16s)<-sample_names(sample_data(read.table(file=file.path(path, "Data", "16s_metadata.csv"), sep=",", header=TRUE, row.names=1)))
                ps16<-SeqDataTable2Phyloseq( SeqDataTablePath=file.path(path, "Data", "PNGFullTest_16s_SeqDataTable.RDS"),
                                            clustering="ESV",
                                            assignment="Idtaxa",
                                            Metadata=metadata_16s)
                #remove water and uncertains
                    ps16<-subset_samples(ps16, !is.na(Site))
                #ensure all samples have x seqs and all taxa have >0 (all 23s already do)
                    ps16<-prune_samples( sample_sums(ps16)>50000, ps16)
                    ps16<-prune_samples( !sample_names(ps16)%in% c("Algae_Con_Dobu_103", "Algae_Low_Dobu_97", "Algae_Med_Dobu_91", "Algae_Con_Illi_219"), ps16)  #removing lower read count duplicate sequencing of same arms plate            
                    ps16<- prune_taxa( taxa_sums(ps16)>0, ps16)

            #23s
                metadata_23s<-read.table(file=file.path(path, "Data", "23s_metadata.csv"), sep=",", header=TRUE, row.names=1) %>% 
                                as_tibble %>%
                                mutate_if(is.character, as.factor) %>%
                                mutate  (pH=factor(pH, order = TRUE,  levels = c("Control pH", "Medium pH", "Low pH")) )  %>%
                                sample_data
                sample_names(metadata_23s)<- metadata_23s$Code



                ps23<-SeqDataTable2Phyloseq( SeqDataTablePath=file.path(path, "Data", "PNGFullTest_23s_SeqDataTable.RDS"),
                                            clustering="ESV",
                                            assignment="Idtaxa",
                                            Metadata=metadata_23s)
                sample_names(ps23)<-paste0("Sample_", sample_names(ps23))
                #ensure all samples have x seqs and all taxa have >0 (all 23s already do)
                    ps23<-prune_samples( sample_sums(ps23)>50000, ps23)
                    ps23<-prune_samples( !sample_names(ps23)%in%c("Sample_219", "Sample_104", "Sample_92", "Sample_97"), ps23)  #removing lower read count duplicate sequencing of same arms plate    
                    ps23<- prune_taxa( taxa_sums(ps23)>0, ps23)

    # 2b. Metabolomic data
        #btab<-read.table(file.path(path, "Data", "Metabolomic Data", "METABOLOMICS-SNETS-V2-9965902b-download_cluster_buckettable-main.tsv"),sep="\t", comment.char = "", check.names=FALSE, header=TRUE, row.names=1)
        btab<-read.table(file.path(path, "Data", "Metabolomic Data", "quantification_table-00000.csv"),sep=",", comment.char = "", check.names=FALSE, header=TRUE, row.names=1)
        x<-strsplit( names(btab), '_')
        names(btab)<-paste0("Vial_", sapply(x, '[', 2))

        metadata<-read.csv(file.path(path, "Data", "Metabolomic Data", "MetabolomeMetadata.csv"), row.names=1)

        annotations<-as.data.frame(data.table::fread(file.path(path, "Data", "Metabolomic Data",  "FEATURE-BASED-MOLECULAR-NETWORKING-f5c65917-view_all_annotations_DB-main.tsv"), sep='\t')      )
        basefordf<-as.data.frame(as.integer(rownames(btab)))
        anot_tab<-left_join(basefordf, annotations, by=c("as.integer(rownames(btab))"='#Scan#'))
        rownames(anot_tab)<-paste0("Compound_", rownames(anot_tab))
        anot_tab$SpectrumID[is.na(anot_tab$SpectrumID)]<-"Unclassified"
        rownames(btab)<-paste0("Compound_", rownames(btab))
        btab<-phyloseq(  otu_table(btab, taxa_are_rows=TRUE),
                             sample_data(metadata), 
                             tax_table(as.matrix(anot_tab)))

        #removes other solvents which only exist for sessile and therefore can't be used for comparisons
        btab<-prune_samples(sample_data(btab)$Solvent=="MeOH", btab)
        btab<- prune_taxa( taxa_sums(btab)>0, btab)
        sort(sample_sums(btab))

saveRDS(ps16, file=file.path(path, "Outputs", "ps16.RDS"))
saveRDS(ps23, file=file.path(path, "Outputs", "ps23.RDS"))
saveRDS(btab, file=file.path(path, "Outputs", "btab.RDS"))


