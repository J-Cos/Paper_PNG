#1. packages, path to project directory and source functions ######
    set.seed(0.1)
    
    library(tidyverse)
    library(ggplot2)
    library(gridExtra)
    library(VennDiagram)
    library(phyloseq)
    library(gplots)
    library(vegan)
    library(labdsv)
    library(boral)
    library(fastDummies)
    library(RColorBrewer)

    library(DECIPHER)

    library(rstanarm)
    library(bayestestR)

    library(reshape2)

    #params
    path<-"/home/j/Dropbox/PNG_Paper" #provide the path to the project folder
    merge_pHs<-"No"
    Water<-"No"
    removeAnomalousSamples<-"Yes"

    #Functions
        function_files<-list.files(file.path(path, "Code","Functions"))
        sapply(file.path(path, "BioinformaticPipeline", "Pipeline", "Functions",function_files),source)

        supportfunction_files<-list.files(file.path(path, "..", "BioinformaticPipeline_Env", "BioinformaticPipeline", "SupportFunctions"))
        sapply(file.path(path, "..", "BioinformaticPipeline_Env", "BioinformaticPipeline", "SupportFunctions",supportfunction_files),source)


# 2 Data prep

    # format as phyloseq
        DadaOutput<-readRDS(file=file.path(path, "Data", "PNGFullTest_16s_DadaOutput.RDS"))
        SeqDatTab<-DadaOutput$SeqDataTable
        Seqs<-DNAStringSet( SeqDatTab$Sequence)
        names(Seqs)<-SeqDatTab$ESV


        SeqTab<-SeqDatTab[,c(3:112)]
        rownames(SeqTab)<-SeqDatTab[,1]
        for (name in 1: length(names(SeqTab))) {

            splits<-strsplit(names(SeqTab), "_")

            name_vec<-as.data.frame(splits[[name]][ -( (length(splits[[name]])-1):length(splits[[name]]) ) ]   )
            name_vec<-name_vec[-1,]
            names(SeqTab)[name] <- paste(name_vec, sep="", collapse="_")

        }

        metadata<-sample_data(read.table(file=file.path(path, "Data", "16s_metadata.csv"), sep=",", header=TRUE))
        rownames(metadata)<-metadata$SampleID

        ps<-phyloseq(   otu_table(SeqTab, taxa_are_rows=TRUE),
                        refseq(Seqs),
                        sample_data(metadata)
                )
    
    #remove non-site samples

        ps<-subset_samples(ps, !is.na(Site))

        if (Water=="No"){
        ps<-subset_samples(ps, pH!="D_Water")
        }

        if (merge_pHs=="Yes") {
            sample_data(ps)$pH<-as.character(sample_data(ps)$pH)
            sample_data(ps)$pH[sample_data(ps)$pH!="A_Control"]<-"B_Acidified"
            sample_data(ps)$pH<-as.factor(sample_data(ps)$pH)
        }

    #fix names0
        ps<-prune_samples( sample_sums(ps)>100, ps)
        ps<- prune_taxa( taxa_sums(ps)>0, ps)

        names(sample_data(ps))[2]<-"Sample"
        for (sam in 1:length(sample_names(ps)) ) {
            splits<- str_split(sample_data(ps)$SampleID, "_")
            if (splits[[sam]][2]=="Sessile") {
                sample_data(ps)$Sample[sam]<-"Sessile"
            } else if (splits[[sam]][2]=="100") {
                sample_data(ps)$Sample[sam]<-"100"
            }
        }

        sample_data(ps)$Holobiont[sample_data(ps)$Sample=="100"] <- "Free-living microbes"
        sample_data(ps)$Holobiont[sample_data(ps)$Sample=="Algae" | sample_data(ps)$Sample=="Sessile"] <- "Holobiont microbes"
        sample_data(ps)$Holobiont[sample_data(ps)$Sample=="RVS" | sample_data(ps)$Sample=="CS"] <- "Single sponge species holobiont microbe only"
        sample_data(ps)$Holobiont<-as.factor(sample_data(ps)$Holobiont)

        sample_data(ps)$pH[sample_data(ps)$pH=="Low"]<-"C_Low"
        sample_data(ps)$pH[sample_data(ps)$pH=="Med"]<-"B_Med"
        sample_data(ps)$pH[sample_data(ps)$pH=="Con"]<-"A_Con"

        sample_data(ps)$pH <- as.factor(sample_data(ps)$pH)

        sample_data(ps)$Site <- as.factor(sample_data(ps)$Site)

        sample_data(ps)$Code <- as.factor(sample_data(ps)$Code)

        sample_data(ps)$Sample <- as.factor(sample_data(ps)$Sample)


    #transform sample read counts
        #In order to conduct further analyses all samples are transformed to have an equal number of reads. 
        #This is as an alternative to rarefying in order to avoid data loss.
        psP = transform_sample_counts(ps, function(x) 10000 * x/sum(x))

    

#3/ analysis


    #3.1 alpha diversity (uses untransformed data for reliable alpha diversity estimates)
        #2.1.1 alpha boxplots
        p<-plot_richness(ps, x="pH", color="Sample", measures=c("Observed"))+
                        geom_boxplot() +
                        ggtitle(paste0(" Alpha diversity by pH level, Observed richness")) +
                        facet_grid(~Holobiont)

        p$layers<-p$layers[-1]      
        p         
                          
        p2<-plot_richness(ps, x="pH", color="Sample", shape="Sample", measures=c("Shannon"))+
                        geom_boxplot() +
                        ggtitle(paste0(" Alpha diversity by pH level, Shannon diversity")) +
                        facet_grid(~Holobiont)

        p2$layers<-p2$layers[-1]      
        p2       


        #by holobiont only
        p3<-plot_richness(ps, x="pH", color="Holobiont", measures=c("Observed", "Shannon"))+
                        geom_boxplot() +
                        ggtitle(paste0(" Alpha diversity by pH level, Observed richness"))

        p3$layers<-p3$layers[-1]      
        p3                  

    #3.2 ordination for the holobiont vs non-holobiont
         psH<-prune_samples( sample_data(ps)$Holobiont!="Sponge Holobiont only", ps)
        psH<- prune_taxa( taxa_sums(psH)>0, psH)

          Figs_H<-list()
            methods<-list(c("PCoA", "unifrac"), c("NMDS", NULL))
            for (method in 1:length(methods) ) {
                    psH_distances<-phyloseq::distance(psH, method="bray") 

                    ord_psH <- ordinate(psH, methods[[method]][1], methods[[method]][2], weighted=FALSE, distance=psH_distances)
                    Figs_H[[paste0(methods[[method]][1])]]<-plot_ordination(psH, ord_psH, type="Sample", color="pH", shape="Sample") +
                                        geom_point(size=3)+
                                        ggtitle(paste0("Holobiont vs Free-living",methods[[method]][1])) +
                                        stat_ellipse(aes(group=pH))#+
                                        #facet_wrap(~Sample,2) 
            }

        grid.arrange(grobs = Figs_H, ncol = 2) ## display plot


    # 3.3 sequencce overlaps
        ps_list<-list()
        for ( Holob in levels(sample_data(psP)$Holobiont) ) {
            for (pH in levels(sample_data(psP)$pH)) {

                        ps_temp<-phyloseq::prune_samples( sample_data(psP)$Holobiont==Holob, psP)
                        ps_temp <-prune_taxa(taxa_sums(ps_temp)>0, ps_temp)

                        ps_temp<-prune_samples(sample_data( ps_temp)$pH==pH, ps_temp)

                        ps_list[[ paste0(Holob, pH) ]]<-prune_taxa(taxa_sums(ps_temp)>0, ps_temp)

            }
        }

                OtuOverlapMat<-matrix( 1:(length(ps_list)^2),
                        nrow=length(ps_list),
                        dimnames=list(paste0("Proportion of ", names(ps_list), " OTUs in:"), names(ps_list)))

                for (i in 1:length(ps_list)) {
                    for (j in 1: length(ps_list)){
                        OtuOverlapMat[i,j]<-signif(sum((taxa_names(ps_list[[i]]) %in% taxa_names(ps_list[[j]])))/length(taxa_names(ps_list[[i]])),2)
                    }
                }

        pheatmap::pheatmap(OtuOverlapMat, display_numbers = T,color = colorRampPalette(c('white', 'red'))(100), fontsize_number=15)





OtuOverlapMat2<-matrix( 1:(length(ps_list)^2),
        nrow=length(ps_list),
        dimnames=list(paste0("Proportion of ", names(ps_list), " OTUs in:"), names(ps_list)))

for (i in 1:length(ps_list)) {
    for (j in 1: length(ps_list)){
        OtuOverlapMat2[i,j]<-signif(sum((taxa_names(ps_list[[i]]) %in% taxa_names(ps_list[[j]]))),4)
    }
}

#con
draw.pairwise.venn( area1= 16740,
                    area2= 15300,
                    cross.area= 4979,
                    category=c("Environment", "Holobiont"),
                    fill=c("Blue", "Red")
                    )
grid.newpage()

#med
draw.pairwise.venn( area1= 9294,
                    area2= 12900,
                    cross.area= 3375,
                    category=c("Environment", "Holobiont"),
                    fill=c("Red", "Blue")
                    )
grid.newpage()

#low
draw.pairwise.venn( area1= 13060,
                    area2= 9661,
                    cross.area= 3545,
                    category=c("Environment", "Holobiont"),
                    fill=c("Blue", "Red")
                    )


    #3.4 holobiont non-holobiont ratios
        #by island
            ps_list<-list()
            for ( Holob in levels(sample_data(psP)$Holobiont) ) {
                for (pH in levels(sample_data(psP)$pH)) {
                    for (Site in levels(sample_data(psP)$Site)) {
                            ps_temp<-phyloseq::prune_samples( sample_data(psP)$Holobiont==Holob, psP)
                            ps_temp <-prune_taxa(taxa_sums(ps_temp)>0, ps_temp)

                            ps_temp<-prune_samples(sample_data( ps_temp)$pH==pH, ps_temp)
                            ps_temp <-prune_taxa(taxa_sums(ps_temp)>0, ps_temp)

                            ps_temp<-prune_samples(sample_data( ps_temp)$Site==Site, ps_temp)

                            ps_list[[ paste0(Holob, pH, Site) ]]<-prune_taxa(taxa_sums(ps_temp)>0, ps_temp)
                    }
                }
            }

            names(ps_list)
            ratiosByIsland<-list(
                "Ratio by island of environmental to holobiont microbial richness"=data.frame(
                    "Control"=c(length(taxa_names(ps_list[[1]])) /length(taxa_names(ps_list[[7]])), length(taxa_names(ps_list[[2]])) /length(taxa_names(ps_list[[8]]))),
                    "Medium"=c(length(taxa_names(ps_list[[3]])) /length(taxa_names(ps_list[[9]])), length(taxa_names(ps_list[[4]])) /length(taxa_names(ps_list[[10]]))),
                    "Low"=c(length(taxa_names(ps_list[[5]])) /length(taxa_names(ps_list[[11]])), length(taxa_names(ps_list[[6]])) /length(taxa_names(ps_list[[12]])))
                )
            )
            rownames(ratiosByIsland[[1]])<-c("Dobu", "Illi")

        #by arms unit
            ps_list<-list()
                ps_s100<-prune_samples( (sample_data(psP)$Sample=="Sessile" | sample_data(psP)$Sample=="100") , psP)
                SharedARMS<-levels(sample_data(ps_s100)$Code) [ table(sample_data(ps_s100)$Code)>1 ]
                samps<-levels(sample_data(ps_s100)$Sample)
                for (ARMS in SharedARMS) {
                    for (samp in samps) {
                            ps_temp<-prune_samples(sample_data(ps_s100)$Code==ARMS, ps_s100)
                            ps_temp <-prune_taxa(taxa_sums(ps_temp)>0, ps_temp)

                            ps_temp<-prune_samples(sample_data( ps_temp)$Sample==samp, ps_temp)
                            pH<- sample_data(ps_temp)$pH
                            ps_list[[ paste0(pH, "_", ARMS, "_", samp) ]] <-prune_taxa(taxa_sums(ps_temp)>0, ps_temp)

                    }

                }

            names(ps_list)

            for (i in (1: length(ps_list)) ) {
                print(sample_data(ps_list[[i]])[,4])
            }


            ratiosByARMS<-list(
                "Ratio of environmental to holobiont microbial richness"=data.frame(
                    "A_Control"=c(length(phyloseq::taxa_names(ps_list[[5]])) /length(taxa_names(ps_list[[6]])), 
                                length(taxa_names(ps_list[[7]])) /length(taxa_names(ps_list[[8]])),
                                length(taxa_names(ps_list[[11]])) /length(taxa_names(ps_list[[12]])),
                                length(taxa_names(ps_list[[13]])) /length(taxa_names(ps_list[[14]])),
                                length(taxa_names(ps_list[[21]])) /length(taxa_names(ps_list[[22]])),
                                length(taxa_names(ps_list[[23]])) /length(taxa_names(ps_list[[24]]))
                                ),
                    "B_Medium"=c(length(taxa_names(ps_list[[9]])) /length(taxa_names(ps_list[[10]])), 
                            length(taxa_names(ps_list[[25]])) /length(taxa_names(ps_list[[26]])),
                            length(taxa_names(ps_list[[27]])) /length(taxa_names(ps_list[[28]])),
                            NA,
                            NA,
                            NA
                            ),
                    "C_Low"=c(length(taxa_names(ps_list[[1]])) /length(taxa_names(ps_list[[2]])), 
                            length(taxa_names(ps_list[[3]])) /length(taxa_names(ps_list[[4]])),
                            length(taxa_names(ps_list[[15]])) /length(taxa_names(ps_list[[16]])),
                            length(taxa_names(ps_list[[17]])) /length(taxa_names(ps_list[[18]])),
                            length(taxa_names(ps_list[[19]])) /length(taxa_names(ps_list[[20]])),
                            NA
                            )
                )
            )
            ratiosByARMS[[1]] %>%
                pivot_longer(everything(), names_to="pH") %>%
                ggplot(aes(x=pH, y=value, color=pH)) + 
                    geom_boxplot() +
                    ggtitle("Ratio of environmental to holobiont microbial richness for each ARMS replicate")



    #pivot longer and append island id
            longdf<-pivot_longer(ratiosByARMS[[1]], everything(), names_to="pH")
            longdf$Island<-c(rep("Dobu", 4), "Illi", rep("Dobu", 2), rep("Illi", 3), NA, rep("Illi", 2), NA, rep("Illi", 2), NA, NA)
            longdf$Island<-as.factor(longdf$Island)
            longdf$pH<-as.factor(longdf$pH)

    #bayes glm
            model_general <- stan_glm(data = longdf, formula=longdf[["value"]] ~ pH+Island, family=gaussian())
            Figs$a_CIplot_general<-plot(model_general, plotfun = "areas", prob = 0.9, pars="beta")
            Figs$a_posteriors_general <- describe_posterior(model_general)



























#supplementary figures
    #ordinations
        #full rdination plots
                    Figs_a<-list()
                    methods<-list(c("PCoA", "unifrac"), c("NMDS", NULL))
                    for (method in 1:length(methods) ) {
                            psP_distances<-phyloseq::distance(psP, method="bray") 

                            ord_psP <- ordinate(psP, methods[[method]][1], methods[[method]][2], weighted=FALSE, distance=psP_distances)
                            Figs_a[[paste0(methods[[method]][1])]]<-plot_ordination(psP, ord_psP, type="Sample", color="pH", shape="Sample") +
                                                geom_point(size=3)+
                                                ggtitle(paste0(methods[[method]][1])) +
                                                stat_ellipse(aes(group=pH))#+
                                                #facet_wrap(~Sample,2) 
                    }

                grid.arrange(grobs = Figs_a, ncol = 2) ## display plot



        # every sample seperately ordinations
                    Figs<-list()
                    methods<-list(c("PCoA", "unifrac"), c("NMDS", NULL))
                    for (method in 1:length(methods) ) {
                        for (sam in c("Algae", "RVS", "CS", "Sessile", "100") ) { # remova and add 23s accordingly
                            psP_sam<-subset_samples(psP, Sample==sam)
                            psP_sam<-prune_taxa(taxa_sums(psP_sam)>0,psP_sam)
                            
                            psP_sam_distances<-phyloseq::distance(psP_sam, method="bray") 

                            ord_psP_sam <- ordinate(psP_sam, methods[[method]][1], methods[[method]][2], weighted=FALSE, distance=psP_sam_distances)
                            Figs[[paste0(sam, methods[[method]][1])]]<-plot_ordination(psP_sam, ord_psP_sam, type="Sample", color="pH", shape="Site") +
                                                geom_point(size=3)+
                                                ggtitle(paste0(sam, methods[[method]][1])) +
                                                stat_ellipse(aes(group=pH))+
                                                facet_wrap(~Sample,2) 
                        }
                    }

                grid.arrange(grobs = Figs, ncol = 5) ## display plot
