# 1. Packages and path to project directory ######
    set.seed(0.1)
    
    library(tidyverse)
    library(ggplot2)
    library(cowplot)
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
    library(ggsci)
    library(indicspecies)
    library(ggsignif)
    library(lme4)
    library(MuMIn)
    library(MsBackendMgf)

    #for mapping
    library(sf)
    library(rnaturalearth)
    library(rnaturalearthdata)
    library(rworldxtra)
    library(ggmap)
    library(ggsn)
    
    #params
        path<-"/home/j/Dropbox/PNG_Paper" #provide the path to the project folder
        theme_set(theme_linedraw())
        Figs<-list()
        Stats<-list()



    #1.1 load functions
    #funcs
    function_files<-list.files(file.path(path, "Functions"))
    sapply(file.path(path, "Functions",function_files),source)
    
    function_files<-list.files(file.path("/home/j/Dropbox/BioinformaticPipeline_Env/BioinformaticPipeline/SupportFunctions"))
    sapply(file.path("/home/j/Dropbox/BioinformaticPipeline_Env/BioinformaticPipeline/SupportFunctions",function_files),source)


# 2. Data prep
    # get data
        ps16<-readRDS(file=file.path(path, "Outputs", "ps16.RDS"))
        ps23<-readRDS(file=file.path(path, "Outputs", "ps23.RDS"))
        btab<-readRDS(file=file.path(path, "Outputs", "btab.RDS"))

    # 2a. Sequencing Data
        #remove additional 16s samples we don't want to use
            #remove water and uncertains
                ps16<-subset_samples(ps16, !is.na(Site))

            #ensure all samples have x seqs and all taxa have >0 (all 23s already do)
                ps16<-prune_samples( sample_sums(ps16)>100, ps16)
                ps16<- prune_taxa( taxa_sums(ps16)>0, ps16)

        
        #wrangle sample data
            #add holobiont column
                sample_data(ps16)$Holobiont[sample_data(ps16)$Sample=="100"] <- "Free-living microbes"
                sample_data(ps16)$Holobiont[sample_data(ps16)$Sample=="Algae" | sample_data(ps16)$Sample=="Sessile"] <- "Holobiont microbes"
                sample_data(ps16)$Holobiont[sample_data(ps16)$Sample=="RVS" | sample_data(ps16)$Sample=="CS"] <- "Single sponge species holobiont microbe only"
                
                sample_data(ps23)$Holobiont[sample_data(ps23)$Sample=="Algae_23s"] <- "Holobiont Macrobes"

            #rename pHs so ordering is automatic
                sample_data(ps16)$pH[sample_data(ps16)$pH=="Low"]<-"C_Low"
                sample_data(ps16)$pH[sample_data(ps16)$pH=="Med"]<-"B_Med"
                sample_data(ps16)$pH[sample_data(ps16)$pH=="Con"]<-"A_Con"

                sample_data(ps23)$pH[sample_data(ps23)$pH=="Low"]<-"C_Low"
                sample_data(ps23)$pH[sample_data(ps23)$pH=="Med"]<-"B_Med"
                sample_data(ps23)$pH[sample_data(ps23)$pH=="Con"]<-"A_Con"

            #remove numeric sample_names on 23s
                sample_names(ps23)<-paste0("Sample_", sample_names(ps23))
            
            #convert all to factors
                cols<-c("Sample", "pH", "Site", "Code", "ARMS", "TechnicalReplicate", "Holobiont")
                sample_data(ps16)[,cols] <- lapply(sample_data(ps16)[,cols], as.factor)
                sample_data(ps23)[,cols] <- lapply(sample_data(ps23)[,cols], as.factor)
                sample_data(ps16)$pH<-relevel(sample_data(ps16)$pH, ref="A_Con")
                sample_data(ps23)$pH<-relevel(sample_data(ps23)$pH, ref="A_Con")

        #create chloroplast only ps
            pschloro<- prune_taxa( as.data.frame(tax_table(ps16) )$Rank_5=="Chloroplast", ps16)
            pschloro_a<-prune_samples(as.data.frame(sample_data(pschloro))$Sample=="Algae", pschloro)
            sample_data(pschloro_a)$Sample<-"Plastids"
            sample_names(pschloro_a)<-paste0(sample_names(pschloro_a), "_Plastids")

        #remove chloroplasts from main ps16
            ps16<- prune_taxa( as.data.frame(tax_table(ps16) )$Rank_5!="Chloroplast", ps16)
            #relabel cyanobacteriia as cyanobacteria
                taxtab<-as.data.frame(tax_table(ps16))
                taxtab$Rank_4[taxtab$Rank_4=="Cyanobacteriia"]<-"Cyanobacteria"
                tax_table(ps16)<-as.matrix(taxtab)

        #create psALL
            ps16ALL<-ps16
            psALL<-merge_phyloseq( merge_phyloseq(ps16ALL, pschloro_a) , ps23)

        #seperate out single holobiont samples
            ps16<-prune_samples( !(sample_data(ps16)$Sample=="RVS" | sample_data(ps16)$Sample=="CS") , ps16)
            ps16<- prune_taxa( taxa_sums(ps16)>0, ps16)
        
        #create algae holobiome ps
            ps16_a<-prune_samples( sample_data(ps16)$Sample=="Algae", ps16)
            ps16_a<- prune_taxa( taxa_sums(ps16_a)>0, ps16_a)
            ps_AlgaeHolobiome<-merge_phyloseq( merge_phyloseq(ps16_a, pschloro_a) , ps23)





        #transform sample read counts
            #In order to conduct further analyses all samples are transformed to have an equal number of reads. 
            #This is as an alternative to rarefying in order to avoid data loss.
            psALLp = transform_sample_counts(psALL, function(x) 10000 * x/sum(x))
            ps16ALLp = transform_sample_counts(ps16ALL, function(x) 10000 * x/sum(x))
            ps16p = transform_sample_counts(ps16, function(x) 10000 * x/sum(x))
            ps23p = transform_sample_counts(ps23, function(x) 10000 * x/sum(x))
            ps_AlgaeHolobiomep = transform_sample_counts(ps_AlgaeHolobiome, function(x) 10000 * x/sum(x))

    # 2b. Metabolomic data
        #btab<-read.table(file.path(path, "Data", "Metabolomic Data", "METABOLOMICS-SNETS-V2-9965902b-download_cluster_buckettable-main.tsv"),sep="\t", comment.char = "", check.names=FALSE, header=TRUE, row.names=1)
        btab<-read.table(file.path(path, "Data", "Metabolomic Data", "quantification_table-00000.csv"),sep=",", comment.char = "", check.names=FALSE, header=TRUE, row.names=1)
        x<-strsplit( names(btab), '_')
        names(btab)<-paste0("Vial_", sapply(x, '[', 2))

        metadata<-read.csv(file.path(path, "Data", "Metabolomic Data", "MetabolomeMetadata.csv"), row.names=1)

        test<-MsBackendMgf::readMgf(f = file.path(path, "Data", "Metabolomic Data", "specs_ms.mgf"))

        annotations<-as.data.frame(data.table::fread(file.path(path, "Data", "Metabolomic Data",  "FEATURE-BASED-MOLECULAR-NETWORKING-f5c65917-view_all_annotations_DB-main.tsv"), sep='\t')      )
        cbind(rownames(btab), )
        annotations[,21]
        basefordf<-as.data.frame(as.integer(rownames(btab)))
        anot_tab<-left_join(basefordf, annotations, by=c("as.integer(rownames(btab))"='#Scan#'))

        rownames(anot_tab)<-paste0("Compound_", rownames(anot_tab))
        rownames(btab)<-paste0("Compound_", rownames(btab))

        btab<-phyloseq(  otu_table(btab, taxa_are_rows=TRUE),
                             sample_data(metadata), 
                             tax_table(as.matrix(anot_tab)))

        sample_data(btab)$Holobiont[sample_data(btab)$Sample=="100-500UM"] <- "Free-living Microbe Associated"
        sample_data(btab)$Holobiont[sample_data(btab)$Sample=="SESSILE"| sample_data(btab)$Sample=="Algae"] <- "Holobiont Associated"
        sample_data(btab)$Holobiont[sample_data(btab)$Sample=="RVS" | sample_data(btab)$Sample=="CS" ] <- "Single sponge species holobiont microbe only"

        sample_data(btab)$pH[sample_data(btab)$pH=="Low"]<-"C_Low"
        sample_data(btab)$pH[sample_data(btab)$pH=="Medium"]<-"B_Med"
        sample_data(btab)$pH[sample_data(btab)$pH=="Control"]<-"A_Con"

        #convert all to factors
            cols<-c("Sample", "pH", "Site", "ARMS", "TechnicalReplicate", "Holobiont")
            sample_data(btab)[,cols] <- lapply(sample_data(btab)[,cols], as.factor)
            sample_data(btab)$pH<-relevel(sample_data(btab)$pH, ref="A_Con")


        #removes other solvents which only exist for sessile and therefore can't be used for comparisons
        btab<-prune_samples(sample_data(btab)$Solvent=="MeOH", btab)
        btab<- prune_taxa( taxa_sums(btab)>0, btab)


        #In order to conduct further analyses all samples are transformed to have an equal number of reads. 
        #This is as an alternative to rarefying in order to avoid data loss.
        btab_p = transform_sample_counts(btab, function(x) 10000 * x/sum(x))

        btab_count = transform_sample_counts(btab, function(x) 100 * x)


# X. Create a custom color scale
    myColors<- c("darkred", "darkturquoise", "chartreuse", "deeppink", "darkgoldenrod1")
    names(myColors) <- c("Holobiont Community Microbiome", "Environmental Microbiome", "Photosynthetic Community Microbiome", "Halisarca Sponge", "Tethya Sponge")
    colScale <- scale_colour_manual(name = "Fraction",values = myColors)
    fillScale <- scale_fill_manual(name = "Fraction",values = myColors)

    myColors_individOnly<-myColors[4:5]
    colScale_individOnly <- scale_colour_manual(name = "Fraction",values = myColors_individOnly)
    fillScale_individOnly <- scale_fill_manual(name = "Fraction",values = myColors_individOnly)


    myColors<- c("darkred", "darkturquoise", "chartreuse", "deeppink", "darkgoldenrod1")
    names(myColors) <- c("Holobiont Community Metabolome", "Environmental Metabolome", "Photosynthetic Community Metabolome", "Halisarca Sponge", "Tethya Sponge")
    ChemcolScale <- scale_colour_manual(name = "Fraction",values = myColors)
    ChemfillScale <- scale_fill_manual(name = "Fraction",values = myColors)

    myColors_individOnly<-myColors[4:5]
    ChemcolScale_individOnly <- scale_colour_manual(name = "Fraction",values = myColors_individOnly)
    ChemfillScale_individOnly <- scale_fill_manual(name = "Fraction",values = myColors_individOnly)

    pHColors<- c("green", "yellow", "orange")
    names(pHColors) <- c("Control pH", "Medium pH", "Low pH")
    pHcolScale <- scale_colour_manual(name = "pH Regime",values = pHColors)
    pHfillScale <- scale_fill_manual(name = "pH Regime",values = pHColors)

# 3. Figure 3 - richness
        #anotations
        anno <- data.frame(x1 = c(1,1), x2 = c(1,3), 
                   y1 = c(0,6000), y2 = c(0,6000), 
                   xstar = c(1,2), ystar = c(0,6500),
                   lab = c("","***"),
                   Sample = c("Environmental", "Algae"))

        # seq violin richness Figure1a
            ps_s100<-prune_samples(sample_data(ps16)$Sample=="Sessile" | sample_data(ps16)$Sample=="100" | sample_data(ps16)$Sample=="Algae", ps16)
            ps_s100<-prune_taxa(taxa_sums(ps_s100)>0, ps_s100)
            levels(sample_data(ps_s100)$Sample)<-c("Environmental", "Algae", "Holobiont")
            Figs[["Figure1a"]]<-plot_richness( ps_s100, x="pH", measures=c("Observed"))+
                    geom_violin(aes(fill=Sample, color=Sample)) +
                    ggtitle(paste0("Figure 3: Multiomic Richness Under Ocean Acidification"), ) +
                    ylab("Sequence Richness") +
                    guides(color="none")+
                    guides(fill="none") +
                    expand_limits(y = c(0, 15000)) +
                    scale_x_discrete(name="", labels=c("", "" ,""))+
                    scale_y_continuous( breaks = c(5000, 10000), expand=c(0,0)) +
                    theme(  strip.background = element_rect(fill = "white"),
                            #panel.border = element_rect(colour = "black", fill = "white"),
                            strip.text = element_text(size=14, color="black", face="bold"),
                            axis.text.x = element_text(face="bold", size=14, angle=0, hjust=0.5),
                            axis.text.y = element_text(face="bold", size=14, angle=0, hjust=0.5),
                            axis.title.y = element_text(face="bold", size=14),
                            plot.title = element_text(size = 20, face = "bold", hjust=0.5),
                            panel.grid.major = element_blank(), 
                            panel.grid.minor = element_blank(),
                            panel.background = element_blank(), 
                            axis.line = element_line(colour = "black")) +
                    geom_text(data = anno[2,], aes(x = xstar,  y = ystar, label = lab), size=7) +
                    geom_segment(data = anno[2,], aes(x = x1, xend = x2, y = y1, yend = y2), size=1,   colour = "black") +
                    facet_wrap(~Sample, ncol=3, labeller = labeller(Sample = 
                                                    c("Environmental" = "Environmental Microbiome",
                                                    "Algae" = "Algae Microbiome",
                                                    "Holobiont" = "Holobiont Microbiome")))+
                    colScale+
                    fillScale

                    
            Figs[["Figure1a"]]$layers<-Figs[["Figure1a"]]$layers[-1]  
            Figs[["Figure1a"]]
        # chemo violin richness Figure1b
            btab_nosponge<-prune_samples( !(sample_data(btab_count)$Sample=="RVS" | sample_data(btab_count)$Sample=="CS" | sample_data(btab_count)$Sample=="TS") , btab_count)
            levels(sample_data(btab_nosponge)$Sample)<-c("Environmental", "Algae", "Holobiont")

            Figs[["Figure1b"]]<-plot_richness(btab_nosponge, x="pH", measures=c("Observed"))+
                            geom_violin(aes(fill=Sample, color=Sample)) +
                            #ggtitle(paste0(" Chemical Richness by pH regime")) +
                            scale_color_npg() + 
                            scale_fill_npg() +
                            ylab("Metabolite Richness") +
                            guides(color="none")+
                            guides(fill="none") +
                            expand_limits(y = c(0, 300)) +
                            theme(  strip.background = element_blank(),
                                    #panel.border = element_blank(),
                                    strip.text = element_blank(),
                                    axis.text.x = element_text(face="bold", size=14, angle=0, hjust=0.5),
                                    axis.text.y = element_text(face="bold", size=14, angle=0, hjust=0.5),
                                    axis.title.y = element_text(face="bold", size=14),
                                    #plot.title = element_text(size = 40, face = "bold", hjust=0.5),
                                    panel.grid.major = element_blank(), 
                                    panel.grid.minor = element_blank(),
                                    panel.background = element_blank(), 
                                    axis.line = element_line(colour = "black")) +
                            scale_x_discrete(name="", labels=c("Control pH", "Medium pH", "Low pH"))+
                            scale_y_continuous(breaks = c(100, 200), expand=c(0,0))+
                            geom_text(data = anno[1,], aes(x = xstar,  y = ystar, label = lab), size=7) +
                            geom_segment(data = anno[1,], aes(x = x1, xend = x2, y = y1, yend = y2), size=1,   colour = "black") +
                            facet_wrap(~Sample, ncol=3, labeller = labeller(Sample = 
                                                            c("Environmental" = "Environmental Microbiome",
                                                            "Algae" = "Algae Microbiome",
                                                            "Holobiont" = "Holobiont Microbiome")))+
                            colScale+
                            fillScale

            Figs[["Figure1b"]]$layers<-Figs[["Figure1b"]]$layers[-1]      
            Figs[["Figure1b"]]

# 4. Figure 4 - distinctness
    # make initial lists
        # seqs 
                psP<-ps16p

            #by arms unit
                ps_list<-list()
                    ps_s100<-prune_samples( (sample_data(psP)$Sample=="Holobiont Community Microbiome" | sample_data(psP)$Sample=="Environmental Microbiome") , psP)
                    SharedARMS<-sort(unique(sample_data(ps_s100)$Code)) [ table(sample_data(ps_s100)$Code)>1 ]
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
            #helper code for correctly extracting desired ps from ps_list
                #        names(ps_list)
                #
                #        for (i in (1: length(ps_list)) ) {
                #            print(sample_data(ps_list[[i]])[,4])
                #        }
        # chems
            #by arms unit
            psb_list<-list()
                psb_s100<-prune_samples( (sample_data(btab_p)$Sample=="Holobiont Community Metabolome" | sample_data(btab_p)$Sample=="Environmental Metabolome") & sample_data(btab_p)$Solvent=="MeOH" , btab_p)
                #sample_data(psb_s100)[sample_data(psb_s100)$ARMS=="29"]
                SharedARMS<-sort(unique(sample_data(psb_s100)$ARMS)) [ table(sample_data(psb_s100)$ARMS)>1 ]
                samps<-unique(sample_data(psb_s100)$Sample)
                for (ARMS in SharedARMS) {
                    for (samp in samps) {
                            psb_temp<-prune_samples(sample_data(psb_s100)$ARMS==ARMS, psb_s100)
                            psb_temp <-prune_taxa(taxa_sums(psb_temp)>0, psb_temp)

                            psb_temp<-prune_samples(sample_data( psb_temp)$Sample==samp, psb_temp)
                            pH<- sample_data(psb_temp)$pH
                            psb_list[[ paste0(pH, "_", ARMS, "_", samp) ]] <-prune_taxa(taxa_sums(psb_temp)>0, psb_temp)

                    }

                }


    # Holobiome/Environmental distinct proportions sequence
        #Sequences distinct
                OverlapByArms<-list(
                    "Proportion of environmental microbial taxa in holobionts"=data.frame(
                        "A_Control"=c( sum((taxa_names(ps_list[[5]]) %in% taxa_names(ps_list[[6]])))/length(taxa_names(ps_list[[5]])),
                                    sum((taxa_names(ps_list[[7]]) %in% taxa_names(ps_list[[8]])))/length(taxa_names(ps_list[[7]])),
                                    sum((taxa_names(ps_list[[11]]) %in% taxa_names(ps_list[[12]])))/length(taxa_names(ps_list[[11]])),
                                    sum((taxa_names(ps_list[[13]]) %in% taxa_names(ps_list[[14]])))/length(taxa_names(ps_list[[13]])),
                                    sum((taxa_names(ps_list[[21]]) %in% taxa_names(ps_list[[22]])))/length(taxa_names(ps_list[[21]])),
                                    sum((taxa_names(ps_list[[23]]) %in% taxa_names(ps_list[[24]])))/length(taxa_names(ps_list[[23]]))
                                    ),
                        "B_Medium"=c(sum((taxa_names(ps_list[[9]]) %in% taxa_names(ps_list[[10]])))/length(taxa_names(ps_list[[9]])), 
                                sum((taxa_names(ps_list[[25]]) %in% taxa_names(ps_list[[26]])))/length(taxa_names(ps_list[[25]])),
                                sum((taxa_names(ps_list[[27]]) %in% taxa_names(ps_list[[28]])))/length(taxa_names(ps_list[[27]])),
                                NA,
                                NA,
                                NA
                                ),
                        "C_Low"=c(sum((taxa_names(ps_list[[1]]) %in% taxa_names(ps_list[[2]])))/length(taxa_names(ps_list[[1]])), 
                                sum((taxa_names(ps_list[[3]]) %in% taxa_names(ps_list[[4]])))/length(taxa_names(ps_list[[3]])),
                                sum((taxa_names(ps_list[[15]]) %in% taxa_names(ps_list[[16]])))/length(taxa_names(ps_list[[15]])),
                                sum((taxa_names(ps_list[[17]]) %in% taxa_names(ps_list[[18]])))/length(taxa_names(ps_list[[17]])),
                                sum((taxa_names(ps_list[[19]]) %in% taxa_names(ps_list[[20]])))/length(taxa_names(ps_list[[19]])),
                                NA
                                )
                    ),
                    "Proportion of holobiont microbial taxa in environment"=data.frame(
                        "A_Control"=c( sum((taxa_names(ps_list[[6]]) %in% taxa_names(ps_list[[5]])))/length(taxa_names(ps_list[[6]])),
                                    sum((taxa_names(ps_list[[8]]) %in% taxa_names(ps_list[[7]])))/length(taxa_names(ps_list[[8]])),
                                    sum((taxa_names(ps_list[[12]]) %in% taxa_names(ps_list[[11]])))/length(taxa_names(ps_list[[12]])),
                                    sum((taxa_names(ps_list[[14]]) %in% taxa_names(ps_list[[13]])))/length(taxa_names(ps_list[[14]])),
                                    sum((taxa_names(ps_list[[22]]) %in% taxa_names(ps_list[[21]])))/length(taxa_names(ps_list[[22]])),
                                    sum((taxa_names(ps_list[[24]]) %in% taxa_names(ps_list[[23]])))/length(taxa_names(ps_list[[24]]))
                                    ),
                        "B_Medium"=c(sum((taxa_names(ps_list[[10]]) %in% taxa_names(ps_list[[9]])))/length(taxa_names(ps_list[[10]])), 
                                sum((taxa_names(ps_list[[26]]) %in% taxa_names(ps_list[[25]])))/length(taxa_names(ps_list[[26]])),
                                sum((taxa_names(ps_list[[28]]) %in% taxa_names(ps_list[[27]])))/length(taxa_names(ps_list[[28]])),
                                NA,
                                NA,
                                NA
                                ),
                        "C_Low"=c(sum((taxa_names(ps_list[[2]]) %in% taxa_names(ps_list[[1]])))/length(taxa_names(ps_list[[2]])), 
                                sum((taxa_names(ps_list[[4]]) %in% taxa_names(ps_list[[3]])))/length(taxa_names(ps_list[[4]])),
                                sum((taxa_names(ps_list[[16]]) %in% taxa_names(ps_list[[15]])))/length(taxa_names(ps_list[[16]])),
                                sum((taxa_names(ps_list[[18]]) %in% taxa_names(ps_list[[17]])))/length(taxa_names(ps_list[[18]])),
                                sum((taxa_names(ps_list[[20]]) %in% taxa_names(ps_list[[19]])))/length(taxa_names(ps_list[[20]])),
                                NA
                                )
                    )
                )


            #pivot longer and append island id
                longdf<-cbind(pivot_longer(OverlapByArms[[1]], everything(), names_to="pH"), pivot_longer(OverlapByArms[[2]], everything(), names_to="pH"))
                longdf<-longdf[c(1,2,4)]
                names(longdf)<- c("pH" ,"Proportion of ESVs from the environment also found in Holobiont",  "Proportion of ESVs from Holobionts also found in the environment")
                longdf$Island<-c(rep("Dobu", 4), "Illi", rep("Dobu", 2), rep("Illi", 3), NA, rep("Illi", 2), NA, rep("Illi", 2), NA, NA)
                longdf$Island<-as.factor(longdf$Island)
                longdf$pH<-as.factor(longdf$pH)


                longerdf<-pivot_longer(longdf, 2:3, names_to="overlapType")
                longerdf$overlapType<-as.factor(longerdf$overlapType)
                levels(longerdf$overlapType)<-c("Holobiont Richness", "Environmental Richness")
                
                longerdf$value<-longerdf$value*100 #make percentage
                longerdf$value<-100-longerdf$value #convert percent shared to percent not shared


                Figs[["ESVOverlap_poster"]] <- ggplot()+
                                            geom_hline(yintercept=50, linetype='dotted')+
                                            geom_violin(data=longerdf[longerdf$overlapType=="Holobiont Richness",], aes(x=pH, y=value), fill="Black") +
                                            #geom_violin(data=longerdf[longerdf$overlapType=="Environmental Richness",], aes(x=pH, y=value, fill=overlapType)) +
                                            geom_boxplot(data=longerdf[longerdf$overlapType=="Holobiont Richness",],aes(x=pH, y=value), width=0.1, color="grey", fill="grey", alpha=0.5) +
                                            #geom_boxplot(data=longerdf[longerdf$overlapType=="Environmental Richness",],aes(x=pH, y=value), width=0.1, color="grey", fill="grey", alpha=0.5) +
                                            labs(   #title="Figure 4: Multiomic Distinctness",
                                                    #subtitle="",
                                                    y="Percentage distinct to holobiont community", 
                                                    #x="pH Regime (control n=6, medium n=3, low n=5)"
                                                    ) +
                                            scale_color_npg() +
                                            scale_fill_lancet() +
                                            annotate("segment", x = 1, xend = 2, y = 30, yend = 30, colour = "black", size=1, alpha=1) +
                                            annotate("text", x = 1.5, y = 28, size=7, label = "**") +
                                            annotate("segment", x = 1, xend = 3, y = 25, yend = 25, colour = "black", size=1, alpha=1) +
                                            annotate("text", x = 2, y = 23, size=7, label = "**") +
                                            expand_limits(y = c(0, 100)) +
                                            theme(  strip.background = element_rect(fill = "white"),
                                                    #panel.border = element_blank(),
                                                    strip.text = element_text(size=20, color="black", face="bold"),
                                                    axis.text.x = element_text(face="bold", size=20, angle=0, hjust=0.5),
                                                    axis.text.y = element_text(face="bold", size=20, angle=0, hjust=0.5),
                                                    axis.title.y = element_text(face="bold", size=20),
                                                    #plot.title = element_text(size = 40, face = "bold", hjust=0.5),
                                                    panel.grid.major.x = element_blank(), 
                                                    panel.grid.minor = element_blank(),
                                                    panel.background = element_blank(), 
                                                    axis.line = element_line(colour = "black")) +
                                            scale_x_discrete(name="", labels=c("Control pH", "Medium pH", "Low pH"))+
                                            scale_y_continuous(breaks = c(0, 100), expand=c(0,0))+
                                            facet_wrap(~overlapType, labeller = labeller(overlapType = 
                                                            c("Holobiont Richness" = "Exact Sequence Variants")))

            #run GLMM
                longerdf$pH<-ordered(longerdf$pH, levels = c("A_Control", "B_Medium", "C_Low"))

                EinH<-longerdf[longerdf$overlapType== "Environmental Richness",]
  
                HinE<-longerdf[longerdf$overlapType=="Holobiont Richness",]

                Stats[[paste0("HinEGLMMpH")]]<-lmer(data = HinE, formula=value ~ pH + (1|Island))
                Stats[[paste0("HinEGLMMnopH")]]<-lmer(data = HinE, formula=value ~(1|Island))

                plot(density(HinE$value, na.rm=TRUE))

                qqnorm(residuals(Stats[["HinEGLMMpH"]]))
                scatter.smooth(residuals(Stats[["HinEGLMMpH"]]) ~ fitted(Stats[["HinEGLMMpH"]]))
                qqnorm(residuals(Stats[["HinEGLMMnopH"]]))
                scatter.smooth(residuals(Stats[["HinEGLMMnopH"]]) ~ fitted(Stats[["HinEGLMMnopH"]]))
                #assumptions look good with poisson distribution (makes sense given count data)

                anova(Stats[["HinEGLMMpH"]], Stats[["HinEGLMMnopH"]])
                    #significant efefct of pH:sample
                confs<-confint(Stats[["HinEGLMMpH"]])
                    # pH confints do not cross zero for all samples
                summary(Stats[["HinEGLMMpH"]])

                r.squaredGLMM(Stats[["HinEGLMMpH"]])

                Stats[[paste0("EinHGLMMpH")]]<-lmer(data = EinH, formula=value ~ pH + (1|Island))
                Stats[[paste0("EinHGLMMnopH")]]<-lmer(data = EinH, formula=value ~(1|Island))

                qqnorm(residuals(Stats[["EinHGLMMpH"]]))
                scatter.smooth(residuals(Stats[["EinHGLMMpH"]]) ~ fitted(Stats[["EinHGLMMpH"]]))
                qqnorm(residuals(Stats[["EinHGLMMnopH"]]))
                scatter.smooth(residuals(Stats[["EinHGLMMnopH"]]) ~ fitted(Stats[["EinHGLMMnopH"]]))
                #assumptions look good with poisson distribution (makes sense given count data)

                anova(Stats[["EinHGLMMpH"]], Stats[["EinHGLMMnopH"]])
                    #no significant efefct of pH:sample
                confs<-confint(Stats[["EinHGLMMpH"]])
                    # pH confints do not cross zero for all samples
                summary(Stats[["EinHGLMMpH"]])

                r.squaredGLMM(Stats[["EinHGLMMpH"]])

    # 5c2. Holobiome/Environmental distinct proportions chemcal 
        #chemo distinct
                OverlapByARMS<-list(

                    "proportion of environmental metabolites also in holobiont"=data.frame(
                        "A_Control"=c(  sum(taxa_names(psb_list[[7]]) %in% taxa_names(psb_list[[8]]) ) /length(taxa_names(psb_list[[7]])),
                                    sum(taxa_names(psb_list[[9]]) %in% taxa_names(psb_list[[10]]) ) /length(taxa_names(psb_list[[9]])),
                                    sum(taxa_names(psb_list[[17]]) %in% taxa_names(psb_list[[18]]) ) /length(taxa_names(psb_list[[17]])),
                                    sum(taxa_names(psb_list[[19]]) %in% taxa_names(psb_list[[20]]) ) /length(taxa_names(psb_list[[19]])),    
                                    sum(taxa_names(psb_list[[29]]) %in% taxa_names(psb_list[[30]]) ) /length(taxa_names(psb_list[[29]])), 
                                    sum(taxa_names(psb_list[[31]]) %in% taxa_names(psb_list[[32]]) ) /length(taxa_names(psb_list[[31]]))
                                 ),
                        "B_Medium"=c(  sum(taxa_names(psb_list[[5]]) %in% taxa_names(psb_list[[6]]) ) /length(taxa_names(psb_list[[5]])),
                                    sum(taxa_names(psb_list[[11]]) %in% taxa_names(psb_list[[12]]) ) /length(taxa_names(psb_list[[11]])),
                                    sum(taxa_names(psb_list[[13]]) %in% taxa_names(psb_list[[14]]) ) /length(taxa_names(psb_list[[13]])),
                                    sum(taxa_names(psb_list[[25]]) %in% taxa_names(psb_list[[26]]) ) /length(taxa_names(psb_list[[25]])),    
                                    sum(taxa_names(psb_list[[33]]) %in% taxa_names(psb_list[[34]]) ) /length(taxa_names(psb_list[[33]])), 
                                    sum(taxa_names(psb_list[[35]]) %in% taxa_names(psb_list[[36]]) ) /length(taxa_names(psb_list[[35]]))
                                 ),
                        "C_Low"=c(  sum(taxa_names(psb_list[[1]]) %in% taxa_names(psb_list[[2]]) ) /length(taxa_names(psb_list[[1]])),
                                    sum(taxa_names(psb_list[[3]]) %in% taxa_names(psb_list[[4]]) ) /length(taxa_names(psb_list[[3]])),
                                    sum(taxa_names(psb_list[[15]]) %in% taxa_names(psb_list[[16]]) ) /length(taxa_names(psb_list[[15]])),
                                    sum(taxa_names(psb_list[[21]]) %in% taxa_names(psb_list[[22]]) ) /length(taxa_names(psb_list[[21]])),    
                                    sum(taxa_names(psb_list[[23]]) %in% taxa_names(psb_list[[24]]) ) /length(taxa_names(psb_list[[23]])), 
                                    sum(taxa_names(psb_list[[27]]) %in% taxa_names(psb_list[[28]]) ) /length(taxa_names(psb_list[[27]]))
                                 )
                    ),
                    "proportion of holobiont metabolites also in environment"=data.frame(
                        "A_Control"=c(  sum(taxa_names(psb_list[[8]]) %in% taxa_names(psb_list[[7]]) ) /length(taxa_names(psb_list[[8]])),
                                    sum(taxa_names(psb_list[[10]]) %in% taxa_names(psb_list[[9]]) ) /length(taxa_names(psb_list[[10]])),
                                    sum(taxa_names(psb_list[[18]]) %in% taxa_names(psb_list[[17]]) ) /length(taxa_names(psb_list[[18]])),
                                    sum(taxa_names(psb_list[[20]]) %in% taxa_names(psb_list[[19]]) ) /length(taxa_names(psb_list[[20]])),    
                                    sum(taxa_names(psb_list[[30]]) %in% taxa_names(psb_list[[29]]) ) /length(taxa_names(psb_list[[30]])), 
                                    sum(taxa_names(psb_list[[32]]) %in% taxa_names(psb_list[[31]]) ) /length(taxa_names(psb_list[[32]]))
                                 ),
                        "B_Medium"=c(  sum(taxa_names(psb_list[[6]]) %in% taxa_names(psb_list[[5]]) ) /length(taxa_names(psb_list[[6]])),
                                    sum(taxa_names(psb_list[[12]]) %in% taxa_names(psb_list[[11]]) ) /length(taxa_names(psb_list[[12]])),
                                    sum(taxa_names(psb_list[[14]]) %in% taxa_names(psb_list[[13]]) ) /length(taxa_names(psb_list[[14]])),
                                    sum(taxa_names(psb_list[[26]]) %in% taxa_names(psb_list[[25]]) ) /length(taxa_names(psb_list[[26]])),    
                                    sum(taxa_names(psb_list[[34]]) %in% taxa_names(psb_list[[33]]) ) /length(taxa_names(psb_list[[34]])), 
                                    sum(taxa_names(psb_list[[36]]) %in% taxa_names(psb_list[[35]]) ) /length(taxa_names(psb_list[[36]]))
                                 ),
                        "C_Low"=c(  sum(taxa_names(psb_list[[2]]) %in% taxa_names(psb_list[[1]]) ) /length(taxa_names(psb_list[[2]])),
                                    sum(taxa_names(psb_list[[4]]) %in% taxa_names(psb_list[[3]]) ) /length(taxa_names(psb_list[[4]])),
                                    sum(taxa_names(psb_list[[16]]) %in% taxa_names(psb_list[[15]]) ) /length(taxa_names(psb_list[[16]])),
                                    sum(taxa_names(psb_list[[22]]) %in% taxa_names(psb_list[[21]]) ) /length(taxa_names(psb_list[[22]])),    
                                    sum(taxa_names(psb_list[[24]]) %in% taxa_names(psb_list[[23]]) ) /length(taxa_names(psb_list[[24]])), 
                                    sum(taxa_names(psb_list[[28]]) %in% taxa_names(psb_list[[27]]) ) /length(taxa_names(psb_list[[28]]))
                                 )
                    )
                )


            #pivot longer and append island id
                longdf<-cbind(pivot_longer(OverlapByARMS[[1]], everything(), names_to="pH"), pivot_longer(OverlapByARMS[[2]], everything(), names_to="pH"))
                longdf<-longdf[c(1,2,4)]
                names(longdf)<- c("pH" ,"Proportion of Metabolites from the environment also found in Holobiont",  "Proportion of Metabolites from Holobionts also found in the environment")
                longdf$Island<-c(rep("Dobu", 9), rep("Illi", 9))
                longdf$Island<-as.factor(longdf$Island)
                longdf$pH<-as.factor(longdf$pH)


                longerdf<-pivot_longer(longdf, 2:3, names_to="overlapType")
                longerdf$overlapType<-as.factor(longerdf$overlapType)
                levels(longerdf$overlapType)<-c("Holobiont Richness", "Environmental Richness")
                
                longerdf$value<-longerdf$value*100 #make percentage
                longerdf$value<-100-longerdf$value #convert percent shared to percent not shared

            #make plot    
                Figs[["MetaboliteOverlap"]] <- ggplot()+
                                            geom_hline(yintercept=50, linetype='dotted')+
                                            geom_violin(data=longerdf[longerdf$overlapType=="Holobiont Richness",], aes(x=pH, y=value),  fill="black") +
                                            #geom_violin(data=longerdf[longerdf$overlapType=="Environmental Richness",], aes(x=pH, y=value, fill=overlapType)) +
                                            geom_boxplot(data=longerdf[longerdf$overlapType=="Holobiont Richness",],aes(x=pH, y=value), width=0.1, color="grey", fill="grey", alpha=0.5) +
                                            #geom_boxplot(data=longerdf[longerdf$overlapType=="Environmental Richness",],aes(x=pH, y=value), width=0.1, color="grey", fill="grey", alpha=0.5) +
                                            labs(  # title="Holobiont Community Metabolome Distinctness - i.e. percentage of holobiont metabolome richness not shared with the environmental metabolome",
                                                    #subtitle="",
                                                    #y="Percentage of Metabolites distinct to holobiont community", 
                                                    x="pH Regime (control n=6, medium n=3, low n=5)") +
                                            scale_color_npg() +
                                            scale_fill_lancet() +
                                            annotate("segment", x = 1, xend = 2, y = 30, yend = 30, colour = "black", size=1, alpha=1) +
                                            annotate("text", x = 1.5, y = 28, size=7, label = "**") +
                                            annotate("segment", x = 1, xend = 3, y = 25, yend = 25, colour = "black", size=1, alpha=1) +
                                            annotate("text", x = 2, y = 23, size=7, label = "**") +
                                            expand_limits(y = c(0, 100)) +
                                            theme(  strip.background = element_rect(fill = "white"),
                                                    #panel.border = element_blank(),
                                                    strip.text = element_text(size=20, color="black", face="bold"),
                                                    axis.text.x = element_text(face="bold", size=20, angle=0, hjust=0.5),
                                                    axis.text.y = element_blank(),
                                                    axis.title.y = element_blank(), 
                                                    #plot.title = element_text(size = 40, face = "bold", hjust=0.5),
                                                    panel.grid.major.x = element_blank(), 
                                                    panel.grid.minor = element_blank(),
                                                    panel.background = element_blank(), 
                                                    axis.line = element_line(colour = "black")) +
                                            scale_x_discrete(name="", labels=c("Control pH", "Medium pH", "Low pH"))+
                                            scale_y_continuous(breaks = c(0, 100), expand=c(0,0))+
                                            facet_wrap(~overlapType, labeller = labeller(overlapType = 
                                                            c("Holobiont Richness" = "Metabolites")))


# 5. Supplementary Figure 1 - 16s NMDS
            ps<-ps16ALLp
            levels(sample_data(ps)$Sample)<-c("Environmental", "Algae", "Tethys Sp.", "Halisarca Sp.", "Holobiont")
            levels(sample_data(ps)$pH)<-c("Control pH", "Medium pH", "Low pH")

            distances<-phyloseq::distance(ps, method="bray") 

                    ord_ps<- ordinate(ps, "NMDS", "unifrac", weighted=FALSE, distance=distances)
                    Figs[["All16sOrd"]]<-plot_ordination(ps, ord_ps, type="Sample", color="Sample", shape="pH") +
                                        geom_point(size=3)+
                                        ggtitle(paste0("Supplementary Figure 1: NMDS of all 16s sequences by fraction (Stress=0.101)")) +
                                        stat_ellipse(aes(group=Sample)) +
                                        colScale +
                                        fillScale+
                                        theme(  strip.background = element_rect(fill = "white"),
                                                #panel.border = element_rect(colour = "black", fill = "white"),
                                                strip.text = element_text(size=14, color="black", face="bold"),
                                                axis.text = element_blank(), 
                                                axis.title = element_blank(), 
                                                plot.title = element_text(size = 20, face = "bold", hjust=0.5),
                                                panel.grid.major = element_blank(), 
                                                panel.grid.minor = element_blank(),
                                                panel.background = element_blank(), 
                                                axis.line = element_line(colour = "black"))


# 6. Supplementary Figure 2 - 23s NMDS
                levels(sample_data(ps_AlgaeHolobiomep)$pH)<-c("Control pH", "Medium pH", "Low pH")
                levels(sample_data(ps_AlgaeHolobiomep)$Site)<-c("Dobu", "Upa-Upasina")

                for (sam in c("Algae_23s", "Algae", "Plastids") ) { # remova and add 23s accordingly
                    pssam<-subset_samples(ps_AlgaeHolobiomep, Sample==sam)
                    pssam<-prune_taxa(taxa_sums(pssam)>0,pssam)
                    
                    ps_sam_distances<-phyloseq::distance(pssam, method="bray") 

                    ord_pssam <- ordinate(pssam, "NMDS", "unifrac", weighted=FALSE, distance=ps_sam_distances, trymax = 50)
                    Figs[[paste0("Ordination_", sam)]]<-plot_ordination(pssam, ord_pssam, type="Sample", color="pH", shape="Site") +
                                        geom_point(size=3)+
                                        ggtitle(paste0(sam, " Fraction")) +
                                        stat_ellipse(aes(group=pH)) +
                                        pHcolScale +
                                        pHfillScale +
                                        theme(  strip.background = element_rect(fill = "white"),
                                                #panel.border = element_rect(colour = "black", fill = "white"),
                                                strip.text = element_text(size=14, color="black", face="bold"),
                                                axis.text = element_blank(), 
                                                axis.title = element_blank(), 
                                                plot.title = element_text(size = 14, face = "bold", hjust=0.5),
                                                panel.grid.major = element_blank(), 
                                                panel.grid.minor = element_blank(),
                                                panel.background = element_blank(), 
                                                axis.line = element_line(colour = "black"))
                }  

# 7. Supplementary Figure 3 - 16s taxa barplots
        # create new ps with taxa grouped at rank 4 and labelled by top 9
            levels(sample_data(ps16p)$Sample)<-c("Environmental", "Algae", "Holobiont")
            levels(sample_data(ps16p)$Site)<-c("Dobu", "Upa-Upasina")


            ps16p4 = tax_glom(ps16p, "Rank_4")
            topotus = names(sort(taxa_sums(ps16p4), TRUE)[1:9])
            toptaxtab = cbind(tax_table(ps16p4), Class = "Other")
            toptaxtab[topotus, "Class"] <- as(tax_table(ps16p4)[topotus, "Rank_4"], "character")
            tax_table(ps16p4) <- (tax_table(toptaxtab))
        #merge samples to simplify complex plot
            sample_data(ps16p4)$mergeid<-as.factor(paste0(sample_data(ps16p4)$pH, "  ", sample_data(ps16p4)$Site, "  ",sample_data(ps16p4)$Sample))
            ps16p4_m<-merge_samples(ps16p4, "mergeid")
            sample_data(ps16p4_m)$mergeid<-levels(sample_data(ps16p4)$mergeid)
            #repair factors
                sample_data(ps16p4_m)$pH<-as.factor(unlist(lapply(strsplit(sample_data(ps16p4_m)$mergeid, "  "), '[', 1)))
                levels( sample_data(ps16p4_m)$pH)<-c("Control pH", "Medium pH", "Low pH")

                sample_data(ps16p4_m)$Site<-as.factor(unlist(lapply(strsplit(sample_data(ps16p4_m)$mergeid, "  "), '[', 2)))
                sample_data(ps16p4_m)$Sample<-as.factor(unlist(lapply(strsplit(sample_data(ps16p4_m)$mergeid, "  "), '[', 3)))
            ps16p4_mp = transform_sample_counts(ps16p4_m, function(x) 10000 * x/sum(x))

        #ensure rank 4 names differentiate "unclassified" and "unclassified bacteria"
            taxtab<-as.data.frame(tax_table(ps16p4_mp))
            taxtab$Class[as.data.frame(tax_table(ps16p4_mp))$Class=="Unclassified"] <-
                    as.data.frame(tax_table(ps16p4_mp))$Rank_3[as.data.frame(tax_table(ps16p4_mp))$Class=="Unclassified"] #rank 3 names where class is unclassified 
                                                                                                                        #(i.e. which differentiate unclassified and unclassified_bacteria) - 
                                                                                                                        #only needed until pipeline adjusted to output this way automatically
            tax_table(ps16p4_mp)<-as.matrix(taxtab)

        #remove single holobiont samples
            ps16p4_mp<-prune_samples( !(sample_data(ps16p4_mp)$Sample=="RVS" | sample_data(ps16p4_mp)$Sample=="CS") , ps16p4_mp)


        #plot
            Figs[["BacterialClasses"]]<-plot_bar(ps16p4_mp, fill="Class", x="pH", facet_grid=sample_Sample~Site) +
                                        geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack")    +
                                        scale_fill_npg() +
                                        scale_color_npg() +
                                        ggtitle("Supplementary Figure 3: Most Abundant Bacteria Classes per Site and pH Combination, Abundance Transformed") +
                                        theme(  strip.background = element_rect(fill = "white"),
                                                #panel.border = element_rect(colour = "black", fill = "white"),
                                                strip.text = element_text(size=14, color="black", face="bold"),
                                                axis.text.y = element_blank(), 
                                                axis.text.x = element_text(face="bold", size=14, angle=0, hjust=0.5), 
                                                axis.title = element_blank(), 
                                                plot.title = element_text(size = 20, face = "bold", hjust=0.5),
                                                panel.grid.major = element_blank(), 
                                                panel.grid.minor = element_blank(),
                                                panel.background = element_blank(), 
                                                axis.line = element_line(colour = "black"))
# 8. Supplementary Figure 4 - 23s taxa barplots
            levels(sample_data(ps23p)$Site)<-c("Dobu", "Upa-Upasina")
            levels( sample_data(ps23p)$pH)<-c("Control pH", "Medium pH", "Low pH")


            ps23p_rand<-subset_samples(ps23p, Code!="122")
            ps23p4 = tax_glom(ps23p_rand, "Rank_6")
            topotus = names(sort(taxa_sums(ps23p4), TRUE)[1:9])
            toptaxtab = cbind(tax_table(ps23p4), Class = "Other")
            toptaxtab[topotus, "Class"] <- as(tax_table(ps23p4)[topotus, "Rank_6"], "character")
            tax_table(ps23p4) <- (tax_table(toptaxtab))
        #plot
            Figs[["AlgalClasses"]]<-plot_bar(ps23p4, fill="Class", x="pH", facet_grid=~Site) +
                                    geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack")    +
                                    scale_fill_aaas() +
                                    scale_color_aaas() +
                                    ggtitle("Supplementary Figure 4: Most Abundant Algal Classes per ARMS, Abundance Transformed") +
                                    theme(  strip.background = element_rect(fill = "white"),
                                            #panel.border = element_rect(colour = "black", fill = "white"),
                                            strip.text = element_text(size=14, color="black", face="bold"),
                                            axis.text.y = element_blank(), 
                                            axis.text.x = element_text(face="bold", size=14, angle=0, hjust=0.5), 
                                            axis.title = element_blank(), 
                                            plot.title = element_text(size = 20, face = "bold", hjust=0.5),
                                            panel.grid.major = element_blank(), 
                                            panel.grid.minor = element_blank(),
                                            panel.background = element_blank(), 
                                            axis.line = element_line(colour = "black"))
                                    
# 9. Supplementary Figure 5 - NMDS holobiont and environmental microbiomes

            levels(sample_data(ps16p)$Site)<-c("Dobu", "Upa-Upasina")
            levels( sample_data(ps16p)$pH)<-c("Control pH", "Medium pH", "Low pH")
            names(sample_data(ps16p))[1]<-"Fraction"


                    ps16p_holos<-prune_samples(sample_data(ps16p)$Fraction=="Holobiont" | sample_data(ps16p)$Fraction=="Environmental", ps16p)
                    ps16p_holos<-prune_taxa(taxa_sums(ps16p_holos)>0,ps16p_holos)

                    ps_holo_distances<-phyloseq::distance(ps16p_holos, method="bray") 

                    ord_holos <- ordinate(ps16p_holos, "NMDS", "unifrac", weighted=FALSE, distance=ps_holo_distances, trymax = 50)
                    Figs[["Holos_Ordination"]]<-plot_ordination(ps16p_holos, ord_holos, type="Fraction", color="pH", shape="Fraction") +
                                        geom_point(size=5)+
                                        ggtitle("Supplementary Figure 5: NMDS showing increasing similarity of holobiont to environmental fraction with OA (Stress=0.14)") +
                                        pHcolScale +
                                        pHfillScale +
                                        theme(  strip.background = element_rect(fill = "white"),
                                                #panel.border = element_rect(colour = "black", fill = "white"),
                                                strip.text = element_text(size=14, color="black", face="bold"),
                                                axis.text = element_blank(), 
                                                axis.title = element_blank(), 
                                                plot.title = element_text(size = 20, face = "bold", hjust=0.5),
                                                panel.grid.major = element_blank(), 
                                                panel.grid.minor = element_blank(),
                                                panel.background = element_blank(), 
                                                axis.line = element_line(colour = "black"))

# 10. Supplementary Figure 6 - Individual microbiome distinctness

   # 5f. holbiont distinctness for specific holobiomes 

               ps_list<-list()
                    ps_holo100<-prune_samples( (sample_data(ps16p)$Sample=="Tethys Sponge Microbiome" | sample_data(ps16p)$Sample=="Halisarca Sponge Microbiome" | sample_data(ps16p)$Sample=="Environmental Microbiome"  )  , ps16p)
                    ps_holo100 <-prune_taxa(taxa_sums(ps_holo100)>0, ps_holo100)

                    ARMSwith100data<-unique(sample_data(ps_holo100)$ARMS [sample_data(ps_holo100)$Sample=="Environmental Microbiome"])
                    samps<-levels(sample_data(ps_holo100)$Sample)
                    for (ARMS in ARMSwith100data) {
                        #get all same arms
                        ps_ARMS<-prune_samples(sample_data(ps_holo100)$ARMS==ARMS, ps_holo100)
                        ps_ARMS <-prune_taxa(taxa_sums(ps_ARMS)>0, ps_ARMS)
                
                        samps<-levels(sample_data(ps_ARMS)$Sample)
                        for (samp in samps) {
                            #get all same sample
                            ps_samp<-prune_samples(sample_data( ps_ARMS)$Sample==samp, ps_ARMS)
                            ps_samp<-prune_taxa(taxa_sums(ps_samp)>0, ps_samp)

                            for (i in (1:length(sample_names(ps_samp)))) { #separate samples
                            #get 1 ps per sample replicate  
                                ps_singlesample<-prune_samples(sample_names( ps_samp)==sample_names( ps_samp)[i], ps_samp)
                                pH<- sample_data(ps_singlesample)$pH

                                ps_list[[ paste0(pH, "_", ARMS, "-", samp, "-", i) ]] <-prune_taxa(taxa_sums(ps_singlesample)>0, ps_singlesample)
                            } 
                        }
                    }

             ARMSnames<-unlist(lapply(str_split(names(ps_list), "-"), '[', 1))
            sample_names<-unlist(lapply(str_split(names(ps_list), "-"), '[', 2))


            ps_df<-data.frame(
                       ps_names= names(ps_list),
                       counter=NA,
                       ARMSnames=ARMSnames,
                       sample_names=sample_names
            )
            numComparisons=length(ARMSnames)-length(unique(ARMSnames)) #as there is 1 100 ARMS for each ARMS names 
            SpecificHolobiontDistinctness<-data.frame(matrix(
                                                            ncol= 5 ,
                                                            nrow=numComparisons,
                                                            dimnames=list(NULL,c("Individual_Holobiont","value", "ARMS", "pH", "site"))
            ))
            ControlSharedESVs<-list()
            MediumSharedESVs<-list()
            LowSharedESVs<-list()

            SharedSeqs<-list()

            counter<-0

             for (ARMSname in unique(ARMSnames)) {
                #ps_df$counter[startsWith(names(ps_list), ARMSname)]<-seq_along(startsWith(names(ps_list), ARMSname))
                ARMSindex<-ps_df$ARMSnames==ARMSname
                
                ARMS_ps_list<-ps_list[ARMSindex]
                ARMS_ps_df<-ps_df[ARMSindex,]
                if( length(ARMS_ps_list)>1 ) {
                    for (i in 2:length(names(ARMS_ps_list))) {
                        counter<-counter+1
                        SpecificHolobiontDistinctness$Individual_Holobiont[counter]<-as.character(sample_data(ARMS_ps_list[[i]])$Sample)
                        SpecificHolobiontDistinctness$value[counter]<-sum((!taxa_names(ARMS_ps_list[[i]]) %in% taxa_names(ARMS_ps_list[[1]])))/length(taxa_names(ARMS_ps_list[[i]]))
                        SpecificHolobiontDistinctness$ARMS[counter]<-as.character(sample_data(ARMS_ps_list[[1]])$ARMS[1])
                        SpecificHolobiontDistinctness$pH[counter]<- as.character(sample_data(ARMS_ps_list[[1]])$pH[1])
                        SpecificHolobiontDistinctness$site[counter]<-as.character(sample_data(ARMS_ps_list[[1]])$Site[1])

                        if (as.character(sample_data(ARMS_ps_list[[1]])$pH[1]) == "A_Con") {
                            ControlSharedESVs[[paste0(counter, as.character(sample_data(ARMS_ps_list[[1]])$pH[1]))]]<-taxa_names(ARMS_ps_list[[i]]) [ taxa_names(ARMS_ps_list[[i]]) %in% taxa_names(ARMS_ps_list[[1]]) ]
                        }
                        if (as.character(sample_data(ARMS_ps_list[[1]])$pH[1]) == "B_Med") {
                            MediumSharedESVs[[paste0(counter, as.character(sample_data(ARMS_ps_list[[1]])$pH[1]))]]<-taxa_names(ARMS_ps_list[[i]]) [ taxa_names(ARMS_ps_list[[i]]) %in% taxa_names(ARMS_ps_list[[1]]) ]
                        }
                        else if (as.character(sample_data(ARMS_ps_list[[1]])$pH[1]) == "C_Low") {
                            LowSharedESVs[[paste0(counter, as.character(sample_data(ARMS_ps_list[[1]])$pH[1]))]]<-taxa_names(ARMS_ps_list[[i]]) [ taxa_names(ARMS_ps_list[[i]]) %in% taxa_names(ARMS_ps_list[[1]]) ]
                        }
                    }
                }
             }

        SpecificHolobiontDistinctness$value<-SpecificHolobiontDistinctness$value*100
        #SpecificHolobiontDistinctness$Individual_Holobiont<-factor( SpecificHolobiontDistinctness$Individual_Holobiont , levels=c("RVS", "CS"))
        #levels(SpecificHolobiontDistinctness$Individual_Holobiont)<-c("Halisarca Sp.", "Tethys Sp.")
        SpecificHolobiontDistinctness$pH<-ordered(SpecificHolobiontDistinctness$pH, levels = c("Control pH", "Medium pH", "Low pH"))


        Figs[["SpecificHolobiontDistinctness_poster"]] <- SpecificHolobiontDistinctness %>% as_tibble %>%
                                                                    mutate_if(is.character, as.factor) %>%
                                                                    mutate(Individual_Holobiont=recode(Individual_Holobiont, 
                         "Tethys Sponge Microbiome"="Tethya Sponge",
                         "Halisarca Sponge Microbiome"="Halisarca Sponge")) %>%
                ggplot(.)+
                                                    geom_hline(yintercept=50, linetype='dotted')+
                                            #geom_violin(aes(x=pH, y=value, fill=Individual_Holobiont)) +
                                            geom_boxplot(aes(x=pH, y=value,fill=Individual_Holobiont))  +
                                            labs(   #title="Supplementary Figure 6: Individual Holobiont Microbiome Distinctness",
                                                    #subtitle="CS=Halisarca sp.; RVS=Tethys sp.",
                                                    y="Percentage ESVs distinct to sponge holobionts", 
                                                    #x="pH Regime (control n=18, medium n=9, low n=16)"
                                                    ) +
                                            colScale_individOnly +
                                            fillScale_individOnly +
                                            annotate("segment", x = 1, xend = 2, y = 30, yend = 30, colour = "black", size=1, alpha=1) +
                                            annotate("text", x = 1.5, y = 28, size=7, label = "**") +
                                            annotate("segment", x = 1, xend = 3, y = 25, yend = 25, colour = "black", size=1, alpha=1) +
                                            annotate("text", x = 2, y = 23, size=7, label = "**") +
                                            expand_limits(y = c(0, 100)) +
                                            theme(  strip.background = element_rect(fill = "white"),
                                                    #panel.border = element_blank(),
                                                    legend.position="none",
                                                    strip.text = element_text(size=14, color="black", face="bold"),
                                                    axis.text.x = element_text(face="bold", size=14, angle=0, hjust=0.5),
                                                    axis.text.y = element_text(face="bold", size=14, angle=0, hjust=0.5),
                                                    #axis.title.y = element_text(face="bold", size=20),
                                                    plot.title = element_text(size = 20, face = "bold", hjust=0.5),
                                                    panel.grid.major.x = element_blank(), 
                                                    panel.grid.minor = element_blank(),
                                                    panel.background = element_blank(), 
                                                    axis.line = element_line(colour = "black")) +
                                            scale_x_discrete(name="", labels=c("Control pH", "Medium pH", "Low pH"))+
                                            scale_y_continuous(breaks = c(0, 100), expand=c(0,0))



# 11. Supplementary Figure 7 - Sponge Microbiome Ordination
            
            ps16_sponges<-prune_samples( (sample_data(ps16ALLp)$Sample=="RVS" | sample_data(ps16ALLp)$Sample=="CS" )  , ps16ALLp)
            ps16_sponges <-prune_taxa(taxa_sums(ps16_sponges)>0, ps16_sponges)
            
            levels(sample_data(ps16_sponges)$pH)<-c("Control pH", "Medium pH", "Low pH")

            ps16_sponge_distances<-phyloseq::distance(ps16_sponges, method="bray") 

            ord_16sponges <- ordinate(ps16_sponges, "NMDS", "bray", weighted=FALSE, distance=ps16_sponge_distances, trymax = 50)
            Figs[["Sponge_Seq_Ordination"]]<-plot_ordination(ps16_sponges, ord_16sponges, type="Sample", color="pH", shape="Sample") +
                                geom_point(size=5)+
                                ggtitle("Supplementary Figure 7: NMDS showing sponge holobiont microbiome disimilarity (Stress = 0.074)") +
                                #stat_ellipse(aes(group=pH)) +
                                facet_wrap(~Sample, labeller = labeller(Sample = 
                                                    c("CS" = "Tethys Sp.",
                                                    "RVS" = "Halisarca Sp."))) +
                                pHcolScale +
                                theme(  strip.background = element_rect(fill = "white"),
                                        #panel.border = element_rect(colour = "black", fill = "white"),
                                        strip.text = element_text(size=14, color="black", face="bold"),
                                        axis.text = element_blank(), 
                                        axis.title = element_blank(), 
                                        plot.title = element_text(size = 20, face = "bold", hjust=0.5),
                                        panel.grid.major = element_blank(), 
                                        panel.grid.minor = element_blank(),
                                        panel.background = element_blank(), 
                                        axis.line = element_line(colour = "black"))

# 12. Supplementary Figure 8 - Individual metabolome distinctness

            btab_holo100<-prune_samples( (sample_data(btab_p)$Sample=="Tethys Sponge Metabolome" | sample_data(btab_p)$Sample=="Halisarca Sponge Metabolome" | sample_data(btab_p)$Sample=="Environmental Metabolome"  )  , btab_p)
            btab_holo100 <-prune_taxa(taxa_sums(btab_holo100)>0, btab_holo100)

            ARMSunits<-unique(sample_data(btab_holo100)$ARMS)

            df<-data.frame(ARMS=rep(NA,100), Fraction=NA, Island=NA, pH=NA, PropDistinct=NA)
            counter<-1

            for (ARMSunit in ARMSunits) {
                temp<-prune_samples (sample_data(btab_holo100)$ARMS==ARMSunit , btab_holo100)
                temp <-prune_taxa(taxa_sums(temp)>0, temp)

                if (1==sum("Environmental Metabolome" == sample_data(temp)$Sample) & dim(sample_data(temp))[1]>1 ) { #i.e. if there is exactly one 100 fraction for an arms unit and at least one other sample
                    EnvSample<-prune_samples (sample_data(temp)$Sample=="Environmental Metabolome" , temp) 
                    EnvSample <-prune_taxa(taxa_sums(EnvSample)>0, EnvSample)

                    HoloSamples<-prune_samples (sample_data(temp)$Sample!="Environmental Metabolome" , temp) 
                    HoloSamples <-prune_taxa(taxa_sums(HoloSamples)>0, HoloSamples)
                    VialNums<-sample_names(HoloSamples)
                    for (vial in VialNums) {
                        OneSample<-prune_samples (sample_names(HoloSamples)==vial , HoloSamples) 
                        OneSample <-prune_taxa(taxa_sums(OneSample)>0, OneSample)

                        df$ARMS[counter]<-ARMSunit
                        df$Fraction[counter]<-as.character(sample_data(OneSample)$Sample)
                        df$Island[counter]<-as.character(sample_data(OneSample)$Site)
                        df$pH[counter]<-as.character(sample_data(OneSample)$pH)
                        df$PropDistinct[counter]<-sum( !taxa_names(OneSample) %in% taxa_names(EnvSample) ) / length(taxa_names(OneSample))*100

                        counter<-counter+1
                    }
                }
            }        


        df<-df[!is.na(df[,1]),] #remove nas
        df$pH<-ordered(df$pH, levels = c("Control pH", "Medium pH", "Low pH"))


        #df$Fraction<-factor( df$Fraction , levels=c("RVS", "CS"))
        #levels(df$Fraction)<-c("Halisarca Sp.", "Tethys Sp.")

            #make plot    
                Figs[["IndividualHolobiontMetabDistinctness"]] <- df %>% as_tibble %>%
                                                                    mutate_if(is.character, as.factor) %>%
                                                                    mutate(Fraction=recode(Fraction, 
                         "Tethys Sponge Metabolome"="Tethya Sponge",
                         "Halisarca Sponge Metabolome"="Halisarca Sponge")) %>%
                ggplot(.)+
                                            geom_hline(yintercept=50, linetype='dotted') +
                                            geom_boxplot( aes(x=pH, y=PropDistinct,fill=Fraction))  +
                                            labs(  # title="Supplementary Figure 8: Individual Holobiont Metabolome Distinctness",
                                                    #subtitle="CS=Halisarca sp.; RVS=Tethys sp.",
                                                    y="Percentage compounds distinct to sponge holobionts", 
                                                    #x="pH Regime (control n=18, medium n=9, low n=16)"
                                                    ) +
                                            ChemcolScale_individOnly +
                                            ChemfillScale_individOnly +
                                            expand_limits(y = c(0, 100)) +
                                            theme(  strip.background = element_rect(fill = "white"),
                                                    #panel.border = element_blank(),
                                                    strip.text = element_text(size=14, color="black", face="bold"),
                                                    axis.text.x = element_text(face="bold", size=14, angle=0, hjust=0.5),
                                                    axis.text.y = element_text(face="bold", size=14, angle=0, hjust=0.5),
                                                    #axis.title.y = element_text(face="bold", size=20),
                                                    plot.title = element_text(size = 20, face = "bold", hjust=0.5),
                                                    panel.grid.major.x = element_blank(), 
                                                    panel.grid.minor = element_blank(),
                                                    panel.background = element_blank(), 
                                                    axis.line = element_line(colour = "black")) +
                                            scale_x_discrete(name="", labels=c("Control pH", "Medium pH", "Low pH"))+
                                            scale_y_continuous(breaks = c(0, 100), expand=c(0,0))



# MAP

    maps<-list()
    pngMAP_df = ggmap::get_map(location = c(140, -12, 155, 2), 
    source = "google", zoom = 7)
    maps[[1]]<-ggmap::ggmap(pngMAP_df) +
                    annotate("rect", xmin=150.5, xmax=151, ymin=-10 , ymax=-9.5, color="red", alpha=0)+ 
                    ggsn::scalebar( x.min = 139.9, x.max = 154.3,  y.min = -11.5, y.max = 0, dist = 200, dist_unit = "km",
                                    transform = TRUE, model = "WGS84") +
                    #ggtitle("Figure 1: Location of study sites, Dobu and Upa-Upasina")+
                    theme(  plot.title =    element_text(size = 20, face = "bold", hjust=0),
                                            axis.title.x = element_blank(),
                                            axis.title.y = element_blank())

            
    
      #  annotate(geom = "text", x = -90, y = 26, label = "Gulf of Mexico", 
     #   fontface = "italic", color = "grey22", size = 6) +

    siteMAP_df = ggmap::get_map(location = c(150.5, -10, 151, -9.5), source = "google", zoom = 12)
    maps[[2]]<-ggmap::ggmap(siteMAP_df)+ 
                    ggsn::scalebar( x.min = 150.6, x.max = 150.98,  y.min = -9.98, y.max = -9.5, dist = 10, dist_unit = "km",
                                transform = TRUE, model = "WGS84")+
                    geom_point(x=150.878843, y=-9.746104, size=5) + 
                    geom_label(label="Dobu", x=150.908843, y=-9.746104, size=5)+
                    geom_point(x=150.827345, y=-9.829434, size=5) +
                    geom_label(label="Upa-Upasina", x=150.885345, y=-9.829434, size=5)+
                    theme(  plot.title =    element_text(size = 20, face = "bold", hjust=0),
                                            axis.title.x = element_blank(),
                                            axis.title.y = element_blank())

            

# 13. Make pdf

pdf(file = file.path(path,"Figs",paste0("PNGPaper_DataGeneratedFigures_",Sys.Date(),".pdf")), width=20, height =8 ) # The height of the plot in inches        

plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
text(5, 8, "Data-generated Figures for Manuscript")
text(5, 7, "Ecosystem stress through the holobiome lens: decline of distinct holobionts under ocean acidification")
text(5, 6, paste0("Jake Williams   ", Sys.Date()))

            grid::grid.newpage()
            grid::grid.draw(cbind( ggplotGrob(maps[[1]]) , ggplotGrob(maps[[2]]) ))

            grid::grid.newpage()
            grid::grid.draw(rbind(ggplotGrob(Figs[["Figure1a"]]), ggplotGrob(Figs[["Figure1b"]])))

            grid.arrange(cbind(ggplotGrob(Figs[["ESVOverlap_poster"]]),  ggplotGrob(Figs[["MetaboliteOverlap"]])), top=textGrob("Figure 4: Multiomic Distinctness", gp=gpar(fontsize = 20, fontface = "bold")))

            Figs[["All16sOrd"]]


            grid.arrange(grobs = list(  Figs[["Ordination_Algae_23s"]]+ theme(legend.position="none"), 
                        Figs[["Ordination_Algae"]]+ theme(legend.position="none"), 
                        Figs[["Ordination_Plastids"]]
                        ), ncol = 3,
                        top=textGrob("Supplementary Figure 2: NMDS of all algae holobiome fractions by pH (Stress=0.12; 0.11; 0.067)", gp=gpar(fontsize = 20, fontface = "bold")))

            Figs[["BacterialClasses"]]

            Figs[["AlgalClasses"]]

            Figs[["Holos_Ordination"]]

            Figs[["SpecificHolobiontDistinctness_poster"]]

            Figs[["Sponge_Seq_Ordination"]]

            Figs[["IndividualHolobiontMetabDistinctness"]]

dev.off()


    CommunityDistinctnessPlot<-egg::ggarrange(Figs[["ESVOverlap_poster"]], Figs[["MetaboliteOverlap"]], nrow=1)# labels=(c("A", "B")),  # from MakePdf.r
                                    #top="Figure 5: Multiomic Distinctness")

saveRDS(CommunityDistinctnessPlot, file=file.path(path, "Outputs", "CommunityDistinctnessPlot.RDS"))

    SpongeDistinctnessPlot<-egg::ggarrange(Figs[["SpecificHolobiontDistinctness_poster"]], Figs[["IndividualHolobiontMetabDistinctness"]], nrow=1)#, labels=(c("A", "B")),  # from MakePdf.r
                                   # top="Figure S5: Sponge Multiomic Distinctness")

saveRDS(SpongeDistinctnessPlot, file=file.path(path, "Outputs", "SpongeDistinctnessPlot.RDS"))
