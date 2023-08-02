#-------------------------
# Distinctness Analysis
#-------------------------
#packages
    library(ggplot2)
    library(phyloseq)
    library(tidyverse)
    library(egg)
    library(grid)

#path
    path<-"/home/j/Dropbox/PNG_Paper" #provide the path to the project folder

#load data and transform to proportions
    ps16<-readRDS(file=file.path(path, "Outputs", "ps16.RDS"))
    btab<-readRDS(file=file.path(path, "Outputs", "btab.RDS"))
    ps16p = transform_sample_counts(ps16, function(x) 10000 * x/sum(x))
    btab_p = transform_sample_counts(btab, function(x) 10000 * x/sum(x))
    #rarefied alternative - effects still significant - no qualtitiative difference
        #ps16p<-rarefy_even_depth(  physeq=ps16, sample.size = 73982, rngseed = 123, 
        #                            replace = TRUE, trimOTUs = TRUE, verbose = TRUE) 
        #btab_p<-rarefy_even_depth(  physeq=btab, sample.size = 100000, rngseed = 123, 
        #                            replace = TRUE, trimOTUs = TRUE, verbose = TRUE) 

#make color palette
    myColors<- c("black", "black", "deeppink", "deeppink", "darkgoldenrod1", "darkgoldenrod1")
    names(myColors) <- c(   "Holobiont Community Microbiome", 
                            "Holobiont Community Metabolome", 
                            "Halisarca Sponge Microbiome", 
                            "Halisarca Sponge Metabolome",
                            "Tethya Sponge Microbiome",
                            "Tethya Sponge Metabolome")
    MyColScale <- scale_colour_manual(name = "Fraction",values = myColors)
    MyFillScale <- scale_fill_manual(name = "Fraction",values = myColors)

#functions
    MakeDistinctnessTibble<-function(ps, envSample, sampletype2, sampletype3="NoSampleType3Selected"){
                ps_list<-list()
                        ps_holo100<-prune_samples( (sample_data(ps)$Sample==envSample | sample_data(ps)$Sample==sampletype2 | sample_data(ps)$Sample==sampletype3  )  , ps)
                        ps_holo100 <-prune_taxa(taxa_sums(ps_holo100)>0, ps_holo100)

                        ARMSwith100data<-unique(sample_data(ps_holo100)$ARMS [sample_data(ps_holo100)$Sample==envSample])
                        samps<-c(envSample, sampletype2, if (sampletype3!="NoSampleType3Selected") {sampletype3})
                        for (ARMS in ARMSwith100data) {
                            #get all same arms
                            ps_ARMS<-prune_samples(sample_data(ps_holo100)$ARMS==ARMS, ps_holo100)
                            ps_ARMS <-prune_taxa(taxa_sums(ps_ARMS)>0, ps_ARMS)

                            if (sum( samps %in% sample_data( ps_ARMS)$Sample )>1 ) {
                                subsamps<-sample_data( ps_ARMS)$Sample 
                                for (subsamp in subsamps) {
                                    #get all same sample
                                    ps_samp<-prune_samples(sample_data( ps_ARMS)$Sample==subsamp, ps_ARMS)
                                    ps_samp<-prune_taxa(taxa_sums(ps_samp)>0, ps_samp)

                                    for (i in (1:length(sample_names(ps_samp)))) { #separate samples
                                    #get 1 ps per sample replicate  
                                        ps_singlesample<-prune_samples(sample_names( ps_samp)==sample_names( ps_samp)[i], ps_samp)
                                        pH<- sample_data(ps_singlesample)$pH

                                        ps_list[[ paste0(pH, "_", ARMS, "-", subsamp, "-", i) ]] <-prune_taxa(taxa_sums(ps_singlesample)>0, ps_singlesample)
                                    } 
                                }
                            }    
                        }

                ps_df<-data.frame(
                        ps_names= names(ps_list),
                        counter=NA,
                        ARMSnames=unlist(lapply(str_split(names(ps_list), "-"), '[', 1)),
                        sample_names=unlist(lapply(str_split(names(ps_list), "-"), '[', 2))
                )
                numComparisons=length(ps_df$ARMSnames)-length(unique(ps_df$ARMSnames)) #as there is 1 100 ARMS for each ARMS names 
                SpecificHolobiontDistinctness<-data.frame(matrix(
                                                                ncol= 7 ,
                                                                nrow=numComparisons,
                                                                dimnames=list(NULL,c("Individual_Holobiont","value", "RichnessHolobiontDistinct", "RichnessHolobiontOverlapping", "ARMS", "pH", "site"))
                ))
                ControlSharedESVs<-list()
                MediumSharedESVs<-list()
                LowSharedESVs<-list()

                SharedSeqs<-list()

                counter<-0

                for (ARMSname in unique(ps_df$ARMSnames)) {
                    #ps_df$counter[startsWith(names(ps_list), ARMSname)]<-seq_along(startsWith(names(ps_list), ARMSname))
                    ARMSindex<-ps_df$ARMSnames==ARMSname
                    
                    ARMS_ps_list<-ps_list[ARMSindex]
                    ARMS_ps_df<-ps_df[ARMSindex,]
                    if( length(ARMS_ps_list)>1 ) {
                        ARMS_ps_list<- c( ARMS_ps_list[grepl(envSample, names(ARMS_ps_list), fixed=TRUE)],
                                        ARMS_ps_list[!grepl(envSample, names(ARMS_ps_list), fixed=TRUE)]) #ensure envSample is first sample
                        for (i in 2:length(names(ARMS_ps_list))) {
                            counter<-counter+1
                            SpecificHolobiontDistinctness$Individual_Holobiont[counter]<-as.character(sample_data(ARMS_ps_list[[i]])$Sample)
                            SpecificHolobiontDistinctness$value[counter]<-sum((!taxa_names(ARMS_ps_list[[i]]) %in% taxa_names(ARMS_ps_list[[1]])))/length(taxa_names(ARMS_ps_list[[i]]))
                            
                            SpecificHolobiontDistinctness$RichnessHolobiontDistinct[counter]<-sum((!taxa_names(ARMS_ps_list[[i]]) %in% taxa_names(ARMS_ps_list[[1]])))
                            SpecificHolobiontDistinctness$RichnessHolobiontOverlapping[counter]<-sum((taxa_names(ARMS_ps_list[[i]]) %in% taxa_names(ARMS_ps_list[[1]])))

                            SpecificHolobiontDistinctness$ARMS[counter]<-as.character(sample_data(ARMS_ps_list[[1]])$ARMS[1])
                            SpecificHolobiontDistinctness$pH[counter]<- as.character(sample_data(ARMS_ps_list[[1]])$pH[1])
                            SpecificHolobiontDistinctness$site[counter]<-as.character(sample_data(ARMS_ps_list[[1]])$Site[1])
                        }
                    }
                }

            

                SpecificHolobiontDistinctness$site<-as.factor(SpecificHolobiontDistinctness$site)
                SpecificHolobiontDistinctness$pH<-as.factor(SpecificHolobiontDistinctness$pH)
                SpecificHolobiontDistinctness$value<-SpecificHolobiontDistinctness$value*100 #make percentage
                SpecificHolobiontDistinctness$pH<-ordered(SpecificHolobiontDistinctness$pH, levels = c("Control pH", "Medium pH", "Low pH"))
            #fix naming of tethya sponge if needed
                if ("Tethys Sponge Microbiome" %in% SpecificHolobiontDistinctness$Individual_Holobiont ) {
                    SpecificHolobiontDistinctness$Individual_Holobiont[SpecificHolobiontDistinctness$Individual_Holobiont=="Tethys Sponge Microbiome"]<-"Tethya Sponge Microbiome"
                }
                if ("Tethys Sponge Metabolome" %in% SpecificHolobiontDistinctness$Individual_Holobiont) {
                    SpecificHolobiontDistinctness$Individual_Holobiont[SpecificHolobiontDistinctness$Individual_Holobiont=="Tethys Sponge Metabolome"]<-"Tethya Sponge Metabolome"
                }
            return(SpecificHolobiontDistinctness)
    }

    PlotDistinctness<-function(DistinctnessTibble, significant, type){
        dodge=position_dodge(width = 0.2)

        p<- DistinctnessTibble %>%
            ggplot(aes(x=pH, y=value, group=interaction(Individual_Holobiont, pH), fill=Individual_Holobiont))+
                geom_hline(yintercept=50, linetype='dotted')+
                #geom_violin(aes(x=pH, y=value, fill=Individual_Holobiont), position=dodge) +
                #geom_violin(data=longerdf[longerdf$overlapType=="Environmental Richness",], aes(x=pH, y=value, fill=overlapType)) +
                geom_boxplot() + #, width=0.1, alpha=0.5,  position = dodge) +
                #geom_boxplot(data=longerdf[longerdf$overlapType=="Environmental Richness",],aes(x=pH, y=value), width=0.1, color="grey", fill="grey", alpha=0.5) +
                MyColScale +
                MyFillScale +
                expand_limits(y = c(0, 100)) +
                theme(  legend.position="none",
                        strip.background = element_rect(fill = "white"),
                        #panel.border = element_blank(),
                        strip.text = element_text(size=20, color="black", face="bold"),
                        axis.text.x = element_text(face="bold", size=16, angle=0, hjust=0.5, color="black"),
                        axis.text.y = element_text(face="bold", size=16, angle=0, hjust=0.5),
                        axis.title.y = element_text(face="bold", size=16),
                        legend.text = element_text(face="bold", size=16),
                        #plot.title = element_text(size = 40, face = "bold", hjust=0.5),
                        panel.grid.major.x = element_blank(), 
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(), 
                        axis.line = element_line(colour = "black")) +
                scale_x_discrete(name="", labels=c("Control pH", "Medium pH", "Low pH"))+
                scale_y_continuous(breaks = c(0, 100), expand=c(0,0))

        if(significant){
            p<-p    +
                annotate("segment", x = 1, xend = 2, y = 30, yend = 30, colour = "black", size=1, alpha=1) +
                annotate("text", x = 1.5, y = 27, size=7, label = "**") +
                annotate("segment", x = 1, xend = 3, y = 25, yend = 25, colour = "black", size=1, alpha=1) +
                annotate("text", x = 2, y = 22, size=7, label = "**") 
        }
        if (type=="Community"){
            p<-p+labs(   #title="Figure 4: Multiomic Distinctness",
                    #subtitle="",
                    y="Percentage distinct to\nbenthic holobiont community", 
                    #x="pH Regime (control n=6, medium n=3, low n=5)"
                    ) 
        } else if (type =="Organism") {
            p<-p+labs(   #title="Figure 4: Multiomic Distinctness",
                    #subtitle="",
                    y="Percentage distinct\nto individual holobiont", 
                    #x="pH Regime (control n=6, medium n=3, low n=5)"
                    )    
        }
        return(p)

    }


#create distinctness tibbles
    DistinctnessTibble_CommunityMicrobiome<-MakeDistinctnessTibble(ps=ps16p, envSample="Environmental Microbiome", sampletype2="Holobiont Community Microbiome", sampletype3="NoSampleType3Selected")
    DistinctnessTibble_CommunityMetabolome<-MakeDistinctnessTibble(ps=btab_p, envSample="Environmental Metabolome", sampletype2="Holobiont Community Metabolome", sampletype3="NoSampleType3Selected")
    DistinctnessTibble_SpongesMicrobiomes<-MakeDistinctnessTibble(ps=ps16p, envSample="Environmental Microbiome", sampletype2="Tethys Sponge Microbiome", sampletype3="Halisarca Sponge Microbiome")
    DistinctnessTibble_SpongesMetabolomes<-MakeDistinctnessTibble(ps=btab_p, envSample="Environmental Metabolome", sampletype2="Tethys Sponge Metabolome", sampletype3="Halisarca Sponge Metabolome")

#make plots
    p1<-PlotDistinctness(DistinctnessTibble_CommunityMicrobiome, type="Community", significant=TRUE) + 
        theme(axis.text.x = element_blank()) + 
        facet_wrap(~Individual_Holobiont, labeller = labeller(Individual_Holobiont = c("Holobiont Community Microbiome" = "Amplicon Sequence Variants")))
    p2<-PlotDistinctness(DistinctnessTibble_CommunityMetabolome, type="Community", significant=TRUE) + 
        theme(axis.text.x = element_blank(), axis.text.y=element_blank() , axis.title.y=element_blank())+
        facet_wrap(~Individual_Holobiont, labeller = labeller(Individual_Holobiont = c("Holobiont Community Metabolome" = "Metabolites")))
    p3<-PlotDistinctness(DistinctnessTibble_SpongesMicrobiomes, type="Organism", significant=TRUE) 
    p4<-PlotDistinctness(DistinctnessTibble_SpongesMetabolomes, type="Organism", significant=FALSE) + 
        theme(axis.text.y = element_blank() , axis.title.y=element_blank(), legend.position = c(0.8, 0.2), legend.title=element_blank())+
        scale_fill_manual(labels=c('Halisarca Sp.', 'Tethya Sp.'), values =c("deeppink", "darkgoldenrod1"))


# run models
        #community sequences distinctness
                DistinctnessTibble_CommunityMicrobiome$pH<-ordered(DistinctnessTibble_CommunityMicrobiome$pH, levels = c("Control pH", "Medium pH", "Low pH"))

                Stats<-lme4::lmer(data = DistinctnessTibble_CommunityMicrobiome, formula=value ~ pH + (1|site))
                Stats_nopH<-lme4::lmer(data = DistinctnessTibble_CommunityMicrobiome, formula=value ~(1|site))

                plot(density(DistinctnessTibble_CommunityMicrobiome$value, na.rm=TRUE))

                qqnorm(residuals(Stats))
                scatter.smooth(residuals(Stats) ~ fitted(Stats))
                qqnorm(residuals(Stats_nopH))
                scatter.smooth(residuals(Stats_nopH) ~ fitted(Stats_nopH))
                #assumptions look good with poisson distribution (makes sense given count data)

                anova(Stats, Stats_nopH)
                    #significant efefct of pH:sample
                confs<-confint(Stats)
                    # pH confints do not cross zero for all samples
                summary(Stats)

                MuMIn::r.squaredGLMM(Stats)
        #sponge seqeunces distinctness
                  #run GLMM
                Stats<-list()
                DistinctnessTibble_SpongesMicrobiomes$pH<-ordered(DistinctnessTibble_SpongesMicrobiomes$pH, levels = c("Control pH", "Medium pH", "Low pH"))


                Stats[[paste0("pH")]]<-lme4::lmer(data = DistinctnessTibble_SpongesMicrobiomes, formula=value ~ pH + (1|site/ARMS) + (1|Individual_Holobiont) )
                Stats[[paste0("nopH")]]<-lme4::lmer(data = DistinctnessTibble_SpongesMicrobiomes, formula=value ~ (1|site/ARMS) + (1|Individual_Holobiont) )

                plot(density(DistinctnessTibble_SpongesMicrobiomes$value, na.rm=TRUE))

                qqnorm(residuals(Stats[["pH"]]))
                scatter.smooth(residuals(Stats[["pH"]]) ~ fitted(Stats[["pH"]]))
                qqnorm(residuals(Stats[["nopH"]]))
                scatter.smooth(residuals(Stats[["nopH"]]) ~ fitted(Stats[["nopH"]]))
                #assumptions look good with poisson distribution (makes sense given count data)

                anova(Stats[["pH"]], Stats[["nopH"]])
                    #significant efefct of pH:sample
                confs<-confint(Stats[["pH"]])
                    # pH confints do not cross zero for all samples
                summary(Stats[["pH"]])

                MuMIn::r.squaredGLMM(Stats[["pH"]])

#additional analysis of absolute quantities
    DistinctnessTibble_CommunityMicrobiome %>%
        ggplot() +
            geom_boxplot(aes(x=pH, y=RichnessHolobiontDistinct, color=site))
    DistinctnessTibble_CommunityMicrobiome %>%
        ggplot() +
            geom_boxplot(aes(x=pH, y=RichnessHolobiontOverlapping))

#save figure
    jpeg(file.path(path,"DistinctnessBoxplot.jpeg"), height = 8.3, width = 11.7, units = 'in', res = 300)
        egg::ggarrange(p1, p2, p3, p4, labels=c("A", "B", "C", "D"), label.args = list(gp=grid::gpar(face="bold"), hjust=-2, vjust=2))
    dev.off()

#end