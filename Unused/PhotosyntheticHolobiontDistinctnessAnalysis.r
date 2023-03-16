# 1. Packages and path to project directory ######
    set.seed(0.1)
    
    library(tidyverse)
    library(ggplot2)
    library(ggsci)
    library(phyloseq)
    library(breakaway)
    library(DivNet)

    library(permute)
    library(vegan)
    library(DESeq2)

    #params
        path<-"/home/j/Dropbox/PNG_Paper" #provide the path to the project folder
        theme_set(theme_linedraw())
        Figs<-list()
        Stats<-list()
        options(scipen=999)


# X. Create a custom color scale
    myColors<- c("darkred", "darkturquoise", "chartreuse", "deeppink", "darkgoldenrod1")
    names(myColors) <- c("Holobiont Community Microbiome", "Environmental Microbiome", "Photosynthetic Community Microbiome", "Halisarca Sponge Microbiome", "Tethys Sponge Microbiome")
    colScale <- scale_colour_manual(name = "Fraction",values = myColors)
    fillScale <- scale_fill_manual(name = "Fraction",values = myColors)

    myColors_individOnly<-myColors[4:5]
    colScale_individOnly <- scale_colour_manual(name = "Fraction",values = myColors_individOnly)
    fillScale_individOnly <- scale_fill_manual(name = "Fraction",values = myColors_individOnly)

    names(myColors_individOnly)<-c("Halisarca Sponge Metabolome", "Tethys Sponge Metabolome")
    chemColScale_individOnly <- scale_colour_manual(name = "Fraction",values = myColors_individOnly)
    chemFillScale_individOnly <- scale_fill_manual(name = "Fraction",values = myColors_individOnly)

    pHColors<- c("green", "yellow", "orange")
    names(pHColors) <- c("Control pH", "Medium pH", "Low pH")
    pHcolScale <- scale_colour_manual(name = "pH Regime",values = pHColors)
    pHfillScale <- scale_fill_manual(name = "pH Regime",values = pHColors)

# 2. Data prep
    # get data
        ps16<-readRDS(file=file.path(path, "Outputs", "ps16.RDS"))
        btab<-readRDS(file=file.path(path, "Outputs", "btab.RDS"))

# 4. Figure 4 - distinctness
    # make initial lists - 
        #all sessile
            #by arms unit
                ps_list<-list()
                    ps_s100<-prune_samples( (sample_data(ps16)$Sample=="Holobiont Community Microbiome" | sample_data(ps16)$Sample=="Environmental Microbiome") , ps16)
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
                    psb_s100<-prune_samples( (sample_data(btab)$Sample=="Holobiont Community Metabolome" | sample_data(btab)$Sample=="Environmental Metabolome") & sample_data(btab)$Solvent=="MeOH" , btab)
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
        #photosynthetic
            #by arms unit
                ps_list_alg<-list()
                    ps_s100<-prune_samples( (sample_data(ps16)$Sample=="Photosynthetic Community Microbiome" | sample_data(ps16)$Sample=="Environmental Microbiome") , ps16)
                    SharedARMS<-sort(unique(sample_data(ps_s100)$ARMS)) [ table(sample_data(ps_s100)$ARMS)>1 ]
                    samps<-levels(sample_data(ps_s100)$Sample)
                    for (ARMS in SharedARMS) {
                        for (samp in samps) {
                                ps_temp<-prune_samples(sample_data(ps_s100)$ARMS==ARMS, ps_s100)
                                ps_temp <-prune_taxa(taxa_sums(ps_temp)>0, ps_temp)

                                ps_temp<-prune_samples(sample_data( ps_temp)$Sample==samp, ps_temp)
                                pH<- sample_data(ps_temp)$pH
                                ps_list_alg[[ paste0(pH, "_", ARMS, "_", samp) ]] <-prune_taxa(taxa_sums(ps_temp)>0, ps_temp)

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
                psb_list_alg<-list()
                    psb_s100<-prune_samples( (sample_data(btab)$Sample=="Photosynthetic Community Metabolome" | sample_data(btab)$Sample=="Environmental Metabolome") & sample_data(btab)$Solvent=="MeOH" , btab)
                    #sample_data(psb_s100)[sample_data(psb_s100)$ARMS=="29"]
                    SharedARMS<-sort(unique(sample_data(psb_s100)$ARMS)) [ table(sample_data(psb_s100)$ARMS)>1 ]
                    samps<-unique(sample_data(psb_s100)$Sample)
                    for (ARMS in SharedARMS) {
                        for (samp in samps) {
                                psb_temp<-prune_samples(sample_data(psb_s100)$ARMS==ARMS, psb_s100)
                                psb_temp <-prune_taxa(taxa_sums(psb_temp)>0, psb_temp)

                                psb_temp<-prune_samples(sample_data( psb_temp)$Sample==samp, psb_temp)
                                pH<- unique(sample_data(psb_temp)$pH)
                                psb_list_alg[[ paste0(pH, "_", ARMS, "_", samp) ]] <-prune_taxa(taxa_sums(psb_temp)>0, psb_temp)

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

    # Photosynthetic/Environmental distinct proportions sequence
        #Sequences distinct
                OverlapByArms<-list(
                    "Proportion of environmental microbial taxa in holobionts"=data.frame(
                        "A_Control"=c( sum((taxa_names(ps_list_alg[[5]]) %in% taxa_names(ps_list_alg[[6]])))/length(taxa_names(ps_list_alg[[5]])),
                                    sum((taxa_names(ps_list_alg[[7]]) %in% taxa_names(ps_list_alg[[8]])))/length(taxa_names(ps_list_alg[[7]])),
                                    sum((taxa_names(ps_list_alg[[11]]) %in% taxa_names(ps_list_alg[[12]])))/length(taxa_names(ps_list_alg[[11]])),
                                    sum((taxa_names(ps_list_alg[[13]]) %in% taxa_names(ps_list_alg[[14]])))/length(taxa_names(ps_list_alg[[13]])),
                                    sum((taxa_names(ps_list_alg[[21]]) %in% taxa_names(ps_list_alg[[22]])))/length(taxa_names(ps_list_alg[[21]])),
                                    sum((taxa_names(ps_list_alg[[23]]) %in% taxa_names(ps_list_alg[[24]])))/length(taxa_names(ps_list_alg[[23]]))
                                    ),
                        "B_Medium"=c(sum((taxa_names(ps_list_alg[[9]]) %in% taxa_names(ps_list_alg[[10]])))/length(taxa_names(ps_list_alg[[9]])), 
                                sum((taxa_names(ps_list_alg[[25]]) %in% taxa_names(ps_list_alg[[26]])))/length(taxa_names(ps_list_alg[[25]])),
                                sum((taxa_names(ps_list_alg[[27]]) %in% taxa_names(ps_list_alg[[28]])))/length(taxa_names(ps_list_alg[[27]])),
                                NA,
                                NA,
                                NA
                                ),
                        "C_Low"=c(sum((taxa_names(ps_list_alg[[1]]) %in% taxa_names(ps_list_alg[[2]])))/length(taxa_names(ps_list_alg[[1]])), 
                                sum((taxa_names(ps_list_alg[[3]]) %in% taxa_names(ps_list_alg[[4]])))/length(taxa_names(ps_list_alg[[3]])),
                                sum((taxa_names(ps_list_alg[[15]]) %in% taxa_names(ps_list_alg[[16]])))/length(taxa_names(ps_list_alg[[15]])),
                                sum((taxa_names(ps_list_alg[[17]]) %in% taxa_names(ps_list_alg[[18]])))/length(taxa_names(ps_list_alg[[17]])),
                                sum((taxa_names(ps_list_alg[[19]]) %in% taxa_names(ps_list_alg[[20]])))/length(taxa_names(ps_list_alg[[19]])),
                                NA
                                )
                    ),
                    "Proportion of holobiont microbial taxa in environment"=data.frame(
                        "A_Control"=c( sum((taxa_names(ps_list_alg[[6]]) %in% taxa_names(ps_list_alg[[5]])))/length(taxa_names(ps_list_alg[[6]])),
                                    sum((taxa_names(ps_list_alg[[8]]) %in% taxa_names(ps_list_alg[[7]])))/length(taxa_names(ps_list_alg[[8]])),
                                    sum((taxa_names(ps_list_alg[[12]]) %in% taxa_names(ps_list_alg[[11]])))/length(taxa_names(ps_list_alg[[12]])),
                                    sum((taxa_names(ps_list_alg[[14]]) %in% taxa_names(ps_list_alg[[13]])))/length(taxa_names(ps_list_alg[[14]])),
                                    sum((taxa_names(ps_list_alg[[22]]) %in% taxa_names(ps_list_alg[[21]])))/length(taxa_names(ps_list_alg[[22]])),
                                    sum((taxa_names(ps_list_alg[[24]]) %in% taxa_names(ps_list_alg[[23]])))/length(taxa_names(ps_list_alg[[24]]))
                                    ),
                        "B_Medium"=c(sum((taxa_names(ps_list_alg[[10]]) %in% taxa_names(ps_list_alg[[9]])))/length(taxa_names(ps_list_alg[[10]])), 
                                sum((taxa_names(ps_list_alg[[26]]) %in% taxa_names(ps_list_alg[[25]])))/length(taxa_names(ps_list_alg[[26]])),
                                sum((taxa_names(ps_list_alg[[28]]) %in% taxa_names(ps_list_alg[[27]])))/length(taxa_names(ps_list_alg[[28]])),
                                NA,
                                NA,
                                NA
                                ),
                        "C_Low"=c(sum((taxa_names(ps_list_alg[[2]]) %in% taxa_names(ps_list_alg[[1]])))/length(taxa_names(ps_list_alg[[2]])), 
                                sum((taxa_names(ps_list_alg[[4]]) %in% taxa_names(ps_list_alg[[3]])))/length(taxa_names(ps_list_alg[[4]])),
                                sum((taxa_names(ps_list_alg[[16]]) %in% taxa_names(ps_list_alg[[15]])))/length(taxa_names(ps_list_alg[[16]])),
                                sum((taxa_names(ps_list_alg[[18]]) %in% taxa_names(ps_list_alg[[17]])))/length(taxa_names(ps_list_alg[[18]])),
                                sum((taxa_names(ps_list_alg[[20]]) %in% taxa_names(ps_list_alg[[19]])))/length(taxa_names(ps_list_alg[[20]])),
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


                Figs[["ESVOverlap_poster_alg"]] <- ggplot()+
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
                                            #annotate("segment", x = 1, xend = 2, y = 30, yend = 30, colour = "black", size=1, alpha=1) +
                                            #annotate("text", x = 1.5, y = 28, size=7, label = "**") +
                                            #annotate("segment", x = 1, xend = 3, y = 25, yend = 25, colour = "black", size=1, alpha=1) +
                                            #annotate("text", x = 2, y = 23, size=7, label = "**") +
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


    #  Holobiome/Environmental distinct proportions chemcal 
        #Photosynthetic distinct
                OverlapByARMS<-list(

                    "proportion of environmental metabolites also in holobiont"=data.frame(
                        "A_Control"=c(  sum(taxa_names(psb_list_alg[[3]]) %in% taxa_names(psb_list_alg[[4]]) ) /length(taxa_names(psb_list_alg[[3]])),
                                    sum(taxa_names(psb_list_alg[[11]]) %in% taxa_names(psb_list_alg[[12]]) ) /length(taxa_names(psb_list_alg[[11]])),
                                    sum(taxa_names(psb_list_alg[[13]]) %in% taxa_names(psb_list_alg[[14]]) ) /length(taxa_names(psb_list_alg[[13]])),
                                    sum(taxa_names(psb_list_alg[[23]]) %in% taxa_names(psb_list_alg[[24]]) ) /length(taxa_names(psb_list_alg[[23]])), 
                                    sum(taxa_names(psb_list_alg[[25]]) %in% taxa_names(psb_list_alg[[26]]) ) /length(taxa_names(psb_list_alg[[25]]))  
                                 ),
                        "B_Medium"=c(  sum(taxa_names(psb_list_alg[[5]]) %in% taxa_names(psb_list_alg[[6]]) ) /length(taxa_names(psb_list_alg[[5]])),
                                    sum(taxa_names(psb_list_alg[[7]]) %in% taxa_names(psb_list_alg[[8]]) ) /length(taxa_names(psb_list_alg[[7]])),
                                    sum(taxa_names(psb_list_alg[[19]]) %in% taxa_names(psb_list_alg[[20]]) ) /length(taxa_names(psb_list_alg[[19]])), 
                                    sum(taxa_names(psb_list_alg[[27]]) %in% taxa_names(psb_list_alg[[28]]) ) /length(taxa_names(psb_list_alg[[27]])),
                                    sum(taxa_names(psb_list_alg[[29]]) %in% taxa_names(psb_list_alg[[30]]) ) /length(taxa_names(psb_list_alg[[29]]))
                                 ),
                        "C_Low"=c(  sum(taxa_names(psb_list_alg[[1]]) %in% taxa_names(psb_list_alg[[2]]) ) /length(taxa_names(psb_list_alg[[1]])),
                                    sum(taxa_names(psb_list_alg[[9]]) %in% taxa_names(psb_list_alg[[10]]) ) /length(taxa_names(psb_list_alg[[9]])),
                                    sum(taxa_names(psb_list_alg[[15]]) %in% taxa_names(psb_list_alg[[16]]) ) /length(taxa_names(psb_list_alg[[15]])),
                                    sum(taxa_names(psb_list_alg[[17]]) %in% taxa_names(psb_list_alg[[18]]) ) /length(taxa_names(psb_list_alg[[17]])),
                                    sum(taxa_names(psb_list_alg[[21]]) %in% taxa_names(psb_list_alg[[22]]) ) /length(taxa_names(psb_list_alg[[21]])) 
                                 )
                    ),
                    "proportion of holobiont metabolites also in environment"=data.frame(
                        "A_Control"=c(  sum(taxa_names(psb_list_alg[[4]]) %in% taxa_names(psb_list_alg[[3]]) ) /length(taxa_names(psb_list_alg[[4]])),
                                    sum(taxa_names(psb_list_alg[[12]]) %in% taxa_names(psb_list_alg[[11]]) ) /length(taxa_names(psb_list_alg[[12]])),
                                    sum(taxa_names(psb_list_alg[[14]]) %in% taxa_names(psb_list_alg[[13]]) ) /length(taxa_names(psb_list_alg[[14]])),
                                    sum(taxa_names(psb_list_alg[[24]]) %in% taxa_names(psb_list_alg[[23]]) ) /length(taxa_names(psb_list_alg[[24]])), 
                                    sum(taxa_names(psb_list_alg[[26]]) %in% taxa_names(psb_list_alg[[25]]) ) /length(taxa_names(psb_list_alg[[26]]))
                                 ),
                        "B_Medium"=c(  sum(taxa_names(psb_list_alg[[6]]) %in% taxa_names(psb_list_alg[[5]]) ) /length(taxa_names(psb_list_alg[[6]])),
                                    sum(taxa_names(psb_list_alg[[8]]) %in% taxa_names(psb_list_alg[[7]]) ) /length(taxa_names(psb_list_alg[[8]])),
                                    sum(taxa_names(psb_list_alg[[20]]) %in% taxa_names(psb_list_alg[[19]]) ) /length(taxa_names(psb_list_alg[[20]])),    
                                    sum(taxa_names(psb_list_alg[[28]]) %in% taxa_names(psb_list_alg[[27]]) ) /length(taxa_names(psb_list_alg[[28]])),
                                    sum(taxa_names(psb_list_alg[[30]]) %in% taxa_names(psb_list_alg[[29]]) ) /length(taxa_names(psb_list_alg[[30]]))
                                 ),
                        "C_Low"=c(  sum(taxa_names(psb_list_alg[[2]]) %in% taxa_names(psb_list_alg[[1]]) ) /length(taxa_names(psb_list_alg[[2]])),
                                    sum(taxa_names(psb_list_alg[[10]]) %in% taxa_names(psb_list_alg[[9]]) ) /length(taxa_names(psb_list_alg[[10]])),
                                    sum(taxa_names(psb_list_alg[[16]]) %in% taxa_names(psb_list_alg[[15]]) ) /length(taxa_names(psb_list_alg[[16]])),
                                    sum(taxa_names(psb_list_alg[[18]]) %in% taxa_names(psb_list_alg[[17]]) ) /length(taxa_names(psb_list_alg[[18]])),
                                    sum(taxa_names(psb_list_alg[[22]]) %in% taxa_names(psb_list_alg[[21]]) ) /length(taxa_names(psb_list_alg[[22]])) 
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
                Figs[["MetaboliteOverlap_alg"]] <- ggplot()+
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
                                            #annotate("segment", x = 1, xend = 2, y = 30, yend = 30, colour = "black", size=1, alpha=1) +
                                            #annotate("text", x = 1.5, y = 28, size=7, label = "**") +
                                            #annotate("segment", x = 1, xend = 3, y = 25, yend = 25, colour = "black", size=1, alpha=1) +
                                            #annotate("text", x = 2, y = 23, size=7, label = "**") +
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


    jpeg(file=file.path(path, "Outputs","PNG_jpegs",  "PhotoSyntheticHolobiontDistinctness.jpeg"), height = 8.3, width = 13, units = 'in', res = 300)
        egg::ggarrange(Figs[["ESVOverlap_poster_alg"]], Figs[["MetaboliteOverlap_alg"]], ncol=2)
    dev.off()                
