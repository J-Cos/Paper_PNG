# 1. Packages and path to project directory ######
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
    library(ggsci)
    library(indicspecies)
    library(ggsignif)
    library(lme4)
    library(MuMIn)
    
    #params
        path<-"/home/j/Dropbox/PNG_Paper" #provide the path to the project folder
        theme_set(theme_linedraw())
        Figs<-list()
        Stats<-list()



# 2. Data prep
    # 2a. Sequencing Data
        # format as phyloseq
            #16s
                metadata_16s<-sample_data(read.table(file=file.path(path, "Data", "16s_metadata.csv"), sep=",", header=TRUE, row.names=1))
                ps16<-SeqDataTable2Phyloseq( SeqDataTablePath=file.path(path, "Data", "PNGFullTest_16s_SeqDataTable.RDS"),
                                            clustering="ESV",
                                            Metadata=metadata_16s)

            #23s
                metadata_23s<-sample_data(read.table(file=file.path(path, "Data", "23s_metadata.csv"), sep=",", header=TRUE, row.names=1))
                ps23<-SeqDataTable2Phyloseq( SeqDataTablePath=file.path(path, "Data", "PNGFullTest_23s_SeqDataTable.RDS"),
                                            clustering="ESV",
                                            Metadata=metadata_23s)
        
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
        btab<-read.table(file.path(path, "Data", "Metabolomic Data", "METABOLOMICS-SNETS-V2-9965902b-download_cluster_buckettable-main.tsv"),sep="\t", comment.char = "", check.names=FALSE, header=TRUE, row.names=1)
        x<-strsplit( names(btab), "_")
        names(btab)<-paste0("Vial_", as.numeric(sapply(x, '[', 2)))

        metadata<-read.csv(file.path(path, "Data", "Metabolomic Data", "MetabolomeMetadata.csv"), row.names=1)

        btab<-otu_table(btab, taxa_are_rows=TRUE)
        btab<-phyloseq(btab, sample_data(metadata))

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



# 3. Compositional Analysis - Ordinations and permanovas for chemical and species richness: MODEL 1 (a, b and c) + SI general ordination plots
    # 3a ordinations - 16 together
            ps<-ps16ALLp
            distances<-phyloseq::distance(ps, method="bray") 

                    ord_ps<- ordinate(ps, "NMDS", "unifrac", weighted=FALSE, distance=distances)
                    Figs[["All16sOrd"]]<-plot_ordination(ps, ord_ps, type="Sample", color="Sample", shape="pH") +
                                        geom_point(size=3)+
                                        ggtitle(paste0("NMDS of all 16s sequences by fraction - Stress=0.101")) +
                                        stat_ellipse(aes(group=Sample)) +
                                        scale_fill_lancet() +
                                        scale_color_lancet()
    
    # 3b. ordinations - sequences - every sample seperately ordinations

            methods<-list(c("PCoA", "unifrac"))
            for (method in 1:length(methods) ) {
                for (sam in c("Sessile", "100", "Algae") ) { # remova and add 23s accordingly
                    ps16p_sam<-subset_samples(ps16p, Sample==sam)
                    ps16p_sam<-prune_taxa(taxa_sums(ps16p_sam)>0,ps16p_sam)
                    
                    ps16p_sam_distances<-phyloseq::distance(ps16p_sam, method="bray") 

                    ord_ps16p_sam <- ordinate(ps16p_sam, methods[[method]][1], methods[[method]][2], weighted=FALSE, distance=ps16p_sam_distances)
                    Figs[[paste0("Ordination_", sam, "_", methods[[method]][1])]]<-plot_ordination(ps16p_sam, ord_ps16p_sam, type="Sample", color="pH", shape="Site") +
                                        geom_point(size=3)+
                                        ggtitle(paste0("Exact Sequence Variants - ", sam, " Fraction")) +
                                        stat_ellipse(aes(group=pH))+
                                        facet_wrap(~Sample,2) +
                                        scale_fill_lancet() +
                                        scale_color_lancet()
                }
                
                ps23p_distances<-phyloseq::distance(ps23p, method="bray") 

                ord_ps23p <- ordinate(ps23p, methods[[method]][1], methods[[method]][2], weighted=FALSE, distance=ps23p_distances)
                Figs[["Ordination_Algae_23s_PCOA"]]<-plot_ordination(ps23p, ord_ps23p, type="Sample", color="pH", shape="Site") +
                                    geom_point(size=3)+
                                    ggtitle(paste0("Exact Sequence Variants - ", "Algae 23s Fraction")) +
                                    stat_ellipse(aes(group=pH))+
                                    facet_wrap(~Sample,2)    +
                                    scale_fill_lancet() +
                                    scale_color_lancet()
                    
            }

    # 3c. ordinations algae holobiome only
                  methods<-list(c("PCoA", "unifrac"))
            for (method in 1:length(methods) ) {
                for (sam in c("Algae_23s", "Algae", "Plastids") ) { # remova and add 23s accordingly
                    pssam<-subset_samples(ps_AlgaeHolobiomep, Sample==sam)
                    pssam<-prune_taxa(taxa_sums(pssam)>0,pssam)
                    
                    ps_sam_distances<-phyloseq::distance(pssam, method="bray") 

                    ord_pssam <- ordinate(pssam, methods[[method]][1], methods[[method]][2], weighted=FALSE, distance=ps_sam_distances)
                    Figs[[paste0("Ordination_", sam, "_", methods[[method]][1])]]<-plot_ordination(pssam, ord_pssam, type="Sample", color="pH", shape="Site") +
                                        geom_point(size=3)+
                                        ggtitle(paste0("Exact Sequence Variants - ", sam, " Fraction")) +
                                        stat_ellipse(aes(group=pH))+
                                        facet_wrap(~Sample,2) +
                                        scale_fill_lancet() +
                                        scale_color_lancet()
                }  
            }
    # 3d. ordinations - chemical - every sample seperately
            methods<-list(c("PCoA", "unifrac"))
            for (method in 1:length(methods) ) {
                for (sam in c("Algae", "SESSILE", "100-500UM") ) {
                    btab_p_sam<-subset_samples(btab_p, Sample==sam)
                    btab_p_sam<-prune_taxa(taxa_sums(btab_p_sam)>0,btab_p_sam)
                    
                    btab_p_sam_distances<-phyloseq::distance(btab_p_sam, method="bray") 

                    ord_btab_p_sam <- ordinate(btab_p_sam, methods[[method]][1], methods[[method]][2], weighted=FALSE, distance=btab_p_sam_distances)
                    Figs[[paste0("Ordination_" , sam, methods[[method]][1])]]<-plot_ordination(btab_p_sam, ord_btab_p_sam, type="Sample", color="pH", shape= "Site") +
                                        geom_point(size=3)+
                                        ggtitle(paste0("Metabolites - ", sam, " Fraction")) +
                                        stat_ellipse(aes(group=pH))+
                                        facet_wrap(~Sample,2)    +
                                        scale_fill_lancet() +
                                        scale_color_lancet()
                }
            }


    # 3e. Permanovas (MODEL 1 (a+b+c))
        #permanova functions
            PermanovaByFraction<-function(ps, fraction) {
                ps_temp<-prune_samples(sample_data(ps)$Sample==fraction, ps)
                ps_temp<-prune_taxa(taxa_sums(ps_temp)>0,ps_temp)     
                df = as(sample_data(ps_temp), "data.frame")
                psP_distances<-phyloseq::distance(ps_temp, method="jaccard")  #bray (abund) or jaccard (p/a) - [results the same so using jaccard]
                output<- adonis(psP_distances ~  Site + pH + Site*pH , df, strata=df$Site)
                return(output)
            }
            PermanovaForCS<-function(ps, fraction="CS") {
                ps_temp<-prune_samples(sample_data(ps)$Sample==fraction, ps)
                ps_temp<-prune_taxa(taxa_sums(ps_temp)>0,ps_temp)     
                df = as(sample_data(ps_temp), "data.frame")
                psP_distances<-phyloseq::distance(ps_temp, method="jaccard") #bray (abund) or jaccard (p/a) - [results the same so using jaccard]
                output<- adonis(psP_distances ~  pH  , df)
                return(output)
            }
        #overall permanova - 16s sequences - are the fractions really different?
            ps<-ps16ALLp
            df = as(sample_data(ps), "data.frame")
            ps_distances<-phyloseq::distance(ps, method="bray")           
            output<- adonis(ps_distances ~  Sample +pH + Site +Sample:pH , df, strata=df$Site)
            #yes compoitional difference between fractions
        #overall permanova -metabolites - are fractions different?   
            ps<-btab_p
            df = as(sample_data(ps), "data.frame")
            ps_distances<-phyloseq::distance(ps, method="bray")           
            output<- adonis(ps_distances ~  Sample +pH + Site +Sample:pH , df, strata=df$Site)
            #yes compoitional difference between fractions
        #permanovas - 16s sequences - every sample seperately - now we know they are, run fraction specific permanovas
            ps<-ps16ALLp
            for (fraction in levels(sample_data(ps)$Sample)[c(1,2,4,5)]) {
               Stats[[paste0("PERMANOVA_Sequence_", fraction)]]<-PermanovaByFraction(ps=ps, fraction=fraction)
            }
            #CS dealt with out of loop as only have this fraction from one site    
                Stats[[paste0("PERMANOVA_Sequence_CS")]]<-PermanovaForCS(ps)
            #seqs
                #100 - Site ph and ph:Site
                #algae - site and pH 
                #RVS - site, ph and ph:site
                #sessile - site, ph and ph:site
                #cs - none

        #permanovas - metabs - every sample seperately - now we know they are, run fraction specific permanovas
            ps<-btab_p
            for (fraction in levels(sample_data(ps)$Sample)[c(1,2,4,5)]) {
               Stats[[paste0("PERMANOVA_Metabolite_", fraction)]]<-PermanovaByFraction(ps=ps, fraction=fraction)
            }
            #CS dealt with out of loop as only have this fraction from one site    
                Stats[[paste0("PERMANOVA_Metabolite_CS")]]<-PermanovaForCS(ps)
            #metabs
                #100 - none
                #algae - site and ph
                #rvs - site and pH
                #sessile - site and ph: site
                #cs none

        #permanova - 23s sequences composiitonal shift in algal macrobes
            ps<-ps23    
            df = as(sample_data(ps), "data.frame")
            ps_distances<-phyloseq::distance(ps, method="bray") 
            output<- adonis(ps_distances ~  Site + pH + Site*pH , df, strata=df$Site)
            #significant effect of pH and Site

        #permanova - plastid sequences composiitonal shift
            ps<-pschloro_a
            df = as(sample_data(ps), "data.frame")
            ps_distances<-phyloseq::distance(ps, method="bray") 
            output<- adonis(ps_distances ~  Site + pH + Site*pH , df, strata=df$Site)
            #significant effect of pH and Site
# 4. Richness Analysis: MODEL 2 and 5 + FIGURE 1 + SI complete richness plot
    # 4a. richness GLMM for seqs - GLMM (model 2)
            ps<-psALL
            alpha_df<-estimate_richness(ps, split = TRUE, measures = NULL) %>%
            merge(sample_data(ps), by=0) 
            alpha_df$pH<-ordered(alpha_df$pH, levels = c("A_Con", "B_Med", "C_Low"))
            Stats[[paste0("FullGLMMpH")]]<-glmer.nb(data = alpha_df, formula=Observed ~ 0+ Sample+ pH:Sample + (1|Site))
            Stats[[paste0("FullGLMMnopH")]]<-glmer.nb(data = alpha_df, formula=Observed ~ 0+ Sample +(1|Site))

            qqnorm(residuals(Stats[["FullGLMMpH"]]))
            scatter.smooth(residuals(Stats[["FullGLMMpH"]]) ~ fitted(Stats[["FullGLMMpH"]]))
            qqnorm(residuals(Stats[["FullGLMMnopH"]]))
            scatter.smooth(residuals(Stats[["FullGLMMnopH"]]) ~ fitted(Stats[["FullGLMMnopH"]]))
            #assumptions look good with poisson distribution (makes sense given count data)

            anova(Stats[["FullGLMMpH"]], Stats[["FullGLMMnopH"]])
                #significant efefct of pH:sample
            confs<-confint(Stats[["FullGLMMpH"]])
                # pH confints do not cross zero for all ALGAL samples
            summary(Stats[["FullGLMMpH"]])
                #significant effect of all lagal holobi0ome seqs (linear in pH)

            r.squaredGLMM(Stats[["FullGLMMpH"]])

    # 4b. richness GLMM for metabolites - GLMM (model 5)
            alpha_dfmetab<-estimate_richness(btab, split = TRUE, measures = NULL) %>%
            merge(sample_data(btab), by=0)
            alpha_dfmetab$pH<-ordered(alpha_dfmetab$pH, levels = c("A_Con", "B_Med", "C_Low"))
            
            Stats[[paste0("MetaboliteGLMMpH")]]<-glmer.nb(data = alpha_dfmetab, formula=Observed ~ 0+ Sample+ pH:Sample + (1|Site))
            Stats[[paste0("MetaboliteGLMMnopH")]]<-glmer.nb(data = alpha_dfmetab, formula=Observed ~ 0+ Sample +(1|Site))

            qqnorm(residuals(Stats[["MetaboliteGLMMpH"]]))
            scatter.smooth(residuals(Stats[["MetaboliteGLMMpH"]]) ~ fitted(Stats[["MetaboliteGLMMpH"]]))
            qqnorm(residuals(Stats[["MetaboliteGLMMnopH"]]))
            scatter.smooth(residuals(Stats[["MetaboliteGLMMnopH"]]) ~ fitted(Stats[["MetaboliteGLMMnopH"]]))
            #assumptions look good with poisson distribution (makes sense given count data)

            anova(Stats[["MetaboliteGLMMpH"]], Stats[["MetaboliteGLMMnopH"]])
                #ph:sample not significant
            confs<-confint(Stats[["MetaboliteGLMMpH"]])
                # pH confints do not cross zero for RVS only
            summary(Stats[["MetaboliteGLMMpH"]])
                #RVS:pH significant (linear in pH)


    # 4c. Violin plots - Species and chemical Richness for key fractions
        # seq violin richness Figure1a
            ps_s100<-prune_samples(sample_data(ps16)$Sample=="Sessile" | sample_data(ps16)$Sample=="100" | sample_data(ps16)$Sample=="Algae", ps16)
            ps_s100<-prune_taxa(taxa_sums(ps_s100)>0, ps_s100)
            levels(sample_data(ps_s100)$Sample)<-c("Environmental", "Algae", "Holobiont")
            Figs[["Figure1a"]]<-plot_richness( ps_s100, x="pH", measures=c("Observed"))+
                    geom_violin(aes(fill=Sample, color=Sample)) +
                    ggtitle(paste0("Community Richness by pH regime")) +
                    scale_color_npg() + 
                    scale_fill_npg() +
                    ylab("Richness") +
                    guides(color="none")+
                    guides(fill="none")
            Figs[["Figure1a"]]$layers<-Figs[["Figure1a"]]$layers[-1]  

        # chemo violin richness Figure1b
            btab_nosponge<-prune_samples( !(sample_data(btab)$Sample=="RVS" | sample_data(btab)$Sample=="CS" | sample_data(btab)$Sample=="TS") , btab)
            levels(sample_data(btab_nosponge)$Sample)<-c("Environmental", "Algae", "Holobiont")

            Figs[["Figure1b"]]<-plot_richness(btab_nosponge, x="pH", measures=c("Observed"))+
                            geom_violin(aes(fill=Sample, color=Sample)) +
                            ggtitle(paste0(" Chemical Richness by pH regime")) +
                            scale_color_npg() + 
                            scale_fill_npg() +
                            ylab("Richness") +
                            guides(color="none")
            Figs[["Figure1b"]]$layers<-Figs[["Figure1b"]]$layers[-1]      

    # 4d. SI box plots - richness for everything
        #seq boxplot      
            Figs[["AllSeqRichness"]]<-plot_richness( psALL, x="pH", measures=c("Observed"))+
                    geom_boxplot(aes(fill=Sample, color=Sample)) +
                    ggtitle(paste0("Community Richness by pH regime")) +
                    scale_color_npg() + 
                    scale_fill_npg() +
                    ylab("Richness") +
                    guides(color="none")
            Figs[["AllSeqRichness"]]$layers<-Figs[["AllSeqRichness"]]$layers[-1]    
        #metab boxplot
            Figs[["AllMetabRichness"]]<-plot_richness(btab, x="pH", measures=c("Observed"))+
                            geom_boxplot(aes(fill=Sample, color=Sample)) +
                            ggtitle(paste0(" Chemical Richness by pH regime")) +
                            scale_color_npg() + 
                            scale_fill_npg() +
                            ylab("Richness") +
                            guides(color="none")
            Figs[["AllMetabRichness"]]$layers<-Figs[["AllMetabRichness"]]$layers[-1]    



    


# 5. Holobiont Distinctness: MODEL 3 AND 4 + FIGURE 2 +SI plot showing increasing compositional similarity
    # BUNK - 5a. Ratio of environmental to holobiont microbial richness (plots and stan_glm)
        psP<-ps16p

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
        #helper code for correctly extracting desired ps from ps_list
        #        names(ps_list)
        #
        #        for (i in (1: length(ps_list)) ) {
        #            print(sample_data(ps_list[[i]])[,4])
        #        }


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
                
             
                Figs[["SequenceRatios"]]<-  ratiosByARMS[[1]] %>%
                                            pivot_longer(everything(), names_to="pH") %>%
                                            ggplot(aes(x=pH, y=value, fill=pH)) + 
                                                geom_hline(yintercept=1) + ylim(0, 2)+
                                                geom_boxplot() +
                                                scale_color_npg() +
                                                scale_fill_npg() +
                                                labs(   title = "Ratio of environmental to holobiont microbial richness for each ARMS",,
                                                        subtitle="",
                                                        y="Ratio of environmental to holobiont microbial richness", x="pH Regime")

        #pivot longer and append island id
                longdf<-pivot_longer(ratiosByARMS[[1]], everything(), names_to="pH")
                longdf$Island<-c(rep("Dobu", 4), "Illi", rep("Dobu", 2), rep("Illi", 3), NA, rep("Illi", 2), NA, rep("Illi", 2), NA, NA)
                longdf$Island<-as.factor(longdf$Island)
                longdf$pH<-as.factor(longdf$pH)

        #bayes glm
                model_general <- stan_glm(data = longdf, formula=longdf[["value"]] ~ pH+Island, family=gaussian())
                Figs$s_CIplot_general<-plot(model_general, plotfun = "areas", prob = 0.9, pars="beta")
                Figs$s_posteriors_general <- describe_posterior(model_general)




    # BUNK - 5b. Same but for chemical ratios (plots and stan_glm)


        #by arms unit
                psb_list<-list()
                    psb_s100<-prune_samples( (sample_data(btab_p)$Sample=="SESSILE" | sample_data(btab_p)$Sample=="100-500UM") & sample_data(btab_p)$Solvent=="MeOH" , btab_p)
                    #sample_data(psb_s100)[sample_data(psb_s100)$ARMS=="29"]
                    SharedARMS<-levels(sample_data(psb_s100)$ARMS) [ table(sample_data(psb_s100)$ARMS)>1 ]
                    samps<-levels(sample_data(psb_s100)$Sample)
                    for (ARMS in SharedARMS) {
                        for (samp in samps) {
                                psb_temp<-prune_samples(sample_data(psb_s100)$ARMS==ARMS, psb_s100)
                                psb_temp <-prune_taxa(taxa_sums(psb_temp)>0, psb_temp)

                                psb_temp<-prune_samples(sample_data( psb_temp)$Sample==samp, psb_temp)
                                pH<- sample_data(psb_temp)$pH
                                psb_list[[ paste0(pH, "_", ARMS, "_", samp) ]] <-prune_taxa(taxa_sums(psb_temp)>0, psb_temp)

                        }

                    }

                names(psb_list)

                for (i in (1: length(psb_list)) ) {
                    print(paste ( i, unlist(sample_data(psb_list[[i]])[,4])))
                }


                ratiosByARMS<-list(
                    "Ratio of environmental to holobiont CHEMICAL richness"=data.frame(
                        "A_Control"=c(length(phyloseq::taxa_names(psb_list[[7]])) /length(taxa_names(psb_list[[8]])), 
                                    length(taxa_names(psb_list[[9]])) /length(taxa_names(psb_list[[10]])),
                                    length(taxa_names(psb_list[[17]])) /length(taxa_names(psb_list[[18]])),
                                    length(taxa_names(psb_list[[19]])) /length(taxa_names(psb_list[[20]])),
                                    length(taxa_names(psb_list[[29]])) /length(taxa_names(psb_list[[30]])),
                                    length(taxa_names(psb_list[[31]])) /length(taxa_names(psb_list[[32]]))
                                    ),
                        "B_Medium"=c(length(taxa_names(psb_list[[5]])) /length(taxa_names(psb_list[[6]])), 
                                    length(taxa_names(psb_list[[11]])) /length(taxa_names(psb_list[[12]])),
                                    length(taxa_names(psb_list[[13]])) /length(taxa_names(psb_list[[14]])),
                                    length(taxa_names(psb_list[[25]])) /length(taxa_names(psb_list[[26]])),
                                    length(taxa_names(psb_list[[33]])) /length(taxa_names(psb_list[[34]])),
                                    length(taxa_names(psb_list[[35]])) /length(taxa_names(psb_list[[36]]))                                    ),
                        "C_Low"=c(length(taxa_names(psb_list[[1]])) /length(taxa_names(psb_list[[2]])), 
                                    length(taxa_names(psb_list[[3]])) /length(taxa_names(psb_list[[4]])),
                                    length(taxa_names(psb_list[[15]])) /length(taxa_names(psb_list[[16]])),
                                    length(taxa_names(psb_list[[21]])) /length(taxa_names(psb_list[[22]])),
                                    length(taxa_names(psb_list[[23]])) /length(taxa_names(psb_list[[24]])),
                                    length(taxa_names(psb_list[[27]])) /length(taxa_names(psb_list[[28]]))                                    )
                    )
                )




        #pivot longer and append island id
                longdf<-pivot_longer(ratiosByARMS[[1]], everything(), names_to="pH")
                longdf$Island<-c(rep("Dobu", 9), rep("Illi", 9))
                longdf$Island<-as.factor(longdf$Island)
                longdf$pH<-as.factor(longdf$pH)


                
                Figs[["ChemRatios"]]<-  ggplot(longdf, aes(x=pH, y=value, fill=pH)) + 
                                                geom_boxplot() +
                                                #facet_wrap(~Island) +
                                                geom_hline(yintercept=1) + ylim(0, 2)+
                                                scale_color_npg() +
                                                scale_fill_npg() +
                                                labs(   title = "Ratio of environmental to holobiont chemical richness for each ARMS",,
                                                        subtitle="",
                                                        y="Ratio of environmental to holobiont chemical richness", x="pH Regime")
         
   
   
        #bayes glm
                model_general <- stan_glm(data = longdf, formula=longdf[["value"]] ~ pH+Island, family=gaussian())
                Figs$c_CIplot_general<-plot(model_general, plotfun = "areas", prob = 0.9, pars="beta")
                Figs$c_posteriors_general <- describe_posterior(model_general)

        grid.arrange(grobs = list(Figs[["SequenceRatios"]], Figs[["ChemRatios"]]), ncol = 2) ## display plot





    # 5c. Holobiome/Environmental distinct proportions
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

            #make plot    
                Figs[["ESVOverlap"]] <- ggplot()+
                                            geom_hline(yintercept=50, linetype='dotted')+
                                            geom_violin(data=longerdf[longerdf$overlapType=="Holobiont Richness",], aes(x=pH, y=value, fill=overlapType)) +
                                            geom_violin(data=longerdf[longerdf$overlapType=="Environmental Richness",], aes(x=pH, y=value, fill=overlapType)) +
                                            geom_boxplot(data=longerdf[longerdf$overlapType=="Holobiont Richness",],aes(x=pH, y=value), width=0.1, color="grey", fill="grey", alpha=0.5) +
                                            geom_boxplot(data=longerdf[longerdf$overlapType=="Environmental Richness",],aes(x=pH, y=value), width=0.1, color="grey", fill="grey", alpha=0.5) +
                                            labs(   title="Microbiome Distinctness - i.e. percentage of community richness not shared between the environmental microbiome and holobiont microbiome",
                                                    subtitle="",
                                                    y="Percentage of community richness distinct to microbiome", x="pH Regime (control n=6, medium n=3, low n=5)") +
                                            scale_color_npg() +
                                            scale_fill_lancet() +
                                            annotate("segment", x = 1, xend = 2, y = 33.5, yend = 33.5, colour = "black", size=1, alpha=1) +
                                            annotate("text", x = 1.5, y = 32.75, size=7, label = "**") +
                                            annotate("segment", x = 1, xend = 3, y = 32.5, yend = 32.5, colour = "black", size=1, alpha=1) +
                                            annotate("text", x = 2, y = 31.75, size=7, label = "**")



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
                    #significant efefct of pH:sample
                confs<-confint(Stats[["EinHGLMMpH"]])
                    # pH confints do not cross zero for all samples
                summary(Stats[["EinHGLMMpH"]])
    
    # 5d. new idea - compositional characteristics of holobiont microbes lost and gained with OA
        ps<-ps16

        ps_holos<-prune_samples(sample_data(ps)$Sample=="Sessile", ps)

        GetComparisonList<-function(Island){
            ps_holos_Dobu<-prune_samples(sample_data(ps_holos)$Site==Island, ps_holos)
            ps_holos_Dobu<-prune_taxa(taxa_sums(ps_holos_Dobu)>0, ps_holos_Dobu)
            ps_list<-list()
            for (i in sample_names(ps_holos_Dobu)) {
                ps_temp<-prune_samples(sample_names(ps_holos_Dobu)==i, ps_holos_Dobu)
                ps_list[[i]]<-prune_taxa(taxa_sums(ps_temp)>0, ps_temp)
            }
            controls<-which(sample_data(ps_holos_Dobu)$pH=="A_Con")
            acidifieds<-which(!sample_data(ps_holos_Dobu)$pH=="A_Con")
            mediums<-which(sample_data(ps_holos_Dobu)$pH=="B_Med")
            lows<-which(sample_data(ps_holos_Dobu)$pH=="B_Low")
            holobiontchange_list<-list()
            for (control in controls) {
                for (acidified in acidifieds) {
                comparison_name<-paste0(sample_names(ps_holos_Dobu)[[control]], "-", sample_names(ps_holos_Dobu)[[acidified]])
                kept_indices<-taxa_names(ps_list[[control]]) %in% taxa_names(ps_list[[acidified]])
                lost_indices<-!taxa_names(ps_list[[control]]) %in% taxa_names(ps_list[[acidified]])  
                holobiontchange_list[[ paste0(comparison_name, ":kept") ]]<- prune_taxa(kept_indices, ps_list[[control]])
                sample_names(holobiontchange_list[[ paste0(comparison_name, ":kept") ]])<-paste0(comparison_name, ":kept")
                holobiontchange_list[[ paste0(comparison_name, ":lost") ]]<- prune_taxa(lost_indices, ps_list[[control]])
                sample_names(holobiontchange_list[[ paste0(comparison_name, ":lost") ]])<-paste0(comparison_name, ":lost")

                gained_indices<-!taxa_names(ps_list[[acidified]]) %in% taxa_names(ps_list[[control]])
                holobiontchange_list[[ paste0(comparison_name, ":gained") ]]<- prune_taxa(gained_indices, ps_list[[acidified]])
                sample_names(holobiontchange_list[[ paste0(comparison_name, ":gained") ]])<-paste0(comparison_name, ":gained")
                }
            }
            return(holobiontchange_list)
        }

        fullList<-c(GetComparisonList("Dobu"),        GetComparisonList("Illi"))

        merged<-fullList[[1]]
        for (i in 2:length(fullList)) {
            merged<-merge_phyloseq(merged, fullList[[i]])
        }
        test<-merge_phyloseq( merged,
                        prune_samples(sample_data(ps)$Sample!="Algae" , ps)#
        )
        sample_data(test)$comparison<-(unlist(lapply(str_split(sample_names(test), ":"),"[", 2)))
        sample_data(test)$comparison[is.na(sample_data(test)$comparison)]<- as.character(sample_data(test)$Sample[is.na(sample_data(test)$comparison)])

        
                    ord <- ordinate(test, "CCA", "bray")
                    Figs[["HoloOrd"]] = plot_ordination(test, ord,color="comparison", shape="Sample") + 
                                        geom_point(size=5) + 
                                        ggtitle("...") +
                                        scale_color_npg() +
                                        stat_ellipse(aes(group=comparison))

        ps<-prune_samples(sample_data(test)$comparison=="gained" | sample_data(test)$comparison=="100", test)
            df = as(sample_data(ps), "data.frame")
            ps_distances<-phyloseq::distance(ps, method="bray") 
            output<- adonis(ps_distances ~  Site + comparison , df, strata=df$Site)
            #significant effect of pH and Site

    # 5e. Ordination showing reduced holobiont distinctness with OA
                    ps16p_holos<-prune_samples(sample_data(ps16p)$Sample=="Sessile" | sample_data(ps16p)$Sample=="100", ps16p)
                    ps16p_holos<-prune_taxa(taxa_sums(ps16p_holos)>0,ps16p_holos)
                    ord <- ordinate(ps16p_holos, "CCA", "bray")
                    Figs[["HoloOrd"]] = plot_ordination(ps16p_holos, ord, type="Sample", color="pH", shape="Sample") + 
                                        geom_point(size=5) + 
                                        ggtitle("Increasing similarity of Sessile to 100um fraction with OA (CCA)") +
                                        scale_color_lancet()
    # BUNK - 5f. METABOLITES distinct proprtion



                ChemOverlapByArms<-list(
                    "Proportion of environmental CHEMICALs in holobionts"=data.frame(
                        "A_Control"=c( sum((taxa_names(psb_list[[7]]) %in% taxa_names(psb_list[[8]])))/length(taxa_names(psb_list[[7]])),
                                    sum((taxa_names(psb_list[[9]]) %in% taxa_names(psb_list[[10]])))/length(taxa_names(psb_list[[9]])),
                                    sum((taxa_names(psb_list[[17]]) %in% taxa_names(psb_list[[18]])))/length(taxa_names(psb_list[[17]])),
                                    sum((taxa_names(psb_list[[19]]) %in% taxa_names(psb_list[[20]])))/length(taxa_names(psb_list[[19]])),
                                    sum((taxa_names(psb_list[[29]]) %in% taxa_names(psb_list[[30]])))/length(taxa_names(psb_list[[29]])),
                                    sum((taxa_names(psb_list[[31]]) %in% taxa_names(psb_list[[32]])))/length(taxa_names(psb_list[[31]]))
                                    ),
                        "B_Medium"=c(sum((taxa_names(psb_list[[5]]) %in% taxa_names(psb_list[[6]])))/length(taxa_names(psb_list[[5]])), 
                                sum((taxa_names(psb_list[[11]]) %in% taxa_names(psb_list[[12]])))/length(taxa_names(psb_list[[11]])),
                                sum((taxa_names(psb_list[[13]]) %in% taxa_names(psb_list[[14]])))/length(taxa_names(psb_list[[13]])),
                                sum((taxa_names(psb_list[[25]]) %in% taxa_names(psb_list[[26]])))/length(taxa_names(psb_list[[25]])), 
                                sum((taxa_names(psb_list[[33]]) %in% taxa_names(psb_list[[34]])))/length(taxa_names(psb_list[[33]])),
                                sum((taxa_names(psb_list[[35]]) %in% taxa_names(psb_list[[36]])))/length(taxa_names(psb_list[[35]]))
                                ),
                        "C_Low"=c(sum((taxa_names(psb_list[[1]]) %in% taxa_names(psb_list[[2]])))/length(taxa_names(psb_list[[1]])), 
                                sum((taxa_names(psb_list[[3]]) %in% taxa_names(psb_list[[4]])))/length(taxa_names(psb_list[[3]])),
                                sum((taxa_names(psb_list[[15]]) %in% taxa_names(psb_list[[16]])))/length(taxa_names(psb_list[[15]])),
                                sum((taxa_names(psb_list[[21]]) %in% taxa_names(psb_list[[22]])))/length(taxa_names(psb_list[[21]])),
                                sum((taxa_names(psb_list[[23]]) %in% taxa_names(psb_list[[24]])))/length(taxa_names(psb_list[[23]])),
                                sum((taxa_names(psb_list[[27]]) %in% taxa_names(psb_list[[28]])))/length(taxa_names(psb_list[[27]]))
                                )
                    ),
                    "Proportion of holobiont CHEMICALs in environment"=data.frame(
                        "A_Control"=c( sum((taxa_names(psb_list[[8]]) %in% taxa_names(psb_list[[7]])))/length(taxa_names(psb_list[[8]])),
                                    sum((taxa_names(psb_list[[10]]) %in% taxa_names(psb_list[[9]])))/length(taxa_names(psb_list[[10]])),
                                    sum((taxa_names(psb_list[[18]]) %in% taxa_names(psb_list[[17]])))/length(taxa_names(psb_list[[18]])),
                                    sum((taxa_names(psb_list[[20]]) %in% taxa_names(psb_list[[19]])))/length(taxa_names(psb_list[[20]])),
                                    sum((taxa_names(psb_list[[30]]) %in% taxa_names(psb_list[[29]])))/length(taxa_names(psb_list[[30]])),
                                    sum((taxa_names(psb_list[[32]]) %in% taxa_names(psb_list[[31]])))/length(taxa_names(psb_list[[32]]))
                                    ),
                        "B_Medium"=c(sum((taxa_names(psb_list[[6]]) %in% taxa_names(psb_list[[5]])))/length(taxa_names(psb_list[[6]])), 
                                sum((taxa_names(psb_list[[12]]) %in% taxa_names(psb_list[[11]])))/length(taxa_names(psb_list[[12]])),
                                sum((taxa_names(psb_list[[14]]) %in% taxa_names(psb_list[[13]])))/length(taxa_names(psb_list[[14]])),
                                sum((taxa_names(psb_list[[26]]) %in% taxa_names(psb_list[[25]])))/length(taxa_names(psb_list[[26]])), 
                                sum((taxa_names(psb_list[[34]]) %in% taxa_names(psb_list[[33]])))/length(taxa_names(psb_list[[34]])),
                                sum((taxa_names(psb_list[[36]]) %in% taxa_names(psb_list[[35]])))/length(taxa_names(psb_list[[36]]))
                                ),
                        "C_Low"=c(sum((taxa_names(psb_list[[2]]) %in% taxa_names(psb_list[[1]])))/length(taxa_names(psb_list[[2]])), 
                                sum((taxa_names(psb_list[[4]]) %in% taxa_names(psb_list[[3]])))/length(taxa_names(psb_list[[4]])),
                                sum((taxa_names(psb_list[[16]]) %in% taxa_names(psb_list[[15]])))/length(taxa_names(psb_list[[16]])),
                                sum((taxa_names(psb_list[[22]]) %in% taxa_names(psb_list[[21]])))/length(taxa_names(psb_list[[22]])),
                                sum((taxa_names(psb_list[[24]]) %in% taxa_names(psb_list[[23]])))/length(taxa_names(psb_list[[24]])),
                                sum((taxa_names(psb_list[[28]]) %in% taxa_names(psb_list[[27]])))/length(taxa_names(psb_list[[28]]))
                                )
                    )
                )


        #pivot longer and append island id
            chemlongdf<-cbind(pivot_longer(ChemOverlapByArms[[1]], everything(), names_to="pH"), pivot_longer(ChemOverlapByArms[[2]], everything(), names_to="pH"))
            chemlongdf<-chemlongdf[c(1,2,4)]
            names(chemlongdf)<- c("pH" ,"Proportion of metabolites from the environment also found in Holobiont",  "Proportion of metabolites from Holobionts also found in the environment")
            chemlongdf$Island<-c(rep("Dobu", 9), rep("Illi", 9))
            chemlongdf$Island<-as.factor(chemlongdf$Island)
            chemlongdf$pH<-as.factor(chemlongdf$pH)

        chemlongerdf<-pivot_longer(chemlongdf, 2:3, names_to="overlapType")
                    
                Figs[["MetaboliteOverlap"]] <- ggplot(chemlongerdf, aes(x=pH, y=value))+
                                            ylim(0.2, 0.7) + geom_hline(yintercept=0.5)+
                                            geom_boxplot(aes(fill=overlapType)) +
                                            labs(   title="Proportion of Metabolites Shared between the Environment and Holobionts",
                                                    subtitle="",
                                                    y="Proprtion of ESVs shared", x="pH Regime") +
                                            scale_fill_npg()


        grid.arrange(grobs = list(Figs[["ESVOverlap"]], Figs[["MetaboliteOverlap"]]), ncol = 2) ## display plot




















    # BUNK - 5g. additional ratios - potentially for SI

                OverlapByArms<-list(
                    "holobiont distinct seqs as proportion of free-living"=data.frame(
                        "A_Control"=c( sum((!taxa_names(ps_list[[6]]) %in% taxa_names(ps_list[[5]])))/length(taxa_names(ps_list[[5]])),
                                    sum((!taxa_names(ps_list[[8]]) %in% taxa_names(ps_list[[7]])))/length(taxa_names(ps_list[[1]])),
                                    sum((!taxa_names(ps_list[[12]]) %in% taxa_names(ps_list[[11]])))/length(taxa_names(ps_list[[11]])),
                                    sum((!taxa_names(ps_list[[14]]) %in% taxa_names(ps_list[[13]])))/length(taxa_names(ps_list[[13]])),
                                    sum((!taxa_names(ps_list[[22]]) %in% taxa_names(ps_list[[21]])))/length(taxa_names(ps_list[[21]])),
                                    sum((!taxa_names(ps_list[[24]]) %in% taxa_names(ps_list[[23]])))/length(taxa_names(ps_list[[23]]))
                                    ),
                        "B_Medium"=c(sum((!taxa_names(ps_list[[10]]) %in% taxa_names(ps_list[[9]])))/length(taxa_names(ps_list[[9]])), 
                                sum((!taxa_names(ps_list[[26]]) %in% taxa_names(ps_list[[25]])))/length(taxa_names(ps_list[[25]])),
                                sum((!taxa_names(ps_list[[28]]) %in% taxa_names(ps_list[[27]])))/length(taxa_names(ps_list[[27]])),
                                NA,
                                NA,
                                NA
                                ),
                        "C_Low"=c(sum((!taxa_names(ps_list[[2]]) %in% taxa_names(ps_list[[1]])))/length(taxa_names(ps_list[[1]])), 
                                sum((!taxa_names(ps_list[[4]]) %in% taxa_names(ps_list[[3]])))/length(taxa_names(ps_list[[3]])),
                                sum((!taxa_names(ps_list[[16]]) %in% taxa_names(ps_list[[15]])))/length(taxa_names(ps_list[[15]])),
                                sum((!taxa_names(ps_list[[18]]) %in% taxa_names(ps_list[[17]])))/length(taxa_names(ps_list[[17]])),
                                sum((!taxa_names(ps_list[[20]]) %in% taxa_names(ps_list[[19]])))/length(taxa_names(ps_list[[19]])),
                                NA
                                )
                    ),
                    "holobiont distinct seqs as proportion of free-living distinct"=data.frame(
                        "A_Control"=c( sum((!taxa_names(ps_list[[6]]) %in% taxa_names(ps_list[[5]])))/sum((!taxa_names(ps_list[[5]]) %in% taxa_names(ps_list[[6]]))),
                                    sum((!taxa_names(ps_list[[8]]) %in% taxa_names(ps_list[[7]])))/sum((!taxa_names(ps_list[[7]]) %in% taxa_names(ps_list[[8]]))),
                                    sum((!taxa_names(ps_list[[12]]) %in% taxa_names(ps_list[[11]])))/sum((!taxa_names(ps_list[[11]]) %in% taxa_names(ps_list[[12]]))),
                                    sum((!taxa_names(ps_list[[14]]) %in% taxa_names(ps_list[[13]])))/sum((!taxa_names(ps_list[[13]]) %in% taxa_names(ps_list[[14]]))),
                                    sum((!taxa_names(ps_list[[22]]) %in% taxa_names(ps_list[[21]])))/sum((!taxa_names(ps_list[[21]]) %in% taxa_names(ps_list[[22]]))),
                                    sum((!taxa_names(ps_list[[24]]) %in% taxa_names(ps_list[[23]])))/sum((!taxa_names(ps_list[[23]]) %in% taxa_names(ps_list[[24]])))
                                    ),
                        "B_Medium"=c(sum((!taxa_names(ps_list[[10]]) %in% taxa_names(ps_list[[9]])))/sum((!taxa_names(ps_list[[9]]) %in% taxa_names(ps_list[[10]]))), 
                                sum((!taxa_names(ps_list[[26]]) %in% taxa_names(ps_list[[25]])))/sum((!taxa_names(ps_list[[25]]) %in% taxa_names(ps_list[[26]]))),
                                sum((!taxa_names(ps_list[[28]]) %in% taxa_names(ps_list[[27]])))/sum((!taxa_names(ps_list[[27]]) %in% taxa_names(ps_list[[28]]))),
                                NA,
                                NA,
                                NA
                                ),
                        "C_Low"=c(sum((!taxa_names(ps_list[[2]]) %in% taxa_names(ps_list[[1]])))/sum((!taxa_names(ps_list[[1]]) %in% taxa_names(ps_list[[2]]))), 
                                sum((!taxa_names(ps_list[[4]]) %in% taxa_names(ps_list[[3]])))/sum((!taxa_names(ps_list[[3]]) %in% taxa_names(ps_list[[4]]))),
                                sum((!taxa_names(ps_list[[16]]) %in% taxa_names(ps_list[[15]])))/sum((!taxa_names(ps_list[[15]]) %in% taxa_names(ps_list[[16]]))),
                                sum((!taxa_names(ps_list[[18]]) %in% taxa_names(ps_list[[17]])))/sum((!taxa_names(ps_list[[17]]) %in% taxa_names(ps_list[[18]]))),
                                sum((!taxa_names(ps_list[[20]]) %in% taxa_names(ps_list[[19]])))/sum((!taxa_names(ps_list[[19]]) %in% taxa_names(ps_list[[20]]))),
                                NA
                                )
                    ),
                    "ratio of env to holobiont distinct richness"=data.frame(
                        "A_Control"=c(length(phyloseq::taxa_names(ps_list[[5]])) /sum(!taxa_names(ps_list[[6]]) %in% taxa_names(ps_list[[5]])) , 
                                    length(taxa_names(ps_list[[7]])) /sum(!taxa_names(ps_list[[8]]) %in% taxa_names(ps_list[[7]])),
                                    length(taxa_names(ps_list[[11]])) /sum(!taxa_names(ps_list[[12]]) %in% taxa_names(ps_list[[11]])),
                                    length(taxa_names(ps_list[[13]])) /sum(!taxa_names(ps_list[[14]]) %in% taxa_names(ps_list[[13]])),
                                    length(taxa_names(ps_list[[21]])) /sum(!taxa_names(ps_list[[22]]) %in% taxa_names(ps_list[[21]])),
                                    length(taxa_names(ps_list[[23]])) /sum(!taxa_names(ps_list[[24]]) %in% taxa_names(ps_list[[23]]))
                                    ),
                        "B_Medium"=c(length(taxa_names(ps_list[[9]])) /sum(!taxa_names(ps_list[[10]]) %in% taxa_names(ps_list[[9]])), 
                                length(taxa_names(ps_list[[25]])) /sum(!taxa_names(ps_list[[26]]) %in% taxa_names(ps_list[[25]])),
                                length(taxa_names(ps_list[[27]])) /sum(!taxa_names(ps_list[[28]]) %in% taxa_names(ps_list[[27]])),
                                NA,
                                NA,
                                NA
                                ),
                        "C_Low"=c(length(taxa_names(ps_list[[1]])) /sum(!taxa_names(ps_list[[2]]) %in% taxa_names(ps_list[[1]])), 
                                length(taxa_names(ps_list[[3]])) /sum(!taxa_names(ps_list[[4]]) %in% taxa_names(ps_list[[3]])),
                                length(taxa_names(ps_list[[15]])) /sum(!taxa_names(ps_list[[16]]) %in% taxa_names(ps_list[[15]])),
                                length(taxa_names(ps_list[[17]])) /sum(!taxa_names(ps_list[[18]]) %in% taxa_names(ps_list[[17]])),
                                length(taxa_names(ps_list[[19]])) /sum(!taxa_names(ps_list[[20]]) %in% taxa_names(ps_list[[19]])),
                                NA
                                )
                    )
                )


        #pivot longer and append island id
            longdf1<-pivot_longer(OverlapByArms[[1]], everything(), names_to="pH")
            longdf2<-pivot_longer(OverlapByArms[[2]], everything(), names_to="pH")
            longdf3<-pivot_longer(OverlapByArms[[3]], everything(), names_to="pH")
            longdf<-cbind(cbind(longdf1, longdf2), longdf3)
            longdf<-longdf[c(1,2,4,6)]
            names(longdf)<- c("pH" ,paste0("Ratio: ",names(OverlapByArms)))
            longdf$Island<-c(rep("Dobu", 4), "Illi", rep("Dobu", 2), rep("Illi", 3), NA, rep("Illi", 2), NA, rep("Illi", 2), NA, NA)
            longdf$Island<-as.factor(longdf$Island)
            longdf$pH<-as.factor(longdf$pH)
            longerdf<-pivot_longer(longdf, cols=starts_with("Ratio:"), names_to="Ratio")
                    
         Figs[["holobDist"]]<-   ggplot(longerdf, aes(x=pH, y=value, fill=Ratio)) + 
                                                geom_hline(yintercept=1) + #ylim(0, 0.5)+
                                                geom_boxplot() +
                                                scale_color_npg() +
                                                scale_fill_npg() +
                                                #scale_y_continuous(trans='log10') +
                                                labs(   title = "Ratio of environmental to Holobiont unqiue microibal richness",,
                                                        subtitle="",
                                                        y="Holobiont Unique richness as proportion of free-living", x="pH Regime")

     #bayes glm
        SelectedRatio<-unique(longerdf$Ratio)
        model<-list()
        for (i in 1:3) {
            dat<-longerdf[longerdf$Ratio==SelectedRatio[i],]
            model[[i]] <- stan_glm(data = dat, formula=dat[["value"]] ~ pH+Island, family=gaussian())

        }
                model_general[[i]] <- stan_glm(data = longerdf[longerdf$Ratio==SelectedRatio[1],], formula=longerdf[longerdf$Ratio==SelectedRatio[1],][["value"]] ~ pH+Island, family=gaussian())
                model_general2 <- stan_glm(data = longerdf[longerdf$Ratio==SelectedRatio[2],], formula=longerdf[longerdf$Ratio==SelectedRatio[2],][["value"]] ~ pH+Island, family=gaussian())
                model_general3 <- stan_glm(data = longerdf[longerdf$Ratio==SelectedRatio[3],], formula=longerdf[longerdf$Ratio==SelectedRatio[3],][["value"]] ~ pH+Island, family=gaussian())

                Figs$s_CIplot_general<-plot(model_general1, plotfun = "areas", prob = 0.9, pars="beta")
                Figs$s_posteriors_general <- describe_posterior(model[[1]])
                 describe_posterior(model[[2]])
                  describe_posterior(model[[3]])

   
   
   
   
    # BUNK - 5h. explorations with algal data - but as subset of ARMS have all three data types this reduces samplesize - nto worth it

        #by arms unit
                ps_list<-list()
                    ps_s100<-prune_samples( (sample_data(psP)$Sample=="Algae" | sample_data(psP)$Sample=="100" )  , psP)
                    SharedARMS<-levels(sample_data(ps_s100)$ARMS) [ table(sample_data(ps_s100)$ARMS)>1 ]
                    samps<-levels(sample_data(ps_s100)$Sample)
                    for (ARMS in SharedARMS) {
                        for (samp in samps) {
                                ps_temp<-prune_samples(sample_data(ps_s100)$ARMS==ARMS, ps_s100)
                                ps_temp <-prune_taxa(taxa_sums(ps_temp)>0, ps_temp)

                                ps_temp<-prune_samples(sample_data( ps_temp)$Sample==samp, ps_temp)
                                nsamples<-dim(sample_data(ps_temp))[1]
                                if (nsamples>1) {
                                        ps_temp<-prune_samples(c(TRUE, rep(FALSE, nsamples-1)), ps_temp)
                                }
                                pH<- sample_data(ps_temp)$pH
                                ps_list[[ paste0(pH, "_", ARMS, "_", samp) ]] <-prune_taxa(taxa_sums(ps_temp)>0, ps_temp)

                        }

                    }
            names(ps_list)
            algaeratios<-list(
                "ratio of env to holobiont distinct richness"=data.frame(
                        "A_Control"=c(length(phyloseq::taxa_names(ps_list[[3]])) /sum(!taxa_names(ps_list[[4]]) %in% taxa_names(ps_list[[3]])) , 
                                    length(taxa_names(ps_list[[11]])) /sum(!taxa_names(ps_list[[12]]) %in% taxa_names(ps_list[[11]])),
                                    length(taxa_names(ps_list[[13]])) /sum(!taxa_names(ps_list[[14]]) %in% taxa_names(ps_list[[13]])),
                                    length(taxa_names(ps_list[[21]])) /sum(!taxa_names(ps_list[[22]]) %in% taxa_names(ps_list[[21]])),
                                    length(taxa_names(ps_list[[23]])) /sum(!taxa_names(ps_list[[24]]) %in% taxa_names(ps_list[[23]])),
                                    NA
                                    ),
                        "B_Medium"=c(length(taxa_names(ps_list[[5]])) /sum(!taxa_names(ps_list[[6]]) %in% taxa_names(ps_list[[5]])), 
                                length(taxa_names(ps_list[[7]])) /sum(!taxa_names(ps_list[[8]]) %in% taxa_names(ps_list[[7]])),
                                length(taxa_names(ps_list[[25]])) /sum(!taxa_names(ps_list[[26]]) %in% taxa_names(ps_list[[25]])),
                                length(taxa_names(ps_list[[27]])) /sum(!taxa_names(ps_list[[28]]) %in% taxa_names(ps_list[[27]])),
                                NA,
                                NA
                                ),
                        "C_Low"=c(length(taxa_names(ps_list[[1]])) /sum(!taxa_names(ps_list[[2]]) %in% taxa_names(ps_list[[1]])), 
                                length(taxa_names(ps_list[[9]])) /sum(!taxa_names(ps_list[[10]]) %in% taxa_names(ps_list[[9]])),
                                length(taxa_names(ps_list[[15]])) /sum(!taxa_names(ps_list[[16]]) %in% taxa_names(ps_list[[15]])),
                                length(taxa_names(ps_list[[17]])) /sum(!taxa_names(ps_list[[18]]) %in% taxa_names(ps_list[[17]])),
                                length(taxa_names(ps_list[[19]])) /sum(!taxa_names(ps_list[[20]]) %in% taxa_names(ps_list[[19]])),
                                NA
                                )
                    ),
                "ratio of env to holobiont richness"=data.frame(
                        "A_Control"=c(length(phyloseq::taxa_names(ps_list[[3]])) /length(taxa_names(ps_list[[4]])) , 
                                    length(taxa_names(ps_list[[11]])) /length(taxa_names(ps_list[[12]]) ),
                                    length(taxa_names(ps_list[[13]])) /length(taxa_names(ps_list[[14]]) ),
                                    length(taxa_names(ps_list[[21]])) /length(taxa_names(ps_list[[22]]) ),
                                    length(taxa_names(ps_list[[23]])) /length(taxa_names(ps_list[[24]]) ),
                                    NA
                                    ),
                        "B_Medium"=c(length(taxa_names(ps_list[[5]])) /length(taxa_names(ps_list[[6]]) ), 
                                length(taxa_names(ps_list[[7]])) /length(taxa_names(ps_list[[8]]) ),
                                length(taxa_names(ps_list[[25]])) /length(taxa_names(ps_list[[26]]) ),
                                length(taxa_names(ps_list[[27]])) /length(taxa_names(ps_list[[28]]) ),
                                NA,
                                NA
                                ),
                        "C_Low"=c(length(taxa_names(ps_list[[1]])) /length(taxa_names(ps_list[[2]]) ), 
                                length(taxa_names(ps_list[[9]])) /length(taxa_names(ps_list[[10]]) ),
                                length(taxa_names(ps_list[[15]])) /length(taxa_names(ps_list[[16]]) ),
                                length(taxa_names(ps_list[[17]])) /length(taxa_names(ps_list[[18]]) ),
                                length(taxa_names(ps_list[[19]])) /length(taxa_names(ps_list[[20]]) ),
                                NA
                                )
                    )
            )

                 #pivot longer and append island id
            longdf1<-pivot_longer(algaeratios[[1]], everything(), names_to="pH")
            longdf2<-pivot_longer(OverlapByArms[[2]], everything(), names_to="pH")
            longdf3<-pivot_longer(OverlapByArms[[3]], everything(), names_to="pH")
            longdf<-cbind(cbind(longdf1, longdf2), longdf3)
            longdf<-longdf[c(1,2,4,6)]
            names(longdf)<- c("pH" ,paste0("Ratio: ",names(OverlapByArms)))
            longdf$Island<-c(rep("Dobu", 4), "Illi", rep("Dobu", 2), rep("Illi", 3), NA, rep("Illi", 2), NA, rep("Illi", 2), NA, NA)
            longdf$Island<-as.factor(longdf$Island)
            longdf$pH<-as.factor(longdf$pH)
            longerdf<-pivot_longer(longdf, cols=starts_with("Ratio:"), names_to="Ratio")
                    
         Figs[["holobDist"]]<-   ggplot(longerdf, aes(x=pH, y=value, fill=Ratio)) + 
                                                geom_hline(yintercept=1) + #ylim(0, 0.5)+
                                                geom_boxplot() +
                                                scale_color_npg() +
                                                scale_fill_npg() +
                                                #scale_y_continuous(trans='log10') +
                                                labs(   title = "Ratio of environmental to Holobiont unqiue microibal richness",,
                                                        subtitle="",
                                                        y="Holobiont Unique richness as proportion of free-living", x="pH Regime")







        # algae sessile comparison
                ps_list<-list()
                    ps_s100<-prune_samples( (sample_data(psP)$Sample=="Algae" | sample_data(psP)$Sample=="Sessile")  , psP)
                    SharedARMS<-levels(sample_data(ps_s100)$ARMS) [ table(sample_data(ps_s100)$ARMS[sample_data(ps_s100)$Sample=="Algae"])>=1 & table(sample_data(ps_s100)$ARMS[sample_data(ps_s100)$Sample=="Sessile"])>=1]
                    samps<-levels(sample_data(ps_s100)$Sample)
                    for (ARMS in SharedARMS) {
                        for (samp in samps) {
                                ps_temp<-prune_samples(sample_data(ps_s100)$ARMS==ARMS, ps_s100)
                                ps_temp <-prune_taxa(taxa_sums(ps_temp)>0, ps_temp)

                                ps_temp<-prune_samples(sample_data( ps_temp)$Sample==samp, ps_temp)
                                nsamples<-dim(sample_data(ps_temp))[1]
                                if (nsamples>1) {
                                        ps_temp<-prune_samples(c(TRUE, rep(FALSE, nsamples-1)), ps_temp)
                                }
                                pH<- sample_data(ps_temp)$pH
                                ps_list[[ paste0(pH, "_", ARMS, "_", samp) ]] <-prune_taxa(taxa_sums(ps_temp)>0, ps_temp)

                        }

                    }
            names(ps_list)
            algaeratios<-list(
                "Sessile non-algal richness"=data.frame(
                        "A_Control"=c(sum(!taxa_names(ps_list[[4]]) %in% taxa_names(ps_list[[3]])) , 
                                    sum(!taxa_names(ps_list[[8]]) %in% taxa_names(ps_list[[7]])),
                                    sum(!taxa_names(ps_list[[10]]) %in% taxa_names(ps_list[[9]])),
                                    sum(!taxa_names(ps_list[[20]]) %in% taxa_names(ps_list[[19]])),
                                    sum(!taxa_names(ps_list[[22]]) %in% taxa_names(ps_list[[21]]))
                                    ),
                        "B_Medium"=c(sum(!taxa_names(ps_list[[6]]) %in% taxa_names(ps_list[[5]])), 
                                sum(!taxa_names(ps_list[[16]]) %in% taxa_names(ps_list[[15]])),
                                sum(!taxa_names(ps_list[[24]]) %in% taxa_names(ps_list[[23]])),
                                sum(!taxa_names(ps_list[[26]]) %in% taxa_names(ps_list[[25]])),
                                NA
                                ),
                        "C_Low"=c(sum(!taxa_names(ps_list[[2]]) %in% taxa_names(ps_list[[1]])), 
                                sum(!taxa_names(ps_list[[12]]) %in% taxa_names(ps_list[[11]])),
                                sum(!taxa_names(ps_list[[14]]) %in% taxa_names(ps_list[[13]])),
                                sum(!taxa_names(ps_list[[18]]) %in% taxa_names(ps_list[[17]])),
                                NA
                                )
                    ))
               
               
    
               ps_list<-list()
                    ps_s100<-prune_samples( (sample_data(psP)$Sample=="Algae" | sample_data(psP)$Sample=="Sessile" | sample_data(psP)$Sample=="100"  )  , psP)
                    SharedARMS<-levels(sample_data(ps_s100)$ARMS) [ table(sample_data(ps_s100)$ARMS[sample_data(ps_s100)$Sample=="Algae"])>=1 & table(sample_data(ps_s100)$ARMS[sample_data(ps_s100)$Sample=="Sessile"])>=1 & table(sample_data(ps_s100)$ARMS[sample_data(ps_s100)$Sample=="100"])>=1 ]
                    samps<-levels(sample_data(ps_s100)$Sample)
                    for (ARMS in SharedARMS) {
                        for (samp in samps) {
                                ps_temp<-prune_samples(sample_data(ps_s100)$ARMS==ARMS, ps_s100)
                                ps_temp <-prune_taxa(taxa_sums(ps_temp)>0, ps_temp)

                                ps_temp<-prune_samples(sample_data( ps_temp)$Sample==samp, ps_temp)
                                nsamples<-dim(sample_data(ps_temp))[1]
                                if (nsamples>1) {
                                        ps_temp<-prune_samples(c(TRUE, rep(FALSE, nsamples-1)), ps_temp)
                                }
                                pH<- sample_data(ps_temp)$pH
                                ps_list[[ paste0(pH, "_", ARMS, "_", samp) ]] <-prune_taxa(taxa_sums(ps_temp)>0, ps_temp)

                        }

                    }
            names(ps_list)


               
            sessileonlyratio<-list(
                "sessile non-algal non-freeliving as proportion of free-living"=data.frame(
                        "A_Control"=c( 
                                        sum(!(taxa_names(ps_list[[6]]) %in% taxa_names(ps_list[[5]]))  & !(taxa_names(ps_list[[6]]) %in%  taxa_names(ps_list[[4]]))) / length(taxa_names(ps_list[[4]])),
                                        sum(!(taxa_names(ps_list[[12]]) %in% taxa_names(ps_list[[11]]))  & !(taxa_names(ps_list[[12]]) %in%  taxa_names(ps_list[[10]]))) / length(taxa_names(ps_list[[10]])),
                                        sum(!(taxa_names(ps_list[[15]]) %in% taxa_names(ps_list[[14]]))  & !(taxa_names(ps_list[[15]]) %in%  taxa_names(ps_list[[13]]))) / length(taxa_names(ps_list[[13]])),
                                        sum(!(taxa_names(ps_list[[27]]) %in% taxa_names(ps_list[[26]]))  & !(taxa_names(ps_list[[27]]) %in%  taxa_names(ps_list[[25]]))) / length(taxa_names(ps_list[[25]])),
                                        sum(!(taxa_names(ps_list[[30]]) %in% taxa_names(ps_list[[29]]))  & !(taxa_names(ps_list[[30]]) %in%  taxa_names(ps_list[[28]]))) / length(taxa_names(ps_list[[28]]))
                                    ),
                        "B_Medium"=c(
                                        sum(!(taxa_names(ps_list[[9]]) %in% taxa_names(ps_list[[8]]))  & !(taxa_names(ps_list[[9]]) %in%  taxa_names(ps_list[[7]]))) / length(taxa_names(ps_list[[7]])),
                                        sum(!(taxa_names(ps_list[[33]]) %in% taxa_names(ps_list[[32]]))  & !(taxa_names(ps_list[[33]]) %in%  taxa_names(ps_list[[31]]))) / length(taxa_names(ps_list[[31]])),
                                        sum(!(taxa_names(ps_list[[36]]) %in% taxa_names(ps_list[[35]]))  & !(taxa_names(ps_list[[36]]) %in%  taxa_names(ps_list[[34]]))) / length(taxa_names(ps_list[[34]])),
                                        NA,
                                        NA
                                    ),
                        "C_Low"=c(
                                        sum(!(taxa_names(ps_list[[3]]) %in% taxa_names(ps_list[[2]]))  & !(taxa_names(ps_list[[3]]) %in%  taxa_names(ps_list[[1]]))) / length(taxa_names(ps_list[[1]])),
                                        sum(!(taxa_names(ps_list[[18]]) %in% taxa_names(ps_list[[17]]))  & !(taxa_names(ps_list[[18]]) %in%  taxa_names(ps_list[[16]]))) / length(taxa_names(ps_list[[16]])),
                                        sum(!(taxa_names(ps_list[[21]]) %in% taxa_names(ps_list[[20]]))  & !(taxa_names(ps_list[[21]]) %in%  taxa_names(ps_list[[19]]))) / length(taxa_names(ps_list[[19]])),
                                        sum(!(taxa_names(ps_list[[24]]) %in% taxa_names(ps_list[[23]]))  & !(taxa_names(ps_list[[24]]) %in%  taxa_names(ps_list[[22]]))) / length(taxa_names(ps_list[[22]])),
                                        NA
                                    )
                    ),
                "sessile non-algal non-freeliving richness"=data.frame(
                        "A_Con"=c( 
                                        sum(!(taxa_names(ps_list[[6]]) %in% taxa_names(ps_list[[5]]))  & !(taxa_names(ps_list[[6]]) %in%  taxa_names(ps_list[[4]]))),
                                        sum(!(taxa_names(ps_list[[12]]) %in% taxa_names(ps_list[[11]]))  & !(taxa_names(ps_list[[12]]) %in%  taxa_names(ps_list[[10]]))) ,
                                        sum(!(taxa_names(ps_list[[15]]) %in% taxa_names(ps_list[[14]]))  & !(taxa_names(ps_list[[15]]) %in%  taxa_names(ps_list[[13]]))) ,
                                        sum(!(taxa_names(ps_list[[27]]) %in% taxa_names(ps_list[[26]]))  & !(taxa_names(ps_list[[27]]) %in%  taxa_names(ps_list[[25]]))),
                                        sum(!(taxa_names(ps_list[[30]]) %in% taxa_names(ps_list[[29]]))  & !(taxa_names(ps_list[[30]]) %in%  taxa_names(ps_list[[28]]))) 
                                    ),
                        "B_Med"=c(
                                        sum(!(taxa_names(ps_list[[9]]) %in% taxa_names(ps_list[[8]]))  & !(taxa_names(ps_list[[9]]) %in%  taxa_names(ps_list[[7]]))) ,
                                        sum(!(taxa_names(ps_list[[33]]) %in% taxa_names(ps_list[[32]]))  & !(taxa_names(ps_list[[33]]) %in%  taxa_names(ps_list[[31]]))) ,
                                        sum(!(taxa_names(ps_list[[36]]) %in% taxa_names(ps_list[[35]]))  & !(taxa_names(ps_list[[36]]) %in%  taxa_names(ps_list[[34]]))) ,
                                        NA,
                                        NA
                                    ),
                        "C_Low"=c(
                                        sum(!(taxa_names(ps_list[[3]]) %in% taxa_names(ps_list[[2]]))  & !(taxa_names(ps_list[[3]]) %in%  taxa_names(ps_list[[1]]))) ,
                                        sum(!(taxa_names(ps_list[[18]]) %in% taxa_names(ps_list[[17]]))  & !(taxa_names(ps_list[[18]]) %in%  taxa_names(ps_list[[16]]))) ,
                                        sum(!(taxa_names(ps_list[[21]]) %in% taxa_names(ps_list[[20]]))  & !(taxa_names(ps_list[[21]]) %in%  taxa_names(ps_list[[19]]))),
                                        sum(!(taxa_names(ps_list[[24]]) %in% taxa_names(ps_list[[23]]))  & !(taxa_names(ps_list[[24]]) %in%  taxa_names(ps_list[[22]]))) ,
                                        NA
                                    )
                    ),
                "Fraction of Holbiont Microbes unique to the Holbiont"=data.frame(
                        "A_Con"=c( 
                                    sum(!(taxa_names(ps_list[[5]]) %in% taxa_names(ps_list[[4]])))/length(taxa_names(ps_list[[5]])),
                                    sum(!(taxa_names(ps_list[[6]]) %in% taxa_names(ps_list[[4]])))/length(taxa_names(ps_list[[6]])),
                                    sum(!(taxa_names(ps_list[[11]]) %in% taxa_names(ps_list[[10]])))/length(taxa_names(ps_list[[11]])),
                                    sum(!(taxa_names(ps_list[[12]]) %in% taxa_names(ps_list[[10]])))/length(taxa_names(ps_list[[12]])),
                                    sum(!(taxa_names(ps_list[[14]]) %in% taxa_names(ps_list[[13]])))/length(taxa_names(ps_list[[14]])),
                                    sum(!(taxa_names(ps_list[[15]]) %in% taxa_names(ps_list[[13]])))/length(taxa_names(ps_list[[15]])),
                                    sum(!(taxa_names(ps_list[[26]]) %in% taxa_names(ps_list[[25]])))/length(taxa_names(ps_list[[26]])),
                                    sum(!(taxa_names(ps_list[[27]]) %in% taxa_names(ps_list[[25]])))/length(taxa_names(ps_list[[27]])),
                                    sum(!(taxa_names(ps_list[[29]]) %in% taxa_names(ps_list[[28]])))/length(taxa_names(ps_list[[29]])),
                                    sum(!(taxa_names(ps_list[[30]]) %in% taxa_names(ps_list[[28]])))/length(taxa_names(ps_list[[30]]))
                                    ),
                        "B_Med"=c(
                                    sum(!(taxa_names(ps_list[[8]]) %in% taxa_names(ps_list[[7]])))/length(taxa_names(ps_list[[8]])),
                                    sum(!(taxa_names(ps_list[[9]]) %in% taxa_names(ps_list[[7]])))/length(taxa_names(ps_list[[9]])),
                                    sum(!(taxa_names(ps_list[[32]]) %in% taxa_names(ps_list[[31]])))/length(taxa_names(ps_list[[32]])),
                                    sum(!(taxa_names(ps_list[[33]]) %in% taxa_names(ps_list[[31]])))/length(taxa_names(ps_list[[33]])),
                                    sum(!(taxa_names(ps_list[[35]]) %in% taxa_names(ps_list[[34]])))/length(taxa_names(ps_list[[35]])),
                                    sum(!(taxa_names(ps_list[[36]]) %in% taxa_names(ps_list[[34]])))/length(taxa_names(ps_list[[36]])),
                                    NA,
                                    NA,
                                    NA,
                                    NA
                                    ),
                        "C_Low"=c(
                                    sum(!(taxa_names(ps_list[[2]]) %in% taxa_names(ps_list[[1]])))/length(taxa_names(ps_list[[2]])),
                                    sum(!(taxa_names(ps_list[[3]]) %in% taxa_names(ps_list[[1]])))/length(taxa_names(ps_list[[3]])),
                                    sum(!(taxa_names(ps_list[[17]]) %in% taxa_names(ps_list[[16]])))/length(taxa_names(ps_list[[17]])),
                                    sum(!(taxa_names(ps_list[[18]]) %in% taxa_names(ps_list[[16]])))/length(taxa_names(ps_list[[18]])),
                                    sum(!(taxa_names(ps_list[[20]]) %in% taxa_names(ps_list[[19]])))/length(taxa_names(ps_list[[20]])),
                                    sum(!(taxa_names(ps_list[[21]]) %in% taxa_names(ps_list[[19]])))/length(taxa_names(ps_list[[21]])),
                                    sum(!(taxa_names(ps_list[[23]]) %in% taxa_names(ps_list[[22]])))/length(taxa_names(ps_list[[23]])),
                                    sum(!(taxa_names(ps_list[[24]]) %in% taxa_names(ps_list[[22]])))/length(taxa_names(ps_list[[24]])),
                                    NA,
                                    NA
                                    )
                    )
            )
            l<-c()
                for (i in (1: length(ps_list)) ) {
                    l[i]<-(sample_data(ps_list[[i]])[,3])
                }


        longdf1<-pivot_longer(sessileonlyratio[[1]], everything(), names_to="pH")
                    longdf1$pH<-as.factor(longdf1$pH)
                    longdf1$Island<-c(rep("Dobu", 4),  rep("Illi", 6), NA, rep("Illi", 2), NA, NA)
                    longdf1$Island<-as.factor(longdf1$Island)
        longdf2<-pivot_longer(sessileonlyratio[[2]], everything(), names_to="pH")
                    longdf2$pH<-as.factor(longdf2$pH)
                    longdf2$Island<-c(rep("Dobu", 4),  rep("Illi", 6), NA, rep("Illi", 2), NA, NA)
                    longdf2$Island<-as.factor(longdf2$Island)
        longdf3<-pivot_longer(sessileonlyratio[[3]], everything(), names_to="pH")
                    longdf3$pH<-as.factor(longdf2$pH)
                    longdf3$AlgaeOrSessile<-rep(c(rep("Algae", 3), rep("Sessile", 3)),5)
                    longdf3$AlgaeOrSessile<-as.factor(longdf3$AlgaeOrSessile)
                    longdf3$Island<-c(  rep("Dobu", 7),  
                                        rep("Illi", 2),
                                        "Dobu",
                                        rep("Illi",9),
                                        NA,
                                        rep("Illi", 2),
                                        NA,
                                        rep("Illi", 2),
                                        rep(NA, 2),
                                        "Illi",
                                        rep(NA, 2)
                            )
            longdf3$Island<-as.factor(longdf3$Island)

        Figs[["sessileuniqueproportion"]]<-   ggplot(longdf1, aes(x=pH, y=value, fill=pH)) + 
                                                geom_violin() +
                                                geom_boxplot(width=0.1, color="grey", fill="grey", alpha=0.5) +
                                                scale_color_npg() +
                                                scale_fill_npg() +
                                                #scale_y_continuous(trans='log10') +
                                                labs(   title = "Holobiont-associated Community Richness as Proportion of Free-living Community Richness",,
                                                        subtitle="(Holobiont-associated community richness = richness in the sessile fraction excluding sequences also found in the algal or 100-500um fractions.)",
                                                        y="Holobiont-associated Community Richness as Proportion of Free-living Community Richness", x="pH Regime")
            
         Figs[["sessileuniquerichness"]]<-   ggplot(longdf2, aes(x=pH, y=value, fill=pH)) + 
                                                geom_violin() +
                                                scale_color_npg() +
                                                scale_fill_npg() +
                                                #scale_y_continuous(trans='log10') +
                                                labs(   title = "sessile uniqe richness",,
                                                        subtitle="",
                                                        y="# unique sessile sequences", x="pH Regime")


        Figs[["FractionHolbiontUnique"]]<-   ggplot(longdf3, aes(x=pH, y=value, fill=pH)) + 
                                                geom_violin() +
                                                geom_boxplot(width=0.1, color="grey", fill="grey", alpha=0.5) +
                                                scale_color_npg() +
                                                scale_fill_npg() +
                                                facet_grid(~AlgaeOrSessile)+
                                                #scale_y_continuous(trans='log10') +
                                                labs(   title = "Fraction of Holbiont Microbes unique to the Holbiont",
                                                        x="pH Regime (control n=5, medium n=3, low n=4)")
            


    # BUNK - 5i. these plots attempting to show increasing compositional similarity between specific holobionts and the overall environemtnal fraction
        #do not show anything
                    ps16p_RVS<-prune_samples(sample_data(ps16ALLp)$Sample=="RVS" | sample_data(ps16ALLp)$Sample=="100", ps16ALLp)
                    ps16p_RVS<-prune_taxa(taxa_sums(ps16p_RVS)>0,ps16p_RVS)
                    ord <- ordinate(ps16p_RVS, "CCA", "bray")
                    Figs[["RVSOrd"]] = plot_ordination(ps16p_RVS, ord, type="Sample", color="pH", shape="Sample") + 
                                        geom_point(size=5) + 
                                        ggtitle("Increasing similarity of Sessile to 100um fraction with OA (CCA)") +
                                        scale_color_lancet()

                    ps16p_CS<-prune_samples(sample_data(ps16ALLp)$Sample=="CS" | sample_data(ps16ALLp)$Sample=="100", ps16ALLp)
                    ps16p_CS<-prune_taxa(taxa_sums(ps16p_CS)>0,ps16p_CS)
                    ord <- ordinate(ps16p_CS, "CCA", "bray")
                    Figs[["CSOrd"]] = plot_ordination(ps16p_CS, ord, type="Sample", color="pH", shape="Sample") + 
                                        geom_point(size=5) + 
                                        ggtitle("Increasing similarity of Sessile to 100um fraction with OA (CCA)") +
                                        scale_color_lancet()







    # NOT STARTED 5j. holbiont distinctness for specific holobiomes -needs further metadata to map ARMS together
               ps_list<-list()
                    ps_s100<-prune_samples( (sample_data(ps16ALLp)$Sample=="RVS" | sample_data(ps16ALLp)$Sample=="CS" | sample_data(ps16ALLp)$Sample=="Sessile"  )  , ps16ALLp)
                    #havent gone past here as we don't have ARMS codes for the songe samples yet
                    #when we get them copy in from examples above
 
# 6. Taxa barplots
    #create informative taxa barplot for 23s using top 7 rank 4s
    #this is per arms unit
        # create new ps with taxa grouped at rank 4 and labelled by top 9
            ps23p4 = tax_glom(ps23p, "Rank_6")
            topotus = names(sort(taxa_sums(ps23p4), TRUE)[1:9])
            toptaxtab = cbind(tax_table(ps23p4), Class = "Other")
            toptaxtab[topotus, "Class"] <- as(tax_table(ps23p4)[topotus, "Rank_6"], "character")
            tax_table(ps23p4) <- (tax_table(toptaxtab))
        #plot
            Figs[["AlgalClasses"]]<-plot_bar(ps23p4, fill="Class", x="TechnicalReplicate", facet_grid=Site~pH) +
            geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack")    +
            scale_fill_npg() +
            scale_color_npg() +
            ggtitle("Most Abundant Algal Classes per ARMS, Abundance Transformed")

    # chloroplast plot
            Figs[["Chloplasts"]]<-plot_bar(pschloro, fill="Rank_5", x="pH", facet_grid=sample_Sample~Site) +
                                        geom_bar(aes(color=Rank_5, fill=Rank_5), stat="identity", position="stack")    +
                                        scale_fill_npg() +
                                        scale_color_npg() +
                                        ggtitle("Chloroplast Abundance Plot")
                                                


        
    # create informative taxa barplot for 16s using only bacteria and top 
        # create new ps with taxa grouped at rank 4 and labelled by top 9
            ps_temp<-prune_samples(sample_data(ps_s100)$Code==ARMS, ps_s100)
            ps_temp <-prune_taxa(taxa_sums(ps_temp)>0, ps_temp)

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
                                        ggtitle("Most Abundant Bacteria Classes per Site and pH Combination, Abundance Transformed")
    # create informative taxa barplot for 16s free-living
        # create new ps with taxa grouped at rank 4 and labelled by top 9
            ps_temp<-prune_samples(sample_data(ps16p)$Sample==100, ps16p)
            ps_temp <-prune_taxa(taxa_sums(ps_temp)>0, ps_temp)

            ps16p_100only = tax_glom(ps_temp, "Rank_4")
            topotus = names(sort(taxa_sums(ps16p_100only), TRUE)[1:9])
            toptaxtab = cbind(tax_table(ps16p_100only), Class = "Other")
            toptaxtab[topotus, "Class"] <- as(tax_table(ps16p_100only)[topotus, "Rank_4"], "character")
            tax_table(ps16p_100only) <- (tax_table(toptaxtab))
        #merge samples to simplify complex plot
            sample_data(ps16p_100only)$mergeid<-as.factor(paste0(sample_data(ps16p_100only)$pH, "  ", sample_data(ps16p_100only)$Site, "  ",sample_data(ps16p_100only)$Sample))
            ps16p_100only_m<-merge_samples(ps16p_100only, "mergeid")
            sample_data(ps16p_100only_m)$mergeid<-levels(sample_data(ps16p_100only)$mergeid)
            #repair factors
                sample_data(ps16p_100only_m)$pH<-as.factor(unlist(lapply(strsplit(sample_data(ps16p_100only_m)$mergeid, "  "), '[', 1)))
                sample_data(ps16p_100only_m)$Site<-as.factor(unlist(lapply(strsplit(sample_data(ps16p_100only_m)$mergeid, "  "), '[', 2)))
                sample_data(ps16p_100only_m)$Sample<-as.factor(unlist(lapply(strsplit(sample_data(ps16p_100only_m)$mergeid, "  "), '[', 3)))
            ps16p_100only_mp = transform_sample_counts(ps16p_100only_m, function(x) 10000 * x/sum(x))

        #ensure rank 4 names differentiate unclassified and unclassified bacteria
            taxtab<-as.data.frame(tax_table(ps16p_100only_mp))
            taxtab$Class[as.data.frame(tax_table(ps16p_100only_mp))$Class=="Unclassified"] <-
                    as.data.frame(tax_table(ps16p_100only_mp))$Rank_3[as.data.frame(tax_table(ps16p_100only_mp))$Class=="Unclassified"] #rank 3 names where class is unclassified 
                                                                                                                        #(i.e. which differentiate unclassified and unclassified_bacteria) - 
                                                                                                                        #only needed until pipeline adjusted to output this way automatically
            tax_table(ps16p_100only_mp)<-as.matrix(taxtab)

        #plot
            Figs[["100-BacterialClasses"]]<-plot_bar(ps16p_100only_mp, fill="Class", x="pH", facet_grid=sample_Sample~Site) +
                                        geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack")    +
                                        scale_fill_npg() +
                                        scale_color_npg() +
                                        ggtitle("100 Fraction - Most Abundant Bacteria Classes per Site and pH Combination, Abundance Transformed")




    ##### gamma deep dive
        # create informative taxa barplot for 16s free-living
        # create new ps with taxa grouped at rank 4 and labelled by top 9
            ps_temp<-prune_samples(sample_data(ps16p)$Sample==100, ps16p)
            ps_temp <-prune_taxa(as.vector(tax_table(ps_temp)[,4]=="Gammaproteobacteria"), ps_temp)
            ps_temp <-prune_taxa(taxa_sums(ps_temp)>0, ps_temp)

            ps16p_100only = tax_glom(ps_temp, "Rank_6")
            topotus = names(sort(taxa_sums(ps16p_100only), TRUE)[1:9])
            toptaxtab = cbind(tax_table(ps16p_100only), Class = "Other")
            toptaxtab[topotus, "Class"] <- as(tax_table(ps16p_100only)[topotus, "Rank_6"], "character")
            tax_table(ps16p_100only) <- (tax_table(toptaxtab))
        #merge samples to simplify complex plot
            sample_data(ps16p_100only)$mergeid<-as.factor(paste0(sample_data(ps16p_100only)$pH, "  ", sample_data(ps16p_100only)$Site, "  ",sample_data(ps16p_100only)$Sample))
            ps16p_100only_m<-merge_samples(ps16p_100only, "mergeid")
            sample_data(ps16p_100only_m)$mergeid<-levels(sample_data(ps16p_100only)$mergeid)
            #repair factors
                sample_data(ps16p_100only_m)$pH<-as.factor(unlist(lapply(strsplit(sample_data(ps16p_100only_m)$mergeid, "  "), '[', 1)))
                sample_data(ps16p_100only_m)$Site<-as.factor(unlist(lapply(strsplit(sample_data(ps16p_100only_m)$mergeid, "  "), '[', 2)))
                sample_data(ps16p_100only_m)$Sample<-as.factor(unlist(lapply(strsplit(sample_data(ps16p_100only_m)$mergeid, "  "), '[', 3)))
            ps16p_100only_mp = transform_sample_counts(ps16p_100only_m, function(x) 10000 * x/sum(x))

        #ensure rank 4 names differentiate unclassified and unclassified bacteria
            taxtab<-as.data.frame(tax_table(ps16p_100only_mp))
            taxtab$Class[as.data.frame(tax_table(ps16p_100only_mp))$Class=="Unclassified"] <-
                    as.data.frame(tax_table(ps16p_100only_mp))$Rank_3[as.data.frame(tax_table(ps16p_100only_mp))$Class=="Unclassified"] #rank 3 names where class is unclassified 
                                                                                                                        #(i.e. which differentiate unclassified and unclassified_bacteria) - 
                                                                                                                        #only needed until pipeline adjusted to output this way automatically
            tax_table(ps16p_100only_mp)<-as.matrix(taxtab)

        #plot
            Figs[["GammaproteobacteriaFamilies"]]<-plot_bar(ps16p_100only_mp, fill="Class", x="pH", facet_grid=sample_Sample~Site) +
                                        geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack")    +
                                        scale_fill_npg() +
                                        scale_color_npg() +
                                        ggtitle("Most Abundant Gammaproteobacteria Families per Site and pH Combination, Abundance Transformed")




    ##### bacteroidia deep dive
     # create informative taxa barplot for 16s free-living
        # create new ps with taxa grouped at rank 4 and labelled by top 9
            ps_temp<-prune_samples(sample_data(ps16p)$Sample==100, ps16p)
            ps_temp <-prune_taxa(as.vector(tax_table(ps_temp)[,4]=="Bacteroidia"), ps_temp)
            ps_temp <-prune_taxa(taxa_sums(ps_temp)>0, ps_temp)

            ps16p_100only = tax_glom(ps_temp, "Rank_6")
            topotus = names(sort(taxa_sums(ps16p_100only), TRUE)[1:9])
            toptaxtab = cbind(tax_table(ps16p_100only), Class = "Other")
            toptaxtab[topotus, "Class"] <- as(tax_table(ps16p_100only)[topotus, "Rank_6"], "character")
            tax_table(ps16p_100only) <- (tax_table(toptaxtab))
        #merge samples to simplify complex plot
            sample_data(ps16p_100only)$mergeid<-as.factor(paste0(sample_data(ps16p_100only)$pH, "  ", sample_data(ps16p_100only)$Site, "  ",sample_data(ps16p_100only)$Sample))
            ps16p_100only_m<-merge_samples(ps16p_100only, "mergeid")
            sample_data(ps16p_100only_m)$mergeid<-levels(sample_data(ps16p_100only)$mergeid)
            #repair factors
                sample_data(ps16p_100only_m)$pH<-as.factor(unlist(lapply(strsplit(sample_data(ps16p_100only_m)$mergeid, "  "), '[', 1)))
                sample_data(ps16p_100only_m)$Site<-as.factor(unlist(lapply(strsplit(sample_data(ps16p_100only_m)$mergeid, "  "), '[', 2)))
                sample_data(ps16p_100only_m)$Sample<-as.factor(unlist(lapply(strsplit(sample_data(ps16p_100only_m)$mergeid, "  "), '[', 3)))
            ps16p_100only_mp = transform_sample_counts(ps16p_100only_m, function(x) 10000 * x/sum(x))

        #ensure rank 4 names differentiate unclassified and unclassified bacteria
            taxtab<-as.data.frame(tax_table(ps16p_100only_mp))
            taxtab$Class[as.data.frame(tax_table(ps16p_100only_mp))$Class=="Unclassified"] <-
                    as.data.frame(tax_table(ps16p_100only_mp))$Rank_3[as.data.frame(tax_table(ps16p_100only_mp))$Class=="Unclassified"] #rank 3 names where class is unclassified 
                                                                                                                        #(i.e. which differentiate unclassified and unclassified_bacteria) - 
                                                                                                                        #only needed until pipeline adjusted to output this way automatically
            tax_table(ps16p_100only_mp)<-as.matrix(taxtab)

        #plot
            Figs[["BacteroidiaFamilies"]]<-plot_bar(ps16p_100only_mp, fill="Class", x="pH", facet_grid=sample_Sample~Site) +
                                        geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack")    +
                                        scale_fill_npg() +
                                        scale_color_npg() +
                                        ggtitle("Most Abundant Bacteroidia Families per Site and pH Combination, Abundance Transformed")




    #indicator species stats
            indval<-multipatt(
                            x=t(otu_table(ps16p)), 
                            cluster=as.numeric(as.factor(paste0(#sample_data(ps16p)$pH#, 
                                                                #sample_data(ps16p)$Site, 
                                                                sample_data(ps16p)$Sample 
                                                                ))),
                            control=how(nperm=999))

            summary(indval,alpha=0.001)

            IndicatorSpecies100Fraction<-indval$sign[which( indval$sign$p.value<0.05 & 
                                                            indval$sign$stat >0.999 &
                                                            indval$sign$s.2 ==0 &
                                                            indval$sign$s.3 ==0
                                                            ),]

            tax_table(ps16p) [ rownames(tax_table(ps16p)) %in% rownames(IndicatorSpecies100Fraction) ]
                #this shows that the following are indicators for the 100 fraction:
                    #ESV_15091 "unclassified_Bacteroidia"         "Unclassified"     "Unclassified" 
                    #ESV_3905  "unclassified_Gammaproteobacteria" "Unclassified"     "Unclassified" 
                    #ESV_4912  "MPNO01"                           "MPNO01"           "MPNO01"       # gtdb seq derived from saline water with plankton
                    #ESV_5810  "Pseudomonadales"                  "Cellvibrionaceae" "Thalassocella"  # known from sea surface waters

        # and for the 16s data grouped into bacterial classes
            indval_class<-multipatt(
                            x=(otu_table(ps16p4_mp)), 
                            cluster=as.numeric(as.factor(paste0(sample_data(ps16p4_mp)$pH #, 
                                                                #sample_data(ps16p4_mp)$Site, 
                                                                #sample_data(ps16p4_mp)$Sample 
                                                                ))),
                            control=how(nperm=999))

            summary(indval_class,alpha=0.05)
            tax_table(ps16p4_mp) [ rownames(tax_table(ps16p4_mp)) == "ESV_26946" ]
                #only one weak (85 indval) indicator class found for pH - a Moduliflexia class bacteria

# 7: chemo:community covariance analysis: SI covariance figures

    
        Rich16<- data.frame(
                        R16=(colSums(otu_table(ps16p)  !=0)),
                        ARMS= as.data.frame(sample_data(ps16p))$ARMS,
                        Sample=sample_data(ps16p)$Sample



        )

        Rich23<-data.frame(
                        R23=colSums(otu_table(ps23p)!=0) ,
                        ARMS=sample_data(ps23p)$ARMS,
                        pH =sample_data(ps23p)$pH, 
                        Site=sample_data(ps23p)$Site,
                        TechnicalReplicate=sample_data(ps23p)$TechnicalReplicate
        )
        
        RichChemo<- data.frame(
                        RC= colSums(otu_table(btab_p)  !=0),
                        ARMS=as.data.frame(sample_data(btab_p))$ARMS,
                        Sample=sample_data(btab_p)$Sample
        )

        Richness_df<-rbind(
                                merge(merge(Rich16[Rich16$Sample=="Algae",], 
                                        Rich23,
                                        by="ARMS"), 
                                            RichChemo[RichChemo$Sample=="Algae",],
                                            by="ARMS"), 
                                merge(merge(Rich16[Rich16$Sample=="100",], 
                                        Rich23,
                                        by="ARMS"), 
                                            RichChemo[RichChemo$Sample=="100-500UM",],
                                            by="ARMS"),
                                merge(merge(Rich16[Rich16$Sample=="Sessile",], 
                                        Rich23,
                                        by="ARMS"), 
                                            RichChemo[RichChemo$Sample=="SESSILE",],
                                            by="ARMS")
                                )

        Figs[["SeqRichnessCovariance"]] <-ggplot(Richness_df, aes(x=R23, y=R16, color= pH, shape =Sample.x, label=ARMS)) +
                                            geom_point(size=12) +
                                            geom_label(size=2.5) +
                                            facet_wrap(~Sample.x) +
                                            labs(   title="Microbe-Algae  Richnesses covaraince ",
                                                    subtitle="Holobiont-associated microbe richness declines alongside declining macrobe richness and with ecosystem perturbation.
                                                                However free-living microbes are not similarly perturbed, when combined with the observed pattern of increasing microbial abundance
                                                                this reveal a major shift away from holobionts under perturbation.",
                                                    y="Microbial (16s) ESV Richness", x="Algae (23s) ESV Richness")
        Figs[["MacChemRichnessCovariance"]] <-ggplot(Richness_df, aes(x=R23, y=RC, color= pH, shape =Site, label=ARMS)) +  
                                            geom_point(size=12) +
                                            geom_label(size=2.5) +
                                            facet_wrap(~Sample.x) +
                                            labs(   title="Metabolite-Algae Richness covariance",
                                                    subtitle="",
                                                    y="Metabolite Richness", x="Algae (23s) ESV Richness")
        Figs[["MicChemRichnessCovariance"]] <-ggplot(Richness_df, aes(x=R16, y=RC, color= pH, shape =Site, label=ARMS)) +  
                                            geom_point(size=12) +
                                            geom_label(size=2.5) +
                                            facet_wrap(~Sample.x) +
                                            labs(   title="Metabolite-Microbe Richness covariance",
                                                    subtitle="",
                                                    y="Metabolite Richness", x="Microbial (16s) ESV Richness")


# 8. plot

    #Mains figures
        #Figure1
            grid.arrange(grobs = list(Figs[["Figure1a"]], Figs[["Figure1b"]]), ncol = 2) ## display plot
        #Figure2
            grid.arrange(grobs = list(Figs[["ESVOverlap"]] ), ncol = 1) ## display plot



    #SI figures
        #ordinations
            grid.arrange(grobs = list(  Figs[[6]], Figs[[4]], Figs[[5]],
                                        Figs[[8]], Figs[[9]], Figs[[10]],
                                        Figs[[7]]
            ), ncol = 3) ## display plot

        # algae ordinations
            grid.arrange(grobs = list(  Figs[["Ordination_Algae_23s_PCoA"]], 
                                        Figs[["Ordination_Algae_PCoA"]], 
                                        Figs[["Ordination_Plastids_PCoA"]]
            ), ncol = 1) ## display plot
        #all richnessses
            grid.arrange(grobs = list(Figs[["AllSeqRichness"]], Figs[["AllMetabRichness"]]), ncol = 2) ## display plot
        #covariances
            grid.arrange(grobs = list(Figs[["SeqRichnessCovariance"]], Figs[["MacChemRichnessCovariance"]], Figs[["MicChemRichnessCovariance"]]  ), ncol = 1) ## display plot
        # holbiont to environment richness ratios
            grid.arrange(grobs = list(Figs[["SequenceRatios"]], Figs[["ChemRatios"]],  Figs[["holobDist"]]), ncol = 2) ## display plot
        #holobiont distinctness ordination
            grid.arrange(grobs = list(Figs[["HoloOrd"]] ), ncol = 1) ## display plot
        #taxabarplots
            grid.arrange(grobs=list(Figs[["AlgalClasses"]],  Figs[["BacterialClasses"]],  Figs[["100-BacterialClasses"]],
                                    Figs[["GammaproteobacteriaFamilies"]], Figs[["BacteroidiaFamilies"]]), 
                        ncol=3)
