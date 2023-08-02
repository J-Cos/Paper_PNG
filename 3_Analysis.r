# 1. Packages and path to project directory ######
    set.seed(0.1)
    
    library(tidyverse)
    library(ggplot2)
    library(phyloseq)
    library(breakaway)
    library(DivNet)
    library(lme4)
    library(lmerTest)

    library(permute)
    library(vegan)
    library(DESeq2)

    #params
        path<-"/home/j/Dropbox/PNG_Paper" #provide the path to the project folder
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



    #1.1 load functions
        #funcs
        function_files<-list.files(file.path(path, "Functions"))
        sapply(file.path(path, "Functions",function_files),source)
        
        function_files<-list.files(file.path("/home/j/Dropbox/BioinformaticPipeline_Env/BioinformaticPipeline/SupportFunctions"))
        sapply(file.path("/home/j/Dropbox/BioinformaticPipeline_Env/BioinformaticPipeline/SupportFunctions",function_files),source)

        GetProportionEsvsAssigned<-function(ps, Rank, UnclassifiedLabels){
            tax_table(ps)%>% 
                as.data.frame %>% 
                select(Rank) %>% 
                filter(! .[[Rank]]  %in% UnclassifiedLabels ) %>%
                dim %>%
                '['(1) %>%
                '/'(length(taxa_names(ps)))*100
        }
        GetProportionReadsAssignedAtPhyla<-function(ps, UnclassifiedLabels){
            ps %>%
                prune_taxa(!tax_table(.)[,3] %in% UnclassifiedLabels, .) %>%
                sample_sums(.) %>%
                sum %>%
                '/' (sum(sample_sums(ps)))*100
        }
        GetProportionCompoundsWithSpectrumID<-function(ps){
            tax_table(ps)%>% 
                as.data.frame %>% 
                select(SpectrumID) %>% 
                filter(! .[["SpectrumID"]]  %in% c("Unclassified") ) %>%
                dim %>%
                '['(1) %>%
                '/'(length(taxa_names(ps)))*100
        }
        GetProportionMoleculesWithSpectrumID<-function(ps){
            ps %>%
                prune_taxa(!tax_table(.)[,2] %in% c("Unclassified"), .) %>%
                sample_sums(.) %>%
                sum %>%
                '/' (sum(sample_sums(ps)))*100
        }
        MakeTableColumnForOneMetabarcodingFraction<-function(ps, UnclassifiedLabels, multiSample=FALSE) {
            #pre calculate results for rate limiting calculations
                ps_phyla<-tax_glom(ps, taxrank="Rank_3")
                alphaEstimates<-ps %>% breakaway %>% summary %>% '['("estimate") %>%unlist  
                if (multiSample){
                    shannonEstimates<-ps_phyla %>% DivNet::divnet(., formula = ~Sample +Sample:pH + Site) %>%'['("shannon") %>% unlist(recursive=FALSE) %>%lapply(., '[', "estimate") %>% unlist %>% as.vector
                } else {
                    shannonEstimates<-ps_phyla %>% DivNet::divnet(., formula = ~pH + Site) %>%'['("shannon") %>% unlist(recursive=FALSE) %>%lapply(., '[', "estimate") %>% unlist %>% as.vector
                }
            #create column
                column<-c(
                            paste0( # num samples and per pH level
                                nsamples(ps), 
                                " ",
                                table(sample_data(ps)$pH)[c(1,3,2)] %>%as.vector %>%paste(collapse=";") %>% paste0("(", ., ")")
                            ),
                            paste0( # average EVSs
                                median(colSums(otu_table(ps)>0))%>% format(big.mark=","),
                                " (",
                                mean(colSums(otu_table(ps)>0))%>% round %>% format(big.mark=","),
                                ")"
                            ),
                            paste0( # average reads
                                median(sample_sums(ps))%>% format(big.mark=","),
                                " (",
                                mean(sample_sums(ps))%>% round %>% format(big.mark=","),
                                ")"
                            ),
                            paste0( # percent identiifieds
                                GetProportionEsvsAssigned(ps=ps, Rank="Rank_3", UnclassifiedLabels=UnclassifiedLabels) %>% signif (3),
                                " (",
                                GetProportionReadsAssignedAtPhyla(ps=ps,  UnclassifiedLabels=UnclassifiedLabels)%>% signif (3),
                                ")"
                            ),
                            paste0( # number phyla occuring
                                ps_phyla%>% otu_table %>% as("matrix") %>% '>'(1) %>% colSums %>% median,
                                " (", 
                                ps_phyla %>% otu_table %>% as("matrix") %>% '>'(1) %>% colSums %>% mean %>% signif (3),
                                ")"
                            ),
                            paste0( # average richenss
                                alphaEstimates %>% median%>% round %>% format(big.mark=","),
                                " (",
                                alphaEstimates %>% mean%>% round %>% format(big.mark=","),
                                ")"
                            ),
                            paste0( # average shannon
                                shannonEstimates %>% median%>% signif (3),
                                " (",
                                shannonEstimates %>% mean%>% signif (3),
                                ")"
                            )
                        )
            return(column)
        }
        MakeTableColumnForOneMetabolomicFraction<-function(ps, multiSample=FALSE) {
            #pre calculate results for rate limiting calculations
                ps<-ps%>% transform_sample_counts(., function(x) floor(x))     
                #not usign breakaway as the structure of metabolomic data is different to genetic and does not make sense to use same statistical estiamtion approach
                alphaEstimates<-estimate_richness(ps, measures=c("Observed")) %>% unlist 
                shannonEstimates<-estimate_richness(ps, measures=c("Shannon")) %>% unlist

            #create column
                column<-c(
                            paste0( # num samples and per pH level
                                nsamples(ps), 
                                " ",
                                table(sample_data(ps)$pH)[c(1,3,2)] %>%as.vector %>%paste(collapse=";") %>% paste0("(", ., ")")
                            ),
                            paste0( # average EVSs
                                median(colSums(otu_table(ps)>0))%>% format(big.mark=","), 
                                " (",
                                mean(colSums(otu_table(ps)>0))%>% round %>% format(big.mark=","),
                                ")"
                            ),
                            paste0( # average reads
                                median(sample_sums(ps))%>%  format(big.mark=","), 
                                " (",
                                mean(sample_sums(ps))%>% round %>% format(big.mark=","),
                                ")"
                            ),
                            paste0( # percent identifieds
                                GetProportionCompoundsWithSpectrumID(ps) %>% signif (3),
                                " (",
                                GetProportionMoleculesWithSpectrumID(ps)%>% signif (3),
                                ")"
                            ),
                            paste0( # number phyla occuring
                                "NA"
                            ),
                            paste0( # average richenss
                                alphaEstimates %>% median%>% round %>% format(big.mark=","),
                                " (",
                                alphaEstimates %>% mean%>% round %>% format(big.mark=","),
                                ")"
                            ),
                            paste0( # average shannon
                                shannonEstimates %>% median%>% signif (3),
                                " (",
                                shannonEstimates %>% mean%>% signif (3),
                                ")"
                            )
                        )
            return(column)
        }

# 2. Data prep
    # get data
        ps16<-readRDS(file=file.path(path, "Outputs", "ps16.RDS")) %>% prune_taxa( as.data.frame(tax_table(.) )$Rank_5!="Chloroplast", .)
        ps23<-readRDS(file=file.path(path, "Outputs", "ps23.RDS"))
        btab<-readRDS(file=file.path(path, "Outputs", "btab.RDS"))

# 3. (descriptive quantities for text)
    #num samples and num ESVs
    ps16
    ps23
    #number ARMS
    length(unique(sample_data(ps16)$ARMS))
    #num samples per site
    table(sample_data(ps16)$Site) + table(sample_data(ps23)$Site)
    #number samples per fraction
    c(table(sample_data(ps16)$Sample),table(sample_data(ps23)$Sample))
    #number samples per pH
    table(sample_data(ps16)$pH) + table(sample_data(ps23)$pH)
    #16s and 23s proportions assigned at phyla
    GetProportionEsvsAssigned(ps=ps16, Rank="Rank_3", UnclassifiedLabels=c("Unclassified","unclassified_Bacteria"))
    GetProportionReadsAssignedAtPhyla(ps=ps16, UnclassifiedLabels=c("Unclassified","unclassified_Bacteria"))
    GetProportionEsvsAssigned(ps=ps23, Rank="Rank_3", UnclassifiedLabels=c("unclassified_Root", "unclassified_Chromista", "unclassified_Plantae"))
    GetProportionReadsAssignedAtPhyla(ps=ps23,  UnclassifiedLabels=c("unclassified_Root", "unclassified_Chromista", "unclassified_Plantae"))
    #number of phyla
    length(unique(tax_table(ps16)[,3]))
    length(unique(tax_table(ps23)[,3]))

    btab
    #number ARMS
    sum(!is.na(unique(sample_data(btab)$ARMS)))
    #num samples per site
    table(sample_data(btab)$Site)
    #number samples per fraction
    table(sample_data(btab)$Sample)
    #number samples per pH
    table(sample_data(btab)$pH) 


    #make table of bacterial class mean abundance and sd per sample type
        newps<-
            ps16 %>%
            prune_samples( sample_sums(.)>100000, .)  %>%
            prune_taxa( taxa_sums(.)>0, .) %>%
            tax_glom(taxrank="Rank_4") %>%
            transform_sample_counts(., function(x) x / sum(x) )

        sampleDat<-newps %>% 
            sample_data %>%
            as_tibble %>%
            mutate(mergeCol=as_factor(paste0(pH, "_", Sample))) %>%
            select(mergeCol)

        table<-otu_table(newps) %>% 
                t 

        table<-cbind(table, sampleDat) %>%
                as_tibble

        meanTable<-table %>%
            dplyr::group_by(mergeCol) %>%
            dplyr::summarise_all(mean) %>%
            mutate(mergeCol=paste0(mergeCol,"_mean"))
        
        sdTable<-table %>%
            dplyr::group_by(mergeCol) %>%
            dplyr::summarise_all(sd)%>%
            mutate(mergeCol=paste0(mergeCol,"_sd"))
        
        combinedTable<-rbind(meanTable, sdTable)


        combinedTable<-column_to_rownames(combinedTable, "mergeCol") 
            
        colnames(combinedTable)<-c( 
                            paste0(tax_table(newps)[,2], "_" ,tax_table(newps)[,3], "_", tax_table(newps)[,4]))
        

        combinedTable<-t(combinedTable)
        write.table(file="Outputs/BacterialRelativeAbundances.csv", row.names=TRUE, col.names=NA, sep=",", combinedTable)

        #ps16<-prune_samples( sample_data(ps16)$ARMS %in% sample_data(ps23)$ARMS, ps16)
        #ps23<-prune_samples( sample_data(ps23)$ARMS %in% sample_data(ps16)$ARMS, ps23)
        ps16<- prune_taxa( taxa_sums(ps16)>0, ps16)
        ps23<- prune_taxa( taxa_sums(ps23)>0, ps23)

        ps23<-rarefy_even_depth (ps23, rngseed=1)
        ps16<-rarefy_even_depth (ps16, rngseed=1)


# 4. Make Table
    #Metabarcoding
        col23<-ps23 %>% 
                MakeTableColumnForOneMetabarcodingFraction(ps=., 
                                                    UnclassifiedLabels=c("unclassified_Root", 
                                                                        "unclassified_Chromista", 
                                                                        "unclassified_Plantae")
                                        )
        colphoto<-ps16 %>% prune_samples(sample_data(.)$Sample %in% c("Algae"), .) %>%  prune_taxa(taxa_sums(.)>1, .) %>%
                MakeTableColumnForOneMetabarcodingFraction(ps=., 
                                                    UnclassifiedLabels=c("Unclassified",
                                                                        "unclassified_Bacteria")
                                        )
        colenv<-ps16 %>% prune_samples(sample_data(.)$Sample %in% c("100"), .) %>%  prune_taxa(taxa_sums(.)>1, .) %>%
                MakeTableColumnForOneMetabarcodingFraction(ps=., 
                                                    UnclassifiedLabels=c("Unclassified",
                                                                        "unclassified_Bacteria")
                                        )
        colholo<-ps16 %>% prune_samples(sample_data(.)$Sample %in% c("Sessile"), .) %>%  prune_taxa(taxa_sums(.)>1, .) %>%
                MakeTableColumnForOneMetabarcodingFraction(ps=., 
                                                    UnclassifiedLabels=c("Unclassified",
                                                                        "unclassified_Bacteria")
                                            )
        colorg<-ps16 %>% prune_samples(sample_data(.)$Sample %in% c("RVS", "CS"), .) %>%  prune_taxa(taxa_sums(.)>1, .) %>%
                MakeTableColumnForOneMetabarcodingFraction(ps=., 
                                                    UnclassifiedLabels=c("Unclassified",
                                                                        "unclassified_Bacteria"),
                                                    multiSample=TRUE
                                            )    
        col16<-ps16 %>% 
                MakeTableColumnForOneMetabarcodingFraction(ps=., 
                                                    UnclassifiedLabels=c("Unclassified",
                                                                        "unclassified_Bacteria"),
                                                    multiSample=TRUE)

    #Metabolomics
        colMetabAlg<-btab%>% prune_samples(sample_data(.)$Sample %in% c("Algae"), .) %>%  prune_taxa(taxa_sums(.)>1, .) %>% MakeTableColumnForOneMetabolomicFraction
        colMetabEnv<-btab%>% prune_samples(sample_data(.)$Sample %in% c("100-500UM"), .) %>%  prune_taxa(taxa_sums(.)>1, .) %>% MakeTableColumnForOneMetabolomicFraction
        colMetabHolo<-btab%>% prune_samples(sample_data(.)$Sample %in% c("SESSILE"), .) %>%  prune_taxa(taxa_sums(.)>1, .) %>% MakeTableColumnForOneMetabolomicFraction
        colMetabOrg<-btab%>% prune_samples(sample_data(.)$Sample %in% c("RVS", "CS"), .) %>%  prune_taxa(taxa_sums(.)>1, .) %>% MakeTableColumnForOneMetabolomicFraction
        colMetabAll<-btab%>% MakeTableColumnForOneMetabolomicFraction

    write.csv(cbind(    col23,
                        colphoto,
                        colenv,
                        colholo,
                        colorg,
                        col16, 
                        colMetabAlg,
                        colMetabEnv,
                        colMetabHolo,
                        colMetabOrg,
                        colMetabAll ), 
                file=file.path(path, "Outputs", "tabledata.csv"))

# 5 . Models
    #5.1 Betta 
        # richness model
            #psmerge<-merge_phyloseq(ps16, ps23)
            psmerge<-ps16  #changed from psmerge to ps16 as droppig 23s data, beta model fraction mutate also adjusted
            meta <- psmerge %>% 
                sample_data %>%
                as_tibble %>%
                mutate("sample_names" = psmerge %>% sample_names )
            richnesspremerge<-psmerge %>% breakaway
            richnessmerge <- meta %>%
                left_join(summary(richnesspremerge),
                            by = "sample_names") %>%
                mutate_if( is.character,as_factor) %>%
                mutate(ARMS=as.factor(ARMS))
        
            bettamod1<-richnessmerge %>%
                mutate  (pH=factor(pH, order = TRUE,  levels = c("Control pH", "Medium pH", "Low pH")) ) %>%
                mutate  (Sample=factor(Sample,  levels = c("Environmental Microbiome", "Holobiont Community Microbiome","Photosynthetic Community Microbiome", "Halisarca Sponge Microbiome", "Tethys Sponge Microbiome"))) %>% #, "Photosynthetic Community Algae"))) %>%
                betta_random(formula = estimate ~ (Sample*pH)|Site , ses = error, data=.)
            
            bettamod1$global

            bettamod1Table<-cbind( bettamod1$table,
                    bettamod1$table[,1]+bettamod1$table[,2]*1.959, 
                    bettamod1$table[,1]-bettamod1$table[,2]*1.959)
            

            write.csv(signif(bettamod1Table,2), file=file.path(path, "Outputs", "MetabarcodingRichnessModel.csv"))
        #shannon model
            #dv23<-ps23 %>% tax_glom(taxrank="Rank_3") %>% DivNet::divnet(., formula = ~pH + Site)
            dv16<-ps16 %>% c%>% DivNet::divnet(., formula = NULL) #~Sample*pH + Site) # changed to null to treat all samples as independent observations

            #shannons <- dv23$shannon %>% summary %>% rbind(dv16$shannon %>% summary )
            shannons <- dv16$shannon %>% summary #changed from merged as 23s removed
            shannonmerge <- meta %>%
                left_join(shannons,
                            by = "sample_names") %>%
                mutate_if( is.character,as_factor) %>%
                mutate(ARMS=as.factor(ARMS))

            bettamod2<-shannonmerge %>%
                mutate  (pH=factor(pH, order = TRUE,  levels = c("Control pH", "Medium pH", "Low pH")) ) %>%
                mutate  (Sample=factor(Sample,  levels = c( "Environmental Microbiome", "Holobiont Community Microbiome", "Photosynthetic Community Microbiome", "Halisarca Sponge Microbiome", "Tethys Sponge Microbiome"))) %>%#, "Photosynthetic Community Algae"))) %>%
                betta_random(formula = estimate ~ Sample*pH|Site, ses = error, data=.)

            bettamod2Table<-cbind( bettamod2$table,
                    bettamod2$table[,1]+bettamod2$table[,2]*1.959, 
                    bettamod2$table[,1]-bettamod2$table[,2]*1.959)
            bettamod2$global

            write.csv(signif(bettamod2Table,2), file=file.path(path, "Outputs", "MetabarcodingShannonModel.csv"))
            

        # non-betta shannon model
            EsvShannon<-#ps23 %>% # ps23 removed
                #estimate_richness(., measures=c("Shannon")) %>%
                #rownames_to_column("sample_names") %>%
                #rbind( 
                    ps16 %>%
                            estimate_richness(., measures=c("Shannon")) %>%
                            rownames_to_column("sample_names")
                #)
                
            ESVShannonMerge <- meta %>%
                left_join(EsvShannon,
                            by = "sample_names") %>%
                mutate_if( is.character,as_factor) %>%
                mutate(ARMS=as.factor(ARMS))


            EsvShannonmod1<-ESVShannonMerge %>%
                #lme4::lmer(data = ., formula=Shannon ~ pH:Sample +Sample +(1|Site))
                nlme::lme(Shannon ~ Sample*pH, data = ., random = ~ 1|Site)

            EsvShannonmodNULL<-ESVShannonMerge %>%
                nlme::lme(Shannon ~ Sample, data = ., random = ~ 1|Site)

                qqnorm(residuals(EsvShannonmod1))
                scatter.smooth(residuals(EsvShannonmod1) ~ fitted(EsvShannonmod1))
                qqnorm(residuals(EsvShannonmodNULL))
                scatter.smooth(residuals(EsvShannonmodNULL)~ fitted(EsvShannonmodNULL))
                #assumptions look good with poisson distribution (makes sense given count data)

                anova(EsvShannonmod1, EsvShannonmodNULL)
                    #no significant efefct of pH:sample
                confs<-nlme::intervals(EsvShannonmod1, which="fixed")
                    # pH confints cross zero for all samples
                summary(EsvShannonmod1)

                MuMIn::r.squaredGLMM(EsvShannonmod1)

            cbind(
                -0.6261659+0.2823125*1.959, 
                -0.6261659-0.2823125*1.959)

    #5.2 chemical richness (non-betta)
            chemRichness<-btab %>%
                transform_sample_counts(., function(x) floor(x))     %>%
                estimate_richness(., measures=c("Observed", "Shannon")) %>%
                rownames_to_column("sample_names")

            chemRichnessMerge <- btab %>%
                    sample_data %>%
                    as_tibble %>%
                    mutate("sample_names" = btab %>% sample_names )%>%
                        left_join(chemRichness,
                                    by = "sample_names") %>%
                        mutate_if( is.character,as_factor) %>%
                        mutate(ARMS=as.factor(ARMS)) %>%
                        mutate  (pH=factor(pH, order = TRUE,  levels = c("Control pH", "Medium pH", "Low pH")) ) %>%
                        mutate  (Sample=factor(Sample,  levels = c("Environmental Metabolome", "Holobiont Community Metabolome", "Photosynthetic Community Metabolome", "Halisarca Sponge Metabolome", "Tethys Sponge Metabolome")))
            #observed
                chemrichnessmod1<-chemRichnessMerge %>%
                        #lme4::lmer(data = ., formula=Observed ~ pH*Sample +(1|Site))
                        nlme::lme(Observed ~ Sample*pH, data = ., random = ~ 1|Site)


                chemrichnessmodNULL<-chemRichnessMerge %>%
                        #lme4::lmer(data = ., formula=Observed ~ Sample +(1|Site))
                        nlme::lme(Observed ~ Sample, data = ., random = ~ 1|Site)
                        qqnorm(residuals(chemrichnessmod1))
                        scatter.smooth(residuals(chemrichnessmod1) ~ fitted(chemrichnessmod1))
                        qqnorm(residuals(chemrichnessmodNULL))
                        scatter.smooth(residuals(chemrichnessmodNULL))
                        #assumptions look good with poisson distribution (makes sense given count data)

                        #anova(chemrichnessmod1, chemrichnessmodNULL)
                            #no significant efefct of pH:sample
                        #confs<-confint(chemrichnessmod1)
                            # pH confints cross zero for all samples
                        summary(chemrichnessmod1)

                        MuMIn::r.squaredGLMM(chemrichnessmod1)
               chemRichnessMerge %>%
                        ggplot() +
                            geom_boxplot(aes(x=pH, y=Observed, fill=Sample))
            #shannon
                chemShannonmod1<-chemRichnessMerge %>%
                    #lme4::lmer(data = ., formula=Shannon ~ pH:Sample +Sample +(1|Site))
                    nlme::lme(Shannon ~ Sample*pH, data = ., random = ~ 1|Site)

                chemShannonmodNULL<-chemRichnessMerge %>%
                    lme4::lmer(data = ., formula=Shannon ~ Sample +(1|Site))

                    qqnorm(residuals(chemShannonmod1))
                    scatter.smooth(residuals(chemShannonmod1) ~ fitted(chemShannonmod1))
                    qqnorm(residuals(chemShannonmodNULL))
                    scatter.smooth(residuals(chemShannonmodNULL)~ fitted(chemShannonmodNULL))
                    #assumptions look good with poisson distribution (makes sense given count data)

                    anova(chemShannonmod1, chemShannonmodNULL)
                        #no significant efefct of pH:sample
                    confs<-confint(chemShannonmod1)
                        # pH confints cross zero for all samples
                    summary(chemShannonmod1)

                    MuMIn::r.squaredGLMM(chemShannonmod1)

                cbind(
                    -0.6261659+0.2823125*1.959, 
                    -0.6261659-0.2823125*1.959)
    #richness figure
    
        sgene<-shannonmerge  %>% filter(Sample!="Photosynthetic Community Algae") %>% select(pH, Sample, estimate) %>% mutate(Type="Metabarcoding", Diversity="Shannon Diversity")%>% 
            mutate  (Sample=factor(Sample,  levels = c( "Environmental Microbiome", "Holobiont Community Microbiome", "Photosynthetic Community Microbiome", #"Photosynthetic Community Algae", 
                                                        "Halisarca Sponge Microbiome", "Tethys Sponge Microbiome"))) %>%    
            mutate(Sample=fct_relevel(Sample,c(#"Photosynthetic Community Algae", 
                                                "Photosynthetic Community Microbiome", "Environmental Microbiome", "Holobiont Community Microbiome", "Halisarca Sponge Microbiome", "Tethys Sponge Microbiome") )) %>%
            ggplot() +
                geom_boxplot(aes(y=estimate, x=pH, fill=pH)) +
                pHfillScale +
                facet_grid(~Sample)+
                ggtitle("Sequence Phylum-level Shannon Diversity") +
                theme(  strip.background = element_rect(fill = "white"),
                        #panel.border = element_rect(colour = "black", fill = "white"),
                        legend.position="none",
                        strip.text = element_blank(), 
                        axis.text.x = element_blank(), 
                        axis.title = element_blank(), 
                        plot.title = element_text(size = 20, face = "bold", hjust=0.5),
                        panel.grid.major = element_blank(), 
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(), 
                        axis.line = element_line(colour = "black"),
                        axis.ticks = element_blank())
        #sponges++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        schem<-chemRichnessMerge  %>% select(pH, Sample,Shannon) %>% mutate(Type="Metabolomics", Diversity="Shannon Diversity") %>% rename(Shannon="estimate")%>% 
            mutate  (Sample=factor(Sample,  levels = c("Environmental Metabolome", "Holobiont Community Metabolome", "Photosynthetic Community Metabolome", "Halisarca Sponge Metabolome", "Tethys Sponge Metabolome")))%>%
            mutate(Sample=fct_relevel(Sample,c("Photosynthetic Community Metabolome", "Environmental Metabolome", "Holobiont Community Metabolome", "Halisarca Sponge Metabolome", "Tethys Sponge Metabolome") )) %>%
            ggplot() +
                geom_boxplot(aes(y=estimate, x=pH, fill=pH)) +
                pHfillScale +
                facet_grid(~Sample)+
                ggtitle("Metabolite Shannon Diversity") +
                theme(  strip.background = element_rect(fill = "white"),
                        #panel.border = element_rect(colour = "black", fill = "white"),
                        legend.title=element_text(size=16, face="bold"),
                        legend.text=element_text(size=16, face="bold"),
                        strip.text = element_blank(), 
                        axis.text.x = element_blank(), 
                        axis.title = element_blank(), 
                        plot.title = element_text(size = 20, face = "bold", hjust=0.5),
                        panel.grid.major = element_blank(), 
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(), 
                        axis.line = element_line(colour = "black"),
                        axis.ticks = element_blank())

        rgene<-richnessmerge %>% filter(Sample!="Photosynthetic Community Algae") %>% select(pH, Sample, estimate) %>% mutate(Type="Metabarcoding", Diversity="Richness")%>% 
            mutate  (Sample=factor(Sample,  levels = c( "Environmental Microbiome", "Holobiont Community Microbiome", "Photosynthetic Community Microbiome", #"Photosynthetic Community Algae", 
                                                        "Halisarca Sponge Microbiome", "Tethys Sponge Microbiome"))) %>%          
            mutate(Sample=fct_relevel(Sample,c(#"Photosynthetic Community Algae", 
                                                "Photosynthetic Community Microbiome", "Environmental Microbiome", "Holobiont Community Microbiome", "Halisarca Sponge Microbiome", "Tethys Sponge Microbiome") )) %>%
            ggplot() +
                geom_boxplot(aes(y=estimate, x=pH, fill=pH)) +
                pHfillScale +
                facet_grid(~Sample, labeller = labeller(Sample = 
                                                    c(  "Environmental Microbiome"= "Environmental\nMicrobiome\n", 
                                                        "Holobiont Community Microbiome" = "Holobiont\nCommunity\nMicrobiome", 
                                                        "Photosynthetic Community Microbiome"="Photosynthetic\nCommunity\nMicrobiome", 
                                                        "Halisarca Sponge Microbiome" ="Halisarca\nSponge\nMicrobiome", 
                                                        "Tethys Sponge Microbiome"="Tethya\nSponge\nMicrobiome"#,
                                                        #"Photosynthetic Community Algae" = "Photosynthetic\nCommunity"
                                                        ))) +
                ggtitle("Sequence Estimated Richness") +
                theme(  strip.background = element_rect(fill = "white"),
                        #panel.border = element_rect(colour = "black", fill = "white"),
                        legend.position="none",
                        strip.text.x = element_text(color="black", face="bold", size=8),
                        axis.text.x = element_blank(), 
                        axis.title = element_blank(), 
                        plot.title = element_text(size = 20, face = "bold", hjust=0.5),
                        panel.grid.major = element_blank(), 
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(), 
                        axis.line = element_line(colour = "black"),
                        axis.ticks = element_blank())

        rchem <-chemRichnessMerge  %>% select(pH, Sample, Observed) %>% mutate(Type="Metabolomics", Diversity="Richness") %>% rename(Observed="estimate")%>% 
            mutate  (Sample=factor(Sample,  levels = c("Environmental Metabolome", "Holobiont Community Metabolome", "Photosynthetic Community Metabolome", "Halisarca Sponge Metabolome", "Tethys Sponge Metabolome")))%>%
            mutate(Sample=fct_relevel(Sample,c("Photosynthetic Community Metabolome", "Environmental Metabolome", "Holobiont Community Metabolome", "Halisarca Sponge Metabolome", "Tethys Sponge Metabolome") )) %>%
            ggplot() +
                geom_boxplot(aes(y=estimate, x=pH, fill=pH)) +
                pHfillScale +
                facet_grid(~Sample, labeller = labeller(Sample = 
                                                    c(  "Environmental Metabolome"= "Environmental\nMetabolome\n", 
                                                        "Holobiont Community Metabolome" = "Holobiont\nCommunity\nMetabolome", 
                                                        "Photosynthetic Community Metabolome"="Photosynthetic\nCommunity\nMetabolome", 
                                                        "Halisarca Sponge Metabolome" ="Halisarca\nSponge\nMetabolome", 
                                                        "Tethys Sponge Metabolome"="Tethya\nSponge\nMetabolome"))) +
                ggtitle("Metabolite Estimated Richness") +
                theme(  strip.background = element_rect(fill = "white"),
                        #panel.border = element_rect(colour = "black", fill = "white"),
                        legend.position="none",
                        strip.text.x = element_text(color="black", face="bold", size=8),
                        axis.text.x = element_blank(), 
                        axis.title = element_blank(), 
                        plot.title = element_text(size = 20, face = "bold", hjust=0.5),
                        panel.grid.major = element_blank(), 
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(), 
                        axis.line = element_line(colour = "black"),
                        axis.ticks = element_blank())

        fullRichnessPlot<-egg::ggarrange(rgene,rchem,sgene, schem)#, labels=c("A", "B", "C", "D"), top="Figure S2: Multiomic richness and Shannon Diversity for all Sample Types")
    #5.3 permanova
        #every fraction individually - revised approach - advised by Emma

        PERMANOVA_inidividual_fractions<-function(ps, fraction, sites=2) {
                Mor<-ps %>%
                    prune_samples(sample_data(.)$Sample %in% fraction, .) %>% 
                    prune_taxa(taxa_sums(.)>1, .) %>% 
                    phyloseq::distance(., method="morisita")  #bray (abund) or jaccard (p/a) - [results the same so using jaccard]

            
            #permanova
            if (sites==2) {
                permMod<-ps %>%
                    prune_samples(sample_data(.)$Sample %in% fraction, .) %>% 
                    prune_taxa(taxa_sums(.)>1, .) %>%
                    sample_data(.) %>%
                    as("data.frame") %>%
                    as_tibble %>%
                    mutate_if(is.character, as_factor) %>%
                    mutate  (pH=factor(pH, order = TRUE,  levels = c("Control pH", "Medium pH", "Low pH")) ) %>%
                    mutate(ARMS=factor(ARMS)) %>%
                    as.data.frame %>%
                    #vegan::adonis2(dist ~   Sample/pH , data=., strata=.$Sample)
                    vegan::adonis2(Mor~ Site + pH, data=., method="morisita",  permutations = how(plots = Plots(strata = .$Site), within = Within(), nperm=9999), by="terms")
            } else {
                permMod<-ps %>%
                    prune_samples(sample_data(.)$Sample %in% fraction, .) %>% 
                    prune_taxa(taxa_sums(.)>1, .) %>%
                    sample_data(.) %>%
                    as("data.frame") %>%
                    as_tibble %>%
                    mutate_if(is.character, as_factor) %>%
                    mutate  (pH=factor(pH, order = TRUE,  levels = c("Control pH", "Medium pH", "Low pH")) ) %>%
                    mutate(ARMS=factor(ARMS)) %>%
                    as.data.frame %>%
                    #vegan::adonis2(dist ~   Sample/pH , data=., strata=.$Sample)
                    vegan::adonis2(Mor~ pH, data=., method="morisita", by="terms")
            } 
                    #pairwiseAdonis::pairwise.adonis2(jacc~ Sample +pH , data=., method="jaccard", strata ='Site')

            #betadisper
                dispermod<-ps %>%
                    prune_samples(sample_data(.)$Sample %in% fraction, .) %>% 
                    prune_taxa(taxa_sums(.)>1, .) %>%
                    sample_data(.) %>%
                    as("data.frame") %>%
                    as_tibble %>%
                    mutate_if(is.character, as_factor) %>%
                    mutate  (pH=factor(pH, order = TRUE,  levels = c("Control pH", "Medium pH", "Low pH")) ) %>%
                    select(pH) %>%
                    as.vector %>%
                    '[['(1) %>%
                    betadisper(d=Mor, group=.) %>%
                    anova

                return(list(dispermod, permMod))
        }

        Teth<-PERMANOVA_inidividual_fractions(ps=ps16, fraction="Tethys Sponge Microbiome")
        Hali<-PERMANOVA_inidividual_fractions(ps=ps16, fraction="Halisarca Sponge Microbiome", sites=1)
        Env<-PERMANOVA_inidividual_fractions(ps=ps16, fraction="Environmental Microbiome")
        Holo<-PERMANOVA_inidividual_fractions(ps=ps16, fraction="Holobiont Community Microbiome")
        Phot<-PERMANOVA_inidividual_fractions(ps=ps16, fraction="Photosynthetic Community Microbiome")
        Phot23<-PERMANOVA_inidividual_fractions(ps=ps23, fraction="Photosynthetic Community Algae")

        Teth_M<-PERMANOVA_inidividual_fractions(ps=btab, fraction="Tethys Sponge Metabolome")
        Hali_M<-PERMANOVA_inidividual_fractions(ps=btab, fraction="Halisarca Sponge Metabolome", sites=1)
        Env_M<-PERMANOVA_inidividual_fractions(ps=btab, fraction="Environmental Metabolome")
        Holo_M<-PERMANOVA_inidividual_fractions(ps=btab, fraction="Holobiont Community Metabolome")
        Phot_M<-PERMANOVA_inidividual_fractions(ps=btab, fraction="Photosynthetic Community Metabolome")

        #whoel arms fractions
            #get distances

                Mor16<-ps16 %>%
                    prune_samples(!sample_data(.)$Sample %in% c("Tethys Sponge Microbiome", "Halisarca Sponge Microbiome"), .) %>% 
                    prune_taxa(taxa_sums(.)>1, .) %>% 
                    phyloseq::distance(., method="morisita")  #bray (abund) or jaccard (p/a) - [results the same so using jaccard]

            #permanova
                ps16 %>%
                    prune_samples(!sample_data(.)$Sample %in% c("Tethys Sponge Microbiome", "Halisarca Sponge Microbiome"), .) %>% 
                    prune_taxa(taxa_sums(.)>1, .) %>%
                    sample_data(.) %>%
                    as("data.frame") %>%
                    as_tibble %>%
                    mutate_if(is.character, as_factor) %>%
                    mutate  (pH=factor(pH, order = TRUE,  levels = c("Control pH", "Medium pH", "Low pH")) ) %>%
                    mutate(ARMS=factor(ARMS)) %>%
                    as.data.frame %>%
                    #vegan::adonis2(dist ~   Sample/pH , data=., strata=.$Sample)
                    vegan::adonis2(Mor16~ Sample + Site + pH + Sample:pH , data=., method="morisita",  permutations = how(plots = Plots(strata = .$Site), within = Within(), nperm=9999), by="terms")

                    #pairwiseAdonis::pairwise.adonis2(jacc~ Sample +pH , data=., method="jaccard", strata ='Site')

            #betadisper
                dispermod16<-ps16 %>%
                    prune_samples(!sample_data(.)$Sample %in% c("Tethys Sponge Microbiome", "Halisarca Sponge Microbiome"), .) %>% 
                    prune_taxa(taxa_sums(.)>1, .) %>%
                    sample_data(.) %>%
                    as("data.frame") %>%
                    as_tibble %>%
                    mutate_if(is.character, as_factor) %>%
                    mutate  (pH=factor(pH, order = TRUE,  levels = c("Control pH", "Medium pH", "Low pH")) ) %>%
                    mutate(ARMS=factor(ARMS)) %>%
                    as.data.frame %>%
                    betadisper(Mor16, groups)
                    

                    ## Perform test
                    anova(mod)

                    ## Permutation test for F
                    permDisper(mod)

                    ## Plot the groups and distances to centroids on the
                    ## first two PCoA axes
                    plot(mod)

                    ## Draw a boxplot of the distances to centroid for each group
                    boxplot(mod)


            #associated nmds
                ord16<-ps16 %>%
                    prune_samples(!sample_data(.)$Sample %in% c("Tethys Sponge Microbiome", "Halisarca Sponge Microbiome"), .) %>% 
                    prune_taxa(taxa_sums(.)>1, .) %>% 
                    ordinate(., method="NMDS", distance=Mor16, trymax = 10000)                

                Figs[["NMDS_WholeARMS"]]<-ps16 %>%
                                            prune_samples(!sample_data(.)$Sample %in% c("Tethys Sponge Microbiome", "Halisarca Sponge Microbiome"), .) %>% 
                                            prune_taxa(taxa_sums(.)>1, .) %>% 
                                            plot_ordination(., ord16, type="Sample", color="pH", shape="Sample") +
                                                geom_point(size=5)+
                                                ggtitle(paste0("Community Microbiome Composition\n(Stress=", signif(ord16$stress,3),")")) +
                                                pHcolScale +
                                                pHfillScale +
                                                scale_shape_manual(values = c(0, 1, 2)) +
                                                theme(  strip.background = element_rect(fill = "white"),
                                                        legend.position="none",
                                                        #panel.border = element_rect(colour = "black", fill = "white"),
                                                        strip.text = element_text(size=14, color="black", face="bold"),
                                                        axis.text = element_blank(), 
                                                        axis.title = element_blank(), 
                                                        axis.ticks=element_blank(),
                                                        plot.title = element_text(size = 12, face = "bold", hjust=0.5),
                                                        panel.grid.major = element_blank(), 
                                                        panel.grid.minor = element_blank(),
                                                        panel.background = element_blank(), 
                                                        axis.line = element_line(colour = "black"))+
                                                guides(color="none")
        #sponges
         #get distances

                MorSponge<-ps16 %>%
                    prune_samples(sample_data(.)$Sample %in% c("Tethys Sponge Microbiome", "Halisarca Sponge Microbiome"), .) %>% 
                    prune_taxa(taxa_sums(.)>1, .) %>% 
                    phyloseq::distance(., method="morisita")  #bray (abund) or jaccard (p/a) - [results the same so using jaccard]

            #permanova
                ps16 %>%
                    prune_samples(sample_data(.)$Sample %in% c("Tethys Sponge Microbiome", "Halisarca Sponge Microbiome"), .) %>% 
                    prune_taxa(taxa_sums(.)>1, .) %>%
                    sample_data(.) %>%
                    as("data.frame") %>%
                    as_tibble %>%
                    mutate_if(is.character, as_factor) %>%
                    mutate  (pH=factor(pH, order = TRUE,  levels = c("Control pH", "Medium pH", "Low pH")) ) %>%
                    mutate(ARMS=factor(ARMS)) %>%
                    as.data.frame %>%
                    #vegan::adonis2(dist ~   Sample/pH , data=., strata=.$Sample)
                    vegan::adonis2(MorSponge~ Sample + Site + pH + Sample:pH , data=., method="morisita",  permutations = how(plots = Plots(strata = .$Site), within = Within(), nperm=9999), by="terms")
                    #pairwiseAdonis::pairwise.adonis2(jacc~ Sample +pH , data=., method="jaccard", strata ='Site')

            #associated nmds
                ordSponge<-ps16 %>%
                    prune_samples(!sample_data(.)$Sample %in% c("Tethys Sponge Microbiome", "Halisarca Sponge Microbiome"), .) %>% 
                    prune_taxa(taxa_sums(.)>1, .) %>% 
                    ordinate(., method="NMDS", distance=MorSponge, trymax = 10000)                

                Figs[["NMDS_Sponge"]]<-ps16 %>%
                                            prune_samples(sample_data(.)$Sample %in% c("Tethys Sponge Microbiome", "Halisarca Sponge Microbiome"), .) %>% 
                                            prune_taxa(taxa_sums(.)>1, .) %>% 
                                            plot_ordination(., ordSponge, type="Sample", color="pH", shape="Sample") +
                                                geom_point(size=5)+
                                                ggtitle(paste0("Sponge Microbiome Composition\n(Stress=", signif(ordSponge$stress,3),")")) +
                                                pHcolScale +
                                                pHfillScale +
                                                scale_shape_manual(values = c(3, 4)) +
                                                theme(  strip.background = element_rect(fill = "white"),
                                                        legend.position="none",
                                                        #panel.border = element_rect(colour = "black", fill = "white"),
                                                        strip.text = element_text(size=14, color="black", face="bold"),
                                                        axis.text = element_blank(), 
                                                        axis.title = element_blank(), 
                                                        axis.ticks=element_blank(),
                                                        plot.title = element_text(size = 12, face = "bold", hjust=0.5),
                                                        panel.grid.major = element_blank(), 
                                                        panel.grid.minor = element_blank(),
                                                        panel.background = element_blank(), 
                                                        axis.line = element_line(colour = "black"))+
                                                guides(color="none")

        #algae gene

                Mor23<-ps23 %>%
                    phyloseq::distance(., method="morisita")  #bray (abund) or jaccard (p/a) - [results the same so using jaccard]

            #permanova
                ps23 %>%
                    sample_data(.) %>%
                    as("data.frame") %>%
                    as_tibble %>%
                    mutate_if(is.character, as_factor) %>%
                    mutate  (pH=factor(pH, order = TRUE,  levels = c("Control pH", "Medium pH", "Low pH")) ) %>%
                    mutate(ARMS=factor(ARMS)) %>%
                    as.data.frame %>%
                    #vegan::adonis2(dist ~   Sample/pH , data=., strata=.$Sample)
                    vegan::adonis2(Mor23~ Site + pH , data=., method="morisita", permutations = how(plots = Plots(strata = .$Site), within = Within(), nperm=9999), by="terms")
                    #pairwiseAdonis::pairwise.adonis2(jacc~ Sample +pH , data=., method="jaccard", strata ='Site')

            #associated nmds
                ord23<-ps23 %>%
                    ordinate(., method="NMDS", distance=Mor23, trymax = 10000)                

                Figs[["NMDS_Algae"]]<-ps23 %>%
                                            plot_ordination(., ord23, type="Sample", color="pH", shape="Sample") +
                                                geom_point(size=5)+
                                                ggtitle(paste0("Photosynthetic Community Composition\n(Stress=", signif(ord23$stress,3),")")) +
                                                pHcolScale +
                                                pHfillScale +
                                                scale_shape_manual(values = c(5)) +
                                                theme(  strip.background = element_rect(fill = "white"),
                                                        legend.position="none",
                                                        #panel.border = element_rect(colour = "black", fill = "white"),
                                                        strip.text = element_text(size=14, color="black", face="bold"),
                                                        axis.text = element_blank(), 
                                                        axis.title = element_blank(), 
                                                        axis.ticks=element_blank(),
                                                        plot.title = element_text(size = 12, face = "bold", hjust=0.5),
                                                        panel.grid.major = element_blank(), 
                                                        panel.grid.minor = element_blank(),
                                                        panel.background = element_blank(), 
                                                        axis.line = element_line(colour = "black"))+
                                                guides(color="none")


        #chemical whoel arms fractions
            #get distances

                Morwachem<-btab %>%
                    prune_samples(!sample_data(.)$Sample %in% c("Tethys Sponge Metabolome", "Halisarca Sponge Metabolome"), .) %>% 
                    prune_taxa(taxa_sums(.)>1, .) %>% 
                    phyloseq::distance(., method="morisita")  #bray (abund) or jaccard (p/a) - [results the same so using jaccard]

            #permanova
                btab %>%
                    prune_samples(!sample_data(.)$Sample %in% c("Tethys Sponge Metabolome", "Halisarca Sponge Metabolome"), .) %>% 
                    prune_taxa(taxa_sums(.)>1, .) %>%
                    sample_data(.) %>%
                    as("data.frame") %>%
                    as_tibble %>%
                    mutate_if(is.character, as_factor) %>%
                    mutate  (pH=factor(pH, order = TRUE,  levels = c("Control pH", "Medium pH", "Low pH")) ) %>%
                    mutate(ARMS=factor(ARMS)) %>%
                    as.data.frame %>%
                    #vegan::adonis2(dist ~   Sample/pH , data=., strata=.$Sample)
                    vegan::adonis2(Morwachem~ Sample + Site + pH + Sample:pH , data=., method="morisita", permutations = how(plots = Plots(strata = .$Site), within = Within(), nperm=9999), by="terms")
                    #pairwiseAdonis::pairwise.adonis2(jacc~ Sample +pH , data=., method="jaccard", strata ='Site')

            #associated nmds
                ordwachem<-btab %>%
                    prune_samples(!sample_data(.)$Sample %in% c("Tethys Sponge Metabolome", "Halisarca Sponge Metabolome"), .) %>% 
                    prune_taxa(taxa_sums(.)>1, .) %>% 
                    ordinate(., method="NMDS", distance=Morwachem, trymax = 10000)                

                Figs[["NMDS_ChemWholeARMS"]]<-btab %>%
                                            prune_samples(!sample_data(.)$Sample %in% c("Tethys Sponge Metabolome", "Halisarca Sponge Metabolome"), .) %>% 
                                            prune_taxa(taxa_sums(.)>1, .) %>% 
                                            plot_ordination(., ordwachem, type="Sample", color="pH", shape="Sample") +
                                                geom_point(size=5)+
                                                ggtitle(paste0("Community Metabolome Composition\n(Stress=", signif(ordwachem$stress,3),")")) +
                                                pHcolScale +
                                                pHfillScale +
                                                scale_shape_manual(labels=c('Environmental', 'Holobiont Community', "Photosynthetic Community"), values = c(0,1,2)) +
                                                scale_color_manual(labels=c('Control pH', 'Medium pH', "Low pH"), values = c("green", "yellow", "orange")) +
                                                theme(  strip.background = element_rect(fill = "white"),
                                                        #legend.position=c(0.875, 0.14),
                                                        legend.key.size = unit(0.5, 'cm'), #change legend key size
                                                        legend.title = element_blank(), #change legend title font size
                                                        legend.text = element_text(size=12, face="bold"), #change legend text font size
                                                        #panel.border = element_rect(colour = "black", fill = "white"),
                                                        strip.text = element_text(size=14, color="black", face="bold"),
                                                        axis.text = element_blank(), 
                                                        axis.title = element_blank(), 
                                                        axis.ticks=element_blank(),
                                                        plot.title = element_text(size = 12, face = "bold", hjust=0.5),
                                                        panel.grid.major = element_blank(), 
                                                        panel.grid.minor = element_blank(),
                                                        panel.background = element_blank(), 
                                                        axis.line = element_line(colour = "black"))


        #chemical sponge fractions
            #get distances

                MorSpongechem<-btab %>%
                    prune_samples(sample_data(.)$Sample %in% c("Tethys Sponge Metabolome", "Halisarca Sponge Metabolome"), .) %>% 
                    prune_taxa(taxa_sums(.)>1, .) %>% 
                    phyloseq::distance(., method="morisita")  #bray (abund) or jaccard (p/a) - [results the same so using jaccard]

            #permanova
                btab %>%
                    prune_samples(sample_data(.)$Sample %in% c("Tethys Sponge Metabolome", "Halisarca Sponge Metabolome"), .) %>% 
                    prune_taxa(taxa_sums(.)>1, .) %>%
                    sample_data(.) %>%
                    as("data.frame") %>%
                    as_tibble %>%
                    mutate_if(is.character, as_factor) %>%
                    mutate  (pH=factor(pH, order = TRUE,  levels = c("Control pH", "Medium pH", "Low pH")) ) %>%
                    mutate(ARMS=factor(ARMS)) %>%
                    as.data.frame %>%
                    #vegan::adonis2(dist ~   Sample/pH , data=., strata=.$Sample)
                    vegan::adonis2(MorSpongechem~ Sample + Site + pH + Sample:pH , data=., method="morisita",  permutations = how(plots = Plots(strata = .$Site), within = Within(), nperm=9999), by="terms")
                    #pairwiseAdonis::pairwise.adonis2(jacc~ Sample +pH , data=., method="jaccard", strata ='Site')

            #associated nmds
                ordSpongechem<-btab %>%
                    prune_samples(sample_data(.)$Sample %in% c("Tethys Sponge Metabolome", "Halisarca Sponge Metabolome"), .) %>% 
                    prune_taxa(taxa_sums(.)>1, .) %>% 
                    ordinate(., method="NMDS", distance=MorSpongechem, maxit=100,trymax = 10000)                

                Figs[["NMDS_ChemSponge"]]<-btab %>%
                                            prune_samples(sample_data(.)$Sample %in% c("Tethys Sponge Metabolome", "Halisarca Sponge Metabolome"), .) %>% 
                                            prune_taxa(taxa_sums(.)>1, .) %>% 
                                            plot_ordination(., ordSpongechem, type="Sample", color="pH", shape="Sample") +
                                                geom_point(size=5)+
                                                ggtitle(paste0("Sponge Metabolome Composition\n(Stress=", signif(ordSpongechem$stress,3),")")) +
                                                pHcolScale +
                                                pHfillScale +
                                                scale_shape_manual(values = c(3,4), labels=c('Halisarca Sp.', 'Tethya Sp.')) +
                                                theme(  strip.background = element_rect(fill = "white"),
                                                        #legend.position=c(0.875, 0.12),
                                                        legend.key.size = unit(0.5, 'cm'), #change legend key size
                                                        legend.title = element_blank(), #change legend title font size
                                                        legend.text = element_text(size=12, face="bold"), #change legend text font size
                                                        #panel.border = element_rect(colour = "black", fill = "white"),
                                                        strip.text = element_text(size=14, color="black", face="bold"),
                                                        axis.text = element_blank(), 
                                                        axis.ticks=element_blank(),
                                                        axis.title = element_blank(), 
                                                        plot.title = element_text(size = 12, face = "bold", hjust=0.5),
                                                        panel.grid.major = element_blank(), 
                                                        panel.grid.minor = element_blank(),
                                                        panel.background = element_blank(), 
                                                        axis.line = element_line(colour = "black"))+
                                                guides(color="none")                                                

        legend<- data.frame( #Sample=as.factor(c("Environmental Microbiome/Metabolome", "Holobiont Community Microbiome/Metabolome", "Photosynthetic Community Microbiome/Metabolome", "Halisarca Sponge Microbiome/Metabolome", "Tethya Sponge Microbiome/Metabolome",
                              #   "Photosynthetic Community Algae")),
                        pH=c("Control pH", "Medium pH", rep("Low pH", 4)),
                        x=rep(1,6), y=rep(1,6)        ) %>%
                #mutate  (Sample=factor(Sample,  
                 #           levels = c("Environmental Microbiome/Metabolome", "Holobiont Community Microbiome/Metabolome", "Photosynthetic Community Microbiome/Metabolome", "Halisarca Sponge Microbiome/Metabolome", "Tethya Sponge Microbiome/Metabolome",
                  #               "Photosynthetic Community Algae"))) %>%
                ggplot(aes(x=x, y=y, #shape=Sample, 
                            color=pH)) +
                    geom_point()+
                    scale_shape_manual(values = 0:5) +
                    pHcolScale +
                    lims(x = c(0,0), y = c(0,0))+
                    theme_void()+
                    theme(  legend.box="horizontal",
                            #legend.margin=margin()
                            legend.position = c(0.5,0.5),
                            #legend.key.size = unit(3, "cm"),
                            legend.text = element_text(size = 12, face = "bold"),
                            legend.title = element_text(size = 12, face = "bold")
                            ) + 
                    guides(colour = guide_legend(override.aes = list(size=10)))


            fullNmdsPlot<-egg::ggarrange(Figs[[1]], Figs[[4]], Figs[[2]], Figs[[5]], ncol=2, nrow=2)#, 
            #labels=(c("A", "", "B", "C", "D", "E")),
             #                               top="Figure S3: NMDS of Multiomic Composition for all Sample Types (Morisita dissimilarity)")
    #5.5 deseq
        #functions
            getPhylaDeseqTable<-function(ps) {
                ps_phylum<-ps %>% tax_glom(taxrank="Rank_3")
                sample_data(ps_phylum)$pH<- factor( sample_data(ps_phylum)$pH, ordered = FALSE )
                sample_data(ps_phylum)$OA<- fct_recode(sample_data(ps_phylum)$pH,"Control"="Control pH", "Acidified"="Medium pH", "Acidified"="Low pH")

                diagdds<-phyloseq_to_deseq2(ps_phylum, ~ Site+OA) %>%
                            DESeq(., test="Wald", fitType="parametric")

                    res = results(diagdds, cooksCutoff = FALSE)
                    alpha = 0.05
                    sigtab = res[which(res$padj < alpha), ]
                    sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps_phylum)[rownames(sigtab), ], "matrix")) %>% arrange(log2FoldChange)
                return(sigtab)
            }

            getFamilyDeseqTable<-function(ps) {
                ps_phylum<-ps %>% tax_glom(taxrank="Rank_7")
                sample_data(ps_phylum)$pH<- factor( sample_data(ps_phylum)$pH, ordered = FALSE )
                sample_data(ps_phylum)$OA<- fct_recode(sample_data(ps_phylum)$pH,"Control"="Control pH", "Acidified"="Medium pH", "Acidified"="Low pH")

                diagdds<-phyloseq_to_deseq2(ps_phylum, ~ Site+OA) %>%
                            DESeq(., test="Wald", fitType="parametric")

                    res = results(diagdds, cooksCutoff = FALSE)
                    alpha = 0.05
                    sigtab = res[which(res$padj < alpha), ]
                    sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps_phylum)[rownames(sigtab), ], "matrix")) %>% arrange(log2FoldChange)
                return(sigtab)
            }

            getPhylaDeseqTable_Sponge<-function(ps) {
                ps_phylum<-ps %>% tax_glom(taxrank="Rank_3")
                sample_data(ps_phylum)$pH<- factor( sample_data(ps_phylum)$pH, ordered = FALSE )
                sample_data(ps_phylum)$OA<- fct_recode(sample_data(ps_phylum)$pH,"Control"="Control pH", "Acidified"="Medium pH", "Acidified"="Low pH")

                diagdds<-phyloseq_to_deseq2(ps_phylum, ~ Sample+Site+OA) %>%
                            DESeq(., test="Wald", fitType="parametric")

                    res = results(diagdds, cooksCutoff = FALSE)
                    alpha = 0.05
                    sigtab = res[which(res$padj < alpha), ]
                    sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps_phylum)[rownames(sigtab), ], "matrix")) %>% arrange(log2FoldChange)
                return(sigtab)
            }

            getCompoundDeseqTable<-function(ps) {
                sample_data(ps)$pH<- factor( sample_data(ps)$pH, ordered = FALSE )
                sample_data(ps)$OA<- fct_recode(sample_data(ps)$pH,"Control"="Control pH", "Acidified"="Medium pH", "Acidified"="Low pH")

                diagdds<-phyloseq_to_deseq2(ps, ~ Site+OA) %>%
                            DESeq(., test="Wald", fitType="parametric")

                    res = results(diagdds, cooksCutoff = FALSE)
                    alpha = 0.05
                    sigtab = res[which(res$padj < alpha), ]
                    sigtab = cbind(as(sigtab, "data.frame"), tax_table(ps)[rownames(sigtab), ][,c(3,42)]  ) %>% arrange(log2FoldChange)
                    head(sigtab, 100)
                return(sigtab)
            }

            getCompoundDeseqTable_Sponge<-function(ps) {
                sample_data(ps)$pH<- factor( sample_data(ps)$pH, ordered = FALSE )
                sample_data(ps)$OA<- fct_recode(sample_data(ps)$pH,"Control"="Control pH", "Acidified"="Medium pH", "Acidified"="Low pH")

                diagdds<-phyloseq_to_deseq2(ps, ~ Sample+Site+OA) %>%
                            DESeq(., test="Wald", fitType="parametric")

                    res = results(diagdds, cooksCutoff = FALSE)
                    alpha = 0.05
                    sigtab = res[which(res$padj < alpha), ]
                    sigtab = cbind(as(sigtab, "data.frame"), tax_table(ps)[rownames(sigtab), ][,c(3,42)]  ) %>% arrange(log2FoldChange)
                return(sigtab)
            }
            #rank_names(ps)
            #tax_table(ps)[!is.na(tax_table(ps)[,42]),42]

        #make sigtabs
            #algae (23s)
                sigtab23<-getPhylaDeseqTable(ps23) %>% cbind(Sample="Photosynthetic Community Algae")
                sigtab23_family<-getFamilyDeseqTable(ps23) %>% cbind(Sample="Photosynthetic Community Algae")

                
                # find percentage of Ochrophyta which are phaeophyceae for text
                    ReadsPerOchEsv<-ps23 %>% 
                        subset_taxa(., Rank_3=="Ochrophyta") %>%
                        taxa_sums()



                    ps23 %>% 
                        subset_taxa(., Rank_3=="Ochrophyta") %>%
                        tax_table %>%
                        cbind(ReadsPerOchEsv) %>%
                        as_tibble %>%
                        mutate(ReadsPerOchEsv= as.integer(ReadsPerOchEsv)) %>%
                        group_by(Rank_4) %>%
                        summarise(ReadsPerOchGroup=sum(ReadsPerOchEsv)) %>%
                        mutate(ProportionReads=ReadsPerOchGroup/sum(ReadsPerOchGroup))
                    #99.7% of reads are phaeophyceae class (brown algae)
                
                    ps23 %>% 
                        subset_taxa(., Rank_3=="Ochrophyta") %>%
                        tax_table %>%
                        cbind(ReadsPerOchEsv) %>%
                        as_tibble %>%
                        mutate(ReadsPerOchEsv= as.integer(ReadsPerOchEsv)) %>%
                        group_by(Rank_7) %>%
                        summarise(ReadsPerOchGroup=sum(ReadsPerOchEsv)) %>%
                        mutate(ProportionReads=ReadsPerOchGroup/sum(ReadsPerOchGroup))
                    # 71.4% are Sargassum genera

            #16s communities and sponges
                sigtabEnv<- ps16 %>% prune_samples(sample_data(.)$Sample %in% c("Environmental Microbiome"), .) %>% prune_taxa(taxa_sums(.)>1, .) %>% getPhylaDeseqTable %>% cbind(Sample="Environmental Microbiome")
                sigtabHolo<- ps16 %>% prune_samples(sample_data(.)$Sample %in% c("Holobiont Community Microbiome"), .) %>% prune_taxa(taxa_sums(.)>1, .) %>% getPhylaDeseqTable %>% cbind(Sample="Holobiont Community Microbiome")
                sigtabPhoto<- ps16 %>% prune_samples(sample_data(.)$Sample %in% c("Photosynthetic Community Microbiome"), .) %>% prune_taxa(taxa_sums(.)>1, .) %>% getPhylaDeseqTable  %>% cbind(Sample="Photosynthetic Community Microbiome")
                sigtabSponge<- ps16 %>% prune_samples(sample_data(.)$Sample %in% c(c("Tethys Sponge Microbiome", "Halisarca Sponge Microbiome")), .) %>% prune_taxa(taxa_sums(.)>1, .) %>% getPhylaDeseqTable_Sponge %>% cbind(Sample="Sponge Microbiomes")

            #chemical commiuities and sponges
                sigtabChemEnv<- btab %>% prune_samples(sample_data(.)$Sample %in% c("Environmental Metabolome"), .) %>% prune_taxa(taxa_sums(.)>1, .) %>% getCompoundDeseqTable %>% cbind(Sample="Environmental Metabolome")
                sigtabChemHolo<- btab %>% prune_samples(sample_data(.)$Sample %in% c("Holobiont Community Metabolome"), .) %>% prune_taxa(taxa_sums(.)>1, .) %>% getCompoundDeseqTable %>% cbind(Sample="Holobiont Community Metabolome")
                sigtabChemPhoto<- btab %>% prune_samples(sample_data(.)$Sample %in% c("Photosynthetic Community Metabolome"), .) %>% prune_taxa(taxa_sums(.)>1, .) %>% getCompoundDeseqTable  %>% cbind(Sample="Photosynthetic Community Metabolome")
                sigtabChemSponge<- btab %>% prune_samples(sample_data(.)$Sample %in% c(c("Tethys Sponge Metabolome", "Halisarca Sponge Metabolome")), .) %>% prune_taxa(taxa_sums(.)>1, .) %>% getCompoundDeseqTable_Sponge %>% cbind(Sample="Sponge Metabolomes")

        #make plots
            bl <- colorRampPalette(c("navy", "royalblue", "lightskyblue"))(200)                      
            re <- colorRampPalette(c("mistyrose", "red2","darkred"))(200)

            AlgGeneChangePlot<-sigtab23 %>% 
                ggplot(., aes(x=Sample, y=Rank_3, fill= log2FoldChange)) + 
                    geom_tile()+
                    scale_fill_gradientn(name = expression(log[2]*"-"*fold~change), colours=c(bl, re), na.value = "grey"  ,limits = c(-30, 30)) +
                    ylab("Phylum")+
                    theme(axis.text = element_text(size=10),
                            axis.title.x=element_blank(),
                            legend.position="none",
                            panel.background=element_rect(fill="white", colour="black"))

            AlgGeneChangePlot<-sigtab23_family %>% 
                ggplot(., aes(x=Sample, y=Rank_7, fill= log2FoldChange)) + 
                    geom_tile(width=1) +
                    scale_y_discrete(position = "right", limits=rev)+
                    scale_x_discrete(labels=c("Benthic\nPhotosynthetic\nCommunity\nAlgae (23S)")) +
                    scale_fill_gradientn(name = expression(atop(log[2]*"-"*fold, change)), colours=c(bl, re), na.value = "grey"  ,limits = c(-30, 30)) +
                    facet_grid(Rank_3~., scales="free_y", space="free_y", switch="y") +
                    theme_bw()+
                    theme(  strip.text.y.left =element_text(angle=0, size=16, face="bold", color="black") ,
                            axis.text.y.right = element_text(size=12),
                            axis.text.x = element_text(size=16, face="bold", color="black"),
                            axis.title=element_blank(),
                            legend.title = element_text(size=14), 
                            legend.text = element_text(size=12),
                            legend.position = "none",
                            legend.background = element_rect(size=0.5, linetype="solid", colour ="black"),
                            strip.placement = "outside",
                            strip.background=element_rect(fill="white", colour="black"),
                            panel.grid=element_blank())
                            #panel.background=element_rect(fill="white", colour="white"),
                            #panel.border=element_rect(fill="white", colour="white"))



            BacGeneChangePlot<-rbind(sigtabEnv, sigtabHolo, sigtabPhoto, sigtabSponge) %>%
                mutate(Sample=fct_relevel(Sample,c( "Environmental Microbiome", "Holobiont Community Microbiome", "Photosynthetic Community Microbiome","Sponge Microbiomes") )) %>%
                ggplot(., aes(x=Sample, y=Rank_3, fill= log2FoldChange)) + 
                    geom_tile(width=1) +
                    scale_y_discrete(position = "left", limits=rev)+
                    scale_x_discrete(labels=c("Sediment\nMicrobiome", "Benthic\nHolobiont\nCommunity\nMicrobiome", "Benthic\nPhotosynthetic\nCommunity\nMicrobiome", "Sponge\nMicrobiome")) +
                    scale_fill_gradientn(name = expression(atop(log[2]*"-"*fold, change)), colours=c(bl, re), na.value = "grey"  ,limits = c(-30, 30)) +
                    #facet_grid(Rank_3~., scales="free_y", space="free_y", switch="y") +
                    theme_bw()+
                    theme(  #strip.text.y.left =element_text(angle=0, size=16, face="bold", color="black") ,
                            axis.text.y.left = element_text(angle=0, size=16, face="bold", color="black") ,
                            axis.text.x = element_text(size=16, face="bold", color="black"),
                            axis.title=element_blank(),
                            legend.title = element_text(size=14), 
                            legend.text = element_text(size=12),
                            legend.position = c(1.15, 0.5),
                            legend.background = element_rect(size=0.5, linetype="solid", colour ="black"),
                            strip.placement = "outside",
                            strip.background=element_rect(fill="white", colour="black"),
                            panel.grid=element_blank())

        #data written to csv and modified based on literature review to group molecules.
            rbind(sigtabChemEnv, sigtabChemHolo, sigtabChemPhoto, sigtabChemSponge) %>%
                '['(complete.cases(.),) %>%
                write.csv(file.path(path, "Outputs", "DeseqChems.csv"))

            ChemChangePlot<-read.csv(file.path(path, "Outputs", "DeseqChems_revised.csv")) %>%
                filter(class!="Contaminants") %>%
                #mutate(Simplified_Compound_Name=fct_reorder(Simplified_Compound_Name, class)) %>%
                #mutate(Code=as.numeric(as.factor(Simplified_Compound_Name)))
                mutate(Figure_Compound_Name=  str_pad( Figure_Compound_Name , width=2, pad=0))   %>% 
                ##mutate(Sample=fct_relevel(Sample,c("Photosynthetic Community Metabolome", "Environmental Metabolome", "Holobiont Community Metabolome", "Sponge Metabolomes") )) %>%
                #arrange((Code)) %>%
                ggplot(., aes(x=Sample, y=Figure_Compound_Name, fill= log2FoldChange)) + 
                    geom_tile(width=1) +
                    scale_y_discrete(position = "right", limits=rev)+
                    scale_x_discrete(labels=c("Sediment\nMetabolome", "Benthic\nHolobiont\nCommunity\nMetabolome", "Benthic\nPhotosynthetic\nCommunity\nMetabolome", "Sponge\nMetabolomes")) +
                    scale_fill_gradientn(name = expression(atop(log[2]*"-"*fold, change)), colours=c(bl, re), na.value = "grey"  ,limits = c(-30, 30)) +
                    facet_grid(class~., scales="free_y", space="free_y", switch="y") +
                    theme_bw()+
                    theme(  strip.text.y.left =element_text(angle=0, size=16, face="bold", color="black") ,
                            axis.text.y.right = element_text(size=12),
                            axis.text.x = element_text(size=16, face="bold", color="black"),
                            axis.title=element_blank(),
                            legend.title = element_text(size=14), 
                            legend.text = element_text(size=12),
                            legend.position = c(1.15, 0.5),
                            legend.background = element_rect(size=0.5, linetype="solid", colour ="black"),
                            strip.placement = "outside",
                            strip.background=element_rect(fill="white", colour="black"),
                            panel.grid=element_blank())
                            #panel.background=element_rect(fill="white", colour="white"),
                            #panel.border=element_rect(fill="white", colour="white"))




        geneDeseqPlot<-egg::ggarrange(AlgGeneChangePlot, BacGeneChangePlot, nrow=2, labels=(c("A", "B")))
                                   # top="Figure S4: Heatmaps of Significant ESV Differential Abundance for all Sample Types")

        chemDeseqPlot<-egg::ggarrange(ChemChangePlot,
                                    top="Figure 4: Heatmaps of Significant Metabolomic Differential Abundance for all Sample Types")








# 6 Non-analytical Figures
    #maps
        maps<-list()
        pngMAP_df =ggmap::get_stamenmap( bbox = c(left = 120, bottom = -30, right = 155, top = 0),
                        zoom = 5, maptype = "terrain-background")
        maps[[1]]<-ggmap::ggmap(pngMAP_df) +
                        annotate("rect", xmin=150.5, xmax=151, ymin=-10 , ymax=-9.5, color="red", fill="red")+ 
                        ggsn::scalebar( x.min = 135.9, x.max = 150,  y.min = -29, y.max = -10, dist = 500, dist_unit = "km",
                                        transform = TRUE, model = "WGS84") +
                        geom_text(label="Papua \n New Guinea", x=146, y=-4, size = 8) +
                        geom_text(label="Indonesia", x=130, y=-3, size = 8)+
                        geom_text(label="Timor Leste", x=129, y=-9, size = 8) +
                        geom_text(label="Australia", x=135, y=-22, size = 8)+

                        #ggtitle("Figure 1: Location of study sites, Dobu and Upa-Upasina")+
                        theme(                  axis.title.x = element_blank(),
                                                axis.title.y = element_blank(),
                                                axis.text= element_text(size = 12, hjust=0))


        #siteMAP_df = ggmap::get_map(location = c(150.5, -10, 151, -9.5), source = "google", zoom = 12)
        siteMAP_df <- ggmap::get_stamenmap( bbox = c(left = 150.5, bottom = -10, right = 151, top = -9.5),
                            zoom = 11, maptype = "terrain-background")
        maps[[2]]<-ggmap::ggmap(siteMAP_df)+ 
                        ggsn::scalebar( x.min = 150.6, x.max = 150.94,  y.min = -9.98, y.max = -9.5, dist = 10, dist_unit = "km",
                                    transform = TRUE, model = "WGS84")+
                        geom_point(x=150.878843, y=-9.746104, size=5) + 
                        geom_text(label="Dobu", x=150.918843, y=-9.746104, size = 8, face = "bold")+
                        geom_point(x=150.827345, y=-9.829434, size=5) +
                        geom_text(label="Upa-Upasina", x=150.735345, y=-9.829434, size = 8, face = "bold")+
                        theme(  plot.title =    element_text(size = 20, face = "bold", hjust=0),
                                                axis.title.x = element_blank(),
                                                axis.title.y = element_blank(),
                                                 axis.text= element_text(size = 12, hjust=0))

            

# 7. Make pdf
    #load figures made else where
    source(file.path(path, "Code", "DescriptiveFigure.R"))
    CommunityDistinctnessPlot<-readRDS(file=file.path(path, "Outputs", "CommunityDistinctnessPlot.RDS")) # MakePdf.R
    SpongeDistinctnessPlot<-readRDS(file=file.path(path, "Outputs", "SpongeDistinctnessPlot.RDS")) # MakePdf.R

save.image(file=file.path(path, 'plottingEnvironment.RData'))
load(file=file.path(path, 'plottingEnvironment.RData'))

    pdf(file = file.path(path,"Figs",paste0("PNGPaper_DataGeneratedFigures_",Sys.Date(),".pdf")), width=20, height =8 ) # The height of the plot in inches        

            plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
            text(5, 8, "Data-generated Figures for Manuscript")
            text(5, 7, "Ecosystem stress through the holobiome lens: decline of distinct holobionts under ocean acidification")
            text(5, 6, paste0("Jake Williams   ", Sys.Date()))

            plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
            text(5, 6, "Main Figures")

            egg::ggarrange(maps[[1]], maps[[2]], nrow=1, #labels=(c("A", "B")),  # from MakePdf.r
                                    top="Figure 1: Study localities, Dobu and Upa-Upasina")
           
            DescriptiveFigure       

            chemDeseqPlot

            CommunityDistinctnessPlot


            plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
            text(5, 6, "Supplementary Figures")


            fullRichnessPlot
            fullNmdsPlot
            geneDeseqPlot
            SpongeDistinctnessPlot

        dev.off()

###############

    jpeg(file=file.path(path, "map"), height = 8.3, width = 16, units = 'in', res = 300)
                egg::ggarrange(maps[[1]], maps[[2]], nrow=1)
    dev.off()

    jpeg(file=file.path(path, "distinctnessBoxplot"), height = 8.3, width = 11.7, units = 'in', res = 300)
        CommunityDistinctnessPlot
    dev.off()

    jpeg(file=file.path(path, "metaboliteChange.jpeg"), height = 8.3, width = 13, units = 'in', res = 300)
        ChemChangePlot
    dev.off()

    jpeg(file=file.path(path, "geneDeseqPlot.jpeg"), height = 8.3, width = 13, units = 'in', res = 300)
        geneDeseqPlot
    dev.off()
#supplementarys
    #richness
    jpeg(file=file.path(path, "FigS1.jpeg"), height = 8.3, width = 13, units = 'in', res = 300)
        fullRichnessPlot
    dev.off()

    jpeg(file=file.path(path, "FigS2.jpeg"), height = 8.3, width = 13, units = 'in', res = 300)
        fullNmdsPlot
    dev.off()

    jpeg(file=file.path(path, "FigS3.jpeg"), height = 8.3, width = 13, units = 'in', res = 300)
        BacGeneChangePlot
    dev.off()

    jpeg(file=file.path(path, "SpongeDistinctnessPlot.jpeg"), height = 8.3, width = 13, units = 'in', res = 300)
        SpongeDistinctnessPlot
    dev.off()