# 1. Packages and path to project directory ######
    set.seed(0.1)
    
    library(tidyverse)
    library(ggplot2)
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
        ps16<-readRDS(file=file.path(path, "Outputs", "ps16.RDS"))
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
            psmerge<-merge_phyloseq(ps16, ps23)
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
                mutate  (Sample=factor(Sample,  levels = c("Environmental Microbiome", "Holobiont Community Microbiome","Photosynthetic Community Microbiome", "Halisarca Sponge Microbiome", "Tethys Sponge Microbiome", "Photosynthetic Community Algae"))) %>%
                betta_random(formula = estimate ~ Sample*pH + (1|Site), ses = error, data=.)
            
            bettamod1$global

            cbind( bettamod1$table,
                    bettamod1$table[,1]+bettamod1$table[,2]*1.959, 
                    bettamod1$table[,1]-bettamod1$table[,2]*1.959)


            write.csv(signif(bettamod1$table,2), file=file.path(path, "Outputs", "MetabarcodingRichnessnModel.csv"))
        #shannon model
            dv23<-ps23 %>% tax_glom(taxrank="Rank_3") %>% DivNet::divnet(., formula = ~pH + Site)
            dv16<-ps16 %>% tax_glom(taxrank="Rank_3") %>% DivNet::divnet(., formula = ~Sample*pH + Site)

            shannons <- dv23$shannon %>% summary %>% rbind(dv16$shannon %>% summary )
            shannonmerge <- meta %>%
                left_join(shannons,
                            by = "sample_names") %>%
                mutate_if( is.character,as_factor) %>%
                mutate(ARMS=as.factor(ARMS))

            bettamod2<-shannonmerge %>%
                mutate  (pH=factor(pH, order = TRUE,  levels = c("Control pH", "Medium pH", "Low pH")) ) %>%
                mutate  (Sample=factor(Sample,  levels = c( "Environmental Microbiome", "Holobiont Community Microbiome", "Photosynthetic Community Microbiome", "Halisarca Sponge Microbiome", "Tethys Sponge Microbiome", "Photosynthetic Community Algae"))) %>%
                betta_random(formula = estimate ~ Sample*pH + (1|Site), ses = error, data=.)

            cbind( bettamod2$table,
                    bettamod2$table[,1]+bettamod2$table[,2]*1.959, 
                    bettamod2$table[,1]-bettamod2$table[,2]*1.959)
            bettamod2$global

            write.csv(signif(bettamod2$table,2), file=file.path(path, "Outputs", "MetabarcodingShannonModel.csv"))
            


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
    
        sgene<-shannonmerge %>% select(pH, Sample, estimate) %>% mutate(Type="Metabarcoding", Diversity="Shannon Diversity")%>% 
            mutate  (Sample=factor(Sample,  levels = c( "Environmental Microbiome", "Holobiont Community Microbiome", "Photosynthetic Community Microbiome", "Halisarca Sponge Microbiome", "Tethys Sponge Microbiome", "Photosynthetic Community Algae"))) %>%            
            ggplot() +
                geom_boxplot(aes(y=estimate, x=pH, fill=pH)) +
                pHfillScale +
                facet_grid(~Sample)+
                ggtitle("Sequence Shannon Diversity") +
                theme(  strip.background = element_rect(fill = "white"),
                        #panel.border = element_rect(colour = "black", fill = "white"),
                        legend.position="none",
                        strip.text = element_text(size=5, color="black", face="bold"),
                        axis.text.x = element_blank(), 
                        axis.title = element_blank(), 
                        plot.title = element_text(size = 20, face = "bold", hjust=0.5),
                        panel.grid.major = element_blank(), 
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(), 
                        axis.line = element_line(colour = "black"))
        #sponges
        schem<-chemRichnessMerge  %>% select(pH, Sample,Shannon) %>% mutate(Type="Metabolomics", Diversity="Shannon Diversity") %>% rename(Shannon="estimate")%>% 
            mutate  (Sample=factor(Sample,  levels = c("Environmental Metabolome", "Holobiont Community Metabolome", "Photosynthetic Community Metabolome", "Halisarca Sponge Metabolome", "Tethys Sponge Metabolome")))%>%
            ggplot() +
                geom_boxplot(aes(y=estimate, x=pH, fill=pH)) +
                pHfillScale +
                facet_grid(~Sample)+
                ggtitle("Metabolite Shannon Diversity") +
                theme(  strip.background = element_rect(fill = "white"),
                        #panel.border = element_rect(colour = "black", fill = "white"),
                        #legend.position="none",
                        strip.text = element_text(size=5, color="black", face="bold"),
                        axis.text.x = element_blank(), 
                        axis.title = element_blank(), 
                        plot.title = element_text(size = 20, face = "bold", hjust=0.5),
                        panel.grid.major = element_blank(), 
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(), 
                        axis.line = element_line(colour = "black"))
        rgene<-richnessmerge %>% select(pH, Sample, estimate) %>% mutate(Type="Metabarcoding", Diversity="Richness")%>% 
            mutate  (Sample=factor(Sample,  levels = c( "Environmental Microbiome", "Holobiont Community Microbiome", "Photosynthetic Community Microbiome", "Halisarca Sponge Microbiome", "Tethys Sponge Microbiome", "Photosynthetic Community Algae"))) %>%          
            ggplot() +
                geom_boxplot(aes(y=estimate, x=pH, fill=pH)) +
                pHfillScale +
                facet_grid(~Sample)+
                ggtitle("Sequence Estimated Richness") +
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
                        axis.line = element_line(colour = "black"))
        rchem <-chemRichnessMerge  %>% select(pH, Sample, Observed) %>% mutate(Type="Metabolomics", Diversity="Richness") %>% rename(Observed="estimate")%>% 
            mutate  (Sample=factor(Sample,  levels = c("Environmental Metabolome", "Holobiont Community Metabolome", "Photosynthetic Community Metabolome", "Halisarca Sponge Metabolome", "Tethys Sponge Metabolome")))%>%
            ggplot() +
                geom_boxplot(aes(y=estimate, x=pH, fill=pH)) +
                pHfillScale +
                facet_grid(~Sample)+
                ggtitle("Metabolite Estimated Richness") +
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
                        axis.line = element_line(colour = "black"))

        fullRichnessPlot<-egg::ggarrange(sgene, schem,rgene,rchem, widths = c(6,5), labels=c("A", "B", "C", "D"), top="Figure S1: Multiomic richness and Shannon Diversity for all Sample Types")
    #5.3 permanova
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
                                                ggtitle(paste0("Community Microbiome Composition (Stress=", signif(ord16$stress,3),")")) +
                                                pHcolScale +
                                                pHfillScale +
                                                scale_shape_manual(values = c(0, 1, 2)) +
                                                theme(  strip.background = element_rect(fill = "white"),
                                                        legend.position="none",
                                                        #panel.border = element_rect(colour = "black", fill = "white"),
                                                        strip.text = element_text(size=14, color="black", face="bold"),
                                                        axis.text = element_blank(), 
                                                        axis.title = element_blank(), 
                                                        plot.title = element_text(size = 10, face = "bold", hjust=0.5),
                                                        panel.grid.major = element_blank(), 
                                                        panel.grid.minor = element_blank(),
                                                        panel.background = element_blank(), 
                                                        axis.line = element_line(colour = "black"))
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
                                                ggtitle(paste0("Sponge Microbiome Composition (Stress=", signif(ordSponge$stress,3),")")) +
                                                pHcolScale +
                                                pHfillScale +
                                                scale_shape_manual(values = c(3, 4)) +
                                                theme(  strip.background = element_rect(fill = "white"),
                                                        legend.position="none",
                                                        #panel.border = element_rect(colour = "black", fill = "white"),
                                                        strip.text = element_text(size=14, color="black", face="bold"),
                                                        axis.text = element_blank(), 
                                                        axis.title = element_blank(), 
                                                        plot.title = element_text(size = 10, face = "bold", hjust=0.5),
                                                        panel.grid.major = element_blank(), 
                                                        panel.grid.minor = element_blank(),
                                                        panel.background = element_blank(), 
                                                        axis.line = element_line(colour = "black"))

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
                                                ggtitle(paste0("Photosynthetic Community Algae Composition (Stress=", signif(ord23$stress,3),")")) +
                                                pHcolScale +
                                                pHfillScale +
                                                scale_shape_manual(values = c(5)) +
                                                theme(  strip.background = element_rect(fill = "white"),
                                                        legend.position="none",
                                                        #panel.border = element_rect(colour = "black", fill = "white"),
                                                        strip.text = element_text(size=14, color="black", face="bold"),
                                                        axis.text = element_blank(), 
                                                        axis.title = element_blank(), 
                                                        plot.title = element_text(size = 10, face = "bold", hjust=0.5),
                                                        panel.grid.major = element_blank(), 
                                                        panel.grid.minor = element_blank(),
                                                        panel.background = element_blank(), 
                                                        axis.line = element_line(colour = "black"))


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
                                                ggtitle(paste0("Community Metabolome Composition (Stress=", signif(ordwachem$stress,3),")")) +
                                                pHcolScale +
                                                pHfillScale +
                                                scale_shape_manual(values = c(6,7,8)) +
                                                theme(  strip.background = element_rect(fill = "white"),
                                                        legend.position="none",
                                                        #panel.border = element_rect(colour = "black", fill = "white"),
                                                        strip.text = element_text(size=14, color="black", face="bold"),
                                                        axis.text = element_blank(), 
                                                        axis.title = element_blank(), 
                                                        plot.title = element_text(size = 10, face = "bold", hjust=0.5),
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
                                                ggtitle(paste0("Sponge Metabolome Composition (Stress=", signif(ordSpongechem$stress,3),")")) +
                                                pHcolScale +
                                                pHfillScale +
                                                scale_shape_manual(values = c(9,10)) +
                                                theme(  strip.background = element_rect(fill = "white"),
                                                        legend.position="none",
                                                        #panel.border = element_rect(colour = "black", fill = "white"),
                                                        strip.text = element_text(size=14, color="black", face="bold"),
                                                        axis.text = element_blank(), 
                                                        axis.title = element_blank(), 
                                                        plot.title = element_text(size = 10, face = "bold", hjust=0.5),
                                                        panel.grid.major = element_blank(), 
                                                        panel.grid.minor = element_blank(),
                                                        panel.background = element_blank(), 
                                                        axis.line = element_line(colour = "black"))

        legend<- data.frame( Sample=as.factor(c("Environmental Metabolome", "Holobiont Community Metabolome", "Photosynthetic Community Metabolome", "Halisarca Sponge Metabolome", "Tethys Sponge Metabolome",
                                "Environmental Microbiome", "Holobiont Community Microbiome", "Photosynthetic Community Microbiome", "Halisarca Sponge Microbiome", "Tethys Sponge Microbiome", "Photosynthetic Community Algae")),
                        pH=c("Control pH", "Medium pH", rep("Low pH", 9)),
                        x=rep(1,11), y=rep(1,11)        ) %>%
                mutate  (Sample=factor(Sample,  levels = c("Environmental Microbiome", "Holobiont Community Microbiome", "Photosynthetic Community Microbiome", "Halisarca Sponge Microbiome", "Tethys Sponge Microbiome", "Photosynthetic Community Algae",
                                                            "Environmental Metabolome", "Holobiont Community Metabolome", "Photosynthetic Community Metabolome", "Halisarca Sponge Metabolome", "Tethys Sponge Metabolome"))) %>%
                ggplot(aes(x=x, y=y, shape=Sample, color=pH)) +
                    geom_point()+
                    scale_shape_manual(values = 0:10) +
                    pHcolScale +
                    lims(x = c(0,0), y = c(0,0))+
                    theme_void()+
                    theme(  legend.box="horizontal",
                            #legend.margin=margin()
                            legend.position = c(0.5,0.5)
                            #legend.key.size = unit(1, "cm"),
                            #legend.text = element_text(size =  12),
                            #legend.title = element_text(size = 15, face = "bold")
                            )
                # guides( #colour = guide_legend(override.aes = list(size=8)),
                    #        fill=guide_legend(nrow=1,byrow=TRUE))

            fullNmdsPlot<-egg::ggarrange(Figs[[1]], Figs[[2]], Figs[[3]], Figs[[4]],Figs[[5]], legend, ncol=2, nrow=3, labels=(c("A", "B", "C", "D", "E", "")),
                                            top="Figure S2: NMDS of Multiomic Composition for all Sample Types (Morisita dissimilarity)")
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
            bl <- colorRampPalette(c("lightskyblue", "royalblue", "navy"))(200)                      
            re <- colorRampPalette(c("darkred", "red2", "mistyrose"))(200)

            AlgGeneChangePlot<-sigtab23 %>% 
                ggplot(., aes(x=Sample, y=Rank_3, fill= log2FoldChange)) + 
                    geom_tile()+
                    scale_fill_gradientn(colours=c(re, bl), na.value = "grey98"  ,
                    limits = c(-30, 30)) +
                    ylab("Phylum")+
                    theme(axis.text = element_text(size=5),
                            axis.title.x=element_blank(),
                            legend.position="none")

            BacGeneChangePlot<-rbind(sigtabEnv, sigtabHolo, sigtabPhoto, sigtabSponge) %>%
                ggplot(., aes(x=Sample, y=Rank_3, fill= log2FoldChange)) + 
                    geom_tile() +
                    scale_fill_gradientn(colours=c(re, bl), na.value = "grey98"  ,
                    limits = c(-30, 30)) +
                    ylab("Phylum")+
                    theme(axis.text = element_text(size=5),
                            axis.title.x=element_blank())


            rbind(sigtabChemEnv, sigtabChemHolo, sigtabChemPhoto, sigtabChemSponge) %>%
                '['(complete.cases(.),) %>%
                write.csv(file.path(path, "Outputs", "DeseqChems.csv"))

            ChemChangePlot<-read.csv(file.path(path, "Outputs", "DeseqChems_revised.csv")) %>%
                filter(class!="Contaminants") %>%
                ggplot(., aes(x=Sample, y=Compound_Name, fill= log2FoldChange)) + 
                    geom_tile() +
                    scale_fill_gradientn(colours=c(re, bl), na.value = "grey98"  ,
                    limits = c(-30, 30)) +
                    theme(  strip.text.y =element_text(angle=0) ,
                            axis.text = element_text(size=5),
                            axis.title.x=element_blank(),
                            legend.position="none")+
                    facet_grid(class~., scales="free_y", space="free_y")





        fullDeseqPlot<-egg::ggarrange(AlgGeneChangePlot, BacGeneChangePlot,ChemChangePlot, nrow=3, labels=(c("A", "B", "C")),
                                    top="Figure S3: Heatmaps of Significant Multiomic Differential Abundance for all Sample Types")







# 6 Non-analytical Figures
    #maps
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

            

# 7. Make pdf
    #load figures made else where
    DescriptiveFigure<-readRDS(file=file.path(path, "Outputs", "DescriptiveFigure.RDS")) #DescriptiveFigure.R
    CommunityDistinctnessPlot<-readRDS(file=file.path(path, "Outputs", "CommunityDistinctnessPlot.RDS")) # MakePdf.R
    SpongeDistinctnessPlot<-readRDS(file=file.path(path, "Outputs", "SpongeDistinctnessPlot.RDS")) # MakePdf.R



    pdf(file = file.path(path,"Figs",paste0("PNGPaper_DataGeneratedFigures_",Sys.Date(),".pdf")), width=20, height =8 ) # The height of the plot in inches        

            plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
            text(5, 8, "Data-generated Figures for Manuscript")
            text(5, 7, "Ecosystem stress through the holobiome lens: decline of distinct holobionts under ocean acidification")
            text(5, 6, paste0("Jake Williams   ", Sys.Date()))

            plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
            text(5, 6, "Main Figures")

            egg::ggarrange(maps[[1]], maps[[2]], nrow=1, labels=(c("A", "B")),  # from MakePdf.r
                                    top="Figure 1: Study localities, Dobu and Upa-Upasina")

            DescriptiveFigure       

            CommunityDistinctnessPlot


            plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
            text(5, 6, "Supplementary Figures")

            fullRichnessPlot
            fullNmdsPlot
            fullDeseqPlot
            SpongeDistinctnessPlot

        dev.off()


