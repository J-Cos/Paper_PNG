# 1. Packages and path to project directory ######
    set.seed(0.1)
    
    library(phyloseq)
    library(tidyverse)
    library(ggplot2)
    library(treemapify)
    library(viridis)


    library(breakaway)
    
    #params
        path<-"/home/j/Dropbox/PNG_Paper" #provide the path to the project folder
        theme_set(theme_linedraw())
        minText<-1 #min size of text in treemaps

    #1.1 load functions
    #funHalisarca Sponge Microbiome
    function_files<-list.files(file.path(path, "Functions"))
    sapply(file.path(path, "Functions",function_files),source)
    
# 2. Load Data
    #16s
    ps16<-readRDS(file=file.path(path, "Outputs", "ps16.RDS")) %>% prune_taxa( as.data.frame(tax_table(.) )$Rank_5!="Chloroplast", .)
    #23s
    ps23<-readRDS(file=file.path(path, "Outputs", "ps23.RDS"))

#2) Mkae visuals of estimate slpha diversity
                    #ensure all samples have x seqs and all taxa have >0 (all 23s already do)
                        ps16<-prune_samples( sample_sums(ps16)>100000, ps16)
                        ps16<-prune_samples(sample_data(ps16)$pH %in% c("Control pH", "Medium pH", "Low pH"), ps16)
                        ps23<-prune_samples( sample_sums(ps23)>100000, ps23)

                        #ps16<-prune_samples( sample_data(ps16)$ARMS %in% sample_data(ps23)$ARMS, ps16)
                        #ps23<-prune_samples( sample_data(ps23)$ARMS %in% sample_data(ps16)$ARMS, ps23)
                        ps16<- prune_taxa( taxa_sums(ps16)>0, ps16)
                        ps23<- prune_taxa( taxa_sums(ps23)>0, ps23)

#Main ARMS richness
            completeARMS<-sample_data(ps16) %>% 
                        as_tibble %>%
                        #mutate_if(is.character, as.factor) %>%
                        filter(!Sample %in% c("Tethys Sponge Microbiome", "Halisarca Sponge Microbiome")) %>%
                        group_by(ARMS, pH, Site) %>%
                        summarise(n=n()) %>%
                        filter(n>2) 
            ps16<-prune_samples(sample_data(ps16)$ARMS %in% completeARMS$ARMS, ps16)
            ps16<- prune_taxa( taxa_sums(ps16)>0, ps16)

            ps16m<-merge_samples(ps16, sample_data(ps16)$ARMS)
            sample_data(ps16m)$Site<-sample_data(completeARMS)$Site
            sample_data(ps16m)$pH<-sample_data(completeARMS)$pH

        GetBreakawayRichness<-function(ps){
            richness <- ps %>% breakaway
            meta <- ps %>%
            sample_data %>%
            as_tibble %>%
            mutate("sample_names" = ps %>% sample_names )

            combined_richness <- meta %>%
            left_join(summary(richness),
                        by = "sample_names") %>%
            mutate_if( is.character,as_factor)

            return(combined_richness)

        }

        Richness16<-GetBreakawayRichness(ps16m)

        boxplotsARMS<-Richness16 %>%
                    mutate  (pH=factor(pH, order = TRUE,  levels = c("Low pH", "Medium pH", "Control pH")) )  %>%
                    ggplot(aes(x=estimate, y=pH)) +
                        geom_violin(fill="black") +
                        geom_boxplot(width=0.1, color="grey", fill="grey", alpha=0.5) +
                        expand_limits(x = c(0, 20000)) +
                        scale_x_continuous( breaks = c(0, 10000, 20000), expand=c(0,0)) +
                    theme(  strip.background = element_rect(fill = "white"),
                            #panel.border = element_blank(),
                            strip.text = element_text( color="black", face="bold"),
                            axis.text.x = element_text(face="bold", size=8, angle=0, hjust=0.5),
                            axis.text.y = element_blank(),
                            axis.title.x = element_text(face="bold", size=8, angle=0, hjust=0.5),
                            axis.title.y = element_blank(),
                            plot.title = element_text(size = 20, face = "bold", hjust=0.5),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.background = element_blank(), 
                            axis.line = element_line(colour = "black"))+
                    labs(x = "Whole ARMS estimated \n bacterial ESV richness")

#23s richness
        Richness23<-GetBreakawayRichness(ps23)
        boxplotsAlg<-Richness23 %>%
                    mutate  (pH=factor(pH, order = TRUE,  levels = c("Low pH", "Medium pH", "Control pH")) )  %>%
                    ggplot(aes(x=estimate, y=pH)) +
                        geom_violin(fill="black") +
                        geom_boxplot(width=0.1, color="grey", fill="grey", alpha=0.5) +
                        expand_limits(x = c(0, 3000)) +
                        scale_x_continuous( breaks = c(0, 1500, 3000), expand=c(0,0)) +
                    theme(  strip.background = element_rect(fill = "white"),
                            #panel.border = element_blank(),
                            strip.text = element_text(color="black", face="bold"),
                            axis.text.x = element_text(face="bold", size=8, angle=0, hjust=0.5),
                            axis.text.y = element_blank(),
                            axis.title.x = element_text(face="bold", size=8, angle=0, hjust=0.5),
                            axis.title.y = element_blank(),
                            plot.title = element_text(size = 20, face = "bold", hjust=0.5),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.background = element_blank(), 
                            axis.line = element_line(colour = "black"))+
                    labs(x = "Photosynthetic community estimated \n algae (23s) sequence richness")


#sponge richness
ps16sponges<-prune_samples(sample_data(ps16)$Sample %in% c("Tethys Sponge Microbiome", "Halisarca Sponge Microbiome"), ps16)
            ps16sponges<- prune_taxa( taxa_sums(ps16sponges)>0, ps16sponges)

        RichnessSponges<-GetBreakawayRichness(ps16sponges)
        boxplotsSponge<-RichnessSponges %>%
                    mutate  (pH=factor(pH, order = TRUE,  levels = c("Low pH", "Medium pH", "Control pH")) )  %>%
                    ggplot(aes(x=estimate, y=pH)) +
                        geom_violin(fill="black") +
                        geom_boxplot(width=0.1, color="grey", fill="grey", alpha=0.5) +
                        expand_limits(x = c(0, 5000)) +
                        scale_x_continuous( breaks = c(0, 2500, 5000), expand=c(0,0)) +
                        scale_y_discrete(labels = c('Low pH','Medium pH','Control pH'), position = "right") +
                    theme(  strip.background = element_rect(fill = "white"),
                            strip.text = element_text(color="black", face="bold"),
                            axis.text.y = element_text(face="bold", size=12, angle=0, hjust=0.5),
                            axis.text.x = element_text(face="bold", size=8, angle=0, hjust=0.5),
                            axis.title.x = element_text(face="bold", size=8, angle=0, hjust=0.5),
                            axis.title.y = element_blank(),
                            plot.title = element_text(size = 20, face = "bold", hjust=0.5),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.background = element_blank(), 
                            axis.line = element_line(colour = "black"))+
                    labs(x = "Sponge estimated \n bacterial ESV richness")

      
    #testing signif        

    Richness16 %>%
        #mutate  (pH=factor(pH, order = TRUE,  levels = c("Control pH", "Medium pH", "Low pH")) ) %>%
        #summarise(estimate=estimate, cumul_score=cumul_score, Site_Longitude=Site_Longitude, Site_Latitude=Site_Latitude, Year=Year, Island=Island, Region=Region , Fraction=Fraction,FIY=as.factor(paste0(Fraction)), error=error) %>%
        betta_random(formula = estimate ~ pH *Sample + (1|Site), ses = error, data=.) %>%
        '['("table")


# 1) Make visual of idealised ARMS plates
                    #ensure all samples have x seqs and all taxa have >0 (all 23s already do)
                        ps16<-prune_samples( sample_sums(ps16)>100000, ps16)
                        ps23<-prune_samples( sample_sums(ps23)>100000, ps23)

                        #ps16<-prune_samples( sample_data(ps16)$ARMS %in% sample_data(ps23)$ARMS, ps16)
                        #ps23<-prune_samples( sample_data(ps23)$ARMS %in% sample_data(ps16)$ARMS, ps23)
                        ps16<- prune_taxa( taxa_sums(ps16)>0, ps16)
                        ps23<- prune_taxa( taxa_sums(ps23)>0, ps23)

        ps23<-rarefy_even_depth (ps23, rngseed=1)
        ps16<-rarefy_even_depth (ps16, rngseed=1)

        GetGroupedSampleData<-function(ps, multisample=FALSE){
            sdm<-sample_data(ps) %>% 
                as_tibble %>%
                mutate_if(is.character, as.factor) %>%
                #group_by(pH, Site, Sample) %>%
                group_by(pH, Sample) %>%
                summarise(n=n()) %>%
                sample_data
            if(multisample) {
                #rownames(sdm)<-paste0(sdm$Site, sdm$pH, sdm$Sample)
                rownames(sdm)<-paste0(sdm$pH, sdm$Sample)
            } else {
                #rownames(sdm)<-paste0(sdm$Site, sdm$pH)
                rownames(sdm)<-paste0(sdm$pH)
            }
            return(sdm)
        }

        ps23m<-merge_samples(ps23, as.character( sample_data(ps23)$pH))
        sample_data(ps23m)<-GetGroupedSampleData(ps23)

        #ps16m<-merge_samples(ps16, as.factor(paste0(sample_data(ps16)$Site, sample_data(ps16)$pH, sample_data(ps16)$Sample)))
        ps16m<-merge_samples(ps16, as.factor(paste0(sample_data(ps16)$pH, sample_data(ps16)$Sample)))
        sample_data(ps16m)<-GetGroupedSampleData(ps16, multisample=TRUE)

        ps23m<-tax_glom(ps23m, "Rank_4")
        ps16m<-tax_glom(ps16m, "Rank_4")

        ps16m_ra<-transform_sample_counts(ps16m, function(x) x/sum(x))
        ps23m_ra<-transform_sample_counts(ps23m, function(x) x/sum(x))

        #p<-plot_bar(ps16m_ra, "Rank_3", fill="Rank_2", facet_grid=pH+Site~sample_Sample)
        #p23<-plot_bar(ps23m_ra, "Rank_3", fill="Rank_2", facet_grid=pH~Site)
        p<-plot_bar(ps16m_ra, "Rank_4", fill="Rank_4", facet_grid=pH~sample_Sample)
        p23<-plot_bar(ps23m_ra, "Rank_4", fill="Rank_4", facet_grid=~pH)


    treemaps16s<-    p$data %>%
            rbind(p23$data) %>% 
            mutate(pH = fct_relevel(pH, c("Control pH", "Medium pH", "Low pH"))) %>%
            filter(sample_Sample %in% c("Environmental Microbiome", "Holobiont Community Microbiome", "Photosynthetic Community Microbiome")) %>%
            mutate(sample_Sample=fct_relevel(sample_Sample,c("Photosynthetic Community Microbiome", "Environmental Microbiome", "Holobiont Community Microbiome") )) %>%
            ggplot(aes(area = Abundance, fill = Rank_3,
                    label = Rank_4, subgroup=Rank_3)) +
                geom_treemap() +
                geom_treemap_text(alpha = 0.25, colour = "black", place = "centre",
                                size = 15, grow = TRUE, min.size=minText) +
                facet_grid(pH~sample_Sample,  labeller = labeller(sample_Sample = 
                                                    c( "Environmental Microbiome"="Environmental\nMicrobiome\n", 
                                                    "Holobiont Community Microbiome"= "Holobiont\nCommunity\nMicrobiome", 
                                                    "Photosynthetic Community Microbiome"="Photosynthetic\nCommunity\nMicrobiome")), 
                                                    switch="y") +
                theme(legend.position="none")+
                geom_treemap_subgroup_text(place = "centre", grow = TRUE,
                                    colour = "white",
                                    fontface = "italic", angle=45, min.size=minText )+
                scale_fill_viridis(option="turbo", discrete=TRUE)+
                theme(  strip.background = element_blank(),
                        strip.text.y.left =  element_text(color="black", face="bold", size=20),
                        strip.text.x = element_text(color="black", face="bold", size=20))

    treemaps23s<-    p$data %>%
            rbind(p23$data) %>% 
            mutate(pH = fct_relevel(pH, c("Control pH", "Medium pH", "Low pH"))) %>%
            filter(sample_Sample %in% c("Photosynthetic Community Algae")) %>%
            ggplot(aes(area = Abundance, fill = Rank_3,
                    label = Rank_4, subgroup=Rank_3)) +
                geom_treemap() +
                geom_treemap_text(alpha = 0.25, colour = "black", place = "centre",
                                size = 15, grow = TRUE, min.size=minText) +
                facet_grid(pH~sample_Sample, labeller = labeller(sample_Sample = 
                                                    c( "Photosynthetic Community Algae"="Photosynthetic\nCommunity\n")),
                                                    switch="y") +
                theme(legend.position="none")+
                geom_treemap_subgroup_text(place = "centre", grow = TRUE,
                                    colour = "white",
                                    fontface = "italic", angle=45, min.size=minText ) +
                scale_fill_viridis(discrete=TRUE) +
                theme(  strip.background = element_blank(),
                        strip.text.y.left =  element_text(color="black", face="bold", size=20),
                        strip.text.x = element_text(color="black", face="bold", size=20))

    treemapsSponges<-    p$data %>%
            rbind(p23$data) %>% 
            mutate(pH = fct_relevel(pH, c("Control pH", "Medium pH", "Low pH"))) %>%
            filter(sample_Sample %in% c("Tethys Sponge Microbiome", "Halisarca Sponge Microbiome")) %>%
            ggplot(aes(area = Abundance, fill = Rank_3,
                    label = Rank_4, subgroup=Rank_3)) +
                geom_treemap() +
                geom_treemap_text(alpha = 0.25, colour = "black", place = "centre",
                                size = 15, grow = TRUE, min.size=minText) +
                facet_grid(pH~sample_Sample, labeller = labeller(sample_Sample = 
                                                    c(  "Tethys Sponge Microbiome" = "Tethya\nSponge\nMicrobiome",
                                                        "Halisarca Sponge Microbiome" = "Halisarca\nSponge\nMicrobiome"))) +
                theme(legend.position="none")+
                geom_treemap_subgroup_text(place = "centre", grow = TRUE,
                                    colour = "white",
                                    fontface = "italic", angle=45, min.size=minText ) +
                scale_fill_viridis(option="turbo", discrete=TRUE)+
                theme(  strip.background = element_blank(),
                        strip.text.y =  element_blank(),
                        strip.text.x = element_text( color="black", face="bold", size=20))


#arranging!!
#DescriptiveFigure<-egg::ggarrange(treemaps16s, boxplotsARMS,treemaps23s,boxplotsAlg, treemapsSponges,boxplotsSponge, widths = c(3,1, 1,1, 2,1),
#                                    top="Figure 3: Descriptive Figure of Metabarcoding Fraction Compositions and Estimated Richnesses")
DescriptiveFigure<-egg::ggarrange(#treemaps23s, 
                                    treemaps16s, treemapsSponges, widths = c(3, 2),
                                    #top="Figure 3: Descriptive Figure of Metabarcoding Fraction Compositions", 
                                    labels=c())
#theme(plot.margin = margin(0.1,0.1,2,0.1, "cm")) 
    jpeg(file=file.path(path, "Fig3.jpeg"), height = 8.3, width = 13, units = 'in', res = 300)
        DescriptiveFigure
    dev.off()

pdf(file = file.path(path,"Figs",paste0("DescriptiveFigure",Sys.Date(),".pdf")), width=20, height =8 ) # The height of the plot in inches        


plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
text(5, 8, "Description of an ARMS Data Set: Papua New Guinea - Ocean Acidification ")
text(5, 7, "Data is aggregated across biological replicates to present the community composition of 'representative ARMS plates' (and analgous communities where a fraction 
            does not map exactly onto a single ARMS plate). A representative ARMS plate is presented for each treatment level, showing Phyla (in white text) and Class (in grey text).")
            #Estimated community richness (Willis, Bunge and Whitman, 2017) is presented for each treatment level, with ARMS plates aggregated to produce richness estimates
            #for whole ARMS.")
text(5, 6, paste0("Jake Williams   ", Sys.Date()))
DescriptiveFigure

dev.off()

