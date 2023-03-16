    #dysbiosis?

            btab_sponges<-prune_samples( (sample_data(btab_p)$Sample=="RVS" | sample_data(btab_p)$Sample=="CS" )  , btab_p)
            btab_sponges <-prune_taxa(taxa_sums(btab_sponges)>0, btab_sponges)
            
            levels(sample_data(btab_sponges)$pH)<-c("Control pH", "Medium pH", "Low pH")

            ps_sponge_distances<-phyloseq::distance(btab_sponges, method="bray") 

            ord_sponges <- ordinate(btab_sponges, "NMDS", "unifrac", weighted=FALSE, distance=ps_sponge_distances, trymax = 50)
            Figs[["Sponge_Metab_Ordination"]]<-plot_ordination(btab_sponges, ord_sponges, type="Sample", color="pH", shape="Sample") +
                                geom_point(size=5)+
                                ggtitle("Supplementary Figure 5: NMDS showing increasing similarity of holobiont to environmental fraction with OA (Stress=0.14)") +
                                stat_ellipse(aes(group=pH)) +
                                facet_wrap(~Sample)
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

            Figs[["Sponge_Metab_Ordination"]]



            ps16_sponges<-prune_samples( (sample_data(ps16ALLp)$Sample=="RVS" | sample_data(ps16ALLp)$Sample=="CS" )  , ps16ALLp)
            ps16_sponges <-prune_taxa(taxa_sums(ps16_sponges)>0, ps16_sponges)
            
            levels(sample_data(ps16_sponges)$pH)<-c("Control pH", "Medium pH", "Low pH")

            ps16_sponge_distances<-phyloseq::distance(ps16_sponges, method="bray") 

            ord_16sponges <- ordinate(ps16_sponges, "NMDS", "bray", weighted=FALSE, distance=ps16_sponge_distances, trymax = 50)
            Figs[["Sponge_Seq_Ordination"]]<-plot_ordination(ps16_sponges, ord_16sponges, type="Sample", color="pH", shape="Sample") +
                                geom_point(size=5)+
                                ggtitle("Supplementary Figure 5: NMDS showing increasing similarity of holobiont to environmental fraction with OA (Stress=0.14)") +
                                stat_ellipse(aes(group=pH)) +
                                facet_wrap(~Sample)
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

            Figs[["Sponge_Seq_Ordination"]]
            
            df_list<-list()


            for (sponge in c("RVS", "CS")) {

                ps16_sponges<-prune_samples( (sample_data(ps16ALLp)$Sample==sponge )  , ps16ALLp)
                levels(sample_data(ps16_sponges)$pH)<-c("Control pH", "Medium pH", "Low pH")


                Con_sponges<-prune_samples(sample_data(ps16_sponges)$pH=="Control pH" , ps16_sponges)
                Con_sponges <-prune_taxa(taxa_sums(Con_sponges)>0, Con_sponges)
                ConDistances<-as.vector(phyloseq::distance(Con_sponges, method="bray") )

                Med_sponges<-prune_samples(sample_data(ps16_sponges)$pH=="Medium pH" , ps16_sponges)
                Med_sponges <-prune_taxa(taxa_sums(Med_sponges)>0, Med_sponges)
                MedDistances<-as.vector(phyloseq::distance(Med_sponges, method="bray") )

                Low_sponges<-prune_samples(sample_data(ps16_sponges)$pH=="Low pH" , ps16_sponges)
                Low_sponges <-prune_taxa(taxa_sums(Low_sponges)>0, Low_sponges)       
                LowDistances<-as.vector(phyloseq::distance(Low_sponges, method="bray") )

                df_list[[sponge]]<-rbind(data.frame(Distances=ConDistances,
                            pH="Control pH",
                            Sponge=sponge),
                        data.frame(Distances=MedDistances,
                            pH="Medium pH",
                            Sponge=sponge),
                        data.frame(Distances=ConDistances,
                            pH="Low pH",
                            Sponge=sponge)
                )
            }

            df<-rbind(df_list[[1]], df_list[[2]])

            #run GLMM
                df$pH<-ordered(df$pH, levels = c("Control pH", "Medium pH", "Low pH"))


                Stats[[paste0("pH")]]<-glmer(data = df, formula=Distances ~ pH + (1|Sponge) )
                Stats[[paste0("nopH")]]<-lmer(data = df,  formula=Distances ~ (1|Sponge) )

                plot(density(df$Distances, na.rm=TRUE))

                qqnorm(residuals(Stats[["pH"]]))
                scatter.smooth(residuals(Stats[["pH"]]) ~ fitted(Stats[["pH"]]))
                qqnorm(residuals(Stats[["nopH"]]))
                scatter.smooth(residuals(Stats[["nopH"]]) ~ fitted(Stats[["nopH"]]))
                #assumptions look good with poisson distribution (makes sense given count data)

                anova(Stats[["pH"]], Stats[["nopH"]])
                    #no significant efefct of pH:sample
                confs<-confint(Stats[["pH"]])
                    # pH confints  cross zero for all samples
                summary(Stats[["pH"]])

                r.squaredGLMM(Stats[["pH"]])           


















                