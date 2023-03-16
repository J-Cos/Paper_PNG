#differential abundance analysis

library("DESeq2")
library("ggplot2")
library(RColorBrewer)


#make color scales

    commonPhyla<-as.data.frame(tax_table(ps16))%>%
        group_by (Rank_3) %>%
        summarise(n=n()) %>%
        arrange(desc(n)) %>%
        head(9) %>%
        '['(,1) %>%
        unlist()

    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


    color = c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999')
    
    myColors<- c(color, "black")
    names(myColors) <- c(commonPhyla, "Other")
    colScale <- scale_colour_manual(name = "Phyla",values = myColors)
    fillScale <- scale_fill_manual(name = "Phyla",values = myColors)



MakePhDifferentialAbundancePlot<-function(ps, Sample, colour, xaxis, title, commonPhyla=NULL) {
            #seperate out single holobiont samples
                ps_func<-prune_samples( (sample_data(ps)$Sample==Sample) , ps)
                ps_func<- prune_taxa( taxa_sums(ps_func)>0, ps_func)

            #taxa agglomeration
            ps_glom = tax_glom(ps_func, colour)
            topotus = names(sort(taxa_sums(ps_glom), TRUE)[1:9])
            toptaxtab = cbind(tax_table(ps_func), GlommedRank = "Other")

            toptax<-as(tax_table(ps_glom)[topotus, colour], "character")
            
            #different behaviour if common phyla specified
            if (is.null(commonPhyla)) {
                toptaxindices<-tax_table(ps_func)[,colour] %in% toptax
            } else { 
                toptaxindices<-tax_table(ps_func)[,colour] %in% commonPhyla
            }

            toptaxtab[toptaxindices, "GlommedRank"]<-toptaxtab[toptaxindices, colour]

            tax_table(ps_func) <- (tax_table(toptaxtab))


        #plot



            diagdds<-phyloseq_to_deseq2(ps_func, ~ pH)
            diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

            res = results(diagdds, cooksCutoff = FALSE)
            alpha = 0.01
            sigtab = res[which(res$padj < alpha), ]
            sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps_func)[rownames(sigtab), ], "matrix"))

            theme_set(theme_bw())
            scale_fill_discrete <- function(palname = "Set1", ...) {
                scale_fill_brewer(palette = palname, ...)
            }
            

            # Phylum order
            x = tapply(sigtab$log2FoldChange, sigtab[["GlommedRank"]], function(x) max(x))
            x = sort(x, TRUE)
            sigtab$Colour = factor(as.character(sigtab[["GlommedRank"]]), levels=names(x))
            # Genus order
            x = tapply(sigtab$log2FoldChange, sigtab[[xaxis]], function(x) mean(x))
            x = sort(x, TRUE)
            sigtab$XAxis = factor(as.character(sigtab[[xaxis]]), levels=names(x))

            PropESVsIncreasing<-c(sum(sigtab$log2FoldChange > 0)/ length(sigtab$log2FoldChange),
                                    length(sigtab$log2FoldChange)
                                )


            plot<-ggplot(sigtab, aes(x=XAxis, y=log2FoldChange, color=Colour)) + geom_point(size=2) + 
            ggtitle(title) +
            theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, size=4))+
            colScale +
            fillScale

            return(list(PropESVsIncreasing, plot ))
    }


#sequences
   a<- MakePhDifferentialAbundancePlot(ps=ps16, Sample="Sessile", colour="Rank_3", xaxis="Rank_6", title="Holobiont", commonPhyla=commonPhyla)
   b<- MakePhDifferentialAbundancePlot(ps=ps16, Sample="100", colour="Rank_3", xaxis="Rank_6", title="Environmental", commonPhyla=commonPhyla)
   c<- MakePhDifferentialAbundancePlot(ps=ps16, Sample="Algae", colour="Rank_3", xaxis="Rank_6", title="Algae", commonPhyla=commonPhyla)
   d<- MakePhDifferentialAbundancePlot(ps=ps16, Sample="CS", colour="Rank_3", xaxis="Rank_6", title="Tethys", commonPhyla=commonPhyla)
   e<- MakePhDifferentialAbundancePlot(ps=ps16, Sample="RVS", colour="Rank_3", xaxis="Rank_6", title="Halisarca", commonPhyla=commonPhyla)
   
   f<- MakePhDifferentialAbundancePlot(ps=ps23, Sample="Algae_23s", colour="Rank_4", xaxis="Rank_4", title="Algae")

            grid.arrange(grobs = list(a[[2]], b[[2]], c[[2]], d[[2]], e[[2]]), ncol = 2,
                        top=textGrob("Sequences", gp=gpar(fontsize = 20, fontface = "bold")))


            data.frame(fraction=c("hol", "env", "alg", "teth", "hali"),
                        propESVincreaing=c(a[[1]][1], b[[1]][1], c[[1]][1], d[[1]][1], e[[1]][1]),
                        totalESVsChanging=c(a[[1]][2], b[[1]][2], c[[1]][2], d[[1]][2], e[[1]][2]),
                        totalIncreasing=c(a[[1]][1], b[[1]][1], c[[1]][1], d[[1]][1], e[[1]][1])*c(a[[1]][2], b[[1]][2], c[[1]][2], d[[1]][2], e[[1]][2])   ,
                        totalDecreasing= (1-c(a[[1]][1], b[[1]][1], c[[1]][1], d[[1]][1], e[[1]][1]))    *c(a[[1]][2], b[[1]][2], c[[1]][2], d[[1]][2], e[[1]][2])
                        )

#  compunds
   a<- MakePhDifferentialAbundancePlot(ps=btab, Sample="SESSILE", colour="class", xaxis="Compound_Name", title="Holo")
   b<- MakePhDifferentialAbundancePlot(ps=btab, Sample="100-500UM", colour="class", xaxis="Compound_Name", title="Env")
   c<- MakePhDifferentialAbundancePlot(ps=btab, Sample="Algae", colour="class", xaxis="Compound_Name", title="Alg")
   d<- MakePhDifferentialAbundancePlot(ps=btab, Sample="CS", colour="class", xaxis="Compound_Name", title="Tethys")
   e<- MakePhDifferentialAbundancePlot(ps=btab, Sample="RVS", colour="class", xaxis="Compound_Name", title="Halisarca")

            grid.arrange(grobs = list( g[[2]], h[[2]], i[[2]], j[[2]], k[[2]]), ncol = 5,
                        top=textGrob("Metabs", gp=gpar(fontsize = 20, fontface = "bold")))

            data.frame(fraction=c("hol", "env", "alg", "teth", "hali"),
                        propESVincreaing=c(a[[1]][1], b[[1]][1], c[[1]][1], d[[1]][1], e[[1]][1]),
                        totalESVsChanging=c(a[[1]][2], b[[1]][2], c[[1]][2], d[[1]][2], e[[1]][2]),
                        totalIncreasing=c(a[[1]][1], b[[1]][1], c[[1]][1], d[[1]][1], e[[1]][1])*c(a[[1]][2], b[[1]][2], c[[1]][2], d[[1]][2], e[[1]][2])   ,
                        totalDecreasing= (1-c(a[[1]][1], b[[1]][1], c[[1]][1], d[[1]][1], e[[1]][1]))    *c(a[[1]][2], b[[1]][2], c[[1]][2], d[[1]][2], e[[1]][2])
                        )