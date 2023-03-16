##################
#Run Sloan Mixture Model on PNG 16s Data
##################


#1. dependencies
    set.seed(0.1)

    #CRAN mirror
    r = getOption("repos")
    r["CRAN"] = "http://cran.us.r-project.org"
    options(repos = r)

    #packages
    if (!require("pacman")) install.packages("pacman")
    pacman::p_load( DECIPHER,
                    Biostrings,
                    tidyverse,
                    ggplot2,
                    gridExtra,
                    phyloseq)

#2. parameters

    path <- "/home/j/Dropbox"
    datapath<- "BioinformaticPipeline_Env/Results/PNGFullTest_16s/PNGFullTest_16s_SeqDataTable.RDS"
    outputpath<-""

    k_param<-0.1


#3. load data and functions
    #funcs
    function_files<-list.files(file.path("Code/Functions"))
    sapply(file.path("Code/Functions",function_files),source)


#4. load data               
                                       
    ps16<-readRDS(file=file.path("Outputs", "ps16.RDS"))
    ps<-rarefy_even_depth(ps16, replace=FALSE)

    ps_condition_list<-list()
                    pHs<-c("Control pH", "Medium pH", "Low pH")
                    samples<-c("Holobiont Community Microbiome", "Environmental Microbiome")
                        for (pH in pHs) {
                            for (sample in samples) {
                                ps_temp<-try(prune_samples(sample_data(ps)$pH==pH, ps))
                                ps_temp<-try(prune_samples(sample_data(ps_temp)$Sample==sample, ps_temp))
                                ps_condition_list[[ paste0(pH) ]][[paste0(sample)]] <-try(prune_taxa(taxa_sums(ps_temp)>0, ps_temp))

                            }
                        }

#run sloan mixture model
m<-list()
for (i in names(ps_condition_list)) {
    spp<-t(otu_table(ps_condition_list[[i]][["Holobiont Community Microbiome"]]))
    pool<-t(otu_table(ps_condition_list[[1]][["Environmental Microbiome"]]))

    m[[i]]<-sncm.fit.new(spp=spp, pool=pool, stats=NULL, taxon=NULL)
}

# mixing weighted migration for the three pH levels
df<-data.frame(
    condition=names(ps_condition_list),
    m.holo=c(m[[1]][[1]]$m.holo*m[[1]][[1]]$mix, m[[2]][[1]]$m.holo*m[[2]][[1]]$mix, m[[3]][[1]]$m.holo*m[[3]][[1]]$mix),
    m.env=c(m[[1]][[1]]$m.env*(1-m[[1]][[1]]$mix), m[[2]][[1]]$m.env*(1-m[[2]][[1]]$mix), m[[3]][[1]]$m.env*(1-m[[3]][[1]]$mix)),
    Rsqr=c(m[[1]][[1]]$Rsqr, m[[2]][[1]]$Rsqr, m[[3]][[1]]$Rsqr)
)
ggplot(df, aes(x=condition))+
    #geom_line(aes(y=m.holo), color="red")+
    #geom_line(aes(y=m.env), color="blue")+
    geom_point(aes(y=m.holo/m.env), color="green")+
    geom_line(aes(y=Rsqr))

