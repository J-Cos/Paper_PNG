library(tidyverse)
library(phyloseq)

dat<-read.table("Data/LPDATA_18S.txt", sep="\t", header=TRUE)

    names(dat)
    otutab<-dat %>% select(-OTUID, -Total, -taxonomy)
    rownames(otutab)<-dat$OTUID

    taxtab<-dat$taxonomy %>%
        str_split(";")
    taxtab<-as.data.frame(do.call(rbind, taxtab))
    taxtab<-tax_table(taxtab)
    rownames(taxtab)<-dat$OTUID

metadat<-read.table("Data/LPMetadata_18S.txt", sep="\t", header=TRUE)
    metadata<-metadat %>% select(-SampleID)
    rownames(metadata)<-metadat$SampleID

ps<-phyloseq(otu_table(otutab, taxa_are_rows=TRUE), tax_table(taxtab), sample_data(metadata))

    proportions_df<-tax_table(ps) %>%
        as.data.frame %>% 
        as_tibble %>%
        cbind(total=taxa_sums(ps)) %>%
        group_by(ta1, ta2, ta3, ta4) %>% 
        summarise(n=sum(total)) %>%
        ungroup %>%
        arrange(desc(n)) %>%
        mutate(proportion=100*n/sum(select(., n)))


proportions_df %>%
        group_by(ta1, ta2) %>% 
        summarise(proportion=sum(proportion)) %>%
        arrange(desc(proportion)) %>% 
        write.csv("Outputs/18sReadAbundanceSummary.csv")

