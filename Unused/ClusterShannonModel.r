#cluster shannon model?
setwd("/rds/general/user/jcw120/home/BioinformaticPipeline_Env") #necessary as it appears different job classes have different WDs.

#CRAN mirror
    r = getOption("repos")
    r["CRAN"] = "http://cran.us.r-project.org"
    options(repos = r)

# 1. Packages and path to project directory ######
#install.packages("remotes")   ## run this line if you do not already have remotes installed
#remotes::install_github("adw96/breakaway")
#remotes::install_github("adw96/DivNet")

    set.seed(0.1)
    
    library(tidyverse)
    library(phyloseq)
    library(DivNet)


    #params
        path<-path <-"../OtherAnalyses"  #provide the path to the project folder


# 2. Data prep
    # get data
        ps16<-readRDS(file=file.path(path, "Data", "PNG","ps16.RDS"))
        ps23<-readRDS(file=file.path(path, "Data", "PNG", "ps23.RDS"))

    #models
        #genetic
            psmerge<-merge_phyloseq(ps16, ps23)

             meta <- psmerge %>%
                sample_data %>%
                as_tibble %>%
                mutate("sample_names" = psmerge %>% sample_names )

            dv16<-ps16  %>% DivNet::divnet(., formula = ~Sample +Sample:pH + Site, base="ESV_2")
            dv23<-ps23  %>% DivNet::divnet(., formula = ~pH + Site)

            shannons <- dv23$shannon %>% summary %>% rbind(dv16$shannon %>% summary )
            shannonmerge <- meta %>%
                left_join(shannons,
                            by = "sample_names") %>%
                mutate_if( is.character,as_factor) %>%
                mutate(ARMS=as.factor(ARMS))%>%
                mutate  (pH=factor(pH, order = TRUE,  levels = c("Con", "Med", "Low")) ) %>%
                mutate  (Sample=factor(Sample,  levels = c("100", "Algae", "Sessile", "CS", "RVS", "Algae_23s")))


#save output

saveRDS(shannonmerge, file=file.path(path, "Outputs", "PNG_EsvShannonDataFrame.RDS"))