#Adjusted Burns code to incorporate two pools and a mixing parameter

#Adam Burns - 2/10/2015
#aburns2@uoregon.edu
#From Burns et al. Contribution of neutral processes to the assembly of the gut microbial communities changes over host development
#Fits the neutral model from Sloan et al. 2006 to an OTU table and returns several fitting statistics. Alternatively, will return predicted occurrence frequencies for each OTU based on their abundance in the metacommunity when stats=FALSE. For use in R.
#spp: A community table for communities of interest with local communities/samples as rows and taxa as columns. All samples must be rarefied to the same depth.
#pool: A community table for defining source community (optional; Default=NULL).
#taxon: A table listing the taxonomic calls for each otu, with OTU ids as row names and taxonomic classifications as columns.
#If stats=TRUE the function will return fitting statistics.
#If stats=FALSE the function will return a table of observed and predicted values for each otu.

#with mixing parameter
sncm.fit.new <- function(spp=NULL, pool=NULL, stats=TRUE, taxon=NULL){
	require(minpack.lm)
	require(Hmisc)
	require(stats4)
	
	options(warn=-1)

	#Calculate the number of individuals per community
	N <- mean(apply(spp, 1, sum, na.rm=TRUE))
	

	#Calculate the average relative abundance of each taxa across communities

		p.m1 <- apply(spp, 2, mean)
		p.m1 <- p.m1[p.m1 != 0]
		p1 <- as.data.frame(p.m1/N)
		p.m2 <- apply(pool, 2, mean)
		p.m2 <- p.m2[p.m2 != 0]
		p2 <- as.data.frame(p.m2/N)

	#Calculate the occurrence frequency of each taxa across communities
	spp.bi <- 1*(spp>0)
	freq <- apply(spp.bi, 2, mean)
	freq <- as.data.frame(freq[freq != 0])

	#Combine
    merge1<-merge(p1, freq, by =0, all=TRUE)
    names(merge1)<-c("ESV","p1", "freq")
    p2<-rownames_to_column(p2, "ESV")
    C<-merge(merge1, p2, by ="ESV")
    names(C)<-c("ESV","p1", "freq", "p2")
	C <- C[order(C[,3]),]
	C <- as.data.frame(C)
	C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),] #Removes rows with any zero (absent in either source pool or local communities)
	p1 <- C.0[,2]
	freq <- C.0[,3]
    p2 <- C.0[,4]
	names(p1) <- C.0[,1]
	names(freq) <- C.0[,1]
	names(p2) <- C.0[,1]

	#Calculate the limit of detection
	d = 1/N

	##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
	m.fit <- nlsLM(freq ~ mix*pbeta(d, N*m.holo*p1, N*m.holo*(1-p1), lower.tail=FALSE) +(1-mix)*pbeta(d, N*m.env*p2, N*m.env*(1-p2), lower.tail=FALSE), start=list(m.holo=0.1,m.env=0.1, mix=0.5))

	m.ci <- confint(m.fit, c('m.holo', 'm.env', 'mix'), level=0.95)
	

	##Calculate goodness-of-fit (R-squared and Root Mean Squared Error)
    estim.m.holo<-coef(m.fit)[1]
	estim.m.env<-coef(m.fit)[2]
    estim.mix<-coef(m.fit)[3]
	freq.pred <- estim.mix*pbeta(d, N*estim.m.holo*p1, N*estim.m.holo*(1-p1), lower.tail=FALSE) +(1-estim.mix)*pbeta(d, N*estim.m.env*p2, N*estim.m.env*(1-p2), lower.tail=FALSE)
	Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
	RMSE <- sqrt(sum((freq-freq.pred)^2)/(length(freq)-1))
	
	pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
	
	##Results
	#if(stats==TRUE){
		fitstats <- data.frame(m.holo=numeric(), m.holo.ci=numeric(), m.env=numeric(), m.env.ci=numeric(), mix=numeric(), mix.ci=numeric(), Rsqr=numeric(), RMSE=numeric(), N=numeric(), Samples=numeric(), Richness=numeric(), Detect=numeric())
		fitstats[1,] <- c(coef(m.fit)[1], coef(m.fit)[1]-m.ci[1], coef(m.fit)[2], coef(m.fit)[2]-m.ci[2], coef(m.fit)[3], coef(m.fit)[3]-m.ci[3], Rsqr, RMSE, N, nrow(spp), length(p1), d)
		#return(fitstats)
	#} else {
		A <- cbind(p1, p2, freq, freq.pred, pred.ci[,2:3])
		A <- as.data.frame(A)
		colnames(A) <- c('p1', "p2", 'freq', 'freq.pred', 'pred.lwr', 'pred.upr')
		if(is.null(taxon)){
			B <- A[order(A[,1]),]
		} else {
			B <- merge(A, taxon, by=0, all=TRUE)
			row.names(B) <- B[,1]
			B <- B[,-1]
			B <- B[order(B[,1]),]
		}
		return(list(fitstats,B))
	#}
}
