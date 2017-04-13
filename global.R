#   ScanGEO: parallel mining of high-throughput gene expression data
#   Copyright (C) <2017>  <Katja Koeppen>

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#    General Public License for more details <http://www.gnu.org/licenses/>.


library(shiny)
library(GEOmetadb)
library(randomForest)

First.wd <- getwd()
con <- dbConnect(SQLite(), "GEOmetadb.sqlite")
Organism <- dbGetQuery(con, "select platform_organism from gds")
Organism <- unlist(strsplit(Organism$platform_organism, ","))
Genus <- strsplit(Organism, " ")
Genus <- unlist(lapply(Genus, function (x){x[1]}))
Genus <- sort(names(head(sort(table(Genus), decreasing = TRUE), 20)))

load("Org.rds")
load("MasterListAll.rds")
load("RF.rds")

ScanGeo <- function(geneList, gdsList, alpha){	
        ScanResult = list()
        ScanPvalues = matrix(nrow = length(geneList), ncol = length(gdsList), 
                             data = NA, dimnames = list(geneList, gdsList))
        for(gds in gdsList){
                incProgress(1/length(gdsList))
                print(gds)
                myF = NULL
                
                G <- try(getGEO(gds, destdir = paste(First.wd, "/softfiles", sep="")))
                if (class(G) == "GDS") {
                
                if(dim(Columns(G))[2] > 3){
                        myF <- apply(apply(Columns(G)[,-c(1,dim(Columns(G))[2])], 2, as.character), 1, 
                                     function(x){paste(x, collapse="_")})
                } else {
                        myF <- as.character(Columns(G)[,-c(1,dim(Columns(G))[2])])
                }
                
                # This requires that all groups have at least two samples
                if(min(table(myF)) > 1){
                        d <- Table(G)
                        dMat <- apply(d[,-c(1,2)], 2, as.numeric)
                        rownames(dMat) <- as.vector(d$ID_REF)
                        for(gene in geneList){
                                print(paste("Finding probes for", gene, gds))
                                
                                probes <-  as.vector(G@dataTable@table$"ID_REF")[as.vector(G@dataTable@table$IDENTIFIER) == gene]
                                
                                # probes <- probes[- which(is.na(probes))]
                                for (p in as.vector(probes)){
                                        if(is.na(p)) {next}
                                        print(paste("Finding values for", gene, gds, p))
                                        if((sum(is.na(dMat[p,])) == 0) & (sum(dMat[p,] > 0 ) == dim(dMat)[2])){
                                                print(paste("ANOVA for", gene, gds, p))
                                                
                                                myP = anova(lm(log2(dMat[p,]) ~ myF))$"Pr(>F)"[1]
                                                
                                                # myP is for a probe, not a gene
                                                # we only update the gene p-value if the probe p value is less than the existing p-value
                                                if (is.nan(myP)) {next}
                                                
                                                if (is.na(ScanPvalues[gene, gds])){
                                                        ScanPvalues[gene, gds] <- myP
                                                } else if (ScanPvalues[gene, gds] > myP) {
                                                        ScanPvalues[gene, gds] <- myP
                                                }
                                                
                                                if (myP < alpha){
                                                        myPDF = sub('/', '_', paste(gene, p, gds, "pdf", sep='.'))
                                                        myCSV = sub('/', '_', paste(gene, p, gds, "csv", sep='.'))
                                                        print(paste("Creating boxplot:", myPDF))
                                                        
                                                        myP = anova(lm(log2(dMat[p,]) ~ myF))$"Pr(>F)"[1]
                                                        CSVmat = matrix(nrow = 1, ncol = length(myF), data = dMat[p,], dimnames = list(gene, myF))
                                                        write.csv(CSVmat, file = myCSV)
                                                        
                                                        pdf(file=myPDF)
                                                        if(nchar(G@header$title) > 60){
                                                                par(cex.main=.6)
                                                        } else {
                                                                par(cex.main=1.2)
                                                        }
                                                        
                                                        if (max(unlist(lapply(myF, nchar)) > 12)){
                                                                par(cex.axis=.5)
                                                                par( oma = c(3, 3, 3, 3))
                                                                par( mar = c( 15, 4.1, 4.1, 2.1))
                                                                
                                                        } else {
                                                                par(cex.axis=1)
                                                                par( oma = c(0, 0, 0, 0))
                                                                par(mar =c( 5.1, 4.1, 4.1, 2.1))
                                                        }
                                                        
                                                        
                                                        ScanResult[[gds]] = c(ScanResult[[gds]], gene)
                                                        
                                                        stripchart(log2(dMat[p,]) ~ factor(myF), 
                                                                   main=paste(G@header$title, paste(paste(gene, p, gds, "pdf", sep='.'), 
                                                                        paste("P =", round(myP, 3)), sep = " : "), sep="\n"),
                                                                   ylab=paste(gene, "RNA Expression (log2)"), las=2, 
                                                                   xlim=c(0.5, length(levels(factor(myF))) + 0.5),
                                                                   vertical = TRUE, pch = 16, col = "black") 
                                                        means <- tapply(log2(dMat[p,]), factor(myF), mean)
                                                        loc <- 1:length(means)
                                                        segments(loc-0.3, means, loc+0.3, means, col="red", lwd=3)
                                                     
                                                        dev.off()
                                                } 
                                                else {
                                                        print (paste("P exceeds target alpha", gene, gds, p, myP))
                                                }
                                        }
                                        else { 
                                                print (paste("Data issues:", gene, gds, p))
                                                dMat[p,]
                                        }
                                }
                                
                        }
                        
                }
        }}
        write.csv(ScanPvalues, file = "03_pValues_summary.csv")
        return(ScanResult)
        
}
