##########################
library(stringr); library(reshape2); library(data.table); library(dplyr);
library(viridis);


memory.limit(size = 24000)

setwd("~/Documents/3way_input/")
outputDirectory='~/Documents/output/' #directory for output files

`%+%` = function(a, b) paste0(a, b)
clipOut = function( x ) { write.table(x,"clipboard-16384", sep="\t", row.names=F, quote=F) }
clipIn = function( ) { read.table("clipboard", header=TRUE, na.strings=c("NA","NaN","?"), sep="\t"); }
numNas = function( dt ) { sapply( pcaPointsForIndiv, function(x){sum(is.na(x))} )   }
cosineDist <- function(x){ as.dist(1 - x%*%t(x)/(sqrt(rowSums(x^2) %*% t(rowSums(x^2))))) }
mbind<-function(...){
    Reduce( function(x,y){cbind(x,y[match(row.names(x),row.names(y)),])}, list(...) )
}


#########################################################
#Model 3 :3-way interactions between two genes (geneA and geneB) with three genomics alterations (geneA_mut:geneA_CNA:geneB_mut) across cancer type 
#WithSecondM = with geneB (the second gene) mutation, WOSecondM = without geneB (the second gene) mutation
#input_file name: (1) Obs_LLR_randomization_AcrossCancer_mutation_Tissue_Loss_Threeway.txt, (2) Obs_LLR_randomization_AcrossCancer_mutation_Tissue_Gain_Threeway.txt

loadInputData = function( fn  ) {
    
    inData = fread(fn, data.table = T)  # glimpse(inData)
    
    inData = inData[, colnames(inData) %>% str_subset( "(Mut|NoMut)_(Del|Loss|WT|Gain|Amp)$") %>% c("Gene","Target", .), with=F ]
    
    inData =
    inData %>% melt( id.vars=c("Gene", "Target"), variable.name="variantType", value.name="N" ) %>%
    mutate( mut = c("Mut"=1, "NoMut"=0)[ str_match( variantType, "([^_]+)_[^_]+" )[,2] ] ) %>%
    mutate( cna = c("Loss"=-1, "WT"=0, "Gain"=1)[ str_match( variantType, "[^_]+_([^_]+)")[,2] ] ) %>%
    select(-variantType)
    
    return(inData)
}


for(tumor in c('ACC', 'DLBC', 'GBM', 'PCPG', 'BLCA', 'UCEC', 'THCA', 'CESC', 'UCS', 'TGCT', 'LIHC', 'CHOL', 'HNSC', 'STAD', 'SKCM', 'UVM', 'PAAD', 'THYM', 'LUSC', 'MESO', 'OV', 'ESCA', 'SARC', 'KIRP', 'LGG', 'LAML', 'COADREAD', 'PRAD', 'LUAD', 'BRCA', 'KIRC', 'KICH')){
    
    Input_Input_direc <- paste("~/3way_input/Obs_LLR_randomization_AcrossCancer_mutation_",tumor,"_Loss_Threeway.txt", sep="")
    Input.dataLong = loadInputData(Input_Input_direc)
    Input.data = Input.dataLong %>% reshape2::dcast( Gene+Target ~ mut+cna, value.var="N")
    
    cat(sprintf( "%s\n", Input_Input_direc ));
    Input.data <- fread(file = Input_Input_direc, data.table = F)
    

    Input.data <- Input.data[,c(1:3,4:7)]
    
    names(Input.data) <- gsub(x = names(Input.data), pattern = "NoMut", replacement = 0) %>%
    gsub(pattern = "Mut", replacement = 1) %>%
    gsub(pattern = "Loss", replacement = -1) %>%
    gsub(pattern = "WT", replacement = 0) %>%
    gsub(pattern = "Gain", replacement = 1) %>%
    gsub(pattern = "WithSecondM", replacement = 1) %>%
    gsub(pattern = "WOSecondM", replacement = 0)
    
    Input.data <- as.data.table(Input.data)
    
    Input.dataLlm = copy(Input.dataLong)
    Input.dataLlm <- as.data.table(Input.dataLlm)
    
    Input.dataLlm[, mut:=as.factor(mut)]
    Input.dataLlm[, cna:=as.ordered(cna)]
    Input_results <- c()
    
    #####################
    cat(sprintf( "Tissue %s\n", tumor ));
    for ( aGene in Input.dataLlm[, unique(Gene)] %>% setdiff("_allCancerGenes") ) {
        
        #cat(sprintf( "Gene %s\n", aGene ));
        
        aLlm <- glm(N ~ mut +cna + Target + cna:mut + cna:Target + mut:Target + cna:mut:Target, family=poisson(link="log"),
        control=glm.control(epsilon = 1e-6,maxit=100), data=Input.dataLlm[Gene == aGene])
        
        aLlCoef = summary(aLlm) %>% coef
        row.names(aLlCoef)[1] <- paste('intercept_',aGene,sep='')
        Input_results <-Input_results %>%rbind(aLlCoef)
    }
    
    Input_results <- as.data.frame(Input_results)
    write.table(Input_results %>% mutate(Tissue=tumor,variable=row.names(Input_results)),
    file= paste(outputDirectory,"Model3_Obs_", tumor,"_Loss1.txt",sep=''),
    quote=FALSE, sep='\t',row.names=F,col.names = T)
}
