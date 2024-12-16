#!/usr/bin/env Rscript
# date:4/10/2023 by Hongjian Jin @ St Jude Children's Research Hospital
######################################################

######################################################
library(optparse)
options(warn=-1)
option_list <- list(
  make_option(c("-L", "--filelist"), type="character",  default=NA, 
    help="character. filenames delimited by text")
  ,make_option(c("-O", "--output"), type="character", default=NA, 
    help="character. output filename")
  ,make_option(c("-s", "--species"), type="character", default=NULL, 
    help="character. use as suffix like hsa or mmu [default %default]")
  )

parser <- OptionParser(usage="%prog -L <filelist> -O <output> [-t <inner|outer>]",
                        description = "Combine sRNAtools count,RPM tables",
                        option_list = option_list,
                        epilogue =paste("Example:\n",
                        " list=`ls *table.result.txt|tr '\\n' ','` ",
                        " sRNAtools_merge_matrix.R -L $list -O hendegrp_287359_miRNAseq ",
                        "\n\n",sep="\n") 
                        )
args<-NA
options(warn=0)
result<-tryCatch({
    args <- parse_args(parser, positional_arguments = TRUE)  # TRUE = c(0, Inf), FALSE,1,  c(1,2)
    }, warning=function(w){
        message(e)
        cat("\n\n")
    }, error=function(e){
        message(e)
        cat("\n\n")
 })

opt <- args$options
if (sum(is.na(opt))>0 ) {
    print_help(parser)
    quit("no")
}

filelist1<- unlist(strsplit(opt$filelist,","))
filelist1 <- filelist1[filelist1!=""]
filelist1 <- filelist1[file.size(filelist1) != 0L]

prefix<- opt$output
species <- opt$species
if (!is.null(species)){
 species <- paste0("_",species)
}
total<-length(filelist1)

if (total ==0 ){
     stop(paste0("No valid sRNAtool outfiles found\n"))
}else{
    cat(paste0("\nInfo:",total, " sRNAtool outfiles found."))
}
#=====================================================
cat("\n\nProcessing inputs:\n")
anno_cols <- c("Tag_name","sncRNA_type","Most_abundant_tag_seq")
dat_cols <- c("Tag_name", "Total_tag_number","Total_RPM")
all_cols <- unique(anno_cols,dat_cols)
anno <- df1 <- df2 <- res <- NULL
ids <- gsub(".table.result.txt","",basename(filelist1))
x <- 0
#suppressPackageStartupMessages(library("data.table"))
for(inFile in filelist1) {
    x <- x+1
    dat <- read.table(inFile, sep="\t",header=TRUE,fill=TRUE,stringsAsFactors = FALSE, quote="",row.names=NULL ,check.names=FALSE,comment.char = "" )
    
    cat("\n\t",x,ids[x],"dim[", nrow(dat),"x",ncol(dat),"]")
    if (!all(all_cols %in% colnames(dat))){
      cat("\nError: invalid sRNA output." )
      next
    }
    df2 <- dat[, dat_cols]
    colnames(df2)[2:3] <- paste(ids[x],".",colnames(df2)[2:3],sep='')
    if (is.null(df1)){
        df1 <- df2
        anno <- dat[,anno_cols]
        if (total == 1){
          res <- df1
        }
        next
    }else if (!is.null(res)){
        df1 <- res
    }
    anno <- rbind(anno, dat[,anno_cols])
    anno <- anno[!duplicated(anno$Tag_name), ]
    res <- merge(x = df1, y = df2, by.x = "Tag_name", by.y = "Tag_name",all = TRUE) #Outer join by colname 
    cat("\n\t merged dim[", nrow(res),"x",ncol(res),"]\n")
}
cat("\n\t annotation rows[", nrow(anno),"]\n")
rownames(anno) <- anno$Tag_name
#print(colnames(res))
#colnames(res) <- gsub("gene_id","geneID",colnames(res))
#=====================================================

cat("\n\nSaving outputs:\n")
final_count <- data.frame(anno[res$Tag_name, ], res[, grep("Total_tag_number",colnames(res))], check.names=FALSE)
colnames(final_count) <- gsub(".Total_tag_number","", colnames(final_count))
outFile1 <- paste(prefix,species,"_all.count.txt",sep="")
final_count[is.na(final_count)] <- 0  # 12/04/2024
write.table(final_count,file=outFile1,  quote=F,sep="\t", row.names=F, col.names=T)
cat("\n\t",outFile1, "[saved", nrow(final_count),"x",ncol(final_count),"]")

final_RPM <- data.frame(anno[res$Tag_name, ], res[, grep("Total_RPM",colnames(res))], check.names=FALSE)
colnames(final_RPM) <- gsub(".Total_RPM","", colnames(final_RPM))
outFile2 <- paste(prefix,species,"_all.RPM.txt",sep="")
write.table(final_RPM,file=outFile2,  quote=F,sep="\t", row.names=F, col.names=T)
cat("\n\t",outFile2, "[saved", nrow(final_RPM),"x",ncol(final_RPM),"]")
if (sum(final_count$sncRNA_type %in% "miRNA")>0){
    miRNA_count <- final_count[final_count$sncRNA_type %in% "miRNA", ]
    miRNA_count[is.na(miRNA_count) | miRNA_count =='NA' ] <- 0

    outFile3 <- paste(prefix,species,"_miRNA.count.txt",sep="")
    write.table(miRNA_count,file=outFile3,  quote=F,sep="\t", row.names=F, col.names=T)
    cat("\n\t",outFile3, "[saved", nrow(miRNA_count),"x",ncol(miRNA_count),"]")

    miRNA_RPM <- final_RPM[final_RPM$sncRNA_type %in% "miRNA", ]
    miRNA_RPM[is.na(miRNA_RPM) | miRNA_RPM=='NA'] <- 0
    outFile4 <- paste(prefix,species,"_miRNA.RPM.txt",sep="")
    write.table(miRNA_RPM,file=outFile4,  quote=F,sep="\t", row.names=F, col.names=T)
    cat("\n\t",outFile4, "[saved", nrow(miRNA_RPM),"x",ncol(miRNA_RPM),"]")
}else{
   cat("\nWarn: no any miRNA detected.\n")
}


cat("\n\nCheers!\n\n")
quit("no")

######################################################
#update log
######################################################
# 4/10/2023, first version 
# 4/20/2023, check file.size, absence of miRNA, single input
