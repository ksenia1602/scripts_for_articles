
options(stringsAsFactors=FALSE)

###libraries
library("DESeq2")
library("tximport")
library("readr")
library("ReportingTools")
library("AnnotationDbi")
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(edgeR)
library(mygene)
library(rtracklayer)
library(reshape2)
library("ggplot2") #Best plots
library("ggrepel") #Avoid overlapping labels
library(dplyr)
library(FactoMineR)
library("ggsci")
library("scales")
library("fgsea")
library(mygene)
library(ReactomePA )
library(clusterProfiler)
library(DT)
require(enrichplot)
require(tidyverse)
require(org.Hs.eg.db)
library(msigdbr)
library(tidyr)
library(dplyr)
############


###input files 
colorr<-rev(pal_simpsons("springfield")(16)[c(1,6,7,13)])

#gff
GTF <- import("C:/Users/User/Dropbox/work_files/genome/gencode.v29.chr_patch_hapl_scaff.annotation.gtf")
gff<-as.data.frame(GTF[,c("gene_id","gene_name","gene_type","transcript_id")])
#gff<-subset(gff,gff$gene_type=="protein_coding")
my_txdf<-as.data.frame(gff[,c("transcript_id","gene_id")])
colnames(my_txdf)<-c("tx_id","gene_id")
gff<-gff[,c("gene_id","gene_name","gene_type")]
gff<-gff[!duplicated(gff),]


options(ggrepel.max.overlaps = Inf)
pca_plot<-function(data,names,type,colorr){
  
  
  
  
  x1<-data[,names]
  x1<-na.omit(x1)
  x1<-apply(x1,1,as.numeric)
  data_cl<-(x1)+1
  data_cl<-log2(data_cl)
  
  
  
  res.pca = PCA(data_cl, ncp=5,  graph=F)
  axis_value<-round(as.numeric(unlist(res.pca$eig[,1]*100/sum(as.numeric(unlist(res.pca$eig[,1])))))[c(1,2)],digits = 1)
  dat<-as.data.frame(res.pca$ind$coord)
  dat$Name<-names
  dat$type<-type
  
  
  
  pca<-ggplot(dat, aes(x =Dim.1,y=Dim.2, 
                       #fill=factor(disease),
                       label=Name)) + 
    scale_color_manual(values = c(colorr)) +
    geom_point(
      aes(colour=factor(type)),
      #colour="red",
      size=3,
      #shape=21, 
      #lwd=2,
      #colour="black",
      alpha = 1)+
    
    
    geom_label_repel(
      aes(Dim.1, Dim.2, label = Name),size = 2.5,
      box.padding = unit(0.2, "lines"),
      point.padding = unit(0.2, "lines"),
      segment.colour = NA)+
    
    theme(panel.background = element_rect(colour = NA))+ 
    theme(legend.justification = 'bottom') + 
    
    labs(color = " ")+
    
    theme_bw() +  
    theme(
      panel.border = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.text.x = element_text(colour="black",size=9,angle=0#,hjust=.9,vjust=.5,face="bold"
      ),
      axis.text.y = element_text(colour="black",size=9,angle=0#,hjust=1,vjust=0,face="bold"
      ),
      axis.title.x = element_text(colour="black",size=9,angle=0#,hjust=.5,vjust=0,face="bold"
      ),
      axis.title.y = element_text(colour="black",size=9,angle=90#,hjust=.5,vjust=.5,face="bold"
      ),
      #panel.background =  element_rect(fill = NA, colour = NA), 
      panel.background = element_rect(fill = "white"),
      #panel.border =      element_rect(fill = NA, colour=NA), 
      panel.grid.major =  element_line(colour =NA , size = 0.2),
      panel.grid.minor =  element_line(colour = NA, size = 0.5),  
      legend.text=element_text(size=9#,face="bold"
      ),
      legend.title=element_text(size=9#,face="bold"
      ),
      plot.title = element_text(size = 9#, face ="bold",hjust = 0.5
      ),
      legend.position="top"
    ) + #ggtitle("Individuals factor map(PCA)")+
    ggtitle(" ")+
    scale_y_continuous(expand = expansion(mult = c(0, .15),
                                          add = c(10, 0)))+
    scale_x_continuous(expand = expansion(mult = c(0, .15),
                                          add = c(10, 0)))+
    ylab(paste0("PC2 (",axis_value[2], "% explained var.)")) + 
    xlab(paste0("PC1 (",axis_value[1], "% explained var.)"))
  
  return(pca)
}
samlmon_output_nonpair_group_tpm<-function(sample_inf,dir,group_comp,my_txdf){
  
  setwd(dir )
  samples<-list.files(path = dir,
                      pattern = "quant.sf", recursive = TRUE,full.names = T,
                      include.dirs = T)
  samples_name<-list.files(path = dir,
                           pattern = "quant.sf", recursive = TRUE,full.names = F,
                           include.dirs = F)
  
  
  names<-unlist(lapply(strsplit(samples_name,"/"), `[[`, 1))
  data_samples<-data.frame(directory=samples,name=names)
  data_samples<-merge(data_samples,sample_inf,by="name")
  # data_samples$or<-as.numeric(data_samples$name)
  #data_samples<-  data_samples[  order(data_samples$or),]
  
  txi <- tximport(data_samples$directory, type="salmon", tx2gene=my_txdf,abundanceCol = "TPM", 
                  countsCol = "NumReads", lengthCol = "Length",ignoreTxVersion=F,
                  countsFromAbundance = "lengthScaledTPM")
  
  ExpDesign <- data.frame(id= data_samples$name,prepost=data_samples$group,clone=data_samples$clone)
  rownames(ExpDesign)<-ExpDesign$id
  dds <- DESeqDataSetFromTximport(txi, ExpDesign, design= ~prepost)
  ddsColl <- collapseReplicates(dds, dds$id, dds$clone)
  ddsColl  <- DESeq( ddsColl )
  design <- model.matrix(~ExpDesign$prepost)
  
  
  ###pca
  
  dds <- DESeq( dds )
  counts_table = counts( dds, normalized=TRUE )
  counts_table<-counts_table[!rowSums(counts_table==0)>=nrow(data_samples)/3, ]
  
  GeneID<-rownames(counts_table)
  filtered_norm_counts<-cbind(counts_table,GeneID)
  data<-as.data.frame(filtered_norm_counts)
  colnames(data)<-data_samples$names_for_col
  
  col<-data_samples$names_for_col
  
  pca<-pca_plot(data, data_samples$names_for_col,
                data_samples$group,
                colorr)
  
  
  ##tpm eddition
  data<-txi$abundance
  data<-as.data.frame(data)
  colnames(data)<-data_samples$names_for_col
  data$GeneID<-rownames(data)
  ########
  res <- results(ddsColl , contrast=c("prepost", group_comp[1],group_comp[2]) , 
                 independentFiltering=TRUE, 
                 pAdjustMethod="fdr")
  res_ordered<-res[order(res$padj),]
  GeneID<-rownames(res_ordered)
  res_ordered<-as.data.frame(res_ordered)
  res_genes<-cbind(res_ordered,GeneID)
  res_genes_merged <- merge(res_genes,data,by=unique("GeneID"))
  res_ordered<-res_genes_merged[order(res_genes_merged$padj),]
  res_ordered<-merge(res_ordered,gff,by.x= "GeneID",by.y="gene_id")
  res_ordered<-res_ordered[order(res_ordered$padj),]
  return(list(res_ordered,pca))
}
enrich_kegg_reactome_go<-function(data){
  
  xx<-as.data.frame(queryMany(data$gene_name, scopes= "symbol", 
                              fields=c("entrezgene"), species="human"))
  xx<-xx[,c("entrezgene","query")]
  xx<-xx[!is.na(xx$entrezgene),]
  xx<-xx[!duplicated(xx$query),]
  xx<-xx[!duplicated(xx$entrezgene),]
  
  kk <- enrichKEGG(xx$entrezgene,
                   organism = "hsa",pAdjustMethod = "fdr",
                   qvalueCutoff=1,
                   pvalueCutoff=1,
                   minGSSize=1,use_internal_data=T 
                   
  )
  kegg<-(as.data.frame(kk))
  
  kk <- enrichPathway(xx$entrezgene,
                      qvalueCutoff=1,
                      pvalueCutoff=1
                      ,
                      minGSSize=1,  readable = T
  )
  reactome<-(as.data.frame(kk))
  
  
  ego <- enrichGO(gene          = xx$entrezgene,
                  
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "fdr",
                  pvalueCutoff  = 1,
                  qvalueCutoff  = 1,
                  readable      = TRUE)
  go<-(as.data.frame(ego))
  return(list(kegg,reactome,go))
}
volc_plot<-function(x1,genesel,col1,col2){
  
  mutateddf <- mutate(x1, sig=ifelse(x1$padj<0.05&abs(x1$log2FoldChange)>1,
                                     "fdr < 0.05 and abs(log2FC)>1", "Not Sig")) 
  
  #convert the rownames to a column
  volc = ggplot(mutateddf, aes(log2FoldChange, -log10(pvalue))) + #volcanoplot with log2Foldchange versus pvalue
    geom_point(aes(col=sig),alpha = 1/2) + #add points colored by significance
    #scale_colour_paletteer_d("ggsci", "blue_material") +
    scale_color_manual(values = c( col1,col2)) +
    #geom_text(data=filter(mutateddf, PValue<0.05), aes(label=gene_name))
    geom_text_repel(data=filter(mutateddf, mutateddf$gene_name%in% genesel), 
                    aes(label=gene_name),size = 5,
                    box.padding = unit(0.35, "lines"),
                    point.padding = unit(0.3, "lines"))+
    theme_bw() +  
    theme(
      panel.border = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.text.x = element_text(colour="black",size=12,angle=0,face="bold"),
      axis.text.y = element_text(colour="black",size=12,angle=0,face="bold"),
      axis.title.x = element_text(colour="black",size=12,angle=0,face="bold"),
      axis.title.y = element_text(colour="black",size=12,angle=90,face="bold"),
      #panel.background =  element_rect(fill = NA, colour = NA), 
      panel.background = element_rect(fill = "white"),
      #panel.border =      element_rect(fill = NA, colour=NA), 
      panel.grid.major =  element_line(colour =NA , size = 0.2),
      panel.grid.minor =  element_line(colour = NA, size = 0.5),  
      legend.text=element_text(size=12,face="bold"),
      legend.title=element_text(size=0,face="bold"),
      plot.title = element_text(size = 0, face ="bold.italic",hjust = 0.5),
      legend.position="top"
    )#adding text for the top 20 genes
  #ggsave("Volcanoplot.jpeg", device="jpeg") #In case you want to easily save to disk
  
  volc<- volc+geom_point(data = filter(mutateddf, mutateddf$gene_name%in% genesel), 
                         aes(log2FoldChange, -log10(pvalue)),
                         col = 'red') 
  
  return(volc)}
gsea_output<-function(data){
  
  genelist<-data$gene_name 
  genelist<- genelist %>%  bitr(fromType = "SYMBOL",
                                toType = c("ENSEMBL", "ENTREZID"),
                                OrgDb = org.Hs.eg.db)
  colnames(genelist)[1]<-"gene_name"
  data<-genelist %>%
    inner_join(data,by='gene_name') %>% 
    dplyr::select(ENTREZID,stat,everything())
  geneList<-data$stat
  names(geneList)<-as.character(data$ENTREZID)
  geneList<-geneList[!duplicated(names(geneList))]
  geneList<-sort(geneList,decreasing=TRUE)
  
  
  m_t2g_h <- msigdbr(species = "Homo sapiens", category = "H") %>% 
    dplyr::select(gs_name, entrez_gene)
  m_t2g_kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>% 
    dplyr::select(gs_name, entrez_gene)
  m_t2g_bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP") %>% 
    dplyr::select(gs_name, entrez_gene)
  m_t2g_reatome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>% 
    dplyr::select(gs_name, entrez_gene)
  
  gsea_kegg<-GSEA(geneList,nPerm = 100000 ,TERM2GENE = m_t2g_kegg,pvalueCutoff = 0.05,verbose = FALSE)
  gsea_halmark<-GSEA(geneList,nPerm = 100000 ,TERM2GENE = m_t2g_h,pvalueCutoff = 0.05,verbose = FALSE)
  gsea_reactome<-GSEA(geneList,nPerm = 100000 ,TERM2GENE = m_t2g_reatome,pvalueCutoff = 0.05,verbose = FALSE)
  gsea_bp<-GSEA(geneList,nPerm = 100000 ,TERM2GENE = m_t2g_bp,pvalueCutoff = 0.05,verbose = FALSE)
  return(list(reactome=gsea_reactome,
              kegg= gsea_kegg,
              halmarks= gsea_halmark,
              gsea_bp=  gsea_bp))
}
kentrez_to_symbol<-function(data,x,y){
  data$col_new<-data[,x]
  data<-data %>% 
    mutate(col_new = strsplit(as.character(col_new), "/")) %>% 
    unnest(col_new)
  genelist<-data$col_new
  genelist<- genelist %>%  bitr(toType = "SYMBOL",
                                fromType = c( "ENTREZID"),
                                OrgDb = org.Hs.eg.db)
  names(data)[names(data) == 'col_new'] <- "ENTREZID"
  
  data<-genelist %>%
    inner_join(data,by="ENTREZID") 
  data<-subset(data,select=c(-ENTREZID))
  data<-aggregate(SYMBOL~., data, 
                  FUN = function(X) paste(unique(X), collapse="/"))
  data<-data[order(data[,y]),]
  return(data)}



sample_inf1<-data.frame(name=c(2:18),
                        group=c(rep("NP normal",8),
                                rep("NP Parkinson",9)),
                        
                        clone=c(rep("PO2_old",2),
                                rep("Olya_old",3),
                                rep("C1_old",3),
                                rep("patientTr5_old",3),
                                rep("patient13_old",3),
                                rep("patient14_old",3)),
                 
                        
                        names_for_col=c(paste0(rep("norma1np",2),"_",c(1,2)),
                                      paste0(rep("norma3np",3),"_",c(1,2,3)),
                                      paste0(rep("norma2np",3),"_",c(1,2,3)),
                                      paste0(rep("PARK2-PD1np",3),"_",c(1,2,3)),
                                      paste0(rep("PARK2-PD2np",3),"_",c(1,2,3)),
                                      paste0(rep("PARK2-PD3np",3),"_",c(1,2,3))),
                        
                        title=c(paste0(rep("norma1 rep",2),c(1,2)," line NP RG2L"),
                                paste0(rep("norma3 rep",2),c(1,2,3)," line NP RFD 3.9 L"),
                                paste0(rep("norma2 rep",2),c(1,2,3)," line NP HD 1.1S"),
                                
                                paste0(rep("PARK2-PD1 rep",2),c(1,2,3)," line NP PDL 1.5L"),
                                paste0(rep("PARK2-PD1 rep",2),c(1,2,3)," line NP PDS13"),
                                paste0(rep("PARK2-PD1 rep",2),c(1,2,3)," line NP PDS14")),
                        source_name="Neural progenitor",
                        organism="Homo Sapiens", 
                        
                        cell_type="neural progenitors",
                        genotype= c(rep("normal",8),
                                    rep("PARK2: 2 exon (del 202-203 AG), 1 intron (splicing mutation IVS1+1G/A)",3),
                                    rep("PARK2: homozygous deletion of 8 exon",3),
                                    rep("PARK2: heterozygous deletion of 2 exon",3)),
                        molecule="RNA",
                        description= c(rep("iPSCs from  healthy donors  were differentiated into uncommitted neural progenitors",8),
                                       rep("iPSCs from   Parkinson's disease patients  were differentiated into uncommitted neural progenitors",9)),
                        
                        
                        processed_data_file="TPM_NP.tsv"
                        
)





sample_inf2<-data.frame(name=c("1-1", "1-2" ,"1-3", "3-1",
                               "3-2", "3-3" ,"5-1", "5-2", "5-3",
                               "6-1", "6-2", "6-3", "7-1", "7-2" ,"7-3" ),
                        
                        group=c(rep("DN normal",6),
                                rep("DN Parkinson",9))
                        ,
                        
                        clone=c(rep("C1",3),
                                rep("PO2",3),
                                rep("patient13",3),
                                rep("patientTr5",3),
                                rep("patient14",3))
                        ,
                        
                     
                        
                        
                        names_for_col=c(paste0(rep("norma2dn",3),"_",c(1,2,3)),
                                      
                                      paste0(rep("norma1dn",3),"_",c(1,2,3)),
                                      paste0(rep("PARK2-PD2dn",3),"_",c(1,2,3)),
                                      paste0(rep("PARK2-PD1dn",3),"_",c(1,2,3)),
                                      paste0(rep("PARK2-PD3dn",3),"_",c(1,2,3))),
                        
                        title=c(
                          paste0(rep("norma2 rep",3),c(1,2,3)," line DN HD 1.1S"),
                          paste0(rep("norma1 rep",3),c(1,2,3)," line DN RG2L"),
                          
                          
                          paste0(rep("PARK2-PD1 rep",3),c(1,2,3)," line DN PDS13"),
                          paste0(rep("PARK2-PD1 rep",3),c(1,2,3)," line DN PDL 1.5L"),
                          
                          paste0(rep("PARK2-PD1 rep",3),c(1,2,3)," line DN PDS14")),
                        
                        source_name="terminally differentiated neurons",
                        organism="Homo Sapiens", 
                        
                        cell_type="neurons presumably dopaminergic",
                        genotype= c(rep("normal",6),
                                    rep("PARK2: homozygous deletion of 8 exon",3),
                                    rep("PARK2: 2 exon (del 202-203 AG), 1 intron (splicing mutation IVS1+1G/A)",3),
                                    
                                    rep("PARK2: heterozygous deletion of 2 exon",3)),
                        molecule="RNA",
                        description= c(rep("iPSCs from  healthy donors  were differentiated into terminally differentiated neurons",6),
                                       rep("iPSCs from   Parkinson's disease patients  were differentiated into terminally differentiated neurons",9)),
                        
                        
                        processed_data_file="TPM_DN.tsv"
                        
)


sample_inf<-rbind(sample_inf1,sample_inf2)

colorr<-rev(pal_simpsons("springfield")(16)[c(1,6,13,7)])
group_comp<-c("DN Parkinson","DN normal")
dir<-"C:/Users/User/Dropbox/work_files/project_files/RNA_seq_data_analysis/IMG_new"
data_neuron<-samlmon_output_nonpair_group_tpm(sample_inf2,dir,group_comp,my_txdf)
pca1<-data_neuron[[2]]
pca1


colorr<-rev(pal_simpsons("springfield")(16)[c(1,6)])
group_comp<-c("NP Parkinson","NP normal")
dir<-"C:/Users/User/Dropbox/work_files/project_files/RNA_seq_data_analysis/IMG_new"
data_np<-samlmon_output_nonpair_group_tpm(sample_inf1,dir,group_comp,my_txdf)
pca2<-data_np[[2]]
pca2

library("cowplot")

plot_grid(pca2,pca1,nrow=1,rel_widths =  c(1,1),labels = c('A', 'B'))

ggsave("C:/Users/User/Dropbox/work_files/project_files/RNA_seq_data_analysis/IMG_new/pca.pdf",
       height = 5,width = 10)



colorr<-rev(pal_simpsons("springfield")(16)[c(1,6,13,7)])
group_comp<-c("NP Parkinson","NP normal")
dir<-"C:/Users/User/Dropbox/work_files/project_files/RNA_seq_data_analysis/IMG_new"
data_all<-samlmon_output_nonpair_group_tpm(sample_inf,dir,group_comp,my_txdf)
pca3<-data_all[[2]]
pca3

ggsave("C:/Users/User/Dropbox/work_files/project_files/RNA_seq_data_analysis/IMG_new/pca_all.pdf",
       height = 5,width = 10)


data_dn<-data_neuron[[1]]
data_dn<-data_dn[,c(1,23,8:22)]

data_np1<-data_np[[1]]
data_np1<-data_np1[,c(1,25,8:24)]


write.table(data_np1,"C:/Users/User/Dropbox/work_files/project_files/RNA_seq_data_analysis/IMG_new/TPM_NP.tsv",row.names = F,quote = F)
write.table(data_dn,"C:/Users/User/Dropbox/work_files/project_files/RNA_seq_data_analysis/IMG_new/TPM_DN.tsv",row.names = F,quote = F)

#library(tools)
#dn<-md5sum("C:/Users/User/Dropbox/work_files/project_files/RNA_seq_data_analysis/IMG_new/TPM_DN.tsv")
#dn
#np<-md5sum("C:/Users/User/Dropbox/work_files/project_files/RNA_seq_data_analysis/IMG_new/TPM_NP.tsv")
#np
