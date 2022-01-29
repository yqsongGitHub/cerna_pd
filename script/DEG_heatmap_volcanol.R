annodata <- function(xdata,ydata,byx="ID",byy="row.names",xselect=c("ID","GeneSymbol")){
  xdata <-  xdata[,xselect]
  data_ann <- merge(xdata,ydata,by.x=byx,by.y=byy)
  data_ann <- data_ann[data_ann$GeneSymbol!='',]
  data_ann[[xselect[2]]]<- gsub("//.*","",data_ann[[xselect[2]]])
  data_ann[[xselect[2]]] <- gsub(" ","", data_ann[[xselect[2]]])
  data_ann <- data_ann[!duplicated(data_ann[[xselect[2]]]),]
  rownames(data_ann) <- data_ann[[xselect[1]]]
  data_ann <- data_ann[,-1]
  return(data_ann)
}
DEGfunc<- function(data,group_list,gpl=NULL,type="array",pvalue=0.05,fc=1.5){
  # type="array" or "RNAseq"
  if(type=="array"){
    group_list <- factor(group_list$grouplist)
    group_list <- relevel(group_list, ref="normal")
    data  <- normalizeBetweenArrays(data)
    design <- model.matrix(~0+factor(group_list))
    colnames(design) <- levels(factor(group_list))
    contrast.matrix <- makeContrasts("disease-normal",
                                     levels = design)
    fit <- lmFit(data,design)
    fit2 <- contrasts.fit(fit, contrast.matrix) 
    fit2 <- eBayes(fit2) 
    DEG <- topTable(fit2, coef=1, n=Inf) 
  }else if(type=="RNAseq"){
    exprSet <- data
    group_list <- c(group_list$grouplist)
    design <- model.matrix(~0+factor(group_list))
    colnames(design) <- levels(factor(group_list))
    rownames(design) <- colnames(exprSet)
    dge <- DGEList(counts=exprSet)
    dge <- calcNormFactors(dge)
    v <- voom(dge,design,plot=FALSE, normalize="quantile")
    fit <- lmFit(v, design)
    cont.matrix <- makeContrasts(contrasts=c('disease-normal'),levels = design)
    fit2 <- contrasts.fit(fit,cont.matrix)
    fit2 <-  eBayes(fit2)
    tempOutput <- topTable(fit2, coef='disease-normal', n=Inf)
    DEG <- na.omit(tempOutput)
  }
  if(!is.null(gpl)){
    gpl <- gpl[,c("ID","GeneSymbol")]
    colnames(gpl)[1] <- "GENE"
    DEG$GENE <- rownames(DEG)
    DEG <- merge(DEG,gpl,by="GENE")
    DEG <- DEG[DEG$GeneSymbol!='',]
    DEG$GeneSymbol <- gsub("//.*","",DEG$GeneSymbol)
    DEG$GeneSymbol <- gsub(" ","",DEG$GeneSymbol)
    DEG <- DEG[!duplicated(DEG$GeneSymbol),]
    rownames(DEG) <- DEG$GENE
    DEG <- DEG[,-1]
  }
  DEGG <- DEG[DEG$P.Value < pvalue&abs(DEG$logFC)>log(fc),]
  re <- list(DEG,DEGG)
  return(re)
  #write.csv(re[[1]],"data-DE.csv")
  #write.csv(re[[2]],"data-DEG.csv")
}
DEGpheatmap <- function(data,group_list,DEG,
                     Scale="row", # "row", "column" and "none"
                     Clusterrows=T,
                     Clustercols=T,
                     Annotation_legend=T,
                     Showrownames=F
                     ){
  heatdata <- data[rownames(data) %in% rownames(DEG),]
  group_list <- factor(group_list$grouplist)
  group_list <- relevel(group_list, ref="normal")
  annotation_col <- data.frame(group_list)
  rownames(annotation_col) <- colnames(heatdata)
  p <- pheatmap(heatdata, 
                cluster_rows = Clusterrows,
                cluster_cols = Clustercols,
                annotation_col =annotation_col, 
                annotation_legend=Annotation_legend, 
                show_rownames = Showrownames,
                scale = Scale, 
                color =colorRampPalette(c("blue", "white","red"))(100),
                #filename = "data-heatmap.pdf",
                #cellwidth = 20, cellheight = 0.1,
                fontsize = 10)
  return(p)
}
DEGvolcano <- function(DEG,pvalue=0.05,fc=1.5,Title="PD vs NC"){
  DEG$Condition=ifelse(DEG$logFC>=fc & DEG$P.Value<=pvalue,"Up",
                       ifelse(DEG$logFC<=fc*(-1) &DEG$P.Value<=pvalue,"Down",
                              "None"))
  DEG$sig=ifelse(DEG$logFC>=fc& DEG$P.Value<=pvalue,"red",
                 ifelse(DEG$logFC<=fc*(-1) &DEG$P.Value<=pvalue,
                        "blue","gray"))
  p <- ggplot(DEG, aes(x = logFC, y = -log10(P.Value), color = Condition)) +
    geom_point(alpha = 0.6, size = 1) +
    scale_colour_manual(values  = c('red2', 'blue2', 'gray'), limits = c('Up', 'Down', 'None')) +
    theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), plot.title = element_text(hjust = 0.5)) +
    theme(legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent'), legend.position = c(0.9, 0.93)) +
    geom_vline(xintercept = c(-log2(fc), log2(fc)), color = 'gray', size = 0.3) +
    geom_hline(yintercept = -log(pvalue, 10), color = 'gray', size = 0.3) +
    xlim(-10, 10) + ylim(0, 6) +
    labs(x = '\nLog2 Fold Change', y = 'Log10 P.Value\n', color = '', title = Title)
  return(p)
}
enrichment_GO_KEGG_Reactome <- function(genesymbol,pvaluecutoff=0.05,qvaluecutoff=1,methods="GO"){
  # methods="GO";"KEGG";"Reactome"
  eg  <-  bitr(genesymbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  if(methods=="GO"){
    all <- enrichGO(eg[,2],'org.Hs.eg.db', keyType = "ENTREZID", ont = "ALL", 
                    pvalueCutoff = pvaluecutoff,qvalueCutoff = qvaluecutoff, readable = TRUE)
    p <- dotplot(all, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
    print(p)
    return(all)
  }else if(methods=="KEGG"){
    KEGG <- enrichKEGG(eg[,2], organism = "hsa", keyType = "kegg", 
                       pvalueCutoff = pvaluecutoff, qvalueCutoff=qvaluecutoff,use_internal_data=T)
    p <- dotplot(KEGG,showCategory=5)
    print(p)
    return(KEGG)
  }else if(methods=="Reactome"){
    rea <- enrichPathway(eg[,2], organism = "human", pvalueCutoff = pvaluecutoff, 
                         qvalueCutoff = qvaluecutoff, readable = TRUE)
    p <-  dotplot(rea, showCategory=5)
    print(p)
    return(rea)
  }
}
DEG_cal_cor <- function(data,gpl,data_lnc,data_mRNA,data_miRNA,method="pearson",corvalue=-0.4,pvalue=0.05){
  data_ann <- annodata(gpl,data)
  rownames(data_ann) <- data_ann$GeneSymbol
  data_ann <- data_ann[,-1]
  list.mRNA <- unique(data_mRNA$GeneSymbol)
  list.miRNA <- unique(data_miRNA$GeneSymbol) 
  list.lnc <- unique(data_lnc$GeneSymbol)
  data_ann1 <- data_ann[row.names(data_ann)%in% c(list.mRNA,list.miRNA,list.lnc),]
  data_ann1 <- t(as.matrix(data_ann1))
  cor.re<- Hmisc::rcorr(data_ann1,type=method)
  cor.value <- cor.re[["r"]]
  cor.pvalue <- cor.re[["P"]]
  cor.value <- data.table::melt(cor.value)
  cor.value1 <-cor.value[cor.value$value<corvalue,]
  cor.value1$pair <- paste0(cor.value1$Var1,"_",cor.value1$Var2)
  cor.pvalue <- data.table::melt(cor.pvalue)
  cor.pvalue1 <- cor.pvalue[cor.pvalue$value<pvalue,]
  cor.pvalue1$pair <- paste0(cor.pvalue1$Var1,"_",cor.pvalue1$Var2)
  # miRNA vs lnc
  cor.value.mi <- cor.value1[cor.value1$Var1 %in% list.miRNA & (cor.value1$Var2 %in%list.lnc),]
  cor.pvalue.mi <- cor.pvalue1[cor.pvalue1$Var1 %in% list.miRNA & (cor.pvalue1$Var2 %in%list.lnc),]
  cor.mi <- merge(cor.value.mi,cor.pvalue.mi,by="pair")
  cor.mi <- cor.mi[,c(2,3,4,7)]
  colnames(cor.mi) <- c("miRNA","lncRNA","cor","pvalue")
  # miRNA vs mRNA
  cor.value.mRNA <- cor.value1[cor.value1$Var1 %in% list.miRNA & (cor.value1$Var2 %in%list.mRNA),]
  cor.pvalue.mRNA <- cor.pvalue1[cor.pvalue1$Var1 %in% list.miRNA & (cor.pvalue1$Var2 %in%list.mRNA),]
  cor.mRNA <- merge(cor.value.mRNA,cor.pvalue.mRNA,by="pair")
  cor.mRNA <- cor.mRNA[,c(2,3,4,7)]
  colnames(cor.mRNA) <- c("miRNA","mRNA","cor","pvalue")
  # all
  cor.mi1 <- cor.mi
  cor.mRNA1 <- cor.mRNA
  colnames(cor.mi1)[2] <- "genename"
  colnames(cor.mRNA1)[2] <- "genename"
  all <- rbind(cor.mi1,cor.mRNA1 )
  re <- list("miRNA-lnc"=cor.mi,
             "miRNA-mRNA"=cor.mRNA,
             "all"=all)
  return(re)
}
calculate_cor <- function(data_ann,list.miRNA,list.mRNA,methods="pearson",pvalue=0.05,corvalue=-0.4){
  ddf0 <- data.frame()
  for (yy in list.miRNA){
    df0 <- data.frame()
    for( xx in list.mRNA){
      val1 <- na.omit(as.numeric(data_ann[rownames(data_ann) %in%yy,,drop=F]))
      val2 <- na.omit(as.numeric(data_ann[rownames(data_ann) %in%xx,,drop=F]))
      if (length(val1) >3 & length(val2) >3){
        cor.mi_m <- cor.test(val1,val2,method=methods)
        df <- data.frame("genelist1"=yy,
                         "genelist2"=xx,
                         "p.value"=cor.mi_m$p.value,
                         "cor"=as.numeric(cor.mi_m$estimate)
        )
        df0 <- rbind(df0,df)
      }
    }
    ddf0 <- rbind(ddf0,df0)
  }
  ddf0 <- ddf0[ddf0$p.value<pvalue & ddf0$cor<corvalue,]
  return(ddf0)
}
pre_database <- function(data,group_list,gpl,RNAmap,miRNA_trans,mi_lnc,mi_m,calculate_DEG_p=0.05,calculate_DEG_fc=1,pre_type="array"){
  DEGG <- DEGfunc(data,group_list,gpl,type=pre_type,pvalue=calculate_DEG_p,fc=calculate_DEG_fc)
  DE <- DEGG[[2]]
  data_all <- merge(DE,RNAmap,by.x="GeneSymbol",by.y="gene_name")
  data_lnc <-data_all[data_all$type=="lncRNA",]
  data_mRNA <- data_all[data_all$type=="mRNA",]
  data_miRNA <- data_all[data_all$type=="miRNA",]
  miRNA_trans1 <- miRNA_trans[miRNA_trans$gene_name %in% data_miRNA$GeneSymbol,]
  data_lnc_mi_conn <- mi_lnc[mi_lnc$miRNA_name %in%miRNA_trans1$miRNA_name & mi_lnc$lnc_name %in% data_lnc$GeneSymbol, ]
  data_lnc_mi_conn1 <- merge(miRNA_trans,data_lnc_mi_conn,by="miRNA_name",all.y= T)
  data_lnc_mi_conn1 <- data_lnc_mi_conn1[,c("miRNA_name","gene_name","lnc_name")]
  colnames(data_lnc_mi_conn1) <- c("miRNA","miRNA_GeneSymbol","lncRNA")
  data_m_mi_conn <- mi_m[mi_m$miRNA_name %in%miRNA_trans1$miRNA_name & mi_m$gene_symbol %in% data_mRNA$GeneSymbol, ]
  data_m_mi_conn1 <- merge(miRNA_trans,data_m_mi_conn,by="miRNA_name",all.y= T)
  colnames(data_m_mi_conn1) <- c("miRNA","miRNA_GeneSymbol","mRNA")
  # all
  data_lnc_mi_conn2 <- data_lnc_mi_conn1
  data_m_mi_conn2 <- data_m_mi_conn1
  colnames(data_lnc_mi_conn2)[3] <- "genename"  
  colnames(data_m_mi_conn2)[3] <- "genename"
  data_mi_lnc_mRNA <- rbind(data_lnc_mi_conn2,data_m_mi_conn2)
  return(list("miRNA-lnc"=data_lnc_mi_conn1,
              "miRNA-mRNA"=data_m_mi_conn1,
              "miRNA-lnc-mRNA"=data_mi_lnc_mRNA
  ))
}
calculate_pre <- function(data,gpl,list.miRNA,list.lnc,list.mRNA,
                          calculate_cor_methods="pearson",
                          calculate_cor_pvalue=0.05,
                          calculate_cor_corvalue=-0.4){
  data_ann <- annodata(gpl,data)
  rownames(data_ann) <- data_ann$GeneSymbol
  data_ann <- data_ann[,-1]
  cor.mi_m <- calculate_cor(data_ann,list.miRNA,list.mRNA,methods = calculate_cor_methods,pvalue = calculate_cor_pvalue,corvalue = calculate_cor_corvalue)
  cor.mi_lnc <- calculate_cor(data_ann,list.miRNA,list.lnc,methods = calculate_cor_methods,pvalue =calculate_cor_pvalue,corvalue = calculate_cor_corvalue)
  mi_m_lnc_merge <- rbind(cor.mi_lnc,cor.mi_m)
  return(mi_m_lnc_merge)
}
calculate_cor_all <- function(data,group_list,
                              gpl,RNAmap,
                              miRNA_trans,
                              mi_lnc,
                              mi_m,
                              calculate_cor_methods="pearson",
                              calculate_cor_DEG_p=0.05,
                              calculate_cor_DEG_fc=1,
                              calculate_cor_pvalue=0.05,
                              calculate_cor_corvalue=-0.4,
                              type="array"
                              ){
  mi_mRNA_lnc <- pre_database(data,group_list,gpl,RNAmap,miRNA_trans,mi_lnc,mi_m,calculate_cor_DEG_p,calculate_cor_DEG_fc,pre_type=type)
  list.miRNA <- intersect(mi_mRNA_lnc[["miRNA-lnc"]]$miRNA_GeneSymbol,mi_mRNA_lnc[["miRNA-mRNA"]]$miRNA_GeneSymbol)
  list.mRNA <- unique(mi_mRNA_lnc[["miRNA-mRNA"]]$mRNA)
  list.lnc <-unique(mi_mRNA_lnc[["miRNA-lnc"]]$lncRNA)
  cor_re <- calculate_pre(data,gpl,list.miRNA,list.lnc,list.mRNA,
                          calculate_cor_methods=calculate_cor_methods,
                          calculate_cor_pvalue=calculate_cor_pvalue,
                          calculate_cor_corvalue=calculate_cor_corvalue)
  return(cor_re)
}

draw_ceRNA <- function(mi_m_lnc,layout_ce=layout.fruchterman.reingold,
                       vertex.size_ce=5,
                       vertex.shape_ce="circle",
                       vertex.label.cex_ce=0.7,
                       vertex.label.color_ce="black",
                       vertex.label.dist_ce=0,
                       edge.width_ce=0.2,
                       edge.color_ce="gray"
){
  # layout_ce:layout.fruchterman.reingold,layout.circle,layout.reingold.tilford,
  #           layout_as_star,layout_with_drl,layout.sphere,layout.random
  # vertex.shape_ce:none,rectangle
  nodedf <- data.frame("RNA"=unique(c(mi_m_lnc[,1],mi_m_lnc[,2],mi_m_lnc[,3])))
  nodedf$group <- ifelse(nodedf$RNA %in%mi_m_lnc[,1],1,ifelse(nodedf$RNA %in%mi_m_lnc[,2],2,3))
  nodedf$ID <- 1:nrow(nodedf)
  colnames(mi_m_lnc) <- rep("col",3)
  linkdf <- do.call(rbind,list(unique(mi_m_lnc[,c(1,2)]),unique(mi_m_lnc[,c(1,3)]),unique(mi_m_lnc[,c(2,1)]),unique(mi_m_lnc[,c(3,1)])))
  colnames(linkdf) <- c("from","to")
  net<-graph.data.frame(linkdf,nodedf,directed = F) 
  colrs<- c("red","green","gray")
  V(net)$color<-colrs[V(net)$group] 
  plot(net,  
       layout=layout_ce,  
       vertex.size=vertex.size_ce,     
       vertex.shape=vertex.shape_ce,     
       vertex.label=NULL,  
       vertex.label.cex=vertex.label.cex_ce,    
       vertex.label.color=vertex.label.color_ce,  
       vertex.label.dist=vertex.label.dist_ce,   
       edge.arrow.size=0,
       edge.width = edge.width_ce, 
       edge.label=NA, 
       edge.curved=0,
       edge.color=edge.color_ce)  
}