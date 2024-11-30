library(DESeq2)
library(sqldf)
library(stringr)
library(gplots)
library(ggplot2)
library(enrichR)
library(patchwork)
library(ggVennDiagram)






################################################################################
#### read mapping wbid - gene symbol ###########################################
################################################################################

map_wbid_genesym_new <- read.csv("finalstrainannotation/genecoord_20230331.csv.gz",sep="\t")
colnames(map_wbid_genesym_new) <- c("wbid","transcript","genestart","geneend","chr","genesym")
map_wbid_genesym <- unique(map_wbid_genesym_new[,c("wbid","genesym")])

####### Read vibrio data
cnt <- read.csv("fromvibrio/cnt_clean.csv", stringsAsFactors = FALSE, sep=" ")
cond <- read.csv("fromvibrio/cond_clean.csv", stringsAsFactors = FALSE, sep=" ")

###########################################################################################
############################### DE analysis - n2 ##########################################
###########################################################################################

write_round <- function(de_wt, fname){
  de_rounded <- de_wt[,c("wbid","genesym","log2FoldChange","pvalue","padj")]
  de_rounded$padj <- signif(de_rounded$padj, digits = 2)
  de_rounded$pvalue <- signif(de_rounded$pvalue, digits = 2)
  de_rounded$log2FoldChange <- round(de_rounded$log2FoldChange, digits = 2)
  write.csv(de_rounded, fname, row.names = FALSE)
}

keep <- cond$strain=="n2"
dds <- DESeqDataSetFromMatrix(countData = cnt[,keep], colData = cond[keep,], design = ~infected)
dds <- DESeq(dds)

de_wt <- as.data.frame(results(dds, contrast = c("infected","inf","uninf")))  

pcutoff <- 1e-20

sum(de_wt$padj<pcutoff & de_wt$log2FoldChange<0, na.rm = TRUE) #407
sum(de_wt$padj<pcutoff & de_wt$log2FoldChange>0, na.rm = TRUE) #551

de_wt$wbid <- rownames(de_wt)

de_wt <- merge(de_wt, map_wbid_genesym)
de_wt$padj[is.na(de_wt$padj)] <- 1
de_wt <- de_wt[order(de_wt$pvalue),]

write.csv(de_wt[,c("wbid","genesym","log2FoldChange","pvalue","padj")], "out/de_wt_inf_vs_uninf.csv", row.names = FALSE)
write_round(de_wt, "out/de_wt_inf_vs_uninf.2digit.csv")

###########################################################################################
############################### DE analysis - ko ##########################################
###########################################################################################

keep <- cond$strain=="g"
dds <- DESeqDataSetFromMatrix(countData = cnt[,keep], colData = cond[keep,], design = ~infected)
dds <- DESeq(dds)

de_ko <- as.data.frame(results(dds, contrast = c("infected","inf","uninf")))  

pcutoff <- 1e-20

sum(de_ko$padj<pcutoff & de_ko$log2FoldChange<0, na.rm = TRUE) #347
sum(de_ko$padj<pcutoff & de_ko$log2FoldChange>0, na.rm = TRUE) #1003

de_ko$wbid <- rownames(de_ko)

de_ko <- merge(de_ko, map_wbid_genesym)
de_ko$padj[is.na(de_ko$padj)] <- 1
de_ko <- de_ko[order(de_ko$pvalue),]
de_ko <- de_ko

write.csv(de_ko[,c("wbid","genesym","log2FoldChange","pvalue","padj")], "out/de_ko_inf_vs_uninf.csv", row.names = FALSE)
write_round(de_ko, "out/de_ko_inf_vs_uninf.2digit.csv")


################################################################################
######################### DE analysis - effect of mutation (inf) ###############
################################################################################

keep <- cond$infected=="inf"
dds <- DESeqDataSetFromMatrix(countData = cnt[,keep], colData = cond[keep,], design = ~strain)
dds <- DESeq(dds)

de_muteffect_inf <- as.data.frame(results(dds, contrast = c("strain","g","n2")))  

pcutoff <- 1e-20

sum(de_muteffect_inf$padj<pcutoff & de_muteffect_inf$log2FoldChange<0, na.rm = TRUE) #
sum(de_muteffect_inf$padj<pcutoff & de_muteffect_inf$log2FoldChange>0, na.rm = TRUE) #

de_muteffect_inf$wbid <- rownames(de_muteffect_inf)

de_muteffect_inf <- merge(de_muteffect_inf, map_wbid_genesym)
de_muteffect_inf$padj[is.na(de_muteffect_inf$padj)] <- 1
de_muteffect_inf <- de_muteffect_inf[order(de_muteffect_inf$pvalue),]

write.csv(de_muteffect_inf[,c("wbid","genesym","log2FoldChange","pvalue","padj")], "out/de_muteffect_inf.csv", row.names = FALSE)
write_round(de_muteffect_inf, "out/de_muteffect_inf.2digit.csv")



################################################################################
##################### DE analysis - effect of mutation (uninf) #################
################################################################################

keep <- cond$infected=="uninf"
dds <- DESeqDataSetFromMatrix(countData = cnt[,keep], colData = cond[keep,], design = ~strain)
dds <- DESeq(dds)

de_muteffect_uninf <- as.data.frame(results(dds, contrast = c("strain","g","n2")))  

pcutoff <- 1e-20

sum(de_muteffect_uninf$padj<pcutoff & de_muteffect_uninf$log2FoldChange<0, na.rm = TRUE) 
sum(de_muteffect_uninf$padj<pcutoff & de_muteffect_uninf$log2FoldChange>0, na.rm = TRUE) 

de_muteffect_uninf$wbid <- rownames(de_muteffect_uninf)

de_muteffect_uninf <- merge(de_muteffect_uninf, map_wbid_genesym)
de_muteffect_uninf$padj[is.na(de_muteffect_uninf$padj)] <- 1
de_muteffect_uninf <- de_muteffect_uninf[order(de_muteffect_uninf$pvalue),]

write.csv(de_muteffect_uninf[,c("wbid","genesym","log2FoldChange","pvalue","padj")], 
          "out/de_muteffect uninf.csv", row.names = FALSE)
write_round(de_muteffect_uninf, "out/de_muteffect uninf.2digit.csv")



################################################################################
############################### Get genes for GO terms #########################
################################################################################


get_sym_for_go <- function(goid){
  library(org.Ce.eg.db)
  retrieved <- AnnotationDbi::select(org.Ce.eg.db, keytype="GOALL", keys=goid, columns="ENSEMBL")  
  unique(de_wt$genesym[de_wt$wbid %in% retrieved$ENSEMBL])
}


#glist_defence <- get_sym_for_go("GO:0050829") #defense response to Gram-negative bacterium
glist_defence_all <- get_sym_for_go("GO:0042742") #defense response to bacterium
glist_stress <- get_sym_for_go("GO:0006950") #Response to stress
glist_immune <- get_sym_for_go("GO:0006955") #Immune response


################################################################################
############# Fig S2 ab - GO analysis for each strain - effect of infection ####
################################################################################


do_signed_go_analysis <- function(list_comp,out_csv,pcutoff){
  
  setEnrichrSite("WormEnrichr")
  listEnrichrDbs()
  dbs <- c("GO_Biological_Process_2018")

  all_go_plots <- NULL
  table_go <- NULL
  all_sgn <- c(-1,1)
  for(one_comp in names(list_comp)){
    for(one_sgn in all_sgn){
      one_de <- list_comp[[one_comp]] 
      
      glist <- as.character(na.omit(one_de$wbid[one_de$padj<pcutoff & one_de$log2FoldChange*one_sgn>0]))  
      print(paste(one_comp,one_sgn,length(glist)))
      
      enriched <- enrichr(map_wbid_genesym$genesym[map_wbid_genesym$wbid %in% glist], dbs)
      
      
      for(cur_db in dbs){
        df <- enriched[[cur_db]]
        df$genotype <- one_comp
        if(one_sgn>0){
          name_genedir <- "UP"
        } else {
          name_genedir <- "DOWN"
        }
        df$genedir <- name_genedir
        df <- df[order(df$Combined.Score, decreasing = TRUE),]
        df$genecount <- unlist(lapply(str_split(df$Genes,";"),length))
        
        sub_df <- df[1:10,]
        sub_df <- sub_df[order(sub_df$P.value, decreasing = TRUE),]
        sub_df$Term <- factor(sub_df$Term,(sub_df$Term))
        onep <- ggplot(sub_df, aes(Term, -log10(Adjusted.P.value))) + 
          coord_flip() +
          ggtitle(paste(one_comp, name_genedir)) +
          xlab("") +
          ylab("-log10(P.adj)")+
          theme_bw() 
        if(one_sgn>0){
          onep <- onep + 
            geom_bar(position = "dodge", stat="identity", fill="darkorange")
        } else {
          onep <- onep + 
            geom_bar(position = "dodge", stat="identity", fill="darkcyan") 
        }
        all_go_plots[[paste(one_comp, one_sgn)]] <- onep
        onep
        
        table_go <- rbind(table_go, df[1:10,])
        
        #ggsave(sprintf("out/go_%s %s.pdf", one_comp, one_sgn))
        #write.csv(df[,c("Term","Overlap","P.value","Adjusted.P.value","genotype","genedir")],
        #          sprintf("out/go_%s %s.csv", one_comp, one_sgn), row.names = FALSE)
      }     
      
    }
  }
  
  write.csv(
    table_go[,c("Term","Overlap","P.value","Adjusted.P.value","genotype","genedir","genecount")],
    out_csv, 
    row.names = FALSE)
  
  ptot <- egg::ggarrange(plots=all_go_plots, ncol=2)
  ptot
}


list_comp <- list(
  ko=de_ko,
  wt=de_wt
)

pcutoff <- 1e-20

ptot <- do_signed_go_analysis(
  list_comp,
  "out/totalgotable.csv",
  pcutoff
)
ggsave(plot=ptot, "out/fig S2abcd all_topgo.pdf", width = 10, height = 5)




################################################################################
#################### fig 6a volcano plot ####################################### effect of mutation #1
################################################################################


plot_volcano_highlight <- function(one_de, highlight_genes, pcutoff) {
  one_de$padj[one_de$padj==0] <- 1e-300 
  
  one_de$color <- "notde"
  one_de$color[one_de$log2FoldChange<0 & one_de$padj<pcutoff] <- "down"
  one_de$color[one_de$log2FoldChange>0 & one_de$padj<pcutoff] <- "up"
  
  ggplot(one_de, aes(log2FoldChange, -log10(padj), color=color, label=genesym)) +
    geom_point() +
    geom_point(data=one_de[one_de$genesym %in% highlight_genes,], color="black", color="black")+
    scale_colour_manual(values=c("#c4c633","lightgray","#acd5db"))+
    theme_bw() +
    ggrepel::geom_text_repel(data=one_de[one_de$genesym %in% highlight_genes,], color="black", max.overlaps = 1000) +
    theme(legend.position = "none")
}

pcutoff <- 1e-10

highlight_genes <- sprintf("catp-%s",1:8)
one_de <- de_muteffect 
plot_volcano_highlight(de_muteffect, highlight_genes, pcutoff) +
  xlab("log2(FC) shk-1 vs WT, no infection")
ggsave("out/fig 6a volcano.pdf", width = 5, height = 5)




################################################################################
#################### fig S4 a volcano plot ##################################### effect of mutation #2
################################################################################

pcutoff <- 1e-10

highlight_genes <- c("zip-8","cav-1")
one_de <- de_muteffect 
plot_volcano_highlight(one_de, highlight_genes, pcutoff)+
  xlab("log2(FC) shk-1 vs WT, no infection") #+
ggsave("out/fig S4 a volcano.pdf", width = 5, height = 5)




################################################################################
############# fig 2ab - GO analysis, effect of mutation UNINF ##################
################################################################################



list_comp <- list(
  mut_uninf=de_muteffect_uninf
)

pcutoff <- 1e-20

ptot <- do_signed_go_analysis(
  list_comp,
  "out/go_mut_uninf.csv",
  pcutoff
)
ggsave(plot=ptot, "out/fig 2ab muteffect uninf go.pdf", width = 20, height = 3)


################################################################################
############# fig 2ab - GO analysis, effect of mutation INF ####################
################################################################################



list_comp <- list(
  mut_inf=de_muteffect_inf
)

pcutoff <- 1e-20

ptot <- do_signed_go_analysis(
  list_comp,
  "out/go_mut_inf.csv",
  pcutoff
)
ggsave(plot=ptot, "out/fig ____ muteffect inf go.pdf", width = 20, height = 3)




###########################################################################################
############################### fig 2e: Heatmap of DE genes ###############################
###########################################################################################


glist <- glist_defence_all

de_comb <- merge(
  data.frame(wt_padj=de_wt$padj, wt_log2FoldChange = de_wt$log2FoldChange, genesym=de_wt$genesym),
  data.frame(ko_padj=de_ko$padj, ko_log2FoldChange = de_ko$log2FoldChange, genesym=de_ko$genesym)
)
de_comb <- de_comb[de_comb$wt_padj < pcutoff | de_comb$ko_padj < pcutoff,]
de_comb <- de_comb[de_comb$genesym %in% glist,]
rownames(de_comb) <- de_comb$genesym
de_comb <- de_comb[,c("wt_log2FoldChange","ko_log2FoldChange")]
colnames(de_comb) <- c("WT","shk-1")
de_comb <- reshape2::melt(as.matrix(de_comb))


ggplot(de_comb, aes(Var1,Var2,fill=value)) + 
  geom_tile() +
  theme_bw() + 
  scale_fill_gradient2(high="red",low="blue", name="Log2 FC") + 
  xlab("") + 
  ylab("") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("out/fig 2e heatmap.pdf", width = 8, height = 1.5)







###########################################################################################
############################### just a number: p.adj, def against bacterium ###############
###########################################################################################


if(FALSE){
  test_def <- function(glist, use_de) {
    pcutoff <- 1e-20
    
    use_de$is_cat <-  use_de$genesym %in% glist
    use_de$is_de  <-  use_de$padj < pcutoff
    
    fisher.test(use_de$is_cat, use_de$is_de)  
  }
  
  test_def(glist_defence_all, de_muteffect)
  test_def(glist_defence_all, de_muteffect_inf) #0.038
  test_def(glist_defence_all, de_ko)
  test_def(glist_defence_all, de_wt)
  
  test_def(glist_immune, de_muteffect)
  test_def(glist_immune, de_muteffect_inf) #0.108
  test_def(glist_immune, de_ko)
  test_def(glist_immune, de_wt)
  

  test_def(glist_stress, de_muteffect) #8e-7  uninf
  test_def(glist_stress, de_muteffect_inf) #0.155
  test_def(glist_stress, de_ko)
  test_def(glist_stress, de_wt)
  
  
}








################################################################################
############# fig s2f - scatter plot, immune response genes #####################
################################################################################

pcutoff <- 1e-20

de_comb <- merge(
  data.frame(wt_padj=de_wt$padj, wt_log2FoldChange = de_wt$log2FoldChange, genesym=de_wt$genesym),
  data.frame(ko_padj=de_ko$padj, ko_log2FoldChange = de_ko$log2FoldChange, genesym=de_ko$genesym)
)
de_comb <- de_comb[de_comb$wt_padj < pcutoff | de_comb$ko_padj < pcutoff,]

ggplot(de_comb[de_comb$genesym %in% glist_immune,], aes(wt_log2FoldChange, ko_log2FoldChange)) + 
  geom_point() +
  theme_bw() +
  xlab("WT log2(FC)") + 
  ylab("shk-1 log2(FC)")
ggsave("out/fig 2f scatter.pdf", width = 5, height = 5)







################################################################################
########### fig S2 g - scatter plot, immune response genes - with highlight ####   
################################################################################

pcutoff <- 1e-20


highlight_genes <- c(
  "clec-67",
  "dod-22",
  "cpr-3",
  "dod-24",
  "clec-4",
  "fil-1",
  "irg-4",
  "clec-66",
  "lys-1",
  "dod-19",
  "clec-85",
  "irg-5",
  "lys-2",
  "fat-6",
  "asm-3",
  "ech-6",
  "fat-7",
  "acdh-1",
  "acdh-2"
)


de_comb <- merge(
  data.frame(wt_padj=de_wt$padj, wt_log2FoldChange = de_wt$log2FoldChange, genesym=de_wt$genesym),
  data.frame(ko_padj=de_ko$padj, ko_log2FoldChange = de_ko$log2FoldChange, genesym=de_ko$genesym)
)
de_comb <- de_comb[de_comb$wt_padj < pcutoff | de_comb$ko_padj < pcutoff,]

de_comb$to_highlight <- de_comb$genesym %in% highlight_genes

de_comb <- de_comb[order(de_comb$to_highlight),]

ggplot(de_comb, aes(wt_log2FoldChange, ko_log2FoldChange, color=to_highlight, label=genesym)) + 
  geom_point() +
  scale_colour_manual(values=c("lightgray","red"), guide = "none") +
  ggrepel::geom_text_repel(data=de_comb[de_comb$genesym %in% highlight_genes,], max.overlaps = 1000, color="black")+
  theme_bw() +
  xlab("WT log2(FC)") + 
  ylab("shk-1 log2(FC)")
ggsave("out/fig S2 g scatter.pdf", width = 5, height = 5)



################################################################################
######################### fig 2c: Venn diagram, all genes ######################
################################################################################

pcutoff <- 1e-20

p1 <- ggvenn::ggvenn(
  list(
    wt=de_wt$genesym[de_wt$padj < pcutoff & de_wt$log2FoldChange<0 ],
    ko=de_ko$genesym[de_ko$padj < pcutoff & de_ko$log2FoldChange<0 ]
  ),
  show_percentage=FALSE
) + ggtitle("GO:0050829, down-regulated")

p2 <- ggvenn::ggvenn(
  list(
    wt=de_wt$genesym[de_wt$padj < pcutoff & de_wt$log2FoldChange>0 ],
    ko=de_ko$genesym[de_ko$padj < pcutoff & de_ko$log2FoldChange>0 ]
  ),
  show_percentage=FALSE
) + ggtitle("GO:0050829, up-regulated")


ptot <- egg::ggarrange(p1,p2, nrow=1)
ggsave(plot=ptot, "out/fig 2c venn.pdf", width = 5, height = 5)



################################################################################
#################### fig 2d: Venn diagram, defence genes #######################
################################################################################

pcutoff <- 1e-20

glist <- glist_defence_all 
ggvenn::ggvenn(
  list(
    wt=de_wt$genesym[de_wt$padj < pcutoff & de_wt$genesym %in% glist],
    ko=de_ko$genesym[de_ko$padj < pcutoff & de_ko$genesym %in% glist]
  ),
  show_percentage=FALSE
) + ggtitle("GO:0050829")

ggsave("out/fig 2d venn.pdf", width = 5, height = 5)






################################################################################
#################### fig xxx: Venn diagram, immune genes #######################
################################################################################

pcutoff <- 1e-20

glist <- glist_immune 
ggvenn::ggvenn(
  list(
    wt=de_wt$genesym[de_wt$padj < pcutoff & de_wt$genesym %in% glist],
    ko=de_ko$genesym[de_ko$padj < pcutoff & de_ko$genesym %in% glist]
  ),
  show_percentage=FALSE
) + ggtitle("GO:0006955")

ggsave("out/fig xxx immune response.pdf", width = 5, height = 5)





















################################################################################
#################### more ###################################################################
################################################################################








## TODO
# take s2e , copy, make stress response
# take s2f , copy, make stress response
# take s2f , copy, make stress response  with background as well

# take s2f , copy, make immune response  with background as well






################################################################################
######## fig xxxxx - scatter plot, immune response genes, with background ######
################################################################################

pcutoff <- 1e-20

de_comb <- merge(
  data.frame(wt_padj=de_wt$padj, wt_log2FoldChange = de_wt$log2FoldChange, genesym=de_wt$genesym),
  data.frame(ko_padj=de_ko$padj, ko_log2FoldChange = de_ko$log2FoldChange, genesym=de_ko$genesym)
)
de_comb <- de_comb[de_comb$wt_padj < pcutoff | de_comb$ko_padj < pcutoff,]
de_comb$inlist <- de_comb$genesym %in% glist_immune
de_comb <- de_comb[order(de_comb$inlist),]

ggplot(de_comb, aes(wt_log2FoldChange, ko_log2FoldChange, color=inlist)) + 
  geom_point() +
  theme_bw() +
  scale_colour_manual(values=c("lightgray","black"))+
  xlab("WT log2(FC)") + 
  ylab("shk-1 log2(FC)") +
  theme(legend.position = "none")
ggsave("out/fig zzzz scatter immune w background.pdf", width = 5, height = 5)









pcutoff <- 1e-20

de_comb <- merge(
  data.frame(wt_padj=de_wt$padj, wt_log2FoldChange = de_wt$log2FoldChange, genesym=de_wt$genesym),
  data.frame(ko_padj=de_ko$padj, ko_log2FoldChange = de_ko$log2FoldChange, genesym=de_ko$genesym)
)
de_comb <- de_comb[de_comb$wt_padj < pcutoff | de_comb$ko_padj < pcutoff,]
de_comb$inlist <- de_comb$genesym %in% glist_stress
de_comb <- de_comb[order(de_comb$inlist),]

ggplot(de_comb, aes(wt_log2FoldChange, ko_log2FoldChange, color=inlist)) + 
  geom_point() +
  theme_bw() +
  scale_colour_manual(values=c("lightgray","black"))+
  xlab("WT log2(FC)") + 
  ylab("shk-1 log2(FC)") +
  theme(legend.position = "none")
ggsave("out/fig zzzz scatter stress w background.pdf", width = 5, height = 5)






pcutoff <- 1e-20

de_comb <- merge(
  data.frame(wt_padj=de_wt$padj, wt_log2FoldChange = de_wt$log2FoldChange, genesym=de_wt$genesym),
  data.frame(ko_padj=de_ko$padj, ko_log2FoldChange = de_ko$log2FoldChange, genesym=de_ko$genesym)
)
de_comb <- de_comb[de_comb$wt_padj < pcutoff | de_comb$ko_padj < pcutoff,]
de_comb$inlist <- de_comb$genesym %in% glist_stress
de_comb <- de_comb[order(de_comb$inlist),]

ggplot(de_comb[de_comb$inlist,], aes(wt_log2FoldChange, ko_log2FoldChange)) + 
  geom_point() +
  theme_bw() +
  scale_colour_manual(values=c("black"))+
  xlab("WT log2(FC)") + 
  ylab("shk-1 log2(FC)") +
  theme(legend.position = "none")
ggsave("out/fig zzzz scatter stress no background.pdf", width = 5, height = 5)














pcutoff <- 1e-20

glist <- glist_stress
ggvenn::ggvenn(
  list(
    wt=de_wt$genesym[de_wt$padj < pcutoff & de_wt$genesym %in% glist],
    ko=de_ko$genesym[de_ko$padj < pcutoff & de_ko$genesym %in% glist]
  ),
  show_percentage=FALSE
) + ggtitle("GO:0050829")

ggsave("out/fig wwwww venn.pdf", width = 5, height = 5)






ggvenn::ggvenn(
  list(
    immune=glist_immune,
    stress=glist_stress
  ),
  show_percentage=FALSE
) 
ggsave("out/fig wwwww venn stress.pdf", width = 5, height = 5)









ggvenn::ggvenn(
  list(
    defence=glist_defence_all,
    stress=glist_stress
  ),
  show_percentage=FALSE
) 


