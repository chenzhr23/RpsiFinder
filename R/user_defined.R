# user_defined
#
# Retain reverse transcription stop information based on user-defined thresholds

#' Retain reverse transcription stop information based on user-defined thresholds
#' @param rocfile ROC input file of single sites information (with suffix _roc_plot.txt) generated from ground_truth() function
#' @param rRNAfile default we recommend hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed
#' @param filtfile File to be filted
#' @param output_dir The path to the output directory
#' @param output_name The name to the output file
#' @export
#'
#'

user_defined <- function(rocfile, rRNAfile, filtfile, output_dir, output_name)
{
  #nohup Rscript $(dirname "$0")/user_defined_evaluation.r -f ${output_path}_roc_plot.txt -r $(dirname "$0")/hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed -s $(dirname "$0")/hg38.psiU.SingleSites.bed -t ${output_path}_user_defined_prediction_total.bed -o ${output_path} > ${output_path}_roc_bestthres.log 2>&1 &

  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  outfile_prefix<-paste(output_dir,"/",output_name,sep="")

  # roc_plot.txt
  ROC_data<-fread(rocfile)
  colnames(ROC_data)<-c("chrom",#1
                        "chromStart",#2
                        "chromEnd",#3
                        "name",#4
                        "foldChange",#5
                        "strand",#6
                        "geneName",#7
                        "geneStart",#8
                        "geneEnd",#9
                        "base",#10
                        "treatPval",#11
                        "ctrlPval",#12
                        "minusPval",#13
                        "treatStopNum",#14
                        "treatStopRPM",#15
                        "treatPreStopNum",#16
                        "treatAfterStopNum",#17
                        "treatReadthroughNum",#18
                        "ctrlStopNum",#19
                        "ctrlStopRPM",#20
                        "ctrlPreStopNum",#21
                        "ctrlAfterStopNum",#22
                        "ctrlReadthroughNum",#23
                        "stopRpmFC",#24
                        "treatPreRpmFold",#25
                        "ctrlPreRpmFold",#26
                        "preRpmFoldRatio",#27
                        "treatAfterRpmFold",#28
                        "ctrlAfterRpmFold",#29
                        "afterRpmFoldRatio",#30
                        "treatStopRatio",#31
                        "ctrlStopRatio",#32
                        "stopRatioFC",#33
                        "treatStopMeanNum",#34
                        "treatStopMeanFold",#35
                        "ctrlStopMeanNum",#36
                        "ctrlStopMeanFold",#37
                        "treatStopMeanFoldRatio",#38
                        "extendSeq",#39
                        "pred_class",#40
                        "real_class")#41
  print("real class:\n")
  print(table(ROC_data$real_class))
  print("prediction class:\n")
  print(table(ROC_data$pred_class))

  tab <- table(Predicted = ROC_data$pred_class,Actual = ROC_data$real_class)
  print(tab)

  #Visualize the prediction effect of the model
  plot_confusion_matrix <- ggplot() +
    geom_confmat(aes(x = ROC_data$real_class, y = ROC_data$pred_class),
                 normalize = TRUE, text.perc = TRUE) +
    labs(x = "Reference", y = "Prediction") +
    scale_fill_gradient2(low = "darkblue", high = "#ec2f2f") +
    theme_bw() +
    theme(plot.margin = unit(c(6, 5, 6, 5), "cm"))

  pdf(paste(outfile_prefix, '_user_defined_confusion_matrix.pdf', sep = ""))
  print(plot_confusion_matrix)
  invisible(dev.off())

  # print(confusionMatrix(factor(ROC_data$pred_class,levels=c("0","1")),factor(ROC_data$real_class,levels=c("0","1"))))
  ROC_data_sel<-ROC_data %>% select(treatPreRpmFold,treatAfterRpmFold,treatStopMeanFold,treatStopRatio,preRpmFoldRatio, afterRpmFoldRatio, stopRatioFC,treatStopMeanFoldRatio,real_class)
  print("total summary (rRNA psi and rRNA non-psi)")
  print(summary(ROC_data_sel))
  print("psi summary (all real rRNA psi)")
  real_rRNA_psi<-ROC_data_sel %>% filter(real_class=="1")
  print(summary(real_rRNA_psi))

  #rRNA-psi-non-psi visualization
  ROC_data_melt<-melt(ROC_data[,c(11:38,41)],id.vars = "real_class")
  ROC_data_melt$real_class<-str_replace(as.character(ROC_data_melt$real_class),"0","non-psi")
  ROC_data_melt$real_class<-str_replace(as.character(ROC_data_melt$real_class),"1","psi")

  data_summary <- function(x) {
    m <- mean(x)
    ymin <- m-sd(x)
    ymax <- m+sd(x)
    return(c(y=m,ymin=ymin,ymax=ymax))
  }

  my_comparisons <- list( c("non-psi", "psi") )

  #treatPreRpmFold
  treatPreRpmFold<-ROC_data_melt%>%filter(str_detect(.$variable,"^treatPreRpmFold$"))
  treatPreRpmFold_bp <- ggplot(treatPreRpmFold, aes(x=real_class, y=log2(value), fill=real_class)) +
    stat_boxplot(geom = "errorbar",
                 width = 0.15) +
    # geom_boxplot(outlier.colour = NA,width=0.2)+
    geom_violin(trim=FALSE,alpha=0.8)+
    stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
    labs(x="group", y = "log2(treatPreRpmFold)")
  treatPreRpmFold_bp<-treatPreRpmFold_bp + scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8),legend.position = "none") + font("xy.text", size = 8)+
    # stat_compare_means(method = "t.test",comparisons = my_comparisons,label = "p.signif")+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(log2(treatPreRpmFold$value))+abs(quantile(log2(treatPreRpmFold$value))[1]))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))

  #treatAfterRpmFold
  treatAfterRpmFold<-ROC_data_melt%>%filter(str_detect(.$variable,"^treatAfterRpmFold$"))
  treatAfterRpmFold_bp <- ggplot(treatAfterRpmFold, aes(x=real_class, y=log2(value), fill=real_class)) +
    stat_boxplot(geom = "errorbar",
                 width = 0.15) +
    geom_violin(trim=FALSE,alpha=0.8)+
    stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
    labs(x="group", y = "log2(treatAfterRpmFold)")
  treatAfterRpmFold_bp<-treatAfterRpmFold_bp + scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8), legend.position = "none") + font("xy.text", size = 8)+
    # stat_compare_means(method = "t.test",comparisons = my_comparisons,label = "p.signif")+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(log2(treatAfterRpmFold$value))+abs(quantile(log2(treatAfterRpmFold$value))[1]))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))

  #preRpmFoldRatio
  preRpmFoldRatio<-ROC_data_melt%>%filter(str_detect(.$variable,"^preRpmFoldRatio$"))
  preRpmFoldRatio_bp <- ggplot(preRpmFoldRatio, aes(x=real_class, y=log2(value), fill=real_class)) +
    stat_boxplot(geom = "errorbar",
                 width = 0.15) +
    geom_violin(trim=FALSE,alpha=0.8)+
    stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
    labs(x="group", y = "log2(preRpmFoldRatio)")
  preRpmFoldRatio_bp<-preRpmFoldRatio_bp + scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    # theme(text=element_text(size=8), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),legend.position = "none") + font("xy.text", size = 8)+
    theme(text=element_text(size=8), legend.position = "none") + font("xy.text", size = 8)+
    # stat_compare_means(method = "t.test",comparisons = my_comparisons,label = "p.signif")+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(log2(preRpmFoldRatio$value))+abs(quantile(log2(preRpmFoldRatio$value))[1]))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))

  #afterRpmFoldRatio
  afterRpmFoldRatio<-ROC_data_melt%>%filter(str_detect(.$variable,"^afterRpmFoldRatio$"))
  afterRpmFoldRatio_bp <- ggplot(afterRpmFoldRatio, aes(x=real_class, y=log2(value), fill=real_class)) +
    stat_boxplot(geom = "errorbar",
                 width = 0.15) +
    geom_violin(trim=FALSE,alpha=0.8)+
    stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
    labs(x="group", y = "log2(afterRpmFoldRatio)")
  afterRpmFoldRatio_bp<-afterRpmFoldRatio_bp + scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8),legend.position = "none") + font("xy.text", size = 8)+
    # stat_compare_means(method = "t.test",comparisons = my_comparisons,label = "p.signif")+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(log2(afterRpmFoldRatio$value))+abs(quantile(log2(afterRpmFoldRatio$value))[1]))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))

  #treatStopRatio
  treatStopRatio<-ROC_data_melt%>%filter(str_detect(.$variable,"^treatStopRatio$"))
  treatStopRatio<-treatStopRatio[treatStopRatio$value!=0,]
  treatStopRatio_bp <- ggplot(treatStopRatio, aes(x=real_class, y=log2(value), fill=real_class)) +
    stat_boxplot(geom = "errorbar",
                 width = 0.15) +
    geom_violin(trim=FALSE,alpha=0.8)+
    stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
    labs(x="group", y = "log2(treatStopRatio)")
  treatStopRatio_bp<-treatStopRatio_bp + scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8),legend.position = "none") + font("xy.text", size = 8)+
    # stat_compare_means(method = "t.test",comparisons = my_comparisons,label = "p.signif")+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(log2(treatStopRatio$value))+abs(quantile(log2(treatStopRatio$value))[1]))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))

  stopRatioFC<-ROC_data_melt%>%filter(str_detect(.$variable,"^stopRatioFC$"))
  stopRatioFC_bp <- ggplot(stopRatioFC, aes(x=real_class, y=log2(value), fill=real_class)) +
    stat_boxplot(geom = "errorbar",
                 width = 0.15) +
    geom_violin(trim=FALSE,alpha=0.8)+
    stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
    labs(x="group", y = "log2(stopRatioFC)")
  stopRatioFC_bp<-stopRatioFC_bp + scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8),legend.position = "none") + font("xy.text", size = 8)+
    # stat_compare_means(method = "t.test",comparisons = my_comparisons,label = "p.signif")+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(log2(stopRatioFC$value))+abs(quantile(log2(stopRatioFC$value))[1]))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))


  pdf(paste(outfile_prefix,"_six_variables_rRNA_violinplot.pdf",sep=""),width=5,height=6)
  print(plot_grid(
    treatPreRpmFold_bp,
    preRpmFoldRatio_bp,
    treatAfterRpmFold_bp,
    afterRpmFoldRatio_bp,
    treatStopRatio_bp,
    stopRatioFC_bp,
    align = "hv",
    labels = c('A','B','C','D','E','F'),ncol=2,nrow=3))
  invisible(dev.off())

  #calculate evaluation indicators
  cat("\n\n=====================Calculate evaluation indicators=====================\n")
  confusion_matrix<-as.data.frame(tab)
  confusion_matrix$Predicted<-str_replace(confusion_matrix$Predicted,"1","psi")
  confusion_matrix$Predicted<-str_replace(confusion_matrix$Predicted,"0","non-psi")
  confusion_matrix$Actual<-str_replace(confusion_matrix$Actual,"1","psi")
  confusion_matrix$Actual<-str_replace(confusion_matrix$Actual,"0","non-psi")
  ud_TP <- confusion_matrix[confusion_matrix$Predicted=="psi"&confusion_matrix$Actual=="psi",]$Freq#True Positives (TP)
  ud_FP <- confusion_matrix[confusion_matrix$Predicted=="psi"&confusion_matrix$Actual=="non-psi",]$Freq#False Positives (FP)
  ud_TN <- confusion_matrix[confusion_matrix$Predicted=="non-psi"&confusion_matrix$Actual=="non-psi",]$Freq#True Negatives (TN)
  ud_FN <- confusion_matrix[confusion_matrix$Predicted=="non-psi"&confusion_matrix$Actual=="psi",]$Freq#False Negatives (FN)
  ud_TPR <- ud_TP / (ud_TP + ud_FN)#sensitivity (true positive rate, TPR)
  ud_TNR <- ud_TN / (ud_TN + ud_FP)#specifity (selectivity or true negative rate, TNR)
  ud_FPR <- 1-ud_TNR#False Positive Rate (FPR) (1 - specificit = FP/​N = FP/(TN + FP), FPR)
  ud_FNR <- 1-ud_TPR#False Negative Rate, FNR)
  ud_Prec <- ud_TP / (ud_TP + ud_FP)#Precision
  ud_Recall <- ud_TP / (ud_TP + ud_FN)#Recall
  ud_ACC <- (ud_TP + ud_TN) / (ud_TP + ud_TN + ud_FP + ud_FN)#accuracy
  ud_F1_score <- (2*ud_Recall*ud_Prec) / (ud_Recall + ud_Prec)#F1_score
  eval<-cbind(ud_TP,ud_FP,ud_TN,ud_FN,ud_TPR,ud_TNR,ud_FPR,ud_FNR,ud_Prec,ud_Recall,ud_ACC,ud_F1_score)
  eval<-round(eval,3)
  print(eval)
  write.table(eval,paste(outfile_prefix, '_ud_eval.txt', sep=""),sep="\t",row.names=FALSE, col.names=TRUE,quote=F)

  #show ud evaluation as pdf table
  ud_eval_t_df<-as.data.frame(t(as.data.frame(eval)))
  colnames(ud_eval_t_df)<-"value_or_percentage"
  tt3 <- ttheme_minimal(
    core=list(bg_params = list(fill = blues9[1:4], col=NA),
              fg_params=list(fontface=3)),
    colhead=list(fg_params=list(col="navyblue", fontface=4L)),
    rowhead=list(fg_params=list(col="orange", fontface=3L)))

  pdf(paste(outfile_prefix, '_user_defined_evaluation.pdf', sep=""), width = 7, height = 7) # Open a new pdf file
  print(grid.arrange(
    tableGrob(ud_eval_t_df, theme=tt3),
    nrow=1))
  invisible(dev.off()) # Close the file

  #filt by best F1 score threshold
  to_filt<-fread(filtfile)
  colnames(to_filt)<-c(
    "chrom",#1
    "chromStart",#2
    "chromEnd",#3
    "name",#4
    "foldChange",#5
    "strand",#6
    "geneName",#7
    "geneStart",#8
    "geneEnd",#9
    "base",#10
    "treatPval",#11
    "ctrlPval",#12
    "minusPval",#13
    "treatStopNum",#14
    "treatStopRPM",#15
    "treatPreStopNum",#16
    "treatAfterStopNum",#17
    "treatReadthroughNum",#18
    "ctrlStopNum",#19
    "ctrlStopRPM",#20
    "ctrlPreStopNum",#21
    "ctrlAfterStopNum",#22
    "ctrlReadthroughNum",#23
    "stopRpmFC",#24
    "treatPreRpmFold",#25
    "ctrlPreRpmFold",#26
    "preRpmFoldRatio",#27
    "treatAfterRpmFold",#28
    "ctrlAfterRpmFold",#29
    "afterRpmFoldRatio",#30
    "treatStopRatio",#31
    "ctrlStopRatio",#32
    "stopRatioFC",#33
    "treatStopMeanNum",#34
    "treatStopMeanFold",#35
    "ctrlStopMeanNum",#36
    "ctrlStopMeanFold",#37
    "treatStopMeanFoldRatio",#38
    "extendSeq",#39
    "pred_class")#40
  to_filt$base<-str_replace_all(to_filt$base, "TRUE", "T")
  to_filt$pred_class<-str_replace(as.character(to_filt$pred_class),"0","non-psi")
  to_filt$pred_class<-str_replace(as.character(to_filt$pred_class),"1","psi")

  print(table(to_filt$pred_class))
  write.table(to_filt,paste(outfile_prefix, '_user_defined_total_prediction.txt', sep=""),sep="\t",row.names=FALSE, col.names=TRUE,quote=F)
  final_pred<-to_filt %>% filter(pred_class=="psi")
  write.table(final_pred,paste(outfile_prefix, '_ud_psi_prediction.txt', sep=""),sep="\t",row.names=FALSE, col.names=TRUE,quote=F)
  write.table(final_pred,paste(outfile_prefix, '_ud_psi_prediction.bed', sep=""),sep="\t",row.names=FALSE, col.names=FALSE,quote=F)

  #read hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed evidence
  evidence<-read.table(rRNAfile,head=F)
  colnames(evidence)<-c("chrom","chromStart","chromEnd","rRNA_anno","score","strand")
  evidence$rRNA_uniq_id<-paste(evidence$chrom,evidence$chromStart,evidence$chromEnd,evidence$strand,sep="_")

  #final_pred miss/hit
  final_pred$rRNA_uniq_id<-paste(final_pred$chrom,final_pred$chromStart,final_pred$chromEnd,final_pred$strand,sep="_")
  final_pred_evidence<-final_pred %>% left_join(evidence,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(!is.na(rRNA_anno))
  write.csv(final_pred_evidence,paste(outfile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_psiFinder_ud_thres_hit.csv",sep=""))
  final_pred_no_evidence<-evidence %>% left_join(final_pred,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(is.na(base))
  write.csv(final_pred_no_evidence,paste(outfile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_psiFinder_ud_thres_miss.csv",sep=""))
  recovery<-paste(round(length(unique(final_pred_evidence$rRNA_uniq_id))/length(unique(evidence$rRNA_uniq_id))*100,2),"%",sep="")
  cat("psiFinder+ud recover (hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed)",recovery,"rRNA psi sites in all known chrom21\n")

  print(paste("User-defined result in",output_dir,sep=" "))
}
