# psiROC
#
# Perform ROC evaluation based on ground truth pseudouridylation data set

#' Perform ROC evaluation based on ground truth pseudouridylation data set
#' @param rocfile ROC input file of single sites information (with suffix _roc_plot.txt) generated from ground_truth() function
#' @param rRNAfile default we recommend hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed
#' @param rRNAfile2 default we recommend hg38.psiU.SingleSites.bed
#' @param filtfile File to be filted
#' @param output_dir The path to the output directory
#' @param output_name The name to the output file
#' @export
#'
#'

#font function from ggpubr
font<-function (object, size = NULL, color = NULL, face = NULL, family = NULL,
          ...)
{
  elmt <- element_text(size = size, color = color, face = face,
                       family = family, ...)
  switch(object, title = theme(plot.title = elmt), subtitle = theme(plot.subtitle = elmt),
         caption = theme(plot.caption = elmt), x = theme(axis.title.x = elmt),
         xlab = theme(axis.title.x = elmt), x.title = theme(axis.title.x = elmt),
         y = theme(axis.title.y = elmt), ylab = theme(axis.title.y = elmt),
         y.title = theme(axis.title.y = elmt), xy = theme(axis.title.x = elmt,
                                                          axis.title.y = elmt), xylab = theme(axis.title.x = elmt,
                                                                                              axis.title.y = elmt), xy.title = theme(axis.title.x = elmt,
                                                                                                                                     axis.title.y = elmt), axis.title = theme(axis.title.x = elmt,
                                                                                                                                                                              axis.title.y = elmt), legendtitle = theme(legend.title = elmt),
         legend.title = theme(legend.title = elmt), legendtext = theme(legend.text = elmt),
         legend.text = theme(legend.text = elmt), x.text = theme(axis.text.x = elmt),
         y.text = theme(axis.text.y = elmt), xy.text = theme(axis.text.x = elmt,
                                                             axis.text.y = elmt), yxtext = theme(axis.text.x = elmt,
                                                                                                 axis.text.y = elmt), axis.text = theme(axis.text.x = elmt,
                                                                                                                                        axis.text.y = elmt), stop("Don't support ", object))
}

#.method_info function from ggpubr
.method_info<-function (method)
{
  if (is.null(method))
    method = "wilcox.test"
  allowed.methods <- list(t = "t.test", t.test = "t.test",
                          student = "t.test", wiloxon = "wilcox.test", wilcox = "wilcox.test",
                          wilcox.test = "wilcox.test", anova = "anova", aov = "anova",
                          kruskal = "kruskal.test", kruskal.test = "kruskal.test")
  method.names <- list(t.test = "T-test", wilcox.test = "Wilcoxon",
                       anova = "Anova", kruskal.test = "Kruskal-Wallis")
  if (!(method %in% names(allowed.methods)))
    stop("Non-supported method specified. Allowed methods are one of: ",
         .collapse(allowed.methods, sep = ", "))
  method <- allowed.methods[[method]]
  method.name <- method.names[[method]]
  list(method = method, name = method.name)
}

#.add_item function from ggpubr
.add_item<-function (.list, ...)
{
  pms <- list(...)
  for (pms.names in names(pms)) {
    .list[[pms.names]] <- pms[[pms.names]]
  }
  .list
}

#.is_p.signif_in_mapping function from ggpubr
.is_p.signif_in_mapping<-function (mapping)
{
  res <- FALSE
  if (!is.null(mapping)) {
    if (!is.null(mapping$label)) {
      .label <- rlang::as_label(mapping$label)
      res <- grepl(pattern = "p\\.signif", .label)
    }
  }
  return(res)
}

#.is_empty function from ggpubr
.is_empty<-function (x)
{
  length(x) == 0
}

#stat_compare_means function from ggpubr
stat_compare_means<-function (mapping = NULL, data = NULL, method = NULL, paired = FALSE,
          method.args = list(), ref.group = NULL, comparisons = NULL,
          hide.ns = FALSE, label.sep = ", ", label = NULL, label.x.npc = "left",
          label.y.npc = "top", label.x = NULL, label.y = NULL, vjust = 0,
          tip.length = 0.03, bracket.size = 0.3, step.increase = 0,
          symnum.args = list(), geom = "text", position = "identity",
          na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, ...)
{
  if (!is.null(comparisons)) {
    method.info <- .method_info(method)
    method <- method.info$method
    method.args <- .add_item(method.args, paired = paired)
    pms <- list(...)
    size <- ifelse(is.null(pms$size), 3.88, pms$size)
    color <- ifelse(is.null(pms$color), "black", pms$color)
    map_signif_level <- FALSE
    if (is.null(label))
      label <- "p.format"
    if (.is_p.signif_in_mapping(mapping) | (label %in% "p.signif")) {
      map_signif_level <- c(`****` = 1e-04, `***` = 0.001,
                            `**` = 0.01, `*` = 0.05, ns = Inf)
      if (hide.ns)
        map_signif_level <- .hide_ns(map_signif_level)
    }
    if (!.is_empty(symnum.args)) {
      symnum.args.isok <- length(symnum.args$cutpoints ==
                                   length(symnum.args$symbols))
      if (!symnum.args.isok)
        stop("Incorrect format detected in symnum.args. ",
             "Check the documentation.")
      map_signif_level <- symnum.args$cutpoints[-1]
      names(map_signif_level) <- symnum.args$symbols
      if (hide.ns)
        map_signif_level <- .hide_ns(map_signif_level)
    }
    if (missing(step.increase)) {
      step.increase <- ifelse(is.null(label.y), 0.12, 0)
    }
    ggsignif::geom_signif(comparisons = comparisons, y_position = label.y,
                          test = method, test.args = method.args, step_increase = step.increase,
                          size = bracket.size, textsize = size, color = color,
                          map_signif_level = map_signif_level, tip_length = tip.length,
                          data = data, vjust = vjust)
  }
  else {
    mapping <- .update_mapping(mapping, label)
    layer(stat = StatCompareMeans, data = data, mapping = mapping,
          geom = geom, position = position, show.legend = show.legend,
          inherit.aes = inherit.aes, params = list(label.x.npc = label.x.npc,
                                                   label.y.npc = label.y.npc, label.x = label.x,
                                                   label.y = label.y, label.sep = label.sep, method = method,
                                                   method.args = method.args, paired = paired, ref.group = ref.group,
                                                   symnum.args = symnum.args, hide.ns = hide.ns,
                                                   na.rm = na.rm, vjust = vjust, ...))
  }
}

psiROC <- function(rocfile, rRNAfile, rRNAfile2, filtfile,output_dir,output_name)
{
  #nohup Rscript $(dirname "$0")/roc.r -f ${outFileName}_roc_plot.txt -t ${outFileName}_roc_filt.bed -r $(dirname "$0")/hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed -s $(dirname "$0")/hg38.psiU.SingleSites.bed -o ${outFileName} > ${outFileName}_roc_bestthres.log 2>&1 &

  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  outfile_prefix<-paste(output_dir,"/",output_name,sep="")

  # load roc_plot.txt
  ROC_data<-as.data.frame(fread(rocfile))
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
                        "class")#40
  # ROC_data$base<-str_replace_all(as.character(ROC_data$base),"TRUE","T")
  ROC_data_sel<-ROC_data %>% select(treatPreRpmFold,treatAfterRpmFold,treatStopMeanFold,treatStopRatio,preRpmFoldRatio, afterRpmFoldRatio, stopRatioFC,treatStopMeanFoldRatio,class)


  #rRNA-psi-non-psi visualization
  ROC_data_melt<-melt(ROC_data[,c(11:38,40)],id.vars = "class")
  ROC_data_melt$class<-str_replace(as.character(ROC_data_melt$class),"0","non-psi")
  ROC_data_melt$class<-str_replace(as.character(ROC_data_melt$class),"1","psi")

  data_summary <- function(x) {
    m <- mean(x)
    ymin <- m-sd(x)
    ymax <- m+sd(x)
    return(c(y=m,ymin=ymin,ymax=ymax))
  }

  my_comparisons <- list( c("non-psi", "psi") )

  #treatPreRpmFold
  treatPreRpmFold<-ROC_data_melt%>%filter(str_detect(.$variable,"^treatPreRpmFold$"))
  treatPreRpmFold_bp <- ggplot(treatPreRpmFold, aes(x=class, y=log2(value), fill=class)) +
    stat_boxplot(geom = "errorbar",
                 width = 0.15) +
    geom_violin(trim=FALSE,alpha=0.8)+
    stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
    labs(x="group", y = "log2(treatPreRpmFold)")
  treatPreRpmFold_bp<-treatPreRpmFold_bp + scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8),legend.position = "none") + font("xy.text", size = 8)+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(log2(treatPreRpmFold$value))+abs(quantile(log2(treatPreRpmFold$value))[1]))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))

  #treatAfterRpmFold
  treatAfterRpmFold<-ROC_data_melt%>%filter(str_detect(.$variable,"^treatAfterRpmFold$"))
  treatAfterRpmFold_bp <- ggplot(treatAfterRpmFold, aes(x=class, y=log2(value), fill=class)) +
    stat_boxplot(geom = "errorbar",
                 width = 0.15) +
    geom_violin(trim=FALSE,alpha=0.8)+
    stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
    labs(x="group", y = "log2(treatAfterRpmFold)")
  treatAfterRpmFold_bp<-treatAfterRpmFold_bp + scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8), legend.position = "none") + font("xy.text", size = 8)+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(log2(treatAfterRpmFold$value))+abs(quantile(log2(treatAfterRpmFold$value))[1]))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))

  #preRpmFoldRatio
  preRpmFoldRatio<-ROC_data_melt%>%filter(str_detect(.$variable,"^preRpmFoldRatio$"))
  preRpmFoldRatio_bp <- ggplot(preRpmFoldRatio, aes(x=class, y=log2(value), fill=class)) +
    stat_boxplot(geom = "errorbar",
                 width = 0.15) +
    geom_violin(trim=FALSE,alpha=0.8)+
    stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
    labs(x="group", y = "log2(preRpmFoldRatio)")
  preRpmFoldRatio_bp<-preRpmFoldRatio_bp + scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8), legend.position = "none") + font("xy.text", size = 8)+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(log2(preRpmFoldRatio$value))+abs(quantile(log2(preRpmFoldRatio$value))[1]))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))

  #afterRpmFoldRatio
  afterRpmFoldRatio<-ROC_data_melt%>%filter(str_detect(.$variable,"^afterRpmFoldRatio$"))
  afterRpmFoldRatio_bp <- ggplot(afterRpmFoldRatio, aes(x=class, y=log2(value), fill=class)) +
    stat_boxplot(geom = "errorbar",
                 width = 0.15) +
    geom_violin(trim=FALSE,alpha=0.8)+
    stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
    labs(x="group", y = "log2(afterRpmFoldRatio)")
  afterRpmFoldRatio_bp<-afterRpmFoldRatio_bp + scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8),legend.position = "none") + font("xy.text", size = 8)+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(log2(afterRpmFoldRatio$value))+abs(quantile(log2(afterRpmFoldRatio$value))[1]))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))

  #treatStopRatio
  treatStopRatio<-ROC_data_melt%>%filter(str_detect(.$variable,"^treatStopRatio$"))
  treatStopRatio<-treatStopRatio[treatStopRatio$value!=0,]
  treatStopRatio_bp <- ggplot(treatStopRatio, aes(x=class, y=log2(value), fill=class)) +
    stat_boxplot(geom = "errorbar",
                 width = 0.15) +
    geom_violin(trim=FALSE,alpha=0.8)+
    stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
    labs(x="group", y = "log2(treatStopRatio)")
  treatStopRatio_bp<-treatStopRatio_bp + scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8),legend.position = "none") + font("xy.text", size = 8)+
    stat_compare_means(method = "t.test",comparisons = my_comparisons,label.y = max(log2(treatStopRatio$value))+abs(quantile(log2(treatStopRatio$value))[1]))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
    scale_fill_manual(values=c("#129a92","#1e95d4"))

  stopRatioFC<-ROC_data_melt%>%filter(str_detect(.$variable,"^stopRatioFC$"))
  stopRatioFC_bp <- ggplot(stopRatioFC, aes(x=class, y=log2(value), fill=class)) +
    stat_boxplot(geom = "errorbar",
                 width = 0.15) +
    geom_violin(trim=FALSE,alpha=0.8)+
    stat_summary(fun.data=data_summary,geom="crossbar", width=0.2,color="black")+
    labs(x="group", y = "log2(stopRatioFC)")
  stopRatioFC_bp<-stopRatioFC_bp + scale_fill_manual(values=c(brewer.pal(9,"Set1")[2],brewer.pal(9,"Set1")[1])) + theme_classic()+
    theme(text=element_text(size=8),legend.position = "none") + font("xy.text", size = 8)+
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


  cat("=====================treatPreRpmFold/controlPreRpmFold=========================\n")
  #treatPreRpmFold controlPreRpmFold
  pre_input=roc(ROC_data$class,ROC_data$ctrlPreRpmFold,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = FALSE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
  pdf(paste(outfile_prefix, '_roc_treatPreRpmFold.pdf', sep=""))
  pre_CMC=roc(ROC_data$class,ROC_data$treatPreRpmFold,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = TRUE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
  print(pre_CMC)
  invisible(dev.off())
  pre_input_label = paste("controlPreRpmFold AUC:",sprintf("%.3f",pre_input$auc))
  pre_CMC_label = paste("treatPreRpmFold AUC:",sprintf("%.3f",pre_CMC$auc))
  preList=list(CMC=pre_CMC, Input=pre_input)
  names(preList) <- c(pre_CMC_label,pre_input_label)
  pre_plot<-ggroc(preList,size=0.8)
  pre_plot<-pre_plot+geom_abline(intercept = 1,slope = 1,col="black",size=0.6,linetype=2)+theme_classic()+theme(legend.position="top",legend.title=element_blank(),panel.background=element_rect(fill="white",color="black"))
  pre_plot<-pre_plot+ guides(colour = guide_legend(nrow = 2))
  pre_plot<-pre_plot+ annotate("text", x = .5, y = .5,
                               label = paste("Method: ",roc.test(pre_input, pre_CMC)$method,"\np.value: ",signif(roc.test(pre_input, pre_CMC)$p.value,3),sep=""),
                               size = 3)
  print(paste(outfile_prefix,"_pre_plot.pdf",sep=""))
  pdf(paste(outfile_prefix,"_pre_plot.pdf",sep=""),width =4 ,height = 3)
  print(pre_plot)
  invisible(dev.off())
  # print(paste(outfile_prefix,"_pre_plot.png",sep=""))
  # png(paste(outfile_prefix,"_pre_plot.png",sep=""))
  # pre_plot
  # invisible(dev.off())

  roc.test(pre_input, pre_CMC)

  cat("=====================treatAfterRpmFold/controlAfterRpmFold=========================\n")
  #treatAfterRpmFold controlAfterRpmFold
  aft_input=roc(ROC_data$class,ROC_data$ctrlAfterRpmFold,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = FALSE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
  pdf(paste(outfile_prefix, '_roc_treatAfterRpmFold.pdf', sep=""))
  aft_CMC=roc(ROC_data$class,ROC_data$treatAfterRpmFold,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = TRUE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
  print(aft_CMC)
  invisible(dev.off())
  aft_input_label = paste("controlAfterRpmFold AUC:",sprintf("%.3f",aft_input$auc))
  aft_CMC_label = paste("treatAfterRpmFold AUC:",sprintf("%.3f",aft_CMC$auc))
  aftList=list(CMC=aft_CMC, Input=aft_input)
  names(aftList) <- c(aft_CMC_label,aft_input_label)
  aft_plot<-ggroc(aftList,size=0.8)
  aft_plot<-aft_plot+geom_abline(intercept = 1,slope = 1,col="black",size=0.6,linetype=2)+theme_classic()+theme(legend.position="top",legend.title=element_blank(),panel.background=element_rect(fill="white",color="black"))
  aft_plot<-aft_plot+ guides(colour = guide_legend(nrow = 2))
  aft_plot<-aft_plot+ annotate("text", x = .5, y = .5,
                               label = paste("Method: ",roc.test(aft_input, aft_CMC)$method,"\np.value: ",signif(roc.test(aft_input, aft_CMC)$p.value,3),sep=""),
                               size = 3)
  print(paste(outfile_prefix,"_aft_plot.pdf",sep=""))
  pdf(paste(outfile_prefix,"_aft_plot.pdf",sep=""),width =4 ,height = 3)
  print(aft_plot)
  invisible(dev.off())
  # print(paste(outfile_prefix,"_aft_plot.png",sep=""))
  # png(paste(outfile_prefix,"_aft_plot.png",sep=""))
  # aft_plot
  # invisible(dev.off())

  roc.test(aft_input, aft_CMC)


  cat("=====================treatStopMeanFold/controlStopMeanFold=========================\n")
  #treatStopMeanFold controlStopMeanFold
  controlStopMeanFold_input=roc(ROC_data$class,ROC_data$ctrlStopMeanFold,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = FALSE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
  pdf(paste(outfile_prefix, '_roc_treatStopMeanFold.pdf', sep=""))
  treatStopMeanFold_CMC=roc(ROC_data$class,ROC_data$treatStopMeanFold,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = TRUE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
  print(treatStopMeanFold_CMC)
  invisible(dev.off())
  controlStopMeanFold_label = paste("controlStopMeanFold AUC:",sprintf("%.3f",controlStopMeanFold_input$auc))
  treatStopMeanFold_label = paste("treatStopMeanFold AUC:",sprintf("%.3f",treatStopMeanFold_CMC$auc))
  StopMeanFoldList=list(CMC=treatStopMeanFold_CMC, Input=controlStopMeanFold_input)
  names(StopMeanFoldList) <- c(treatStopMeanFold_label,controlStopMeanFold_label)
  StopMeanFold_plot<-ggroc(StopMeanFoldList,size=0.8)
  StopMeanFold_plot<-StopMeanFold_plot+geom_abline(intercept = 1,slope = 1,col="black",size=0.6,linetype=2)+theme_classic()+theme(legend.position="top",legend.title=element_blank(),panel.background=element_rect(fill="white",color="black"))
  StopMeanFold_plot<-StopMeanFold_plot+ guides(colour = guide_legend(nrow = 2))
  StopMeanFold_plot<-StopMeanFold_plot+ annotate("text", x = .5, y = .5,
                                                 label = paste("Method: ",roc.test(controlStopMeanFold_input, treatStopMeanFold_CMC)$method,"\np.value: ",signif(roc.test(controlStopMeanFold_input, treatStopMeanFold_CMC)$p.value,3),sep=""),
                                                 size = 3)
  print(paste(outfile_prefix,"_StopMeanFold_plot.pdf",sep=""))
  pdf(paste(outfile_prefix,"_StopMeanFold_plot.pdf",sep=""),width =4 ,height = 3)
  print(StopMeanFold_plot)
  invisible(dev.off())

  # print(paste(outfile_prefix,"_StopMeanFold_plot.png",sep=""))
  # png(paste(outfile_prefix,"_StopMeanFold_plot.png",sep=""))
  # StopMeanFold_plot
  # invisible(dev.off())

  roc.test(controlStopMeanFold_input, treatStopMeanFold_CMC)


  cat("=====================treatStopRatio/controlStopRatio=========================\n")
  #treatStopRatio controlStopRatio
  stopratio_input=roc(ROC_data$class,ROC_data$ctrlStopRatio,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = FALSE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
  pdf(paste(outfile_prefix, '_roc_treatStopRatio.pdf', sep=""))
  stopratio_CMC=roc(ROC_data$class,ROC_data$treatStopRatio,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = TRUE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
  print(stopratio_CMC)
  invisible(dev.off())
  stopratio_input_label = paste("controlStopRatio AUC:",sprintf("%.3f",stopratio_input$auc))
  stopratio_CMC_label = paste("treatStopRatio AUC:",sprintf("%.3f",stopratio_CMC$auc))
  stopratioList=list(CMC=stopratio_CMC, Input=stopratio_input)
  names(stopratioList) <- c(stopratio_CMC_label,stopratio_input_label)
  stopratio_plot<-ggroc(stopratioList,size=0.8)
  stopratio_plot<-stopratio_plot+geom_abline(intercept = 1,slope = 1,col="black",size=0.6,linetype=2)+theme_classic()+theme(legend.position="top",legend.title=element_blank(),panel.background=element_rect(fill="white",color="black"))
  stopratio_plot<-stopratio_plot+ guides(colour = guide_legend(nrow = 2))
  stopratio_plot<-stopratio_plot+ annotate("text", x = .5, y = .5,
                                           label = paste("Method: ",roc.test(stopratio_input, stopratio_CMC)$method,"\np.value: ",signif(roc.test(stopratio_input, stopratio_CMC)$p.value,3),sep=""),
                                           size = 3)
  print(paste(outfile_prefix,"_stopratio_plot.pdf",sep=""))
  pdf(paste(outfile_prefix,"_stopratio_plot.pdf",sep=""),width =4 ,height = 3)
  print(stopratio_plot)
  invisible(dev.off())
  # print(paste(outfile_prefix,"_stopratio_plot.png",sep=""))
  # png(paste(outfile_prefix,"_stopratio_plot.png",sep=""))
  # stopratio_plot
  # invisible(dev.off())

  roc.test(stopratio_input, stopratio_CMC)


  cat("=====================Ratio: preRpmFoldRatio, afterRpmFoldRatio, stopRatioFC,treatStopMeanFoldRatio=========================\n")

  mycol = brewer.pal(8, "Set2")[c(1:4)]

  #ratio: preRpmFoldRatio, afterRpmFoldRatio, stopRatioFC,treatStopMeanFoldRatio
  pdf(paste(outfile_prefix, '_roc_preRpmFoldRatio.pdf', sep=""))
  preRpmFoldRatio=roc(ROC_data$class,ROC_data$preRpmFoldRatio,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = TRUE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
  print(preRpmFoldRatio)
  invisible(dev.off())
  pdf(paste(outfile_prefix, '_roc_afterRpmFoldRatio.pdf', sep=""))
  afterRpmFoldRatio=roc(ROC_data$class,ROC_data$afterRpmFoldRatio,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = TRUE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
  print(afterRpmFoldRatio)
  invisible(dev.off())
  pdf(paste(outfile_prefix, '_roc_stopRatioFC.pdf', sep=""))
  stopRatioFC=roc(ROC_data$class,ROC_data$stopRatioFC,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = TRUE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
  print(stopRatioFC)
  invisible(dev.off())
  pdf(paste(outfile_prefix, '_roc_treatStopMeanFoldRatio.pdf', sep=""))
  treatStopMeanFoldRatio=roc(ROC_data$class,ROC_data$treatStopMeanFoldRatio,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = TRUE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
  print(treatStopMeanFoldRatio)
  invisible(dev.off())
  preRpmFoldRatio_label = paste("preRpmFoldRatio AUC:",sprintf("%.3f",preRpmFoldRatio$auc))
  afterRpmFoldRatio_label = paste("afterRpmFoldRatio AUC:",sprintf("%.3f",afterRpmFoldRatio$auc))
  stopRatioFC_label = paste("stopRatioFC AUC:",sprintf("%.3f",stopRatioFC$auc))
  treatStopMeanFoldRatio_label = paste("treatStopMeanFoldRatio AUC:",sprintf("%.3f",treatStopMeanFoldRatio$auc))
  ratioList=list(preRpmFoldRatio = preRpmFoldRatio, afterRpmFoldRatio = afterRpmFoldRatio,  stopRatioFC = stopRatioFC, treatStopMeanFoldRatio = treatStopMeanFoldRatio)
  names(ratioList) <- c(preRpmFoldRatio_label, afterRpmFoldRatio_label, stopRatioFC_label, treatStopMeanFoldRatio_label)
  ratio_plot<-ggroc(ratioList,size=0.8)
  ratio_plot<-ratio_plot+geom_abline(intercept = 1,slope = 1,col="black",size=0.6,linetype=2)+theme_classic()+theme(legend.position="top",legend.title=element_blank(),panel.background=element_rect(fill="white",color="black"))
  ratio_plot<-ratio_plot+ guides(colour = guide_legend(nrow = 2, ncol = 2))+scale_color_manual(values = mycol)
  print(paste(outfile_prefix,"_ratio_plot.pdf",sep=""))
  pdf(paste(outfile_prefix,"_ratio_plot.pdf",sep=""),width =9 ,height = 9)
  ratio_plot
  print(ratio_plot)
  invisible(dev.off())
  # print(paste(outfile_prefix,"_ratio_plot.png",sep=""))
  # png(paste(outfile_prefix,"_ratio_plot.png",sep=""))
  # ratio_plot
  # invisible(dev.off())


  cat("=====================six variables: treatPreRpmFold, preRpmFoldRatio, treatAfterRpmFold, afterRpmFoldRatio, stopRatioFC,treatStopRatio=========================\n")

  mycol = brewer.pal(8, "Dark2")[c(1:6)]

  #ratio: preRpmFoldRatio, afterRpmFoldRatio, stopRatioFC,treatStopRatio
  pdf(paste(outfile_prefix, '_roc_preRpmFoldRatio.pdf', sep=""))
  preRpmFoldRatio=roc(ROC_data$class,ROC_data$preRpmFoldRatio,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = TRUE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
  print(preRpmFoldRatio)
  invisible(dev.off())
  pdf(paste(outfile_prefix, '_roc_afterRpmFoldRatio.pdf', sep=""))
  afterRpmFoldRatio=roc(ROC_data$class,ROC_data$afterRpmFoldRatio,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = TRUE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
  print(afterRpmFoldRatio)
  invisible(dev.off())
  pdf(paste(outfile_prefix, '_roc_stopRatioFC.pdf', sep=""))
  stopRatioFC=roc(ROC_data$class,ROC_data$stopRatioFC,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = TRUE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
  print(stopRatioFC)
  invisible(dev.off())
  pdf(paste(outfile_prefix, '_roc_treatStopRatio.pdf', sep=""))
  treatStopRatio=roc(ROC_data$class,ROC_data$treatStopRatio,smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = TRUE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T)
  print(treatStopRatio)
  invisible(dev.off())
  preRpmFoldRatio_label = paste("preRpmFoldRatio AUC:",sprintf("%.3f",preRpmFoldRatio$auc))
  afterRpmFoldRatio_label = paste("afterRpmFoldRatio AUC:",sprintf("%.3f",afterRpmFoldRatio$auc))
  stopRatioFC_label = paste("stopRatioFC AUC:",sprintf("%.3f",stopRatioFC$auc))
  treatStopRatio_label = paste("treatStopRatio AUC:",sprintf("%.3f",treatStopRatio$auc))
  ratioList=list(treatPreRpmFold = pre_CMC, preRpmFoldRatio = preRpmFoldRatio, treatAfterRpmFold= aft_CMC, afterRpmFoldRatio = afterRpmFoldRatio,  stopRatioFC = stopRatioFC, treatStopRatio = treatStopRatio)
  names(ratioList) <- c(pre_CMC_label, preRpmFoldRatio_label, aft_CMC_label, afterRpmFoldRatio_label, stopRatioFC_label, treatStopRatio_label)
  six_variables_ratio_plot<-ggroc(ratioList,size=0.8)
  six_variables_ratio_plot<-six_variables_ratio_plot+geom_abline(intercept = 1,slope = 1,col="black",size=0.6,linetype=2)+
    theme_classic()+
    theme(legend.position="top",
          legend.key.width= unit(0.3, 'cm'),
          legend.text = element_text(size=7),
          legend.title=element_blank(),
          panel.background=element_rect(fill="white",color="black"))
  six_variables_ratio_plot<-six_variables_ratio_plot+ guides(colour = guide_legend(nrow = 2, ncol = 3))+scale_color_manual(values = mycol)
  print(paste(outfile_prefix,"_six_variables_ratio_plot.pdf",sep=""))
  pdf(paste(outfile_prefix,"_six_variables_plot.pdf",sep=""),width =6 ,height = 6)
  six_variables_ratio_plot
  print(six_variables_ratio_plot)
  invisible(dev.off())

  pdf(paste(outfile_prefix,"_roc_summary.pdf",sep=""),width=16,height=8)
  print(plot_grid(pre_plot, aft_plot, StopMeanFold_plot, stopratio_plot, ratio_plot,six_variables_ratio_plot, align = "hv",labels = c('A', 'B','C','D','E','F')))
  invisible(dev.off())

  ROC_data_sel<-as.data.frame(ROC_data_sel)
  invisible(lapply(seq_along(ROC_data_sel[,-length(ROC_data_sel)]),function(i){
    assign(colnames(ROC_data_sel)[i],roc(ROC_data$class,ROC_data_sel[,i],smooth=FALSE, print.auc=TRUE, col="#e41a1c", plot = FALSE, print.thres="best", print.thres.best.method="youden", levels = c(0,1), direction='<',auc=T, ci=T))
    assign(paste(colnames(ROC_data_sel)[i],"_thres",sep=""),coords(get(colnames(ROC_data_sel)[i]), "best",transpose = TRUE)[1],envir = .GlobalEnv)
  }))

  cat("=====================pre_CMC=========================\n")
  # treatPreRpmFold_thres<-coords(pre_CMC, "best",transpose = TRUE)[1]
  print(coords(pre_CMC, "best", ret = "all", transpose = TRUE))

  cat("=====================aft_CMC=========================\n")
  # treatAfterRpmFold_thres<-coords(aft_CMC, "best",transpose = TRUE)[1]
  print(coords(aft_CMC, "best", ret = "all", transpose = TRUE))

  cat("=====================treatStopMeanFold_CMC=========================\n")
  # treatStopMeanFold_thres<-coords(treatStopMeanFold_CMC, "best",transpose = TRUE)[1]
  print(coords(treatStopMeanFold_CMC, "best", ret = "all", transpose = TRUE))

  cat("=====================preRpmFoldRatio=========================\n")
  # preRpmFoldRatio_thres<-coords(preRpmFoldRatio, "best",transpose = TRUE)[1]
  print(coords(preRpmFoldRatio, "best", ret = "all", transpose = TRUE))

  cat("=====================afterRpmFoldRatio=========================\n")
  # afterRpmFoldRatio_thres<-coords(afterRpmFoldRatio, "best",transpose = TRUE)[1]
  print(coords(afterRpmFoldRatio, "best", ret = "all", transpose = TRUE))

  cat("=====================treatStopRatio=========================\n")
  # treatStopRatio_thres<-coords(treatStopRatio, "best",transpose = TRUE)[1]
  print(coords(treatStopRatio, "best", ret = "all", transpose = TRUE))

  cat("=====================stopRatioFC=========================\n")
  # stopRatioFC_thres<-coords(stopRatioFC, "best",transpose = TRUE)[1]
  print(coords(stopRatioFC, "best", ret = "all", transpose = TRUE))

  cat("=====================treatStopMeanFoldRatio=========================\n")
  # treatStopMeanFoldRatio_thres<-coords(treatStopMeanFoldRatio, "best",transpose = TRUE)[1]
  print(coords(treatStopMeanFoldRatio, "best", ret = "all", transpose = TRUE))


  thres<-cbind(treatPreRpmFold_thres,treatAfterRpmFold_thres,treatStopMeanFold_thres,treatStopRatio_thres,preRpmFoldRatio_thres, afterRpmFoldRatio_thres, stopRatioFC_thres,treatStopMeanFoldRatio_thres)
  thres
  write.table(thres,paste(outfile_prefix, '_roc_bestthres.txt', sep=""),sep="\t",row.names=FALSE, col.names=FALSE,quote=F)
  write.table(thres,paste(outfile_prefix, '_roc_bestthres_colname.txt', sep=""),sep="\t",row.names=FALSE, col.names=TRUE,quote=F)

  # Values of all the confusion matrix terms were calculated at the optimal threshold
  invisible(lapply(seq_along(ROC_data_sel[,-length(ROC_data_sel)]),function(x){
    assign(paste(colnames(ROC_data_sel)[x],"_TP",sep=""),dim(ROC_data_sel[as.numeric(ROC_data_sel$class)==1 & ROC_data_sel[,x] > get(paste(colnames(ROC_data_sel)[x],"_thres",sep="")), ])[1],envir = .GlobalEnv)#True Positives (TP)
    assign(paste(colnames(ROC_data_sel)[x],"_FP",sep=""),dim(ROC_data_sel[as.numeric(ROC_data_sel$class)==0 & ROC_data_sel[,x] > get(paste(colnames(ROC_data_sel)[x],"_thres",sep="")), ])[1],envir = .GlobalEnv)#False Positives (FP)
    assign(paste(colnames(ROC_data_sel)[x],"_TN",sep=""),dim(ROC_data_sel[as.numeric(ROC_data_sel$class)==0 & ROC_data_sel[,x] <= get(paste(colnames(ROC_data_sel)[x],"_thres",sep="")), ])[1],envir = .GlobalEnv)#True Negatives (TN)
    assign(paste(colnames(ROC_data_sel)[x],"_FN",sep=""),dim(ROC_data_sel[as.numeric(ROC_data_sel$class)==1 & ROC_data_sel[,x] <= get(paste(colnames(ROC_data_sel)[x],"_thres",sep="")), ])[1],envir = .GlobalEnv)#False Negatives (FN)
    assign(paste(colnames(ROC_data_sel)[x],"_TPR",sep=""),round(get(paste(colnames(ROC_data_sel)[x],"_TP",sep=""))/(get(paste(colnames(ROC_data_sel)[x],"_TP",sep=""))+get(paste(colnames(ROC_data_sel)[x],"_FN",sep=""))),3),envir = .GlobalEnv)#sensitivity (true positive rate, TPR)
    assign(paste(colnames(ROC_data_sel)[x],"_TNR",sep=""),round(get(paste(colnames(ROC_data_sel)[x],"_TN",sep=""))/(get(paste(colnames(ROC_data_sel)[x],"_TN",sep=""))+get(paste(colnames(ROC_data_sel)[x],"_FP",sep=""))),3),envir = .GlobalEnv)#specifity (selectivity or true negative rate, TNR)
    assign(paste(colnames(ROC_data_sel)[x],"_FPR",sep=""),round(1-get(paste(colnames(ROC_data_sel)[x],"_TNR",sep="")),3),envir = .GlobalEnv)#False Positive Rate (FPR) (1 - specificit = FP/​N = FP/(TN + FP), FPR)
    assign(paste(colnames(ROC_data_sel)[x],"_FNR",sep=""),round(1-get(paste(colnames(ROC_data_sel)[x],"_TPR",sep="")),3),envir = .GlobalEnv)#False Negative Rate, FNR)
    assign(paste(colnames(ROC_data_sel)[x],"_Prec",sep=""),round(get(paste(colnames(ROC_data_sel)[x],"_TP",sep=""))/(get(paste(colnames(ROC_data_sel)[x],"_TP",sep=""))+get(paste(colnames(ROC_data_sel)[x],"_FP",sep=""))),3),envir = .GlobalEnv)#Precision
    assign(paste(colnames(ROC_data_sel)[x],"_Recall",sep=""),round(get(paste(colnames(ROC_data_sel)[x],"_TP",sep=""))/(get(paste(colnames(ROC_data_sel)[x],"_TP",sep=""))+get(paste(colnames(ROC_data_sel)[x],"_FN",sep=""))),3),envir = .GlobalEnv)#Recall
    assign(paste(colnames(ROC_data_sel)[x],"_ACC",sep=""),round((get(paste(colnames(ROC_data_sel)[x],"_TP",sep=""))+get(paste(colnames(ROC_data_sel)[x],"_TN",sep="")))/(get(paste(colnames(ROC_data_sel)[x],"_TP",sep=""))+get(paste(colnames(ROC_data_sel)[x],"_TN",sep=""))+get(paste(colnames(ROC_data_sel)[x],"_FP",sep=""))+get(paste(colnames(ROC_data_sel)[x],"_FN",sep=""))),3),envir = .GlobalEnv)#accuracy
    assign(paste(colnames(ROC_data_sel)[x],"_F1_score",sep=""),round((2*get(paste(colnames(ROC_data_sel)[x],"_Recall",sep=""))*get(paste(colnames(ROC_data_sel)[x],"_Prec",sep="")))/(get(paste(colnames(ROC_data_sel)[x],"_Recall",sep=""))+get(paste(colnames(ROC_data_sel)[x],"_Prec",sep=""))),3),envir = .GlobalEnv)#F1_score
    assign(paste(colnames(ROC_data_sel)[x],"_eval",sep=""),cbind(get(paste(colnames(ROC_data_sel)[x],"_TP",sep="")),get(paste(colnames(ROC_data_sel)[x],"_FP",sep="")),get(paste(colnames(ROC_data_sel)[x],"_TN",sep="")),get(paste(colnames(ROC_data_sel)[x],"_FN",sep="")),get(paste(colnames(ROC_data_sel)[x],"_TPR",sep="")),get(paste(colnames(ROC_data_sel)[x],"_TNR",sep="")),get(paste(colnames(ROC_data_sel)[x],"_FPR",sep="")),get(paste(colnames(ROC_data_sel)[x],"_FNR",sep="")),get(paste(colnames(ROC_data_sel)[x],"_Prec",sep="")),get(paste(colnames(ROC_data_sel)[x],"_Recall",sep="")),get(paste(colnames(ROC_data_sel)[x],"_ACC",sep="")),get(paste(colnames(ROC_data_sel)[x],"_F1_score",sep=""))),envir = .GlobalEnv)
    tmp<-as.data.frame(get(paste(colnames(ROC_data_sel)[x],"_eval",sep="")))
    names(tmp)<-c(
      paste(colnames(ROC_data_sel)[x],"_TP",sep=""),
      paste(colnames(ROC_data_sel)[x],"_FP",sep=""),
      paste(colnames(ROC_data_sel)[x],"_TN",sep=""),
      paste(colnames(ROC_data_sel)[x],"_FN",sep=""),
      paste(colnames(ROC_data_sel)[x],"_TPR",sep=""),
      paste(colnames(ROC_data_sel)[x],"_TNR",sep=""),
      paste(colnames(ROC_data_sel)[x],"_FPR",sep=""),
      paste(colnames(ROC_data_sel)[x],"_FNR",sep=""),
      paste(colnames(ROC_data_sel)[x],"_Prec",sep=""),
      paste(colnames(ROC_data_sel)[x],"_Recall",sep=""),
      paste(colnames(ROC_data_sel)[x],"_ACC",sep=""),
      paste(colnames(ROC_data_sel)[x],"_F1_score",sep="")
    )
    assign(paste(colnames(ROC_data_sel)[x],"_eval",sep=""), tmp, envir = .GlobalEnv)
  }))


  confusion_matrix_and_indicators<-as.data.frame(rbind(
    t(as.data.frame(treatPreRpmFold_eval)),
    t(as.data.frame(treatAfterRpmFold_eval)),
    t(as.data.frame(treatStopMeanFold_eval)),
    t(as.data.frame(preRpmFoldRatio_eval)),
    t(as.data.frame(afterRpmFoldRatio_eval)),
    t(as.data.frame(treatStopRatio_eval)),
    t(as.data.frame(stopRatioFC_eval)),
    t(as.data.frame(treatStopMeanFoldRatio_eval))))
  colnames(confusion_matrix_and_indicators)<-"value_or_percentage"
  confusion_matrix_and_indicators_arrange<-arrange(confusion_matrix_and_indicators,desc(value_or_percentage))
  write.table(confusion_matrix_and_indicators,paste(outfile_prefix, '_roc_confusion_matrix_and_indicators.txt', sep=""),sep="\t",row.names=TRUE, col.names=TRUE,quote=F)
  write.table(confusion_matrix_and_indicators_arrange,paste(outfile_prefix, '_roc_confusion_matrix_and_indicators_arrange.txt', sep=""),sep="\t",row.names=TRUE, col.names=TRUE,quote=F)


  best_eval<-rownames(confusion_matrix_and_indicators_arrange)[str_detect(rownames(confusion_matrix_and_indicators_arrange),"_F1_score")][1]
  assign("best_eval_t_df",as.data.frame(t(as.data.frame(get(str_replace(best_eval,"_F1_score","_eval"))))))
  colnames(best_eval_t_df)<-"value_or_percentage"

  #show afterRpmFoldRatio_eval_t_df as pdf table
  tt3 <- ttheme_minimal(
    core=list(bg_params = list(fill = blues9[1:4], col=NA),
              fg_params=list(fontface=3)),
    colhead=list(fg_params=list(col="navyblue", fontface=4L)),
    rowhead=list(fg_params=list(col="orange", fontface=3L)))

  pdf(paste(outfile_prefix, '_roc_best_evaluation.pdf', sep=""), width = 7, height = 7) # Open a new pdf file
  grid.arrange(
    tableGrob(best_eval_t_df, theme=tt3),
    nrow=1)
  invisible(dev.off()) # Close the file

  #read hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed evidence
  evidence<-read.table(rRNAfile,head=F)#"hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed"
  colnames(evidence)<-c("chrom","chromStart","chromEnd","rRNA_anno","score","strand")
  evidence$rRNA_uniq_id<-paste(evidence$chrom,evidence$chromStart,evidence$chromEnd,evidence$strand,sep="_")

  #read hg38.psiU.SingleSites.bed
  evidence2<-read.table(rRNAfile2,head=F)#"hg38.psiU.SingleSites.bed"
  colnames(evidence2)<-c("chrom","chromStart","chromEnd","rRNA_anno","score","strand")
  evidence2$rRNA_uniq_id<-paste(evidence2$chrom,evidence2$chromStart,evidence2$chromEnd,evidence2$strand,sep="_")

  #filt by stopRatioFC_thres
  to_filt<-as.data.frame(fread(filtfile,skip=1))
  colnames(to_filt)<-c("chrom",#1
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
                       "extendSeq")#39
  to_filt<-to_filt %>% filter(grepl("^chr[0-9|a-z|A-Z]*$",chrom),base=="T",treatStopNum>10)

  to_filt$roc_class<-as.numeric(to_filt[,str_replace(best_eval,"_F1_score","")]>get(str_replace(best_eval,"_F1_score","_thres")))
  to_filt$roc_class<-str_replace(as.character(to_filt$roc_class),"0","non-psi")
  to_filt$roc_class<-str_replace(as.character(to_filt$roc_class),"1","psi")
  table(to_filt$roc_class)
  write.table(to_filt,paste(outfile_prefix, '_roc_total_prediction.txt', sep=""),sep="\t",row.names=FALSE, col.names=TRUE,quote=F)
  final_pred<-to_filt[to_filt[,str_replace(best_eval,"_F1_score","")]>get(str_replace(best_eval,"_F1_score","_thres")),]
  write.table(final_pred,paste(outfile_prefix, '_roc_psi_prediction.txt', sep=""),sep="\t",row.names=FALSE, col.names=TRUE,quote=F)
  write.table(final_pred,paste(outfile_prefix, '_roc_psi_prediction.bed', sep=""),sep="\t",row.names=FALSE, col.names=FALSE,quote=F)

  #known_data miss/hit hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed
  ROC_data$rRNA_uniq_id<-paste(ROC_data$chrom,ROC_data$chromStart,ROC_data$chromEnd,ROC_data$strand,sep="_")
  ROC_data_evidence<-ROC_data %>% left_join(evidence,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% as_tibble() %>% filter(!is.na(rRNA_anno))
  write.csv(ROC_data_evidence,paste(outfile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_psiFinder_hit.csv",sep=""))
  ROC_data_no_evidence<-evidence %>% left_join(ROC_data,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(is.na(class))
  write.csv(ROC_data_no_evidence,paste(outfile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_psiFinder_miss.csv",sep=""))
  recovery<-paste(round(length(unique(ROC_data_evidence$rRNA_uniq_id))/length(unique(evidence$rRNA_uniq_id))*100,2),"%",sep="")
  cat("psiFinder recover (hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed)",recovery,"rRNA psi sites in all known chrom21\n")

  #known_data miss/hit hg38.psiU.SingleSites.bed
  ROC_data$rRNA_uniq_id<-paste(ROC_data$chrom,ROC_data$chromStart,ROC_data$chromEnd,ROC_data$strand,sep="_")
  ROC_data_evidence2<-ROC_data %>% left_join(evidence2,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% as_tibble() %>% filter(!is.na(rRNA_anno))
  write.csv(ROC_data_evidence2,paste(outfile_prefix,"_hg38.psiU.SingleSites.bed_psiFinder_hit.csv",sep=""))
  ROC_data_no_evidence2<-evidence2 %>% left_join(ROC_data,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(is.na(class))
  write.csv(ROC_data_no_evidence2,paste(outfile_prefix,"_hg38.psiU.SingleSites.bed_psiFinder_miss.csv",sep=""))
  recovery<-paste(round(length(unique(ROC_data_evidence2$rRNA_uniq_id))/length(unique(evidence2$rRNA_uniq_id))*100,2),"%",sep="")
  cat("psiFinder recover (hg38.psiU.SingleSites.bed.bed)",recovery,"rRNA psi sites in all known chrom21\n")

  #final_pred miss/hit
  final_pred$rRNA_uniq_id<-paste(final_pred$chrom,final_pred$chromStart,final_pred$chromEnd,final_pred$strand,sep="_")
  final_pred_evidence<-final_pred %>% left_join(evidence,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% as_tibble() %>% filter(!is.na(rRNA_anno))
  write.csv(final_pred_evidence,paste(outfile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_psiFinder_roc_stopRatioFC_thres_hit.csv",sep=""))
  final_pred_no_evidence<-evidence %>% left_join(final_pred,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(is.na(base))
  write.csv(final_pred_no_evidence,paste(outfile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_psiFinder_roc_stopRatioFC_thres_miss.csv",sep=""))
  recovery<-paste(round(length(unique(final_pred_evidence$rRNA_uniq_id))/length(unique(evidence$rRNA_uniq_id))*100,2),"%",sep="")
  cat("psiFinder+ROC recover (hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed)",recovery,"rRNA psi sites in all known chrom21\n")

  print(paste("ROC evaluation in",output_dir,sep=" "))
}
