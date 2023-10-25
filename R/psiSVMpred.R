# psiSVMpred
#
# Build SVM model based on ground truth pseudouridylation data set

#' Build SVM model based on ground truth pseudouridylation data set
#' @param filtfile File to be filted (File ready to be predicted)
#' @param svmmodelfile SVM model file
#' @param rRNAfile default we recommend hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed
#' @param output_dir The path to the output directory
#' @param output_name The name to the output file
#' @export
#'
#'

psiSVMpred <- function(filtfile, svmmodelfile, rRNAfile, output_dir,output_name)
{
  #nohup Rscript $(dirname "$0")/svm_predict_totalRNA.r -f ${input%.bed}_svm_filt.bed -k $SVM_model -r $(dirname "$0")/hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed -s $(dirname "$0")/hg38.psiU.SingleSites.bed  -o ${input%.bed} > ${input%.bed}.log

  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  outfile_prefix<-paste(output_dir,"/",output_name,sep="")

  # mymodel<-readRDS(file = svmmodelfile)#"my-svm.RData"
  load(file = svmmodelfile)#"my-svm.RData"

  #filt by svm best model
  cat("\n\n=====================Filt by svm best model=====================\n")
  cat("below is your input data ready to be predicted...\n")
  # get prediction
  to_pred<-as.data.frame(fread(filtfile, skip=1))
  colnames(to_pred)<-c("chrom",#1
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
  to_pred_var<-to_pred %>% select(treatPreRpmFold, preRpmFoldRatio, treatAfterRpmFold, afterRpmFoldRatio, treatStopRatio, stopRatioFC)
  print(str(to_pred_var))

  pred <- predict(mymodel,to_pred_var,decision.values=TRUE,probability=TRUE)

  #get attr
  pred.decision.values<-as.vector(attr(pred, "decision.values"))
  pred.probabilities<-attr(pred, "probabilities")
  pred.probabilities<-as.data.frame(pred.probabilities)

  #add attr to SVM_pred_data
  SVM_pred_data<-to_pred
  SVM_pred_data$pred.decision.values<-pred.decision.values
  SVM_pred_data$svm_pred_desc_class<-ifelse(SVM_pred_data$pred.decision.values>0,"psi","non-psi")
  SVM_pred_data<-cbind(SVM_pred_data,pred.probabilities,pred)
  SVM_pred_data$svm_pred_prob_class<-ifelse(SVM_pred_data$psi>=0.5,"psi","non-psi")
  print("Psi predicted by SVM model:")
  print(table(SVM_pred_data$svm_pred_prob_class))

  write.xlsx(SVM_pred_data,paste(outfile_prefix, '_svm_total_prediction.xlsx', sep=""),overwrite = TRUE)
  write.table(SVM_pred_data,paste(outfile_prefix, '_svm_total_prediction.txt', sep=""),sep="\t",row.names=FALSE, col.names=TRUE,quote=F)
  write.table(SVM_pred_data,paste(outfile_prefix, '_svm_total_prediction.bed', sep=""),sep="\t",row.names=FALSE, col.names=FALSE,quote=F)

  #known evidence
  evidence<-read.table(rRNAfile,head=F)#"hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed"
  colnames(evidence)<-c("chrom","chromStart","chromEnd","rRNA_anno","score","strand")
  evidence$rRNA_uniq_id<-paste(evidence$chrom,evidence$chromStart,evidence$chromEnd,evidence$strand,sep="_")

  final_pred<-SVM_pred_data[SVM_pred_data$svm_pred_prob_class=="psi",]
  final_pred <- final_pred[final_pred$foldChange > 2,]
  final_pred<-final_pred %>% arrange(desc(pred.decision.values))
  print("Final psi prediction (stopRatioFC >2):")
  print(table(final_pred$svm_pred_prob_class))

  write.xlsx(final_pred, paste(outfile_prefix, '_svm_psi_prediction.xlsx', sep = ""), overwrite = TRUE)
  write.table(final_pred,paste(outfile_prefix, '_svm_psi_prediction.txt', sep=""),sep="\t",row.names=FALSE, col.names=TRUE,quote=F)
  write.table(final_pred,paste(outfile_prefix, '_svm_psi_prediction.bed', sep=""),sep="\t",row.names=FALSE, col.names=FALSE,quote=F)

  final_pred$rRNA_uniq_id<-paste(final_pred$chrom,final_pred$chromStart,final_pred$chromEnd,final_pred$strand,sep="_")
  final_pred_evidence<-final_pred %>% left_join(evidence,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% as_tibble() %>% filter(!is.na(rRNA_anno))
  write.csv(final_pred_evidence,paste(outfile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_psiFinder_svm_hit.csv",sep=""))
  final_pred_no_evidence<-evidence %>% left_join(final_pred,by=c("rRNA_uniq_id"="rRNA_uniq_id")) %>% filter(is.na(base))
  write.csv(final_pred_no_evidence,paste(outfile_prefix,"_hg38_human_chr21_rRNA_known_pseudoU_SingleSites_psiFinder_svm_miss.csv",sep=""))
  recovery<-paste(round(length(unique(final_pred_evidence$rRNA_uniq_id))/length(unique(evidence$rRNA_uniq_id))*100,2),"%",sep="")
  cat("psiFinder+SVM recover (hg38_human_chr21_rRNA_known_pseudoU_SingleSites.bed)",recovery,"rRNA psi sites in all known chrom21\n")

  print(paste("SVM prediction in",output_dir,sep=" "))
}
