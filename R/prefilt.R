# prefilt
#
# Filt out undesired psiFinder information

#' Filt out undesired psiFinder information
#' @param RpsiFinder_res The R result of RpsiFinder function
#' @param target_chrom The target chromsome to be retained
#' @param target_base The target base to be retained
#' @param treatStopNum_thres The filtering threshold for treatStopNum
#' @param output_dir The path to the output directory
#' @param output_name The output file name
#' @export
#'
prefilt <- function(RpsiFinder_res, target_chrom="^chr[0-9|a-z|A-Z]*$", target_base="T",treatStopNum_thres=10,output_dir,output_name)
{
  #awk 'FS=OFS="\t" {if($1~/^chr[0-9|a-z|A-Z]*$/ && $10=="T" && $14>10){print $0}}' ${outFileName}_roc_all.bed > ${outFileName}_roc_filt.bed
  prefilt_res<-RpsiFinder_res %>% filter(grepl(target_chrom,`#chrom`),base==target_base,treatStopNum>treatStopNum_thres)

  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  output_dir_output_name<-paste(output_dir,"/",output_name,sep="")

  fwrite(x=prefilt_res,file=paste(output_dir_output_name,".txt",sep=""),sep="\t")
  write.table(prefilt_res,paste(output_dir_output_name,".bed",sep=""),row.names=F,col.names=F,quote=F,sep="\t")
  return(prefilt_res)
}
