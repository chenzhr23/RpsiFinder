# RpsiFinder
#
# This function run psiFinder and return overall reverse transcription stop information

#' This function run psiFinder and return overall reverse transcription stop information
#' @param treat_bam The path to the treat bam file
#' @param input_bam The path to the input bam file
#' @param gene_bed12 The path to the gene annotation bed12 file
#' @param output_dir The path to the output directory
#' @param output_name The name to the output file
#' @param options Additional options to pass to the program
#' @export
#'
RpsiFinder <- function(genome_fa,
                       genome_fai,
                       treat_bam,
                       input_bam,
                       gene_bed12,
                       output_dir,
                       output_name,
                       options = "-p 1.5 -t 5 -r 0.05 -M 1 -f 1 -m 0 -s -n -w 20")
  {
  psiFinder_path <- system.file("program", "psiFinder", package = "RpsiFinder")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  output_dir_output_name<-paste(output_dir,"/",output_name,sep="")
  cmd <- paste(psiFinder_path,
               "--fa", genome_fa,
               "--fai", genome_fai,
               "--treat", treat_bam,
               "--input", input_bam,
               options,
               "--gene", gene_bed12,
               "-o", paste(output_dir_output_name,".txt",sep=""),
               paste("2>", output_dir, "/", output_name,".log",sep=""),
               sep=" ")

  print(cmd)
  system(cmd)
  assign(output_name,fread(output_dir_output_name),envir = .GlobalEnv)
  }
