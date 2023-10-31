\name{udfilt}
\alias{udfilt}
\title{Filt out undesired psiFinder information}
\usage{
udfilt <- function(RpsiFinder_res_file,
                   treatPreRpmFold_thres,
                   treatAfterRpmFold_thres,
                   treatStopRatio_thres,
                   preRpmFoldRatio_thres,
                   afterRpmFoldRatio_thres,
                   stopRatioFC_thres,
                   output_dir,
                   output_name)
}
\arguments{
  \item{RpsiFinder_res_file}{The file path of psiFidner result (RpsiFinder output file).}
  \item{treatPreRpmFold_thres}{The filtering threshold for treatPreRpmFold.}
  \item{treatAfterRpmFold_thres}{The filtering threshold for treatAfterRpmFold.}
  \item{treatStopRatio_thres}{The filtering threshold for treatStopRatio.}
  \item{preRpmFoldRatio_thres}{The filtering threshold for preRpmFoldRatio.}
  \item{afterRpmFoldRatio_thres}{The filtering threshold for afterRpmFoldRatio.}
  \item{stopRatioFC_thres}{The filtering threshold for stopRatioFC.}
  \item{output_dir}{Directory of the output file.}
  \item{output_name}{Name of the output file.}
}
\description{This function will retain psiFinder information with desired identification metrics greater than thresholds (treatPreRpmFold, treatAfterRpmFold, treatStopRatio, preRpmFoldRatio, afterRpmFoldRatio, stopRatioFC).}
\examples{
RpsiFinder_example()#when finished, view your global environment, and you will detect Data named 'RpsiFinder_example'

prefilt_out <- prefilt(RpsiFinder_res=RpsiFinder_example,
                        target_chrom="^chr[0-9|a-z|A-Z]*$",
                        target_base="T",
                        treatStopNum_thres=10,
                        output_dir="../test_out/RpsiFinder_out",
                        output_name="test_prefilt")

udfilt_out <- udfilt(RpsiFinder_res_file="../test_out/RpsiFinder_out/test_prefilt.txt",
                   treatPreRpmFold_thres=1.5,
                   treatAfterRpmFold_thres=1.5,
                   treatStopRatio_thres=0.04,
                   preRpmFoldRatio_thres=1.5,
                   afterRpmFoldRatio_thres=1.5,
                   stopRatioFC_thres=1.5,
                   output_dir="../test_out/RpsiFinder_out",
                   output_name="test_udfilt")

}
\value{Retained desired psiFinder information (by user-defined argument)}
