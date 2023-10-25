\name{prefilt}
\alias{prefilt}
\title{Filt out undesired psiFinder information}
\usage{
prefilt <- function(RpsiFinder_res,
                    target_chrom="^chr[0-9|a-z|A-Z]*$",
                    target_base="T",
                    treatStopNum_thres=10,
                    output_dir,
                    output_name)
}
\arguments{
  \item{RpsiFinder_res}{The R result of RpsiFinder function.}
  \item{target_chrom}{The target chromsome to be retained, default is "^chr[0-9|a-z|A-Z]*$", meaning retain those psiFinder information with normal chromsome.}
  \item{target_base}{The target base to be retained, default is "T", meaning retain those psiFinder information with T base (where pseudouridylation take place).}
  \item{treatStopNum_thres}{The filtering threshold for treatStopNum (CMC-treated reverse transcription stop number of target_base).}
  \item{output_dir}{Directory of the output file.}
  \item{output_name}{Name of the output file.}
}
\description{This function will retain psiFinder information with desired target chrom/target base and treatStopNum > treatStopNum_thres.}
\examples{
RpsiFinder_example()#when finished, view your global environment, and you will detect Data named 'RpsiFinder_example'

prefilt_out <- prefilt(RpsiFinder_res=RpsiFinder_example,
                        target_chrom="^chr[0-9|a-z|A-Z]*$",
                        target_base="T",
                        treatStopNum_thres=10,
                        output_dir="/home/chenzhr23/test_out/RpsiFinder_out",
                        output_name="test_prefilt")

}
\value{Retained desired psiFinder information (by pre-filt argument)}
