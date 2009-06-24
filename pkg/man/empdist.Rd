\name{empdist}
\Rdversion{1.1}
\alias{empdist}
\title{
Function that uses permutation to obtain empirical null distributions that can in turn be used to identify which parameters differ between populations.
}
\description{
This function identifies the best-fit parameters for a pair of samples and an a specified joint model or models, then permutes data points between the samples to obtain empirical null distributions for those models in order to determine whether the respective parameters of the two samples are different from each other. This function can also use resampling with replacement to obtain the distributions of fitted parameters, MLEs, and any other statistic calculated by \code{findpars}.
}
\usage{
empdist(x, y, label, nil = 0, bnil = 0, pf = mean, rounds = 5000, type = "p", models = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A vector of ages at death (integers, usually representing days).
}
  \item{y}{
A second vector of ages at death (integers, usually representing days). In contrast with \code{findpars} and \code{simdist} this is a required argument.
}
  \item{label}{
A character string prefix based on which the permutation output files are named. You may specify a file-path in the prefix, but do not specify a file extension. The output file will be saved as 'LABEL.txt' where LABEL is replaced by the value of this argument. The permuted/resampled data will likewise be saved as either 'LABEL_perms.txt' or 'LABEL_resam.txt' depending on whether the \code{type} argument is 'p' (permutation), or 'r' (sampling with replacement). This is a required argument.
}
  \item{nil}{
Passed to \code{findpars}, see the help for that function. Normally does not need to be set.
}
  \item{bnil}{
Passed to \code{findpars}, see the help for that function. Normally does not need to be set.
}
  \item{pf}{
A function to call when two different parameters are constrained and need to produce a single starting value. In addition to \code{mean} (default), other valid functions include \code{median}, \code{gmean}, \code{max}, and \code{min}.
}
  \item{rounds}{
An integer indicating how many permutations or resamplings to perform. Default is 5000.
}
  \item{type}{
A character which is either 'p' (default), indicating that data should be permuted between the two samples or 'r', indicating that data should be resampled within samples.
}
  \item{models}{
A character vector of models to try which can be any combination of: 'w','g','gm','l', or 'lm'. In some cases, models that are not listed in this argument are fitted anyway because they are needed to obtain starting parameters for the models that are listed.
}
}
\details{
This function has two practical uses at the moment. The first is obtaining empirical null distributions for the log-ratios of unconstrained model fits versus model fits constrained in each parameter in turn. This is accomplished by running \code{findpars} on two samples equal in length to the original data but with the datapoints randomly reassigned between them and repeating this process \code{rounds} times. The above is what happens when the \code{type} argument is equal to 'p', which it is by default. If the \code{type} argument is instead set to 'r', then the samples are resampled with replacement \code{rounds} times from the original and analyzed by \code{findpars}. The purpose of this is to obtain distributions for estimates of individual parameters, MLEs, or any other statistic that \code{findpars} calculates.

In either case, two tab-delimited text files are saved: one containing the permuted or resampled data, ending in '_perms.txt' or '_resam.txt' respectively, and the other having a name identical to the \code{label} argument.
}

\value{
This function does not return anything to the console, but instead saves its output as two files. The first is named 'LABEL.txt' and has identical columns to the output from \code{findpars}, but each value in the 'id' column corresponds the column (from the left) in the permuted/resampled data file that contains the data which produced that given row. 

The second file is named either 'LABEL_perms.txt' or 'LABEL_resam.txt' depending on whether \code{type} was set to 'p' (or omitted) or set to 'r'. In the former case, the file is just a grid of zeroes and ones, which indicate a datapoint's membership in the control or experimental group, respectively, for that permutation. In the latter case, the file is a grid of resampled data points, where the first N are in the control group and the rest are in the experimental group (the value of N and the total sample size is obtained from the \code{x} and \code{y} arguments and not stored in the output file).
}
\references{
Pletcher,S.D., Khazaeli,A.A., and Curtsinger,J.W. (2000). Why do life spans differ? Partitioning mean longevity differences in terms of age-specific mortality parameters. Journals of Gerontology Series A-Biological Sciences and Medical Sciences 55, B381-B389
}
\note{
Even though the examples use 'rounds=100' for demonstration purposes, for publication it is recommended to leave the 'rounds' argument at its default of 5000 or give it a larger value. This may take several hours or even days depending on the dataset. This is why the output is saved to a text file rather than an object within R-- the user can monitor the progress of this function by accessing the data file with \code{\link{constrshow}} from another R session.

Currently this function spams the console with periods and abbreviated model names. This behavior is intentional, for troubleshooting purposes and to indicate that R has not crashed or hung.
}
\author{
Alex F. Bokov
}
\seealso{
\code{\link{findpars}}, \code{\link{constrshow}}
}
\Warning{
These examples may take 10 minutes or more to finish running.
}
\examples{
data(demo1);
data(demo2);
empdist(demo1,demo2,'test',rounds=100);
empdist(demo1,demo2,'resample_test',type='r',rounds=100);
test<-read.table('test.txt',header=1,sep='\t',check.names=F);
test[1:10,];
resample_test<-read.table('resample_test.txt',header=1,sep='\t',check.names=F);
resample_test[1:10,];
}
\keyword{ survival }
