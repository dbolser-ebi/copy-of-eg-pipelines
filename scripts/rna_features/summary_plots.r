args = commandArgs(trailingOnly = TRUE)
for(i in 1:length(args))
{
  if (args[i] == "-l" || args[i] == "--library") {
    lib=args[i+1]
  }
}
library(optparse, lib.loc=lib)
library(plyr, lib.loc=lib)
library(reshape, lib.loc=lib)

option_list = list(
  make_option(c("-i", "--inputfile"), type="character", help="Input file summarising the alignments"),
  make_option(c("-e", "--evalue"), type="double", default=1e-06, help="E-value cut-off [default %default]"),
  make_option(c("-b", "--biotypesfile"), type="character", default="biotypes.svg", help="Output file for boxplot of total aligments per biotype [default %default]"),
  make_option(c("-d", "--distinctfile"), type="character", default="distinct.svg", help="Output file for bar chart of distinct models per species [default %default]"),
  make_option(c("-c", "--plotcolour"), type="character", default="forestgreen", help="Colour for the plots [default %default]"),
  make_option(c("-l", "--library"), type="character", default="R_lib", help="R library location")
)
parser = OptionParser(usage="%prog [options]", option_list=option_list)
options = parse_args(parser, args=args)

if (is.null(options$inputfile)) {
  print("--inputfile is required")
  quit()
}

aln_summary = read.delim(options$inputfile, header=FALSE)
col_names = c('Evalue', 'Species', 'Biotype', 'Name')
names(aln_summary) = col_names

aln_filtered = subset(aln_summary, Evalue <= options$evalue | Evalue == 0)
aln_filtered$Evalue = NULL

# Biotypes boxplot
biotypes = cast(aln_filtered, Species ~ Biotype)
biotypes$Species = NULL
biotypes[biotypes==0] = NA

svg(options$biotypesfile, width=ncol(biotypes)*0.75, height=6)
par(las=2)
par(mar=c(10, 5, 2, 2)+ 0.1)
boxplot(biotypes, log="y", ylab="Total alignments", col=c(options$plotcolour))
dev.off()


# Model frequency bar chart
distinct_names = count(aln_filtered, vars = c("Species", "Name"))
distinct_names$freq = NULL
name_freqs = count(distinct_names, vars = c("Species"))

svg(options$distinctfile, width=nrow(name_freqs)/3, height=6)
par(las=2)
par(mar=c(12, 5, 2, 2)+ 0.1)
max_y = (ceiling(max(name_freqs$freq)/50))*50
barplot(name_freqs$freq, ylab="Distinct Rfam models", ylim=c(0, max_y), names.arg=name_freqs$Species, col=c(options$plotcolour))
dev.off()
