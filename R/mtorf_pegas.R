data <- read.dna("~/Downloads/poc_mtorf_alignment.fas", format = "fasta")
dataAli<-clustal(data)
checkAlignment(data)
data2 <- data[-7,]
checkAlignment(data2)
dataHaplo <- haplotype(data2)
dataHaplo
#
dataHaplo<-sort(dataHaplo, what = "labels")
dataNet<-haploNet(dataHaplo)

countHap <- function(hap = h, dna = x){
  with(
    stack(setNames(attr(hap, "index"), rownames(hap))),
    table(hap = ind, pop = attr(dna, "dimnames")[[1]][values])
  )
}
write.table(countHap (dataHaplo, data2),file="FILE_NAME.txt", sep="\t", quote=FALSE)
#
pdf(file="~/Desktop/mtORF_haplotypes.pdf", width = 8, height = 15, pointsize = 10)
plot(dataNet, size=attr(dataNet, "freq"), scale.ratio=0.2, pie=countHap(dataHaplo, data2), show.mutation=3)
legend("bottomleft", colnames(countHap(dataHaplo, data2)), col=rainbow(ncol(countHap(dataHaplo, data2))), pch=19, ncol=2)
dev.off()
#
pdf(file="~/Desktop/mtORF_haplotypes_sites.pdf", width = 8, height = 15, pointsize = 10)
mydata <- as.data.frame(countHap(dataHaplo, data2))
good <- mydata[mydata$Freq == 1,]
IDs <- strsplit(as.character(good$pop), "_")
IDs <- sapply(IDs, "[[", 3)
new.hap <- table(good$hap, IDs)
plot(dataNet, size=attr(dataNet, "freq"), scale.ratio=0.2, pie=new.hap, show.mutation=3)
legend("bottomright", colnames(new.hap), col=rainbow(ncol(new.hap)), pch=19, ncol=2)
dev.off()