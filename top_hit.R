library(tidyverse)


#read in blast results
blast <- read_table2("~/blast/blastoutput.txt", col_names = F)  
colnames(blast) <- c("qseqid","sseqid","pident","qlen","length","qstart","qend","sstart","send","evalue")
blast$cover <- blast$length/blast$qlen

#keep hit with highest cover, pident, then evalue. what to do when cover=pident=evalue?
blast <- group_by(blast, qseqid) %>% top_n(1, cover) %>% top_n(1, pident) %>% top_n(-1, evalue)

write.table(blast, "~/blast/blastfiltered.txt", sep="\t",row.names=FALSE, col.names =FALSE, quote=FALSE)
