library(GenomicFeatures)
library(rtracklayer)

gtf <- makeTxDbFromGFF("/Users/lukas/CSE/inquilineGenomics/ORevo/data/McKenzie/SupplementaryFile4.AcepOrs.final.clean.gff") #change me!
exons <- exonsBy(gtf, by="gene")

#make introns
exons <- reduce(exons)
exons <- exons[sapply(exons, length) > 1]

introns <- lapply(exons, function(x) {
    #Make a "gene" GRange object
    gr = GRanges(seqnames=seqnames(x)[1], ranges=IRanges(start=min(start(x)),
        end=max(end(x))), 
        strand=strand(x)[1])
    db = disjoin(c(x, gr))
    ints = db[countOverlaps(db, x) == 0]
    #Add an ID
    if(as.character(strand(ints)[1]) == "-") {
        ints$exon_id = c(length(ints):1)
    } else {
        ints$exon_id = c(1:length(ints))
    }
    ints
})
introns <- GRangesList(introns)