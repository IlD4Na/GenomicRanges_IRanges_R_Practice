install.packages("BiocManager")
library(BiocManager)


BiocManager::install("IRanges")
BiocManager::install("GenomicRanges", dependencies=TRUE)


library(IRanges) # Gestione di Intervalli genomici 
library(GenomicRanges)


BiocManager::install("Biostrings")

library(Biostrings)

"""
Il pacchetto Biostrings di R contiene funzioni e classi per la 
rappresentazione di stringhe biologiche, come DNA, RNA e 
amino acidi. Contiene anche funzioni per l'allineamento:
sia locale (SW) che globale (NW).

Ha inoltre funzioni per leggere e scrivere FASTA file
"""


"""
Possiamo assegnare al DNA i vari nucleotidi o altri caratteri dello 
IUPAC_CODE_MAP
"""


dna1 <- DNAString("ACGT") #Contiene dentro una sequenza di DNA

dna1

IUPAC_CODE_MAP 
# Praticamente qui vediamo V che sta a significare che la base può essere o A o C o G,
# oppure N che va a indicare che potrebbe essere qualunque. 

dna1 <- DNAString("ACGTMRWSYKVHDBN")

dna1

# Il set di stringhe ricorda un pò la struttura dei IRanges, qui vediamo 3 seq
# di lunghezza differente e vediamo che a sinistra ci fa vedere il loro spessore

dna2 <- DNAStringSet(c("ACGT","CGTN","NNNGGCGCGTTAA"))
dna2

names(dna2) <- c("seq1","seq2","seq3")


#C'è tutta una serie di funzioni anche per l'RNA che servono a creare seq di RNA

rna1 <- RNAString("UUUACCCGGGN")
rna1
rna2 <- RNAStringSet(c("UUUUAAAACCC","NUGCAUGGGA"))
rna2

# Per gli AA usiamo AAString e AAStringset 

aa1 <- AAString("ARNXZZZZ")
aa1
AA_STANDARD

# Confronto tra le seq in dna1 e dna2
dna1==dna2
unique(dna2)
sort(dna2)

rev(dna2) # applicata ad un set restituisce le sequenze in ordine invertito

rev(dna1) # Invece se applicato ad una sola sequenza mi da il contrario, me la ribalta


## 2 Parte ## 

dna_1<- DNAString("ACTGGC")
dna_1
translate(dna_1) #Traduce in aa
rna_1<-RNAString("AUG")
rna_1
translate(rna_1) #Traduce in aa

#Ottenere la sequenza complementare
complement(dna_1)

#Ottenere la sequenza reverse complementare
reverseComplement(dna_1)

letterFrequency(dna1, "GCN") # Quante volte è presente G e C e N

letterFrequency(dna1, "N") # Quante volte è presente N

dinucleotideFrequency(dna1) #Mi dice tutte le combinazioni di dinucleotidi

trinucleotideFrequency(dna1) # Controlla tutte le frequenze di tutti i trinucleotidi


oligonucleotideFrequency(dna1, width = 4)

"""
consensusMatrix prende un insieme di sequenze tutte della stessa 
lunghezza (allineamento) e costruisce una matrice che indica quante volte
compare ogni nucleotide o aa in ogni posizione. 
"""
consensusMatrix(dna2, as.prob=T) 


dna1_lunga <- DNAString("ACTACAGTTCAAATGGTAAACCGGTTTACCCAAA")

"""
Scorre la sequenza usando una finestra (sliding window) di lunghezza 3 
in ogni finestra conta quante A sono presenti.
"""
letterFrequencyInSlidingView(dna1_lunga,view.width = 3,letters = "A")


letterFrequencyInSlidingView(dna1_lunga,view.width = 3,letters = "AC") # conta in ogni sliding window sia quante che A che quante C



## GenomicRanges  GRanges 


"""
Estende l'IRanges per includere annotazioni come nome del cromosoma, 
nome del gene o il tipo di tessuto.
"""


gr <- GRanges(seqnames = "chr1", strand = c("+","-","*"),
              ranges = IRanges(start = c(1,5,7),width = 3))

strand(gr)
seqnames(gr)
ranges(gr)
start(gr)
end(gr)
width(gr)
seqinfo(gr) #ci dice varie informazioni utili


seqlengths(gr) <- 3
genome(gr) <- "hg19"
isCircular(gr) <- T

seqinfo(gr) # Adesso ricontrollando vedo più cose.

"""
Se io faccio class(seqnames(gr)) mi dirà che appartiene alla
classe Rle che è una speciale classe di Bioconductor, 
rappresenta un vettore compresso ovvero un 
Run-Lenght-Encoding ovvero invece di memorizzare
gli elementi ripetuti, memorizza il valore e il numero di 
ripetizioni. 
"""


