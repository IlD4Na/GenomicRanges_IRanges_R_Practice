install.packages("BiocManager")
library(BiocManager)


BiocManager::install("IRanges")
BiocManager::install("GenomicRanges", dependencies=TRUE)


library(IRanges) # Gestione di Intervalli genomici 
library(GenomicRanges)


## USO OGGETTO IRanges ## 


### Iranges è basato su intervalli, ogni intervallo è chiamato range


# Un oggetto con 1 range e 0 colonne di metadata
sin_interval <- IRanges::IRanges(start=1,end=5)

# Definizione di un range dove impostiamo lo start e lo spessore (width)
# width = Lunghezza dell'intervallo comprendendo l'inizio e la fine 
# width = end - start + 1

sin_interval_0 <- IRanges::IRanges(star=1,width = 6)
# calcola lui la posizione di end 

# Multiple ranges

ir <- IRanges::IRanges(start = c(1,2,5),end = c(3,5,7))

ir2 <- IRanges::IRanges(start = c(1,2,5), width = 3)

ir3 <- IRanges::IRanges(start = c(1,2,5),width = c(3,4,3))

all.equal(ir,ir2) # Non sono gli stessi oggetti perchè qui ho impostato l'intervallo uguale in tutti e 3 i casi 

all.equal(ir,ir3) # Adesso sono uguali perchè ho specificato l'intervallo di ogni range 


# Funzioni che servono per estrapolare informazioni dai nostri ranges
# Faccio sempre riferimento al pacchetto con IRanges::


IRanges::start(ir)
IRanges::end(ir)
IRanges::width(ir)

# Assegnazione di valori tramite start,end e width

IRanges::start(ir) <- c(1,2,4) # Assegno gli inizi come 1,2,4

IRanges::start(ir)[2] <- 4 # Cambio il 2 a 4

IRanges::width(ir)[c(1,3)] <- c(1,1) # Assegno la lunghezza uguale a 1 al primo e terzo range 


# Assegnazione di nomi a ogni range 

# Ho dato un nome ad ogni range 
names(ir) <- paste("Range -",1:length(ir)) 

length(ir) # È il numero di ranges 

# Subsettings

ir["Range - 1"] #Seleziono l'elemento che si chiama "Range - 1"

ir[1] #Seleziono il 1 elemento

ir[c("Range - 1","Range - 3")] #Seleziono gli elementi che hanno quel nome


# Unione di vari subranges 

ir5 <- c(ir,ir2,ir3)

names(ir5) <- paste("Range:",1:length(ir5))

## NORMAL IRANGES ##

'''  
Gli intervalli normalizzati sono unici e non si sovrappongono, se io ho 3 
intervalli di questo tipo [1–5], [3–7], [6–10] posso usare la funzione 
reduce(ir) e vado ad ottenere gli intervalli normalizzati che rappresentano la
fusione dei gruppi sovrapposti, in questo caso sarebbe [1-10].
''' 
## Uso delle funzioni reduce(), disjoin(), gaps()


# Funzione che stampa un grafico di Iranges con i vari intervalli
plot.IRanges <- 
  function(x,
           xlim = NULL,
           ylim = NULL,
           height=0.8,
           title = NULL,
           col= "grey",
           border = NULL,
           add = FALSE,
           xlab=NULL,
           ylab=NULL,
           add.labels=FALSE,
           col.labels = "black",
           ...) {
    if (is.null(xlab))
      xlab <- ""
    if (is.null(ylab))
      ylab <-""
    xrng <- range(x)
    group <- GenomicRanges::disjointBins(x)
    if (is.null(xlim))
      xlim <- c(IRanges::start(xrng),IRanges::end(xrng))
    if (is.null(ylim))
      ylim <- c(1-height / 2, max(group) + height/2)
    
    if (!add) {
      graphics::plot.default(
        NULL,
        xlim=xlim,
        ylim=ylim,
        main =title,
        axes = FALSE,
        xlab = xlab,
        ylab= ylab
      )
      abline(h = unique(group), col = "grey90")
      axis(1, at=pretty(xlim), tick = TRUE)
    }
    
    rect(
      xleft = start(x),
      ybottom = group - height / 2,
      xright = end(x),
      ytop=group + height / 2,
      col=col,
      border=border
    )
    if (add.labels & !is.null(names(x))) {
      text(
        x = mid(x),
        y=group,
        labels = names(x),
        col = col.labels
      )
    }
  }


ir <- IRanges::IRanges(start=c(2,4,8,9),end = c(5,5,9,11))


plot.IRanges(ir, title = "Il mio primo Iranges")


# Uso della funzione di reduce 

red_ir <- IRanges::reduce(ir) # I 4 ranges sono stati sintetizzati in due range



par(mfrow=c(2,1)) # Con questo parametro stiamo separando la finestra del plot in 2, 2 righe e 1 colonna

plot.IRanges(ir, title = "Classic Iranges")
plot.IRanges(red_ir, title = "Normalized Iranges")

#Il Disjoin che suddivide in intervalli non sovrapposti 

disj_ir <- disjoin(ir)

par(mfrow=c(3,1))

plot.IRanges(ir, title = "Classic Iranges")
plot.IRanges(disj_ir, title = "Disjoint Iranges")
plot.IRanges(red_ir, title = "Normalized Iranges")


#Uso di gaps 

"""

gaps() restituisce gli intervalli che non sono presenti nell'IRanges originale 
ovvero dove nessun intervallo originale arriva 

"""

par(mfrow=c(2,1))

ir_withgap <- IRanges(start = c(1,6,13),end = c(4,9,16))

plot.IRanges(ir_withgap, title = "IR with gaps", xlim = c(0,17))

fill_ir <- gaps(ir_withgap)

plot.IRanges(fill_ir, title = "Gaps(ir)", xlim = c(0,17))


"""

Quelle fatte fino ad adesso sono trasformazioni Inter ranges perchè trasformano 
tutti gli intervalli insieme per produrre un nuovo set di ranges. 

"""

## MANIPOLAZIONE INTRA RANGE ##

"""
Prendono singolarmente ogni intervallo dentro range e lo trasformano
individualmente senza considerare gli altri. 

Sono trasformazioni dove ogni range viene mappato in un nuovo range ed alcuni
esempi sono: Shift, Narrow, Flank, Resize, Restrict
"""


ir <- IRanges::IRanges(start= c(1,5,7), end= c(5,9,12))

# Shift ovvero shifta di tot posizioni l'inizio e la fine e lascia invariata la lunghezza dell'intervallo

shift(ir, shift = 3) # shift in avanti

shift(ir, shift = -2) # shift indietro

# Narrow restringe i nostri intervalli, diamo un inizio e fissa la fine

narrow(ir,start = 4)

narrow(ir, end=4) #significa che la 5 posizione dopo lo start sarà la loro ultima
                  #posizione

# Resize che ridimensiona i nostri ranges in modo tale che abbiano uno width specifico

resize(ir, width = 4,fix = "start") #viene mantenuto fisso lo start

resize(ir, width = 4,fix = "center") #viene mantenuto fisso il centro dell'intervallo

resize(ir, width = 4,fix = "end") #viene mantenuto fisso la fine

# Definiamo un valore preciso di spessore in questo caso 4


#flank permette di generare delle regioni fiancheggianti 

flank(ir,width = 3, start = T) # genera delle porzioni fiancheggianti rispetto
                               # alle posizioni di start

flank(ir,width = 3, start = F) #Qui genera le posizioni fiancheggianti 


#restrict serve a 


restrict(ir,start = 5, end = 10)

"""
Con restrict noi andiamo a impostare un valore minimo di start di 5, quindi nell'esempio
vediamo che nel primo caso da 1 è diventato 5 nel secondo è rimasto invariato nel terzo
pure lo start è rimasto invariato però l'end è cambiata perche abbiamo un valore massimo 
di 10 impostato come end. 
"""




