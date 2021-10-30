library("seqinr")
setwd("C:/Users/Jonathan Quirino/Desktop/Itesm_SEM2/Periodo 2/Biologica Comput/evidencia_1")

AY508724 <- read.fasta("AY508724.fasta")
class(AY508724)

AY485277 <- read.fasta("AY485277.fasta")
class(AY485277)

AY390556 <- read.fasta("AY390556.fasta")
class(AY390556)

AY278489 <- read.fasta("AY278489.fasta")
class(AY278489)

MN908947 <- read.fasta("MN908947.fasta")
class(MN908947)

MN985325 <- read.fasta("MN985325.fasta")
class(MN985325)

MT292571 <- read.fasta("MT292571.fasta")
class(MT292571)

MW000351 <- read.fasta("MW000351.fasta")
class(MW000351)

MT873893 <- read.fasta("MT873893.fasta")
class(MT873893)

JX869059 <- read.fasta("JX869059.fasta")
class(JX869059)

#Obtener la longitud de genomas virales

length(AY508724[[1]])

length(AY485277[[1]])

length(AY390556[[1]])

length(AY278489[[1]])

length(MN908947[[1]])

length(MN985325[[1]])

length(MT292571[[1]])

length(MW000351[[1]])

length(MT873893[[1]])

length(JX869059[[1]])

#Ejercicio 4

secLen <- function(sec) {
  return(length(sec))
}

basePorcentaje <- function(sec) {
  a <- 0
  c <- 0
  t <- 0
  g <- 0
  sec_len <- secLen(sec)
  
  for (i in 1:sec_len) {
    if (sec[i] == "a") {
      a <- a + 1
    }
    if (sec[i] == "c") {
      c <- c + 1
    }
    if (sec[i] == "t") {
      t <- t + 1
    }
    if (sec[i] == "g") {
      g <- g + 1
    }
  }
  
  Nombre <- c('Adenina', 'Timina', 'Citosina', 'Guanina')
  Cantidad <- c(a, c, t, g)
  Porcentaje <- c((a/sec_len)*100, (c/sec_len)*100, (t/sec_len)*100, (g/sec_len)*100)
  
  df <- data.frame(Nombre, Cantidad, Porcentaje); return(df)
}

vector_AY508724 <- AY508724[[1]]
class(vector_AY508724)

vector_AY485277 <- AY485277[[1]]
class(vector_AY485277)

vector_AY390556 <- AY390556[[1]]
class(vector_AY390556)

vector_AY278489 <- AY278489[[1]]
class(vector_AY278489)

vector_MN908947 <- MN908947[[1]]
class(vector_MN908947)


vector_MN985325 <- MN985325[[1]]
class(vector_MN985325)

vector_MT292571 <- MT292571[[1]]
class(vector_MT292571)

vector_MW000351 <- MW000351[[1]]
class(vector_MW000351)

vector_MT873893 <- MT873893[[1]]
class(vector_MT873893)

vector_JX869059 <- JX869059[[1]]
class(vector_JX869059)

basePorcentaje(vector_AY508724)
basePorcentaje(vector_AY485277)
basePorcentaje(vector_AY390556)
basePorcentaje(vector_AY278489)
basePorcentaje(vector_MN908947)

basePorcentaje(vector_MN985325)
basePorcentaje(vector_MT292571)
basePorcentaje(vector_MW000351)
basePorcentaje(vector_MT873893)
basePorcentaje(vector_JX869059)

invert <- function(dna){
  return(rev(dna))
}



