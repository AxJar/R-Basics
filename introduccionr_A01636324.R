# Axel Jarquin Morga / A01636324

# 1. Crear una funcion que crea una seciuencia aleatoria de nucleotidos de DNA, size n
secuencia <- function(n) {
  dna <- c('A', 'C', 'T', 'G')
  return(sample(dna, n, replace = T))
}

sec <- secuencia(30); sec

# 2. Crea una funcion que calcule el tamano de la secuencia
secLen <- function(sec) {
  return(length(sec))
}

sec_len <- secLen(sec); sec_len

# 3. Crea una funcion que calcule el porcentaje de cada base de una secuencia 
basePorcentaje <- function(sec) {
  a <- 0
  c <- 0
  t <- 0
  g <- 0
  sec_len <- secLen(sec)
  
  for (i in 1:sec_len) {
    if (sec[i] == "A") {
      a <- a + 1
    }
    if (sec[i] == "C") {
      c <- c + 1
    }
    if (sec[i] == "T") {
      t <- t + 1
    }
    if (sec[i] == "G") {
      g <- g + 1
    }
  }
  
  Nombre <- c('Adenina', 'Timina', 'Citosina', 'Guanina')
  Cantidad <- c(a, c, t, g)
  Porcentaje <- c((a/sec_len)*100, (c/sec_len)*100, (t/sec_len)*100, (g/sec_len)*100)
  
  df <- data.frame(Nombre, Cantidad, Porcentaje); return(df)
}

basePorcentaje(sec)

# 4. Crea una funcion que cuando recibe una hebra directa, regrese una hebra inversa
hebraInv <- function(sec) {
  return(rev(sec))
}

hebra_inversa <- hebraInv(sec); hebra_inversa

# 5. Crea una funcion que cuando recibe una hebra directa, regrese una hebra complementaria
complemento <- function(sec) {
  cdna <- c()
  sec_len <- secLen(sec)
  
  for(i in 1:sec_len) {
    if (sec[i] == "A") {
      cdna[i] <- "T"
    }
    if (sec[i] == "C") {
      cdna[i] <- "G"
    }
    if (sec[i] == "T") {
      cdna[i] <- "A"
    }
    if (sec[i] == "G") {
      cdna[i] <- "C"
    }
  }
  return(cdna)
}

sec_complemento <- complemento(sec); sec_complemento

# 6. Crea una funcion que transcribe DNA a RNA
transcripcion <- function(sec) {
  mrna <- c()
  sec_len <- secLen(sec)
  
  for(i in 1:sec_len) {
    if (sec[i] == "A") {
      mrna[i] <- "U"
    }
    if (sec[i] == "C") {
      mrna[i] <- "G"
    }
    if (sec[i] == "T") {
      mrna[i] <- "A"
    }
    if (sec[i] == "G") {
      mrna[i] <- "C"
    }
  }
  return(mrna)
}

mrna <- transcripcion(sec); mrna

# 7. Crea una funcion que crea codones de inicio y de terminacion de una secuencia
rnaToCodones <- function(mrna) {
  mrna_len <- secLen(mrna)
  codones <- c()
  i <- 1
  j <- 1
  
  while (i < mrna_len) {
    codones[j] <- sprintf("%s%s%s", mrna[i], mrna[i+1], mrna[i+2])
    j <- j + 1
    i <- i + 3
  }
  
  codones_resultado <- c()
  
  for (i in 1:length(codones)) {
    if (codones[i] == 'AUG') {
      codones_resultado[i] <- "START"
    } else if (codones[i] == 'UAA' || codones[i] == 'UAG' || codones[i] == 'UGA') {
      codones_resultado[i] <- "STOP"
      break
    } else if (codones[i] != 'AUG' || codones[i] != 'UAA' || codones[i] != 'UAG' || codones[i] != 'UGA') {
      codones_resultado[i] <- codones[i]
    }
  }
  
  return(codones_resultado)
}

codones <- rnaToCodones(mrna); codones

# 8. Crea una funcion que traduce una secuencia de RNA a una secuencia de aminoacidos
rnaToAmino <- function(mrna) {
  mrna_len <- secLen(mrna)
  amino_acido <- c()
  
  for (i in 1:mrna_len) {
    if (mrna[i] == "UUU" || mrna[i] == "UUC") {
      amino_acido[i] <- "Phe"
    }
    if (mrna[i] == "UUA" || mrna[i] == "UUG" || mrna[i] == "CUU" || mrna[i] == "CUC" || mrna[i] == "CUA" || mrna[i] == "CUG") {
      amino_acido[i] <- "Leu"
    }
    if (mrna[i] == "UCU" || mrna[i] == "UCC" || mrna[i] == "UCA" || mrna[i] == "UCG" || mrna[i] == "AGU" || mrna[i] == "AGC") {
      amino_acido[i] <- "Ser"
    }
    if (mrna[i] == "UAU" || mrna[i] == "UAC") {
      amino_acido[i] <- "Tyr"
    }
    if (mrna[i] == "UGU" || mrna[i] == "UGC") {
      amino_acido[i] <- "Cys"
    } 
    if (mrna[i] == "UGG") {
      amino_acido[i] <- "Trp"
    }
    if (mrna[i] == "CCU" || mrna[i] == "CCC" || mrna[i] == "CCA" || mrna[i] == "CCG") {
      amino_acido[i] <- "Pro"
    }
    if (mrna[i] == "CAU" || mrna[i] == "CAC") {
      amino_acido[i] <- "His"
    }
    if (mrna[i] == "CAA" || mrna[i] == "CAG") {
      amino_acido[i] <- "Gln"
    }
    if (mrna[i] == "CGU" || mrna[i] == "CGC" || mrna[i] == "CGA" || mrna[i] == "CGG" || mrna[i] == "AGA" || mrna[i] == "AGG") {
      amino_acido[i] <- "Arg"
    }
    if (mrna[i] == "GUU" || mrna[i] == "GUC" || mrna[i] == "GUA" || mrna[i] == "GUG") {
      amino_acido[i] <- "Val"
    }
    if (mrna[i] == "GCU" || mrna[i] == "GCC" || mrna[i] == "GCA" || mrna[i] == "GCG") {
      amino_acido[i] <- "Ala"
    }
    if (mrna[i] == "GAU" || mrna[i] == "GAC") {
      amino_acido[i] <- "Asp"
    }
    if (mrna[i] == "GAA" || mrna[i] == "GAG") {
      amino_acido[i] <- "Glu"
    }
    if (mrna[i] == "GGU" || mrna[i] == "GGC" || mrna[i] == "GGA" || mrna[i] == "GGG") {
      amino_acido[i] <- "Gly"
    }
    if (mrna[i] == "UAA" || mrna[i] == "UAG" || mrna[i] == "UGA") {
      amino_acido[i] <- "STOP"
      break
    }
    if (mrna[i] == "AUU" || mrna[i] == "AUC" || mrna[i] == "AUA") {
      amino_acido[i] <- "Ile"
    }
    if (mrna[i] == "ACU" || mrna[i] == "ACC" || mrna[i] == "ACA" || mrna[i] == "ACG") {
      amino_acido[i] <- "Thr"
    }
    if (mrna[i] == "AAU" || mrna[i] == "AAC") {
      amino_acido[i] <- "Asn"
    }
    if (mrna[i] == "AAA" || mrna[i] == "AAG") {
      amino_acido[i] <- "Lys"
    }
    if (mrna[i] == "AUG") {
      amino_acido[i] <- "Met"
    }
    
  }
  
  return(amino_acido)
}

rnaToAmino(codones)

