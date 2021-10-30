library("seqinr")
setwd("C:/Users/axela/OneDrive - Instituto Tecnologico y de Estudios Superiores de Monterrey/Materias 2do Semestre/Biologia computacional")

sarsCovid <- read.fasta("SarsCovid_virus.fasta")
class(sarsCovid)

wuhanVirus <- read.fasta("Wuhan_virus.fasta")
class(wuhanVirus)

zikaVirus <- read.fasta("Zika_virus.fasta")
class(zikaVirus)

middleVirus <- read.fasta("MiddleEast_virus.fasta")
class(middleVirus)

dengueVirus <- read.fasta("Dengue_virus.fasta")
class(dengueVirus)

#Obtener la longitud de genomas virales
length(sarsCovid[[1]])

length(wuhanVirus[[1]])

length(zikaVirus[[1]])

length(middleVirus[[1]])

length(dengueVirus[[1]])

#Obtener el contenido de nucleotidos
count(seq = sarsCovid[[1]],wordsize = 1)

count(seq = wuhanVirus[[1]],wordsize = 1)

count(seq = zikaVirus[[1]],wordsize = 1)

count(seq = middleVirus[[1]],wordsize = 1)

count(seq = dengueVirus[[1]],wordsize = 1)

#Obtener la secuencia complementaria
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

# 2. Crea una funcion que calcule el tamano de la secuencia
secLen <- function(sec) {
  return(length(sec))
}

sec_len <- secLen(sec); sec_len


comp(sarsCovid[[1]])

comp(wuhanVirus[[1]])

comp(zikaVirus[[1]])

comp(middleVirus[[1]])

comp(dengueVirus[[1]])

#obtencion de secuencias de acidos nucleicos o proteinas en formato fasta
vector_sarsCovid <- sarsCovid[[1]]
class(vector_sarsCovid)
typeof(vector_sarsCovid)

vector_wuhanVirus <- wuhanVirus[[1]]
class(vector_wuhanVirus)

vector_zikaVirus <- zikaVirus[[1]]
class(vector_zikaVirus)

vector_middleVirus <- middleVirus[[1]]
class(vector_middleVirus)

vector_dengueVirus <- dengueVirus[[1]]
class(vector_dengueVirus)

#Calculo del porcentaje de nucleotidos
length(vector_sarsCovid)
length(vector_dengueVirus)
length(vector_middleVirus)
length(vector_zikaVirus)
length(vector_wuhanVirus)

basePorcentaje(vector_sarsCovid)
basePorcentaje(vector_dengueVirus)
basePorcentaje(vector_middleVirus)
basePorcentaje(vector_zikaVirus)
basePorcentaje(vector_wuhanVirus)

base.content <- function(dna) {
  a <- 0 
  t <- 0
  g <- 0
  c <- 0
  for (i in 1:length(dna)) {
    if(dna[i] == "a") {
      a <- a + 1
    } else if(dna[i] == "t"){
      t <- t + 1
    } else if(dna[i] == "c"){
      c <- c + 1
    } else if(dna[i] == "g"){
      g <- g + 1
    }
  }
  print("Contenido de Adenina: ")
  print(a)
  print("Contenido de Timina: ")
  print(t)
  print("Contenido de Citosina: ")
  print(c)
  print("Contenido de Guanina: ")
  print(g)
}

base.content(vector_sarsCovid)
base.content(vector_dengueVirus)
base.content(vector_middleVirus)
base.content(vector_wuhanVirus)
base.content(vector_zikaVirus)



length(sarsCovid[[1]])

