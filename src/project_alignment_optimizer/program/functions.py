# coding=utf-8

import variables as Var
import logging as log
from Bio import SeqIO
import os
import sys

# -------------------
# Funciones Principales
# ---------------------

def configureVariables():
  # Pido por consola todas las variables a utilizar por el sistema
  printAndLog("Solicito Variables")

def loadFile(filename):
  # Validar path
  # Validar formato
  # Obtengo el archivo con el alineamiento inicial
  checkIsValidPath(filename)
  with open(filename, "r") as handle:
    fasta = SeqIO.parse(handle, "fasta")
    if (any(fasta)):
      printAndLog("Getting alignments")
      sequences=[]
      for f in fasta:
        seq = str(f.seq)
        sequences.append(seq)
      printAndLog(str(len(sequences))+" alignments were obtained correctly")
      return sequences
    else:
        print("Invalid file format")

def checkIsValidPath(filename):
   if not os.path.isfile(filename):
      print("Invalid file path")
      sys.exit()

def calculateScore(alignment):
  # Sacar el score del alineamiento obtenido por CLUSTAL
  finalScore = len(alignment)
  return finalScore

def filterAlignment(alignment):
  # Saco la secuencia que tiene mas gaps
  # Si hay mas de una, elimino la que tenga menor score con respecto a la secuencia query, usando alineamiento de a pares
  # (Mas adelante, en esta funcion tambien usamos lo de las secuencias homologas)
  cantGaps = 0
  for seq in range(0,len(alignment)):
    if (seq != 0):#var.querySequence()):
      gaps = alignment[seq].count("-")
      if (gaps > cantGaps):
        minGapSeq = seq
        cantGaps = gaps
      elif (gaps == cantGaps):
        scoreMinGapSeq = scorePar(alignment[minGapSeq], alignment[0])#alignment[var.querySequence()])
        scoreNewSeq = scorePar(alignment[seq], alignment[0])#alignment[var.querySequence()])
        if (scoreNewSeq < scoreMinGapSeq):
          minGapSeq = seq
          cantGaps = gaps # El valor es gual al anterior que estaba guardado
  response = alignment
  printAndLog("Se procede a eliminar del alineamiento la secuencia: " + alignment[minGapSeq])
  response.pop(minGapSeq)
  return response

def getOriginalSequences(alignment):
  # Obtengo del alineamiento pasado por parametro las secuencias originales
  # Para eso elimino todos los gaps
  originalSequences = alignment
  return originalSequences

def generateAlignment(originalSequences):
  # Genero el nuevo alineamiento por medio de CLUSTAL
  newAlignment = originalSequences
  return newAlignment

def generateTree(alignment):
  # Genero el árbol filogenético
  tree = alignment
  return tree


# ---------------------
# Funciones Auxiliares
# ---------------------

def calculateScorePar(sequences):
  score = 0
  for u in range(0,len(sequences)):
    for v in range(u+1,len(sequences)):
      score += scorePar(sequences[u],sequences[v])
  return score

def scorePar(u,v):
  score = 0
  for i in range(0,len(u)):
    score += compareCost(u[i],v[i])
  return score

def compareCost(u,v):
    match = 1#var.match()
    mismatch = -1#var.mismatch()
    gap = -1#var.gapPenalty()
    if (u == v):
        return match
    else:
        if(u == "-" or v == "-"):
            return gap
        return mismatch

def printAndLog(msg):
    print(msg)
    log.info(msg)