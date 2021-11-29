# coding=utf-8

from project_alignment_optimizer.program.variables import Var
import logging as log
import pathlib
from Bio import SeqIO
from Bio import AlignIO
from Bio import Align
from Bio.Seq import Seq
from Bio.Align.Applications import ClustalwCommandline
import re
import os
import sys

from Bio import Entrez,pairwise2
from Bio.SeqRecord import SeqRecord

# -------------------
# Funciones Principales
# ---------------------


def configureVariables():
    # Pido por consola todas las variables a utilizar por el sistema
    printAndLog("Solicito Variables")


def loadFile(filename):
    # Validar path
    # Validar formato
    # Validar que sea la cantidad de secuencias sea mayor al valor minimo
    # Obtengo el archivo con el alineamiento inicial
    checkIsValidPath(filename)
    with open(filename, "r") as handle:
        originalAlignment = AlignIO.read(handle, "fasta")
        if (any(originalAlignment)):
            printAndLog("Getting alignments")
            printAndLog(str(len(originalAlignment)) +
                        " aligned sequences were obtained correctly")
            return originalAlignment
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


def filterAlignment(anAlignment):
    # Saco la secuencia que tiene mas gaps
    # Se puede agregar la comparacion con la query usando pairwise si hay dos con la misma cantidad de gaps
    # (Mas adelante, en esta funcion tambien usamos lo de las secuencias homologas)
    mostGappedSeq = anAlignment[0]
    for sequence in anAlignment:
        mostGappedSeq = mostGapped(mostGappedSeq, sequence,anAlignment[0])
    printAndLog("Sequence {id} has been removed from the alignment.".format(
        id=mostGappedSeq.id))

    sequences = list(filter(lambda seq: seq.id != mostGappedSeq.id, anAlignment))
    printAndLog(str(len(sequences)) + " sequences left.")
    return sequences


def mostGapped(aSequence, anotherSequence, aQuerySequence):
    # Obtengo la secuencia con mayor cantidad de gaps
    if gaps(aSequence) > gaps(anotherSequence):
        return aSequence
    elif gaps(aSequence) == gaps(anotherSequence):
        aligner = Align.PairwiseAligner()
        alignment1 = aligner.align(aSequence.seq.ungap(), aQuerySequence.seq.ungap())
        alignment2 = aligner.align(anotherSequence.seq.ungap(),aQuerySequence.seq.ungap())
        if alignment1.score > alignment2.score:
            return aSequence
        else:
            return anotherSequence

    else:
        return anotherSequence


def gaps(aSequence):
    return len(aSequence.seq) - len(aSequence.seq.ungap())


def getungappedSequences(anAlignment):
    # Obtengo del alineamiento pasado por parametro las secuencias originales
    # Para eso elimino todos los gaps
    baseSequences = []
    for sequence in anAlignment:
        sequence.seq = sequence.seq.ungap()
        baseSequences.append(sequence)
    return baseSequences


def parseScore(aClustalOutputString):
    return re.search("(?<=Alignment Score )\d*", aClustalOutputString).group()


def generateAlignmentAndCalculateScore(originalSequences):
    # Genero el nuevo alineamiento por medio de CLUSTAL
    tempDir = str(pathlib.Path(__file__).parent.resolve())
    SeqIO.write(originalSequences, (tempDir + "/seqs.fasta"), "fasta")
    command = ClustalwCommandline(
        Var().clustalWPath(), infile=(tempDir + "/seqs.fasta"))
    clusalAlignmentOutput = command()
    score = parseScore(clusalAlignmentOutput[0])
    printAndLog("New alignment Score: " + score)
    return score


def loadCurrentAlignment():
    tempDir = str(pathlib.Path(__file__).parent.resolve())
    return AlignIO.read(tempDir + "/seqs.aln", "clustal")


def generateTree(alignment):
    # Genero el árbol filogenético
    tree = alignment
    return tree

def getHomologuesSequencesOrderedByMaxScore(seqQuery):
    #Se pasa como parametro un SeqRecord
    #y devuelve una lista de SeqRecord 
    #ordenada de mayor a menor por score con la seqQuery
    #Habria que ver si solo traer, por ejemplo, un N definido de seq homologas
    result = []
    idsProteins = getIdsHomologuesSequences(seqQuery.id)
    if idsProteins:
        idString = ""
        for idProtein in idsProteins:
            idString = idString + idProtein + ","
        Entrez.email = "A.N.Other@example.com"
        idString = idString[0:len(idString)-1]
        handle = Entrez.efetch(db="protein", id=idString, rettype="gb", retmode="xml")
        output = Entrez.read(handle)
        result = getSequencesOrderedByMaxScore(seqQuery, output)
    return result

# ---------------------
# Funciones Auxiliares
# ---------------------

def calculateScorePar(sequences):
    score = 0
    for u in range(0, len(sequences)):
        for v in range(u+1, len(sequences)):
            score += scorePar(sequences[u], sequences[v])
    return score


def scorePar(u, v):
    score = 0
    for i in range(0, len(u)):
        score += compareCost(u[i], v[i])
    return score


def compareCost(u, v):
    match = 1  # var.match()
    mismatch = -1  # var.mismatch()
    gap = -1  # var.gapPenalty()
    if (u == v):
        return match
    else:
        if(u == "-" or v == "-"):
            return gap
        return mismatch

def printAndLog(msg):
    print(msg)
    log.info(msg)

def getIdsHomologuesSequences(idProtein):    
    Entrez.email = "A.N.Other@example.com"
    handle = Entrez.efetch(db="protein", id=idProtein,rettype='ipg', retmode='xml')
    output = Entrez.read(handle)
    proteinList = output['IPGReport']['ProteinList']
    result = []
    if proteinList:
        for indx in range(0,len(proteinList)):
            result.append(proteinList[indx].attributes['accver'])

    return result

def getSequencesOrderedByMaxScore(seqQuery,output):
    result = []
    if output:
        for protein in output:
            if protein:
                descriptor = protein["GBSeq_locus"]
                seqHomologue = SeqRecord(Seq(protein["GBSeq_sequence"].upper()),id =descriptor,name=descriptor,description=descriptor)
                result.append((seqHomologue,pairwise2.align.globalxx(seqQuery.seq, seqHomologue.seq,score_only=True)))
        result.sort(key=takeSecond,reverse=True)
    return [i[0] for i in result] 

def takeSecond(elem):
    return elem[1]

