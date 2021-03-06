# coding=utf-8

from urllib.request import Request
import Bio
from project_alignment_optimizer.program.constants import CLUSTALW_PATH, TEMP_DIR, GAPEXT, GAPOPEN, MATRIX, N_HOMOLOGOUS_SEQUENCES, ADMIT_HOMOLOGOUS, MATCH, MIN_SEQUENCES, MISMATCH, NOT_VALID_QUERY_NO, VALID_QUERY_YES, DB_HOMOLOGOUS_SEQUENCES, PURIFY_START, PURIFY_END
from Bio import SeqIO, AlignIO, Align, Entrez, pairwise2, Phylo
from Bio.Seq import Seq
from Bio.Align.Applications import ClustalwCommandline
from Bio.SeqRecord import SeqRecord
import logging as log
import pathlib
import re
import os
import sys
import copy as c
import datetime


# ---------------------
# Funciones Principales
# ---------------------

def find_alignment_by_header(currentAlignment, query_sequence_header):
    for alignment in currentAlignment:
        if alignment.id == query_sequence_header:
            return alignment

    printAndLogCritical(f"The header {query_sequence_header} has not been found.")
    sys.exit()


def executeFirstAlgorithm(lastAlignment, querySeq, homologousSequences, lastScore, env_variables):
    # Saco la secuencia que tenga mas aminoacidos en las columnas donde la query tenga gaps
    # El tipo de filtrado es el 0
    alignmentFiltered = filterSequence(
        0, lastAlignment, querySeq, homologousSequences, env_variables)
    # Ejecuto el algoritmo
    return executeAlgorithm(alignmentFiltered, lastAlignment, homologousSequences, lastScore, env_variables)


def excuteSecondAlgorithm(lastAlignment, querySeq, homologousSequences, lastScore, env_variables):
    # Saco la secuencia que tenga mas gaps de todo el alineamiento
    # El tipo de filtrado es el 1
    alignmentFiltered = filterSequence(
        1, lastAlignment, querySeq, homologousSequences, env_variables)
    # Ejecuto el algoritmo
    return executeAlgorithm(alignmentFiltered, lastAlignment, homologousSequences, lastScore, env_variables)


def executeAlgorithm(alignmentFiltered, lastAlignment, homologousSequences, lastScore, env_variables):
    min_sequences = env_variables[MIN_SEQUENCES]
    if len(alignmentFiltered) < min_sequences:
        alignmentWithHomologousSequence, homologousSequence = addHomologousSequence(
            alignmentFiltered, homologousSequences)
        newAlignment, newScore = generateNewAlignmentAndScore(
            alignmentWithHomologousSequence, env_variables)
        if(newScore > lastScore):
            return True, newAlignment, newScore
        else:
            homologousSequences.insert(0, homologousSequence)
            return False, lastAlignment, lastScore
    else:
        # Calculo el nuevo score y alineamiento
        copyAl = c.deepcopy(lastAlignment)
        newAlignment, newScore = generateNewAlignmentAndScore(
            alignmentFiltered, env_variables)
        # Compruebo si mejoro el alineamiento

        if(newScore > lastScore):
            return True, newAlignment, newScore
        else:
            printAndLogInfo("The alignment score didn't improve ????")
            # Como no mejoro, si acepta secuiencias homologas y todavia quedan por agregar, agrego una de ellas
            admit_homologous = env_variables[ADMIT_HOMOLOGOUS]
            if (admit_homologous == 1) and (len(homologousSequences) > 0):
                alignmentWithHomologousSequence, homologousSequence = addHomologousSequence(
                    newAlignment, homologousSequences)
                newAlignment, newScore = generateNewAlignmentAndScore(
                    alignmentWithHomologousSequence, env_variables)
                if(newScore > lastScore):
                    return True, newAlignment, newScore
                else:
                    homologousSequences.insert(0, homologousSequence)
                    return False, copyAl, lastScore
            else:
                return False, copyAl, lastScore

def align(args, env_variables):

    printAndLogInfo("---------------------------------------")
    fileName = args.file
    query_sequence_header = args.query_sequence_header

    # Cargo el archivo con el alineamiento inicial que me pasa el usuario
    try:
        currentAlignment = loadFile(fileName)
    except BaseException as err:
        if not isinstance(err, KeyboardInterrupt):
            printAndLogCritical("ERROR: "+str(err))
        else:
            printAndLogCritical("Keyboard interrupt")
        sys.exit()
    printAndLogInfo("---------------------------------------")

    outputDir = createOutputDir(fileName, query_sequence_header)

    if(alignmentHasNMinSequences(currentAlignment, env_variables)):

        # Recorto, cuando corresponda, las secuencias purificadoras
        currentAlignment = trimPurifyingSequences(currentAlignment, env_variables)

        # Obtengo las Secuencia Query
        querySeq = find_alignment_by_header(currentAlignment, query_sequence_header)
        printAndLogInfo("Query Sequence: " + querySeq.id)
        printAndLogInfo("---------------------------------------")

        # Obtengo las secuencias homologas
        homologousSequences = getHomologousSequences(
            querySeq, currentAlignment, env_variables, args.homologous_sequences_path)

        # Obtengo las secuencias originales
        seqsAux = c.deepcopy(currentAlignment)
        ungappedSequences = getUngappedSequences(seqsAux)

        # Genero el alineamiento para obtener el score perteneciente al alineamiento inicial pasado por el usuario
        # Este alineamiento lo descarto, ya que no me sirve
        currentScore = generateAlignmentAndCalculateScore(
            ungappedSequences, env_variables)

        # Genero nuevos alineamientos y sus scores correspondientes
        # Mientras aumente el score actual sigo aplicando los filtros
        bestAlignment = currentAlignment
        bestScore = currentScore
        better = True
        printAndLogInfo("Current Alignment: " + str(len(currentAlignment)))
        while(better):
            lastAlignment = currentAlignment
            copyLastAlignment = c.deepcopy(lastAlignment)
            lastScore = currentScore
            printAndLogInfo("---------------------------------------")
            # Hago el primer filtrado (Saco la secuencia que tenga mas aminoacidos en las columnas donde la query tenga gaps)
            printAndLogInfo(
                "ALGORITHM 1 - Filter sequence that introduces most gaps to the Query Sequence:")

            improve, currentAlignment, currentScore = executeFirstAlgorithm(
                lastAlignment, querySeq, homologousSequences, lastScore, env_variables)
            querySeq = find_alignment_by_header(currentAlignment, query_sequence_header)
            if(improve):
                bestAlignment = currentAlignment
                bestScore = currentScore
                lastScore = currentScore
                lastAlignment = copyLastAlignment
                printAndLogInfo("The alignment score improved ????")
            else:
                currentAlignment = copyLastAlignment
                printAndLogInfo("The alignment score didn't improve ????")
                printAndLogInfo("---------------------------------------")
                # Hago el segundo filtrado (Saco la secuencia que tenga mas gaps de todo el alineamiento)
                printAndLogInfo("ALGORITHM 2 - Filter sequence with most gaps:")
                # aca se caga, currentAligmnet no es un alineamiento sino secuencias sin gaps.
                improve, currentAlignment, currentScore = excuteSecondAlgorithm(
                    copyLastAlignment, querySeq, homologousSequences, currentScore, env_variables)
                querySeq = find_alignment_by_header(
                    currentAlignment, query_sequence_header)

                if(improve):
                    bestAlignment = currentAlignment
                    bestScore = currentScore
                    lastScore = currentScore
                    lastAlignment = currentAlignment
                    printAndLogInfo("The alignment score improved ????")
                else:
                    # Como no mejoro mas con ninguno de los filtrados termino con la busqueda
                    better = False
                    printAndLogInfo("The alignment score didn't improve ????")
                    printAndLogInfo("---------------------------------------")

        printAndLogInfo("The best alignment obtained contains " +
                        str(len(bestAlignment)) + " sequences")
        printAndLogInfo("The final score is " + str(bestScore))
        printAndLogInfo("---------------------------------------")
        # Genero el ??rbol filogen??tico
        generateTree(bestAlignment, outputDir)
        printAndLogInfo("---------------------------------------")
        # Exportamos el alineamiento final y el ??rbol 
        exportFinalAlignment(bestAlignment, outputDir)
        printAndLogInfo("---------------------------------------")
        #Movemos el log
        os.replace(TEMP_DIR + "/alignment_optimizer.log", outputDir + "/alignment_optimizer.log")
        #Limpiamos los archivos temporales
        [f.unlink() for f in pathlib.Path(TEMP_DIR).glob("*") if f.is_file()] 

    else:
        printAndLogInfo('Current amount of sequences provided is ' + str(len(currentAlignment)) +
                        ' and it is less than the minimum of sequences established for the alignment')


def loadFile(filename):
    # Validar path
    # Validar formato
    # Validar que sea la cantidad de secuencias sea mayor al valor minimo
    # Obtengo el archivo con el alineamiento inicial
    checkIsValidPath(filename)
    printAndLogInfo("Getting Alignments")
    if filename.lower().endswith(('.fasta')):
        with open(filename, "r") as handle:
            originalAlignment = AlignIO.read(handle, "fasta")
            if (any(originalAlignment)):
                printAndLogInfo(str(len(originalAlignment)) +
                                " aligned sequences were obtained correctly")
                return originalAlignment
            else:
                raise Exception("Invalid file format")
    else:
        raise Exception("Invalid extension file")


def checkIsValidPath(filename):
    if not os.path.isfile(filename):
        raise Exception("Invalid file path")


def trimPurifyingSequences(anAlignment, env_variables):
    removeStart = env_variables[PURIFY_START]
    removeEnd = env_variables[PURIFY_END]
    if removeStart > 0 or removeEnd > 0:
        lenSeq = len(anAlignment[0].seq)
        # Primero tengo que revisar que el largo de la cadena es mayor a las secuencias que voy a cortar del incio y el final
        if lenSeq > (removeStart + removeEnd):
            printAndLogInfo("Cutting the purifying sequences")
            # Corto las secuencias
            for i in range(0, len(anAlignment)):
                seq = anAlignment[i].seq
                start = removeStart
                end = len(anAlignment[i].seq) - removeEnd
                seq = seq[start:end]
                anAlignment[i].seq = Seq(seq)
            printAndLogInfo("The alignment was cut " + str(removeStart) +
                            " positions at the beginning and " + str(removeEnd) + " positions at the end")
            printAndLogInfo("---------------------------------------")
        else:
            printAndLogCritical(
                "The length of the alignment is shorter than the total positions to be cut from the purifiers.")
            sys.exit()
    return anAlignment


def getHomologousSequences(querySeq, sequences, env_variables, homologous_path):
    admit_homologous = env_variables[ADMIT_HOMOLOGOUS]
    if admit_homologous == 1:
        printAndLogInfo("Getting Homologous Sequences")
        db_hs = env_variables[DB_HOMOLOGOUS_SEQUENCES]
        if db_hs == 0:
            homologousSequences = getHomologousSequencesOrderedByMaxScore(
                querySeq, env_variables)
        elif db_hs == 1:
            path = None

            if homologous_path is not None:
                path = homologous_path

            if path:
                checkIsValidPath(path)
                homologousSequences = getHomologousSequencesForFastaOrderByMaxScore(
                    querySeq, path, env_variables)
            else:
                printAndLogCritical("Invalid homologous database path")
                sys.exit()

        response = list(filter(lambda seq: seq.id !=
                        sequences[0].id, homologousSequences))
        for ind in range(1, len(sequences)):
            response = list(filter(lambda seq: seq.id != sequences[ind].id, response))

        printAndLogInfo(str(len(response)) +
                        " homologous sequences were obtained correctly")
        printAndLogInfo("---------------------------------------")
        return response
    else:
        return []


def getHomologousSequencesForFastaOrderByMaxScore(querySeq, path, env_variables):
    try:
        sequences = loadFile(path)
        result = []
        nHomologous = env_variables[N_HOMOLOGOUS_SEQUENCES]
        if len(sequences) < nHomologous:
            nHomologous = len(sequences)
        tmp = sequences[:nHomologous]
        for seq in tmp:
            result.append((seq, pairwise2.align.globalxx(
                querySeq.seq, seq.seq, score_only=True)))
        result.sort(key=takeSecond, reverse=True)
        return [i[0] for i in result]
    except BaseException as err:
        if not isinstance(err, KeyboardInterrupt):
            printAndLogCritical("ERROR: "+str(err))
        else:
            printAndLogCritical("Keyboard interrupt")
        sys.exit()


def filterSequence(filterType, anAlignment, querySeq, homologousSequences, env_variables):
    min_sequences = env_variables[MIN_SEQUENCES]
    admit_homologous = env_variables[ADMIT_HOMOLOGOUS]
    condition1 = len(anAlignment) > min_sequences
    condition2 = len(anAlignment) == min_sequences
    condition3 = admit_homologous == 1
    condition4 = len(homologousSequences) > 0
    if condition1 or (condition2 and condition3 and condition4):

        # Realizo el filtrado correspondiente
        if filterType == 0:
            # Sacola secuencia que mas gap genera en la secuencia query
            seqToRemove = sequenceProvidesMostGaps(anAlignment, querySeq)
        else:
            # Saco la secuencia que tiene mas gaps
            seqToRemove = sequenceWithMostGaps(anAlignment, querySeq, env_variables)

        printAndLogInfo(
            "Sequence {id} has been removed from the alignment.".format(id=seqToRemove.id))
        sequences = list(filter(lambda seq: seq.id != seqToRemove.id, anAlignment))
        printAndLogInfo(str(len(sequences)) + " sequences left.")
        printAndLogInfo("Current Alignment: " + str(len(sequences)))
        return sequences
    else:
        if condition2 and (not condition3):
            # La longitud del alineamiento es igual al minimo de secuencias y no acepta homologas
            printAndLogInfo(
                "The alignment cannot be filtered because the minimum has already been reached")
        elif condition2 and condition3 and (not condition4):
            # La longitud del alineamiento es igual al minimo de secuencias, acepta secuencias homologas, pero no hay mas para agregar
            printAndLogInfo(
                "The alignment cannot be filtered because there are no homologous sequences left to add and the minimum has already been reached")
        return anAlignment


def sequenceProvidesMostGaps(anAlignment, querySeq):
    # Me fijo en que posiciones hay gaps y guardo los indices en una lista
    indices = []
    for i in range(0, len(querySeq.seq)):
        if querySeq.seq[i] == '-':
            indices.append(i)
    # Recorro todas las secuencias y me voy quedando con los indices que esten en la lista
    countSeqToRemove = -1
    for x in range(0, len(anAlignment)):
        seq = anAlignment[x].seq
        count = 0
        for y in indices:
            # Una vez obtenido esas secuencias filtradas, busco la que tenga mayor numero de aminoaciodos y esa es la que voy a eliminar
            if seq[y] != '-':
                count += 1
        if count > countSeqToRemove and querySeq.id != anAlignment[x].id:
            countSeqToRemove = count
            seqToRemove = anAlignment[x]
    return seqToRemove


def sequenceWithMostGaps(anAlignment, querySeq, env_variables):
    if querySeq.id != anAlignment[0].id:
        mostGappedSeq = anAlignment[0]
        start = 1
    else:
        mostGappedSeq = anAlignment[1]
        start = 2
    for ind in range(start, len(anAlignment)):
        if querySeq.id != anAlignment[ind].id:
            mostGappedSeq = mostGapped(
                mostGappedSeq, anAlignment[ind], querySeq, env_variables)
    return mostGappedSeq


def homologousSequenceToAdd(homologousSequences):
    return homologousSequences.pop(0)


def addHomologousSequence(sequences, homologousSequences):
    # Agrega una secuencia que no este y que no haya sido agregada previamente
    seqToAdd = homologousSequenceToAdd(homologousSequences)
    printAndLogInfo(seqToAdd.id + " add to the Alignment")
    nNewSeq = len(seqToAdd.seq)
    nSeqAlignment = len(sequences[0].seq)
    nTotal = nSeqAlignment - nNewSeq
    seqFinal = seqToAdd + Seq("-"*nTotal)
    responseSequences = list(filter(lambda seq: True, sequences))
    responseSequences.append(seqFinal)
    printAndLogInfo("Current Alignment: " + str(len(responseSequences)))
    return responseSequences, seqToAdd


def alignmentHasNMinSequences(sequences, env_variables):
    return len(sequences) >= env_variables[MIN_SEQUENCES]


def mostGapped(aSequence, anotherSequence, aQuerySequence, env_variables):
    # Obtengo la secuencia con mayor cantidad de gaps
    if gaps(aSequence) > gaps(anotherSequence):
        return aSequence
    elif gaps(aSequence) == gaps(anotherSequence):
        aligner = Align.PairwiseAligner()
        aligner.match_score = env_variables[MATCH]
        aligner.mismatch_score = env_variables[MISMATCH]
        alignment1 = aligner.align(aSequence.seq.ungap(), aQuerySequence.seq.ungap())
        alignment2 = aligner.align(anotherSequence.seq.ungap(),
                                   aQuerySequence.seq.ungap())
        if alignment1.score < alignment2.score:
            return aSequence
        else:
            return anotherSequence
    else:
        return anotherSequence


def gaps(aSequence):
    return len(aSequence.seq) - len(aSequence.seq.ungap())


def getUngappedSequences(anAlignment):
    # Obtengo del alineamiento pasado por parametro las secuencias originales
    # Para eso elimino todos los gaps
    baseSequences = []
    for sequence in anAlignment:
        sequence.seq = sequence.seq.ungap()
        baseSequences.append(sequence)
    return baseSequences


def parseScore(aClustalOutputString):
    score = re.search("(?<=Alignment Score )-?\d*", aClustalOutputString).group()
    return int(score)


def generateAlignmentAndCalculateScore(originalSequences, env_variables):
    # Genero el nuevo alineamiento por medio de CLUSTAL
    SeqIO.write(originalSequences, (TEMP_DIR + "/seqs.fasta"), "fasta")
    try:
        command = ClustalwCommandline(
            CLUSTALW_PATH, infile=(TEMP_DIR + "/seqs.fasta"), gapopen=env_variables[GAPOPEN], gapext=env_variables[GAPEXT], matrix=env_variables[MATRIX])
        clusalAlignmentOutput = command()
        score = parseScore(clusalAlignmentOutput[0])
        printAndLogInfo(f"New alignment Score: {score}")
        return score
    except BaseException as err:
        printAndLogCritical("ERROR: " + str(err))
        sys.exit()


def loadCurrentAlignment():

    return AlignIO.read(TEMP_DIR + "/seqs.aln", "clustal")


def generateTree(alignment, outputDir):
    # Genero el ??rbol filogen??tico
    printAndLogInfo("Generating Phylogenetic Tree")
    SeqIO.write(alignment, (TEMP_DIR + "/finalAlignment.fasta"), "fasta")
    command = ClustalwCommandline(
        CLUSTALW_PATH, infile=(TEMP_DIR + "/finalAlignment.fasta"))
    command()
    tree = Phylo.read(TEMP_DIR + "/finalAlignment.dnd", "newick")
    outputTree = outputDir + "/finalAlignmentTree.dnd"
    Phylo.write(tree, outputTree, "newick")
    print("Phylogenetic Tree:")
    Phylo.draw_ascii(tree)
    printAndLogInfo("---------------------------------------")
    printAndLogInfo(f"Phylogenetic Tree exported as {outputTree}")


def getHomologousSequencesOrderedByMaxScore(seqQuery, env_variables):
    # Se pasa como parametro un SeqRecord
    # y devuelve una lista de SeqRecord
    # ordenada de mayor a menor por score con la seqQuery
    result = []
    idsProteins = getIdsHomologousSequences(seqQuery.id)
    if idsProteins:
        idString = ""
        for idProtein in idsProteins:
            idString = idString + idProtein + ","
        Entrez.email = "A.N.Other@example.com"
        idString = idString[0:len(idString)-1]
        try:
            handle = Entrez.efetch(db="protein", id=idString, rettype="gb",
                                   retmode="xml", retmax=env_variables[N_HOMOLOGOUS_SEQUENCES])
            output = Entrez.read(handle)
            if output:
                result = getSequencesOrderedByMaxScore(seqQuery, output)
            return result
        except Exception as err:
            if not isinstance(err, KeyboardInterrupt):
                if isinstance(err,OSError) and err.reason.errno == -3:
                    printAndLogCritical("ERROR: Internet connection error ???????????")
                else:
                    printAndLogCritical("ERROR: "+str(err))
            else:
                printAndLogCritical("Keyboard interrupt")
            sys.exit()
    else:
        return result


def getIdsHomologousSequences(idProtein):
    Entrez.email = "A.N.Other@example.com"
    try:
        handle = Entrez.efetch(db="protein", id=idProtein, rettype='ipg', retmode='xml')
        output = Entrez.read(handle, validate=False)
        result = []
        if 'IPGReport' in output:
            ipgRep = output['IPGReport']
            if 'ProteinList' in ipgRep:
                proteinList = ipgRep['ProteinList']
                if proteinList:
                    for indx in range(0, len(proteinList)):
                        result.append(proteinList[indx].attributes['accver'])
            else:
                printAndLogInfo('Entrez did not return a list of homologous sequences')
        else:
                raise Exception("The header of the query sequence was not found in Entrez")
        return result
    except Exception as err:
        if not isinstance(err, KeyboardInterrupt):
            if isinstance(err,OSError) and err.reason.errno == -3:
                printAndLogCritical("ERROR: Internet connection error ???????????")
            else:
                printAndLogCritical("ERROR: "+str(err))
        else:
            printAndLogCritical("Keyboard interrupt")
        sys.exit()


def getSequencesOrderedByMaxScore(seqQuery, output):
    result = []
    if output:
        for protein in output:
            if protein:
                descriptor = protein["GBSeq_locus"]
                seqHomologous = SeqRecord(Seq(protein["GBSeq_sequence"].upper(
                )), id=descriptor, name=descriptor, description=descriptor)
                result.append((seqHomologous, pairwise2.align.globalxx(
                    seqQuery.seq, seqHomologous.seq, score_only=True)))
        result.sort(key=takeSecond, reverse=True)
    return [i[0] for i in result]


def takeSecond(elem):
    return elem[1]


def generateNewAlignmentAndScore(alignmentFiltered, env_variables):
    ungappedSequences = getUngappedSequences(alignmentFiltered)
    currentScore = generateAlignmentAndCalculateScore(ungappedSequences, env_variables)
    currentAlignment = loadCurrentAlignment()
    return currentAlignment, currentScore


def exportFinalAlignment(bestAlignment, outputDir):
    dir = str(pathlib.Path(outputDir).resolve())
    name = "alignment.fasta"
    filePath = dir + "/" + name
    AlignIO.write(bestAlignment, filePath, "fasta")
    printAndLogInfo(f"Output Alignment exported as {filePath}")


def createOutputDir(inputPath, querySeq):
    parent = str(pathlib.Path(__file__).parent.parent.parent.parent.resolve())
    OUTPUT_PATH = "/OUTPUTS"

    isExist = os.path.exists(parent + OUTPUT_PATH)
    if not isExist:
        os.makedirs(parent + OUTPUT_PATH)

    # Limpiamos el nombre en el caso de que tenga '/' ej: MPVIND2009-KOL/1945/09
    clean_querySeq = querySeq.replace("/", "_")

    inputName = str(pathlib.Path(inputPath).name).split(".")[0]
    timestamp = str(datetime.datetime.now())[:-7].replace(" ", "_")
    folderName = "[" + timestamp + "]" + inputName + "-" + clean_querySeq
    folderPath = parent + OUTPUT_PATH + "/" + folderName
    pathlib.Path(folderPath).mkdir()

    return folderPath


def query_yes_no(question):

    print(f"{question} (y/n):")
    input_value = input()

    if input_value is None:
        raise ValueError(f"invalid default answer: '{input_value}'")
    if input_value.lower() in VALID_QUERY_YES:
        return True
    elif input_value.lower() in NOT_VALID_QUERY_NO:
        return False
    else:
        print("Please respond with 'yes' or 'no' " "(or 'y' or 'n').\n")
        query_yes_no(question)


def create_temp_file():
    parent = str(pathlib.Path(__file__).parent.resolve())
    TEMP = "/temp"

    isExist = os.path.exists(parent + TEMP)
    if not isExist:
        os.makedirs(parent + TEMP)

def printAndLogInfo(msg):
    print(msg)
    log.info(msg)


def printAndLogCritical(msg):
    print(msg)
    log.critical(msg)


def printAndLogError(msg):
    print(msg)
    log.error(msg)


def printAndLogWarning(msg):
    print(msg)
    log.warning(msg)
