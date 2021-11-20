from variables import Var
import functions as func
import logging as log


# ---------------------
# Logs
# ---------------------
log.basicConfig(filename='alignment_optimizer.log', format='%(asctime)s - %(levelname)s: %(message)s', level=log.DEBUG)
log.debug('Probando debug')
log.info('Probando info')
log.warning('Probando warning')
log.error('Probando error')


# ---------------------
# Programa Principal
# ---------------------

def main():

    # Configuro las variables
    # --- Todo esto se lo tenemos que pedir al usuario y comentar cuales son los valores por defecto ---
    func.configureVariables()

    # Cargo el archivo con el alineamiento inicial que me pasa el usuario
    lastAlignment = func.loadFile()

    # Obtengo las secuencias originales
    originalSequences = func.getOriginalSequences(lastAlignment)

    # Genero el alineamiento para obtener el score perteneciente al alineamiento inicial pasado por el usuario
    # Este alineamiento lo descarto, ya que no me sirve
    lastAlignmentAux = func.generateAlignment(originalSequences)
    lastScore = func.calculateScore(lastAlignmentAux)

    # Realizo el filtrado de secuencias del alineamiento inicial pasado por el usuario
    aligmentFiltered = func.filterAlignment(lastAlignment)

    # Obtengo las secuencias originales
    originalSequences = func.getOriginalSequences(aligmentFiltered)

    # Genero el nuevo alineamiento y calculo su score
    currentAlignment = func.generateAlignment(originalSequences)
    currentScore = func.calculateScore(currentAlignment)

    # Genero nuevos alineamientos y sus scores correspondientes 
    # mientras aumente el score actual o llegue al minimo de secuencias
    while(currentScore > lastScore or len(currentAlignment) >= Var().nMinSequences()):
        lastAlignment = currentAlignment
        lastScore = currentScore
        aligmentFiltered = func.filterAlignment(lastAlignment)
        originalSequences = func.getOriginalSequences(aligmentFiltered)
        currentAlignment = func.generateAlignment(originalSequences)
        currentScore = func.calculateScore(currentAlignment)

    # Genero el árbol filogenético y lo retorno
    tree = func.generateTree(currentAlignment)
    return tree


if __name__ == '__main__':
    main()

# Para correrlo:
# python src/project_alignment_optimizer/program/main.py