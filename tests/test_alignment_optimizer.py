import pytest

from src.project_alignment_optimizer.program import functions


def test_filtrado():

    alignment = ["AA--C","A---C","AABDC","---DC","-ABD-"]
    # querySequence = 0 por defecto, es decir "AA--C"

    response = functions.filterAlignment(alignment)
    # Saca la secuencia "---DC" por ser la de mayor cantidad de gaps y menor score respecto de la query
    assert response == ["AA--C","A---C","AABDC","-ABD-"]

    response = functions.filterAlignment(alignment)
    # Saca la secuencia "A---C" por ser la de mayor cantidad de gaps
    assert response == ["AA--C","AABDC","-ABD-"]
