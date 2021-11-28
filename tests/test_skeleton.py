import pytest
#from project_alignment_optimizer.program.variables import Var
#import project_alignment_optimizer.program.functions as func
import os

from project_alignment_optimizer.main import main
from project_alignment_optimizer.fib_model import fib
from project_alignment_optimizer import fib_model
from project_alignment_optimizer.program import functions,variables

__author__ = "ItuFede"
__copyright__ = "ItuFede"
__license__ = "MIT"


def test_fib():
    """API Tests"""
    assert fib(1) == 1
    assert fib(2) == 1
    assert fib(7) == 13
    with pytest.raises(AssertionError):
        fib(-10)


def test_main(capsys):
    """CLI Tests"""
    # capsys is a pytest fixture that allows asserts agains stdout/stderr
    # https://docs.pytest.org/en/stable/capture.html
    main(["7"])
    captured = capsys.readouterr()
    assert "The 7-th Fibonacci number is 13" in captured.out

def test_file_with_minimun_amount_of_seqs_doesnt_optimize():
    #minimo de 35 secuencias. el archivo tiene 35. no tiene sacar ninguna secuencia y tiene que quedar = el score
    
    functions.configureVariables()
    #ver como tomar el fasta de tests 
    file = os.read('/alintest1.fasta')
    # Cargo el archivo con el alineamiento inicial que me pasa el usuario
    lastAlignment = functions.loadFile(file)
    #verifico que tengo un score 
    #alineo
    #optimizo
    #verifico que no cambie el score 
    assert functions.calculateScore(lastAlignment) == 3

def test_aligment_get_optimized():

    #se importa un archivo y se mejora el socre
    file = os.read('/alintest2.fasta')