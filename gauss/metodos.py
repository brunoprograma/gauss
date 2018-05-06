from .utils import *


def jacobi(matriz_completa, vetor_xk=None, E=0.01, maxiter=1000):
    """

    :param matriz_completa:
    :param vetor_xk:
    :param E:
    :param maxiter:
    :return:
    """
    print("Método Gauss-Jacobi")
    matriz_completa_ordenada = ordena_matriz(matriz_completa)
    matriz_A = []
    vetor_b = []

    # separa a matriz completa Ax = b em uma matriz A e um vetor b e verifica se a matriz tem tamanho válido,
    # ou seja, todas as suas linhas tem o mesmo tamanho
    tam = None
    valid = True
    for linha in matriz_completa_ordenada:
        vetor_b.append(linha[-1])
        matriz_A.append(linha[:-1])

        if tam is not None and tam != len(linha):
            valid = False
            break

        tam = len(linha)

    if not vetor_xk:
        vetor_xk = [0 for i in range(len(matriz_A[0]))] # vetor de zeros
    elif valid and (len(vetor_xk) != tam): # tamanho do vetor inicial inválido
        valid = False

    if valid:
        conv = linhas(matriz_A)
        vetor_xk_1 = vetor_xk[:]

        if conv:
            print("De acordo com o critério das linhas o método DEVERÁ CONVERGIR.")
        else:
            print("ATENÇÃO! De acordo com o critério das linhas o método PODERÁ NÃO CONVERGIR!")

        vetor_x = []
        iteracoes = 0
        while not parada_jacobi(vetor_xk, vetor_xk_1, E):
            iteracoes += 1
            vetor_xk_1 = vetor_xk[:]  # cópia de vetores
            novo_vetor_xk = []

            for i in range(len(matriz_A)):
                if len(vetor_x) < i+1:
                    linha = matriz_A[i][:]  # copiamos a linha
                    pivo = linha[i]
                    vetor_x.append(list(map(lambda x: -x/pivo, linha)))  # invertemos os sinais
                    vetor_x[i][i] = 0
                    vetor_x[i].append(vetor_b[i]/pivo)  # valor de b

                x = vetor_x[i][-1] # valor de b
                for j in range(len(vetor_x[i])-1):
                    x += vetor_x[i][j] * vetor_xk[j]
                novo_vetor_xk.append(x)

            vetor_xk = novo_vetor_xk[:]

            if iteracoes >= maxiter:
                print("O sistema não convergiu após {} iterações".format(maxiter))
                return

        print("Solução {} =".format(["x{}".format(i+1) for i in range(len(vetor_xk))]), vetor_xk, "com {} iterações e".format(iteracoes), "Epsilon (E) = {}".format(E))
    else:
        print("As entradas informadas são inválidas!")


def seidel(matriz_completa, vetor_xk=None, E=0.01, maxiter=1000):
    """

    :param matriz_completa:
    :param vetor_xk:
    :param E:
    :param maxiter:
    :return:
    """
    print("Método Gauss-Seidel")
    matriz_completa_ordenada = ordena_matriz(matriz_completa)
    matriz_A = []
    vetor_b = []

    # separa a matriz completa Ax = b em uma matriz A e um vetor b e verifica se a matriz tem tamanho válido,
    # ou seja, todas as suas linhas tem o mesmo tamanho
    tam = None
    valid = True
    for linha in matriz_completa_ordenada:
        vetor_b.append(linha[-1])
        matriz_A.append(linha[:-1])

        if tam is not None and tam != len(linha):
            valid = False
            break

        tam = len(linha)

    if not vetor_xk:
        vetor_xk = [0 for i in range(len(matriz_A[0]))] # vetor de zeros
    elif valid and (len(vetor_xk) != tam): # tamanho do vetor inicial inválido
        valid = False

    if valid:
        conv = sassenfeld(matriz_A)
        vetor_xk_1 = vetor_xk[:]

        if conv:
            print("De acordo com o critério de Sassenfeld o método DEVERÁ CONVERGIR.")
        else:
            print("ATENÇÃO! De acordo com o critério de Sassenfeld o método PODERÁ NÃO CONVERGIR!")

        vetor_x = []
        iteracoes = 0
        while not parada_jacobi(vetor_xk, vetor_xk_1, E):
            iteracoes += 1
            vetor_xk_1 = vetor_xk[:]  # cópia de vetores

            for i in range(len(matriz_A)):
                if len(vetor_x) < i+1:
                    linha = matriz_A[i][:]  # copiamos a linha
                    pivo = linha[i]
                    vetor_x.append(list(map(lambda x: -x/pivo, linha)))  # invertemos os sinais
                    vetor_x[i][i] = 0
                    vetor_x[i].append(vetor_b[i]/pivo)  # valor de b

                x = vetor_x[i][-1] # valor de b
                for j in range(len(vetor_x[i])-1):
                    x += vetor_x[i][j] * vetor_xk[j]
                vetor_xk[i] = x

            if iteracoes >= maxiter:
                print("O sistema não convergiu após {} iterações".format(maxiter))
                return

        print("Solução {} =".format(["x{}".format(i+1) for i in range(len(vetor_xk))]), vetor_xk, "com {} iterações e".format(iteracoes), "Epsilon (E) = {}".format(E))
    else:
        print("As entradas informadas são inválidas!")