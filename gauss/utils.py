###########################################################
# Funções auxiliares para testar a convergência dos métodos
###########################################################


def linhas(matriz):
    """
    Verifica se a matriz atende ao critério de convergência das linhas
    :param matriz: list of lists: matriz de qualquer tamanho
    :return: bool: True quando atende o critério das linhas
    """
    for i in range(len(matriz)):
        linha = list(map(abs, matriz[i]))
        pivo = linha[i]
        soma = sum(linha)-pivo
        if soma >= pivo:
            return False

    return True


def sassenfeld(matriz):
    """
    Verifica se a matriz atende ao critério de convergência de Sassenfeld
    :param matriz: list of lists: matriz de qualquer tamanho
    :return: True quando atende o critério de Sassenfeld
    """
    lista_b = []

    for i in range(len(matriz)):
        linha = matriz[i][:]
        pivo = abs(linha[i])
        linha[i] = 0 # remove pivo da lista
        linha_abs = list(map(abs, linha)) # valores absolutos
        soma = sum(linha_abs)

        # multiplica os valores da linha pelos valores de B já encontrados
        for j in range(len(lista_b)):
            soma -= linha_abs[j]
            soma += linha_abs[j]*lista_b[j]

        soma = soma/pivo
        lista_b.append(soma)

        if soma > 1:
            return False

    return True


class itemgetterabs:
    """
    Return a callable object that fetches the given item(s) from its operand.
    After f = itemgetter(2), the call f(r) returns r[2].
    After g = itemgetter(2, 5, 3), the call g(r) returns (r[2], r[5], r[3])
    """
    __slots__ = ('_items', '_call')

    def __init__(self, item, *items):
        item = abs(item)
        items = list(map(abs, items))
        if not items:
            self._items = (item,)
            def func(obj):
                return obj[item]
            self._call = func
        else:
            self._items = items = (item,) + items
            def func(obj):
                return tuple(obj[i] for i in items)
            self._call = func

    def __call__(self, obj):
        return self._call(obj)

    def __repr__(self):
        return '%s.%s(%s)' % (self.__class__.__module__,
                              self.__class__.__name__,
                              ', '.join(map(repr, self._items)))

    def __reduce__(self):
        return self.__class__, self._items


def ordena_matriz(matriz):
    matriz_resultado = []
    matriz_iter = matriz[:]
    for i in range(len(matriz)):
        matriz_iter = sorted(matriz_iter, key=itemgetterabs(i), reverse=True)
        matriz_resultado.append(matriz_iter[0])
        matriz_iter = matriz_iter[1:]

    return matriz_resultado



def parada_jacobi(vetor_xk, vetor_xk_1, E):
    """
    Critério de parada do algoritmo de Jacobi, calcula a distância relativa entre duas iterações
    :param vetor_xk: list: vetor da iteração atual
    :param vetor_xk_1: list: vetor da iteracao anterior
    :param E: float: Épsilon
    :return: bool: True se o algoritmo deve parar
    """
    vetor_dk = []
    vetor_xki = []

    for i in range(len(vetor_xk)):
        vetor_dk.append(abs(vetor_xk[i] - vetor_xk_1[i]))
        vetor_xki.append(abs(vetor_xk[i]))

    vetor_dk.sort(reverse=True)
    vetor_xki.sort(reverse=True)

    try:
        resultado = vetor_dk[0]/vetor_xki[0] < E
    except ZeroDivisionError:
        return False

    return resultado