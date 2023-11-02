import numpy as np
from math import log
from math import exp

def PrimeraDerivadaAdelante2(f,df,a,h):
    h1=h
    D11 = (f(a + h1) - f(a)) / h1
    h2=h/2
    D12 = (f(a + h2) - f(a)) / h2
    Error1= abs(D11-df(a))
    Error2= abs(D12-df(a))

    #Tabla de Resultados
    print(' \n   h        D2 Aprox        Error')
    print('{:4.2e}    {:8.8f}     {:4.4e}'.format(h1,D11,Error1))
    print('{:4.2e}    {:8.8f}     {:4.4e}'.format(h2,D12,Error2))
    print(' ')

    #Calculo del Orden de Convergencia
    p = log(Error1 / Error2) / log(h1 / h2)
    print('1ra Derivada Aprox Adelante 2 pts.')
    print('El orden de convergencia es:', p)
    print(' ')

    return D11

def PrimeraDerivadaAtras2(f,df,a,h):
    h1=h
    D11 = (f(a) - f(a-h1)) / h1
    h2=h/2
    D12 = (f(a) - f(a-h2)) / h2
    Error1= abs(D11-df(a))
    Error2= abs(D12-df(a))

    #Tabla de Resultados
    print(' \n   h        D2 Aprox        Error')
    print('{:4.2e}    {:8.8f}     {:4.4e}'.format(h1,D11,Error1))
    print('{:4.2e}    {:8.8f}     {:4.4e}'.format(h2,D12,Error2))
    print(' ')

    #Calculo del Orden de Convergencia
    p = log(Error1 / Error2) / log(h1 / h2)
    print('1ra Derivada Aprox Atrás 2 pts.')
    print('El orden de convergencia es:', p)
    print(' ')

    return D11

def PrimeraDerivadaCentrada3(f,df,a,h):
    h1=h
    D11 = (f(a + h1) - f(a - h1)) / (2*h1)
    h2=h/2
    D12 = (f(a + h2) - f(a - h2)) / (2*h2)
    Error1= abs(D11-df(a))
    Error2= abs(D12-df(a))

    #Tabla de Resultados
    print(' \n   h        D2 Aprox        Error')
    print('{:4.2e}    {:8.8f}     {:4.4e}'.format(h1,D11,Error1))
    print('{:4.2e}    {:8.8f}     {:4.4e}'.format(h2,D12,Error2))
    print(' ')

    #Calculo del Orden de Convergencia
    p = log(Error1 / Error2) / log(h1 / h2)
    print('1ra Derivada Aprox Centrada 3 pts.')
    print('El orden de convergencia es:', p)
    print(' ')

    return D11

def SegundaDerivadaCentrada3(f,d2f,a,h):
    h1=h
    D21 = (f(a + h1) - 2 * f(a) + f(a - h1)) / h1 ** 2
    h2=h/10
    D22 = (f(a + h2) - 2 * f(a) + f(a - h2)) / h2 ** 2
    Error1= abs(D21-d2f(a))
    Error2= abs(D22-d2f(a))

    #Tabla de Resultados
    print(' \n   h        D2 Aprox        Error')
    print('{:4.2e}    {:8.8f}     {:4.4e}'.format(h1,D21,Error1))
    print('{:4.2e}    {:8.8f}     {:4.4e}'.format(h2,D22,Error2))
    print(' ')

    #Calculo del Orden de Convergencia
    p = log(Error1 / Error2) / log(h1 / h2)
    print('2da Derivada Aprox Centrada 3 pts.')
    print('El orden de convergencia es:', p)
    print(' ')

    return D22

def SegundaDerivadaCentrada5(f,d2f,a,h):
    h1=h
    D21 = (-f(a+2*h1)+16*f(a+h1)-30*f(a)+16*f(a-h1)-f(a-2*h1))/(12*h1**2)
    h2=h/10
    D22 = (-f(a+2*h2)+16*f(a+h2)-30*f(a)+16*f(a-h2)-f(a-2*h2))/(12*h2**2)
    Error1= abs(D21-d2f(a))
    Error2= abs(D22-d2f(a))

    #Tabla de Resultados
    print(' \n   h        D2 Aprox        Error')
    print('{:4.2e}    {:8.8f}     {:4.4e}'.format(h1,D21,Error1))
    print('{:4.2e}    {:8.8f}     {:4.4e}'.format(h2,D22,Error2))
    print(' ')

    #Calculo del Orden de Convergencia
    p = log(Error1 / Error2) / log(h1 / h2)
    print('2da Derivada Aprox Centrada 5 pts.')
    print('El orden de convergencia es:', p)
    print(' ')

    return D22

def ReglaTrapecio(g,a,b):

    h = b - a
    IT = 0.5 * h * (g(a) + g(b))

    return IT

def ReglaSimpson(g,a,b):

    h = (b - a)/2
    IS =  (h/3) * (g(a) +4*g(a+h) + g(a+2*h))

    return IS

def ReglaTrapecioCompuesta(g,a,b,m=10):

    h = (b-a)/m
    xi = np.linspace(a,b,m+1)

    Suma = 0  #Inicialice suma a cero
    for k in range(0,m):
        Si = ReglaTrapecio(g,xi[k],xi[k+1])
        Suma = Suma + Si  #Actualice la suma

    return Suma

def ReglaSimpsonCompuesta(g,a,b,m=10):

    h = (b-a)/(2*m)
    xi = np.linspace(a,b,2*m+1)

    Suma = 0  #Inicialice suma a cero
    for k in range(0,m):
        Si = ReglaSimpson(g,xi[2*k],xi[2*k+2])
        Suma = Suma + Si  #Actualice la suma

    return Suma

def SimpsonCompuesta(f,a,b,m=4):
    h = (b-a)/(2*m)
    x = np.linspace(a,b,2*m+1)

    S = f(a)+f(b)  #Sume los valores en los extremos del intervalo

    for i in range(1,m+1):  #Los nodos impares tienen un peso de 4
        S=S+4*f(x[2*i-1][0])

    for i in range(1,m):    #Los nodos pares tienen un peso de 2
        S=S+2*f(x[2*i][0])



    IS = h*S/3

    #print('Regla de Simpson Compuesta')
    #print('Subintervalos  = {:3d}   h = {:3.3e}'.format(m,h))
    #print(IS, '\n')

    return IS,h

def TablaReglasCompuesta(f,a,b,ms,F,R):
    Is = np.zeros(np.shape(ms))  #Vector de aproximaciones
    hs = np.zeros(np.shape(ms))  #Vector de hs
    Es = np.zeros(np.shape(ms))  #Vector de errores
    p  = np.zeros(np.shape(ms))  #Vector de ordenes de convergencia

    Iexac= F(b)-F(a)    #valor exacto de la integral

    m = ms[0][0]    #Integral aproximada 1er valor de m
    if R=='T':
        [Is[0],hs[0]] = ReglaTrapecioCompuesta(f, a, b, m)
    else:
        [Is[0], hs[0]] = ReglaSimpsonCompuesta(f, a, b, m)
    Es[0]=abs(Is[0]-Iexac)

    #Encabezado de la Tabla de Resultados
    print('h             Iaprox           Error       p ')
    print('{:2.4e}   {:3.6f}    {:2.4e} '.format(hs[0][0],Is[0][0],Es[0][0]))
    for i in range(1,len(ms)):
        m= ms[i][0]
        if R=='T':
            [Is[i],hs[i]]=ReglaTrapecioCompuesta(f, a, b, m)
        else:
            [Is[i], hs[i]]=ReglaSimpsonCompuesta(f, a, b, m)
        Es[i]=abs(Is[i]-Iexac)
        p[i] =log(Es[i]/Es[i-1])/log(hs[i]/hs[i-1])
        #Imprima el resto de la tabla
        print('{:2.4e}   {:3.6f}   {:2.4e}    {:2.3f}'.format(hs[i][0], Is[i][0], Es[i][0],p[i][0]))

    print(' ')
    #Devuelve los valores aproximados y los tamaños de pasos
    return Is, hs



