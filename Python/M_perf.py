#!/usr/bin/env python3
import numpy
import scipy.io
from Funcion_Qr_1 import Funcion_Qr_1
from Funcion_Qp_1 import Funcion_Qp_1
from time import perf_counter
from statistics import mean

def Inversa(x,y):
    p=numpy.array([[x],[y]])
    Qp=Funcion_Qp_1(p)
    Qr=Funcion_Qr_1(p)
    Posq=1023-(Qp*7.815144409)
    Posq=int(numpy.round(Posq))
    Posr=512-(Qr*195.3786102)
    Posr=int(numpy.round(Posr))
    return  Posr,Posq

ti=scipy.io.loadmat('Walk_50.mat')['xy']
#ti=scipy.io.loadmat('Circle.mat')['xy']
#ti=scipy.io.loadmat('Ellipse.mat')['xy']
x=ti[0]
y=ti[1]
c=0
tiempo=[0]*10
while c<10:
    a=perf_counter()
    for idx in range(0,50):
        q1,q2=Inversa(x[idx],y[idx])
    b=perf_counter()
    tiempo[c]=b-a
    c=c+1
print(mean(tiempo))


