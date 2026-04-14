#!/usr/bin/env python

from numpy import *
from numpy.linalg import inv, eig, det
import copy

Degree = pi / 180
Identity = array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

# ==============================================================================
def roundC (z):

    return round(z.real, 6) + round(z.imag, 6) * 1j

# ==============================================================================


def length (v):

# Calculate the length of a vector expressed in orthogonal coordinates.


    return sqrt (add.reduce (v * v))

# ==============================================================================

def normalize (v):

# Normalize a vector expressed in orthogonal coordinates.

    v = array (v)
    l = length (v)
    if l != 0: v = divide (v, l)
    else: v = [0., 0., 0.]

    return v

# ==============================================================================

def cross (v1, v2):

    return array ([ v1[1] * v2[2] - v1[2] * v2[1],
            v1[2] * v2[0] - v1[0] * v2[2],
            v1[0] * v2[1] - v1[1] * v2[0]])

# ==============================================================================

def rationalize (x, tolerance=1.0e-4):
# Return two integers so that abs (n/d - x) <= tolerance.

    if "complex" in str(type(x)):
        if abs(x.imag)<tolerance: x = x.real
        else: return 0,0  
    xk	 = x
    k	 = 0
    n	 = (1, 0)	# (n(k), n(k-1))
    d	 = (0, 1)	# (d(k), d(k-1))	
    while 1:
            k = k + 1
            ak = int (xk)
            xk = xk - ak
            n = ak * n[0] + n[1], n[0]
            d = ak * d[0] + d[1], d[0]
            if abs (float (n[0]) / float (d[0]) - x) <= tolerance: break
            xk = 1 / xk
    return n[0], d[0]

# ==============================================================================
def gcd (a, b):

# Return the greatest common divisor of integers a and b.
# By convention gcd (0,0) = 0, and g is always a positive integer.

    if (a, b) != (int (a), int(b)):
        raise TypeError, "gcd needs integer arguments"
    else: a, b = int (a), int(b)

    u = abs (a)
    v = abs (b)

    while v != 0:
        q = u  / v 	# !!! integer division !!!
        t = u - q * v
        u = v
        v = t

    return u

# ==============================================================================

def lcm (a, b):

# Return the least common multiple of integer a and b

    if (a, b) != (int (a), int(b)):
        raise TypeError, "lcm needs integer arguments"
    else: a, b = int (a), int(b)
    return a * b / gcd (a, b)

# ==============================================================================
def gcd3D (a, b, c):

    return gcd(gcd(a,b),c)

# ==============================================================================

def Bezout (a, b):
# Return the Bezout integers u,v such that au+bv = pgcd(a,b)

    if b==0: return (1,0)
    elif a==0: return (0,1)
    else:
        (u,v) = Bezout(b,a%b)
        uf,vf = v,u-(a/b)*v
        if a*uf+b*vf < 0:
            uf,vf = -uf,-vf
        return (uf,vf) # to write a//b in Python3
    
# ==============================================================================

def Bezout3D (a, b, c):
# Return the Bezout integers u,v,w such that au+bv+cw = pgcd(a,b,c)
# and two vectors U, V such that the set of solution is [u,v,w] + nU + mV
    if a==0 and b== 0: sol,U,V = [0,0,1],[1,0,0],[0,1,0]
    elif b==0 and c== 0: sol,U,V = [1,0,0],[0,1,0],[0,0,1]
    elif a==0 and c== 0: sol,U,V =[0,1,0],[1,0,0],[0,0,1]
    else:
        y,z = Bezout(b,c)
        bc = y*b+z*c # = gcd(b,c)
        b1,c1 = b/bc,c/bc 
        y0,z0 = y,z
        x1,Y1 = Bezout(a,bc)
        y1,z1 = Y1*y,Y1*z
        sol = [x1,y1,z1]
        U = [bc,-a*y0,-a*z0]
        V = [0,c1,-b1]

    return array(sol),array(U),array(V)
 
# ==============================================================================   

def ReduceAxe(axe):
    # to be use if axe is of type [nu,nvmnw]
    gcd = gcd3D(axe[0],axe[1],axe[2])
    if gcd: naxe = [axe[0]/gcd,axe[1]/gcd,axe[2]/gcd]
    else: naxe = axe
    return naxe

# ==============================================================================   

def LeftInv(m):
    # left inverse of a non-square matrix
    t = transpose(m)
    im = dot(inv(dot(t,m)),t)
    return im

# ==============================================================================

def IntegriseAxe (v, precision=0.01):
    # en entree une liste et en sortie une liste

    if len(v)>0:
        verif = map (lambda x:int(isinstance(x, int)), v)
        if verif == [1, 1, 1]:
            vres = v
        elif abs(v[0]) < 0.0001 and abs(v[1]) < 0.0001 and abs(v[2]) < 0.0001 :
            vres = [0, 0, 0]
        else:
            vsign = map(sign, v)
            vReal = v
            tolerance = precision

            maximum = max (map (abs, vReal))
            newvReal = vReal / maximum

            while 1:
                vRational = []
                lSquare	 = 0
                for number in newvReal:
                    n, d = rationalize (number, tolerance)
                    vRational.append ([n, d])
                    lSquare = lSquare + pow ((float (n) / float (d)), 2)
                if lSquare == 0:
                    tolerance = tolerance / 2.0
                else: break

            dLCM = lcm (vRational[0][1], vRational[1][1])
            dLCM = lcm (dLCM, vRational[2][1])
            vInteger = array ([dLCM * vRational[0][0] / vRational[0][1],
                                dLCM * vRational[1][0] / vRational[1][1],
                                dLCM * vRational[2][0] / vRational[2][1]])
            cGCD = gcd (vInteger[0], vInteger[1])
            cGCD = gcd (cGCD, vInteger[2])
            vInteger = vInteger / cGCD

            vInteger = array(map (abs, vInteger.tolist ()))
            vres = vsign * vInteger
            try: vres = map(int, vres)
            except:vres = map(roundC, vres)
    else:
        vres = []

    return vres

# QUELQUES BASES POUR CALCULER DES SOUS-GROUPES EN PIOCHANT DES ELEMENTS DANS DES COSETS (BRUTAL)
#========================================================================================================
def IndginG(g0, G):
    # G = groupe de matrice
    res = -1
    i = 0
    for g in G:
        if allclose(g, g0):
            res = i
            break
        i += 1
    return res
#========================================================================================================        

def ConjugateSubGroups(G_H):
    H = G_H[0]
    listConj = []
    for gH in G_H:
        Conj = []
        g = gH[0]
        for h in H:
            ghinvg = dot(g, dot(h, inv(g)))
            Conj.append(ghinvg)
        listConj.append(Conj)
    return listConj

#========================================================================================================
def ListVectTrans(G_H):
    listConj = ConjugateSubGroups(G_H)

    listVect = []
    for conj in listConj:
        for m in conj:
            res = calculeRot(m)
            if not res: res = calculeRot(dot(-1, m))
            if res[0]:
                vect = res[1] # vecteur commun a phase mere et phase fille
                listVect.append(vect)
                break
    return listVect

#========================================================================================================

def TabCompoGroup(G):

    TabCompo = [[-1 for i in range(len(G))] for j in range(len(G))]
    i = 0
    for gi in G:
        j = 0
        for gj in G:
            gij = dot(gi, gj)
            k = IndginG(gij, G)
            TabCompo[i][j] = k
            j += 1
        i += 1

    return TabCompo

#========================================================================================================
def CompoList(ll1, l2):
    # retourne la liste de listes [l1,el] for el in l2
    res = []
    for l1 in ll1:
        for el in l2:
            l = copy.deepcopy(l1)
            l.append(el)
            res.append(l)
    return res
#========================================================================================================

def  AllList(n, m):
    # generate all the possible n liste with integer between 0 and m-1
    lm = range(m)
    l = [[]]
    for i in range(n):
        l = CompoList(l, lm)
    return l

#========================================================================================================

def IsGroup(K, tabCompo):
    flag = 1
    for ligne in tabCompo:
        if 0 not in ligne:
            flag = 0
            break
        for el in ligne:
            if el not in K:
                flag = 0
                break
    return flag

 #========================================================================================================

def FindSubGroups(G_H):
    # G_H = G/H ensemble quotient
    # cherche tous les sous-groupes en prenant un element par coset

    G = []
    for gH in G_H:
        for g in gH: G.append(g)
    orderG = len(G)
    H = G_H[0]
    orderH = len(H)
    TabCompo = TabCompoGroup(G)

    lH = range(orderH)
    lG_H = range(len(G_H))
    alllist = AllList(len(G_H), orderH)
    listSubGroups = []

    p = 0
    for l in alllist:
        setel = [l[i] + i * orderH for i in lG_H]
        tabcompo = [[TabCompo[i][j] for i in setel] for j in setel]
        if IsGroup(setel, tabcompo): listSubGroups.append(setel) # grosse methode brutale parce que je suis nul en math
        p += 1

    listSubGroupsMat = []
    for subgind in listSubGroups:
        listSubGroupsMat.append([G[ind] for ind in subgind])

    return listSubGroupsMat

# ==============================================================================
#   LES FORMULES DE ROTATIONS NE DEPENDENT PAS DU SYST. CRISTALLIN
# ==============================================================================

def calculeRot_axe_brute (mat):
    # modifie le 28/08/09
    # seulement qd angle de rotation est pi

    evalues, evectors = eig(mat)
    compt = 0
    vect = []
    for i in range(len(evalues)):
        if abs(evalues[i] - 1) <= 0.0001 :
                compt = compt + 1
                indice = i
    if compt == 1 :
        vect = evectors[:, indice] # il faut prendre la colonne ici, avec Numeric cetait ligne
        #vectrat = IntegriseAxe(vect)

    return vect

# ==============================================================================

def calculeRot_notsoold (mat):
    # used in the version 1.8 and gave up because of zirconia
    # for all crystals
    # return the angle in Degree

    if det(mat) < 0:  # pensez a le retirer plus tard pour acceler les calculs
        res = []

    else:
        x = (mat[2][1] - mat[1][2])
        y = (mat[0][2] - mat[2][0])
        z = (mat[1][0] - mat[0][1])
        t = trace(mat)
##        if t>3: t=3
##        if t<-1 : t=-1
        ang = arccos(float(t - 1) / 2) # le methode en arctan2(t,r-1) ne marche que pour les cubiques
        ang = round(ang / Degree, 5)
        vect = array([x,y,z])
        if abs(ang-180) < 0.0001:
            lev = eig(mat)[0]
            ind1=0
            for ind in range(len(lev)):
                ev=lev[ind]
                if abs(ev-1)<0.0000001:
                    ind1= ind
                    break
            vect = real(transpose(eig(mat)[1])[ind1])
        if abs(ang) > 0.0001:
            # pour determiner le signe de ang, si negatif on inv l axe de rot
            #-------
            u = array([1, 0, 0])
            v = dot(mat, u) - u
            if (length(v) <= 0.000000001):
                u = array([0, 1, 0])
                v = dot(mat, u) - u
                if (length(v) <= 0.00000001) :
                    u = array([0, 0, 1])
                    v = dot(mat, u) - u
            if det(array([vect, u, v])) < 0 : vect = -vect
            #-------
        elif ang < 0.0001:
            ang = 0
            vect = [0,0,1]

        res = [ang, vect]
        return res
    
# ==============================================================================

def calculeRot (mat):
    # for all crystals
    # return the angle in Degree

    if det(mat) < 0:  # pensez a le retirer plus tard pour acceler les calculs
        res = []

    else:
        x = (mat[2][1] - mat[1][2])
        y = (mat[0][2] - mat[2][0])
        z = (mat[1][0] - mat[0][1])
        ang = arccos(float(trace(mat) - 1) / 2) # le methode en arctan2(t,r-1) ne marche que pour les cubiques
        ang = round(ang / Degree, 5)
        if abs(x) + abs(y) + abs(z) > 0.00001 and abs(ang - 180) > 0.000001:
            vectrat = IntegriseAxe([x, y, z])
            # pour determiner le signe de ang, si negatif on inv l axe de rot
            #-------
            u = array([1, 0, 0])
            v = dot(mat, u) - u
            if (length(v) <= 0.00001):
                u = array([0, 1, 0])
                v = dot(mat, u) - u
                if (length(v) <= 0.00001) :
                    u = array([0, 0, 1])
                    v = dot(mat, u) - u
            if det(array([vectrat, u, v])) < 0 : vectrat = map(lambda x:-x, vectrat)
            #-------
        else:
            if allclose(mat, Identity):
                ang = 0
                vectrat = [1, 0, 0]
            else:
                ang = 180
                vectrat=IntegriseAxe(calculeRot_axe_brute(mat))

        res = [ang, vectrat]
        return res


# ==============================================================================

def calculeMatRot(angD, axe):
    # ang en degree et axe dans le repere du cristal

    x, y, z = normalize(axe)
    ang = angD * Degree
    cc = (1 - cos(ang))
    s = sin(ang)
    MRot = array([[1 + cc * (x ** 2 - 1), -z * s + cc * x * y, y * s + cc * x * z],
                      [z * s + cc * x * y, 1 + cc * (y ** 2 - 1), -x * s + cc * y * z],
                      [-y * s + cc * x * z, x * s + cc * y * z, 1 + cc * (z ** 2 - 1)]])

    return MRot

    
# ==============================================================================

def Polar(F):
    # polar decomposition F = QU with Q rotation and U diagonal
    
    F2 = dot(transpose(F),F)
    v1,v2,v3= eig(F2)[0]
    T = transpose(eig(F2)[1])
    e1,e2,e3 = T[0],T[1],T[2]
    # print "valeurs propres au carre:", v1,v2,v3
    U = array([[sqrt(round(v1,3)),0,0],[0,sqrt(round(v2,3)),0],[0,0,sqrt(round(v3,3))]])
    Q = dot(F,inv(U))
    return Q, U



# ==============================================================================
# ==============================================================================

if __name__ == '__main__':
    print "Bezout"
    print Bezout(25,71)
    print "Bezout3D"
    print Bezout3D(65,297,-224)
    print
    axe = array([2.983647, 5.289, 3.12458])
    axe = array([ 1.01954006, -0.97924159, -0.04029847])
    print axe
    print IntegriseAxe(axe, 0.1)
    print IntegriseAxe(axe, 0.01)
    print IntegriseAxe(axe, 0.001)
    axe = array([0, 2, -1])
    print axe
    print IntegriseAxe(axe)
    va = array([1, 1, 0])
    vb = array([1, -1, 0])
    m1 = array([[-0., 1., -0.],
       [ 1., -0., -0.],
       [-0., 0., -1.]])
    m2 = array([[ 0., -0., 1.],
       [ 0., -1., 0.],
       [ 1., -0., -0.]])
    m3 = array([[-0., 0., -1.],
       [ 0., -1., 0.],
       [-1., 0., 0.]])
    m4 = array([[ 0., -1., 0.],
       [-1., 0., 0.],
       [-0., 0., -1.]])
    print eig(m1)
    print eig(m2)
    print eig(m3)
    print eig(m4)
    print calculeRot(m1)
    print calculeRot(m2)
    print calculeRot(m3)
    print calculeRot(m4)
    m = array([[ 1, -1, 0], [ 1, 0, 0], [ 0, 0, 1]])
    print calculeRot(m)
    m = ([[ 0, 0, -1], [ 0, 1, 0], [ 1, 0, 0]])
    print calculeRot(m)
    m = array([[  1.00000000e+00, 0.00000000e+00, -2.31674668e-17],
       [  1.00000000e+00, -1.00000000e+00, -7.28967355e-18],
       [  0.00000000e+00, 0.00000000e+00, -1.00000000e+00]])
    print m.round(5)
    print det(m)
    print calculeRot(m)
    m = calculeMatRot(10, array([1, -1, 0]))
    print m
    print IntegriseAxe(dot(m, array([1, 0, 0])))
    print IntegriseAxe(dot(m, array([0, 1, 0])))
    print IntegriseAxe(array([0.356802, -0.244023, 0.901746]),0.05)









