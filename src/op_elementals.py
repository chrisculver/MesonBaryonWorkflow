from WickContractions.ops.operator import *
from WickContractions.ops.elemental import *
from WickContractions.ops.quarks import *
from WickContractions.ops.commuting import *

cIdx = lambda c : 'c_{' + str(c) + '}'
sIdx = lambda s : 's_{' + str(s) + '}'
epsilonTensor = lambda c0 : EpsilonTensor([cIdx(c0), cIdx(c0+1), cIdx(c0+2)])
gammaMatrix = lambda name, s0: SpinMatrix(name, [sIdx(s0), sIdx(s0+1)])
baryonSpinTensor = lambda name, s0: SpinMatrix(name, [sIdx(s0), sIdx(s0+1), sIdx(s0+2)])
deltaTensor = lambda c0 : IndexedObject('\\delta', [cIdx(c0), cIdx(c0+1)])

nextSpinIdx = 0
nextColorIdx= 0

def quark(anti, f, t, x):
    global nextSpinIdx, nextColorIdx
    q=Quark(anti, f, sIdx(nextSpinIdx), cIdx(nextColorIdx), t, x)
    nextSpinIdx+=1
    nextColorIdx+=1
    return q

def baryonSource(terms, x, g):
    elementals=[]
    for flavors,coef in terms.items():
        fs=[f for f in flavors]
        elementals.append(baryonSourceElemental(coef, fs, x, g))
    return Operator(elementals)

def baryonSourceElemental(coef, flavors, x, gammaName):
    eTensor = epsilonTensor(nextColorIdx)
    sTensor = baryonSpinTensor(gammaName, nextSpinIdx)

    q0 = quark(True, flavors[0], 't_i', x)
    q1 = quark(True, flavors[1], 't_i', x)
    q2 = quark(True, flavors[2], 't_i', x)

    return ElementalOperator(coef,
        [eTensor, sTensor],
        [q0,q1,q2]
    )

def baryonSink(terms, x, g):
    elementals=[]
    for flavors,coef in terms.items():
        fs=[f for f in flavors]
        elementals.append(baryonSinkElemental(coef, fs, x, g))
    return Operator(elementals)


def baryonSinkElemental(coef, flavors, x, gammaName):
    eTensor = epsilonTensor(nextColorIdx)
    sTensor = baryonSpinTensor(gammaName, nextSpinIdx)

    q0 = quark(False, flavors[0], 't_f', x)
    q1 = quark(False, flavors[1], 't_f', x)
    q2 = quark(False, flavors[2], 't_f', x)

    return ElementalOperator(coef,
        [eTensor, sTensor],
        [q0,q1,q2]
    )

def mesonSource(terms, x, g):
    elementals=[]
    for flavors,coef in terms.items():
        fs=[f for f in flavors]
        elementals.append(mesonSourceElemental(coef, fs, x, g))
    return Operator(elementals)

def mesonSourceElemental(coef, flavor, x, gammaName):
    dTensor = deltaTensor(nextColorIdx)
    sTensor = gammaMatrix(gammaName, nextSpinIdx)

    q0=quark(True, flavor[0], 't_i', x)
    q1=quark(False, flavor[1], 't_i', x)

    return ElementalOperator(
        coef,
        [dTensor, sTensor],
        [q0,q1]
    )

def mesonSink(terms, x, g):
    elementals=[]
    for flavors,coef in terms.items():
        fs=[f for f in flavors]
        elementals.append(mesonSinkElemental(coef, fs, x, g))
    return Operator(elementals)

def mesonSinkElemental(coef, flavor, x, gammaName):
    dTensor = deltaTensor(nextColorIdx)
    sTensor = gammaMatrix(gammaName, nextSpinIdx)

    q0=quark(True, flavor[0], 't_f', x)
    q1=quark(False, flavor[1], 't_f', x)

    return ElementalOperator(
        coef,
        [dTensor, sTensor],
        [q0,q1]
    )

def mesonBaryonSource(mesonData, baryonData):
    mesonIsospin=mesonData['isospin']
    xm=mesonData['x']
    gm=mesonData['g']
    baryonIsospin=baryonData['isospin']
    xb=baryonData['x']
    gb=baryonData['g']

    elementals=[]
    for mesonFs, mesonC in mesonIsospin.items():
        for baryonFs, baryonC in baryonIsospin.items():
            mfs=[f for f in mesonFs]
            bfs=[f for f in baryonFs]
            elementals.append(MesonBaryonSourceElemental(mesonC*baryonC, 
                mfs, bfs, xm, xb, gm, gb
            ))
    return Operator(elementals)

def MesonBaryonSourceElemental(coef, mesonFs, baryonFs, xm, xb, gm, gb):
    eTensor = epsilonTensor(nextColorIdx)
    sTensorB = baryonSpinTensor(gb, nextSpinIdx)

    q0 = quark(True, baryonFs[0], 't_i', xb)
    q1 = quark(True, baryonFs[1], 't_i', xb)
    q2 = quark(True, baryonFs[2], 't_i', xb)

    dTensor = deltaTensor(nextColorIdx)
    sTensorM = gammaMatrix(gm, nextSpinIdx)

    q3=quark(True, mesonFs[0], 't_i', xm)
    q4=quark(False, mesonFs[1], 't_i', xm)


    return ElementalOperator(coef,
        [eTensor, dTensor, sTensorB, sTensorM],
        [q0,q1,q2,q3,q4]
    )

def mesonBaryonSink(mesonData, baryonData):
    mesonIsospin=mesonData['isospin']
    xm=mesonData['x']
    gm=mesonData['g']
    baryonIsospin=baryonData['isospin']
    xb=baryonData['x']
    gb=baryonData['g']

    elementals=[]
    for mesonFs, mesonC in mesonIsospin.items():
        for baryonFs, baryonC in baryonIsospin.items():
            mfs=[f for f in mesonFs]
            bfs=[f for f in baryonFs]
            elementals.append(MesonBaryonSinkElemental(mesonC*baryonC, 
                mfs, bfs, xm, xb, gm, gb
            ))
    return Operator(elementals)

def MesonBaryonSinkElemental(coef, mesonFs, baryonFs, xm, xb, gm, gb):
    eTensor = epsilonTensor(nextColorIdx)
    sTensorB = baryonSpinTensor(gb, nextSpinIdx)

    q0 = quark(False, baryonFs[0], 't_f', xb)
    q1 = quark(False, baryonFs[1], 't_f', xb)
    q2 = quark(False, baryonFs[2], 't_f', xb)

    dTensor = deltaTensor(nextColorIdx)
    sTensorM = gammaMatrix(gm, nextSpinIdx)

    q3=quark(True, mesonFs[0], 't_f', xm)
    q4=quark(False, mesonFs[1], 't_f', xm)


    return ElementalOperator(coef,
        [eTensor, dTensor, sTensorB, sTensorM],
        [q0,q1,q2,q3,q4]
    )