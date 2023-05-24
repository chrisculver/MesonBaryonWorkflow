from WickContractions.ops.operator import *
from WickContractions.ops.elemental import *
from WickContractions.ops.quarks import *
from WickContractions.ops.commuting import *

cIdx = lambda c : 'c_{' + str(c) + '}'
sIdx = lambda s : 's_{' + str(s) + '}'
epsilonTensor = lambda c0 : EpsilonTensor([cIdx(c0), cIdx(c0+1), cIdx(c0+2)])
spinTensor = lambda name, s0: SpinMatrix(name, [sIdx(s0), sIdx(s0+1)])
deltaTensor = lambda c0 : IndexedObject('\\delta', [cIdx(c0), cIdx(c0+1)])

nextSpinIdx = 0
nextColorIdx= 0

def quark(anti, f, t, x):
    global nextSpinIdx, nextColorIdx
    q=Quark(anti, f, sIdx(nextSpinIdx), cIdx(nextColorIdx), t, x)
    nextSpinIdx+=1
    nextColorIdx+=1
    return q

def baryonSource(coef, flavors, x, gammaName):
    eTensor = epsilonTensor(nextColorIdx)
    sTensor = spinTensor(gammaName, nextSpinIdx)

    q0 = quark(True, flavors[0], 't_i', x)
    q1 = quark(True, flavors[1], 't_i', x)
    q2 = quark(True, flavors[2], 't_i', x)

    return ElementalOperator(coef,
        [eTensor, sTensor],
        [q0,q1,q2]
    )

def baryonSink(coef, flavors, x, gammaName):
    eTensor = epsilonTensor(nextColorIdx)
    sTensor = spinTensor(gammaName, nextSpinIdx)

    q0 = quark(False, flavors[0], 't_f', x)
    q1 = quark(False, flavors[1], 't_f', x)
    q2 = quark(False, flavors[2], 't_f', x)

    return ElementalOperator(coef,
        [eTensor, sTensor],
        [q0,q1,q2]
    )

def mesonSource(coef, flavor, x, gammaName):
    dTensor = deltaTensor(nextColorIdx)
    sTensor = spinTensor(gammaName, nextSpinIdx)

    q0=quark(True, flavor[0], 't_i', x)
    q1=quark(False, flavor[1], 't_i', x)

    return ElementalOperator(
        coef,
        [dTensor, sTensor],
        [q0,q1]
    )

def mesonSink(coef, flavor, x, gammaName):
    dTensor = deltaTensor(nextColorIdx)
    sTensor = spinTensor(gammaName, nextSpinIdx)

    q0=quark(True, flavor[0], 't_i', x)
    q1=quark(False, flavor[1], 't_i', x)

    return ElementalOperator(
        coef,
        [dTensor, sTensor],
        [q0,q1]
    )

def MesonBaryonSource(coef, baryonData, mesonData):
    eTensor = epsilonTensor(nextColorIdx)
    sTensorB = spinTensor(baryonData['gamma'], nextSpinIdx)

    q0 = quark(True, baryonData['flavor'][0], 't_i', baryonData['x'])
    q1 = quark(True, baryonData['flavor'][1], 't_i', baryonData['x'])
    q2 = quark(True, baryonData['flavor'][2], 't_i', baryonData['x'])

    dTensor = deltaTensor(nextColorIdx)
    sTensorM = spinTensor(mesonData['gamma'], nextSpinIdx)

    q3=quark(True, mesonData['flavor'][0], 't_i', mesonData['x'])
    q4=quark(False, mesonData['flavor'][1], 't_i', mesonData['x'])


    return ElementalOperator(coef,
        [eTensor, dTensor, sTensorB, sTensorM],
        [q0,q1,q2,q3,q4]
    )

def MesonBaryonSink(coef, baryonData, mesonData):
    eTensor = epsilonTensor(nextColorIdx)
    sTensorB = spinTensor(baryonData['gamma'], nextSpinIdx)

    q0 = quark(False, baryonData['flavor'][0], 't_i', baryonData['x'])
    q1 = quark(False, baryonData['flavor'][1], 't_i', baryonData['x'])
    q2 = quark(False, baryonData['flavor'][2], 't_i', baryonData['x'])

    dTensor = deltaTensor(nextColorIdx)
    sTensorM = spinTensor(mesonData['gamma'], nextSpinIdx)

    q3=quark(True, mesonData['flavor'][0], 't_i', mesonData['x'])
    q4=quark(False, mesonData['flavor'][1], 't_i', mesonData['x'])


    return ElementalOperator(coef,
        [eTensor, dTensor, sTensorB, sTensorM],
        [q0,q1,q2,q3,q4]
    )