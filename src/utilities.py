import copy

def simplify_laphDiagrams(laphDiagrams):
    for d in laphDiagrams:
        d.create_baryon_blocks()
        d.combine_indices()

    # the below has no guarantee to be NC agnostic nor work for any multi meson&baryon operator
    for d in laphDiagrams:
        d.create_baryon_source()


def find_all_baryon_tensors(laphDiagrams):
    allBaryonTensors = {}
    idx=0
    for d in laphDiagrams:
        for tensor in d.commuting:
            if tensor.id() not in allBaryonTensors:
                allBaryonTensors[tensor.id()]=idx
                idx+=1
    return allBaryonTensors


def diagram_as_graph(d, allBaryonTensors):
    tst=copy.deepcopy(d)

    graph = {}

    for b in reversed(tst.commuting):
        tst.commuting.remove(b)
        bIdx=allBaryonTensors[b.id()]
    
        for idx in reversed(b.indices):
            for e in tst.commuting:
                if idx in e.indices:
                    fIdx=allBaryonTensors[e.id()]
                    contraction=str(bIdx)+","+str(fIdx)
                    if contraction in graph:
                        graph[contraction].append([b.indices.index(idx),e.indices.index(idx)])
                    else:
                        graph[contraction]=[[b.indices.index(idx),e.indices.index(idx)]]
                    
                    break

    return graph