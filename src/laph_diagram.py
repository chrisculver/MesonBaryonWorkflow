import WickContractions.corrs.diagram
import copy
from WickContractions.ops.indexed import IndexedFunction

# TODO: Users should probably make this themselves?

class LDiagram(WickContractions.corrs.diagram.Diagram):
    def __init__(self, diag):
        if(type(diag) is not WickContractions.corrs.diagram.Diagram):
            raise ValueError("Must make Laph diagram out of src.diags.Diagram")
        self.coef = diag.coef
        self.commuting = diag.commuting
        self.props = diag.props
        self.Tblocks = {}
        self.Bblocks= {}
        
        lCnt=0
        for p in self.props:
            p.name='\\tau'
            for i,indices in enumerate([p.left_indices,p.right_indices]):
                # swap color for eigenvector idx
                color_idx = indices.c
                if color_idx[0:2] == 'c_':
                    idx='l_'+str(color_idx.split('_')[1])
                else:
                    idx='l_{'+str(lCnt)+'}'
                lCnt+=1
                indices.c=idx
                args=[]
                if i==0:
                    args.append(p.xi)
                    args.append(p.ti)
                else:
                    args.append(p.xf)
                    args.append(p.tf)
                if i==0: 
                    self.commuting.append(IndexedFunction('V^*',[color_idx, idx],args))
                else:
                    self.commuting.append(IndexedFunction('V',[color_idx, idx], args))

    def get_first_idx_ends_of(self,name):
        """ generates a list of the end indices on a named object
        
        :param: name: A string
        
        :return: list of strings that are the end of the index or an empty list
                 if no object of that name is in the list.
        
        """
        for elem in self.commuting:
            if(elem.name==name):
                idx_ends=[]
                for idx in elem.indices:
                    idx_ends.append(idx[2:])
                return idx_ends
        return []
    
    
    def create_T_blocks(self):
        color_obj_name = '\\epsilon' 
        epsilons=[]
        for elem in reversed(self.commuting):
            if elem.name==color_obj_name:
                epsilons.append(elem)
                self.commuting.remove(elem)
        
        for e in epsilons:
            # start a T function with no indices or arguments
            Tblock = IndexedFunction('T',[],[])
            TblockEquals = [e]
            for idx in e.indices:
                for elem in self.commuting:
                    if elem.name not in ['T','T^*']:
                        if idx in elem.indices:
                            TblockEquals.append(elem)
                            new_indices = list(elem.indices)
                            new_indices.remove(idx)
                            for new_idx in new_indices:
                                if new_idx not in Tblock.indices:
                                    Tblock.indices.append(new_idx)
                                else:
                                    print("Warning this index is already in TBlock...\n DIdn't expect this")
                            if type(elem) == IndexedFunction:
                                for new_arg in elem.arguments:
                                    if new_arg not in Tblock.arguments:
                                        Tblock.arguments.append(new_arg)
                                    # here I expect repeats of arguments
                            
                            self.commuting.remove(elem)
                            
            rhsString=''
            otherNames=[]
            for elem in TblockEquals:
                rhsString+=str(elem)
                if elem.name != '\\epsilon' and elem.name not in otherNames:
                    otherNames.append(elem.name)
            
            if len(otherNames)>1:
                print("T has different kinds of V's in it...")
                Tblock.name+='^?'
            elif otherNames[0]=='V^*':
                Tblock.name+='^*'
            
            self.Tblocks[str(Tblock)]=rhsString
            self.commuting.append(Tblock)
               