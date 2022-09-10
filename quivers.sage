r"""
Code about quivers
Need to write a function that tests whether a quiver, represented by a matrix, is a fork or not
Use that to spit out the forkless mutation graph of a quiver.

Copyright stuff is here because I hosted it on github as well, and they insist on something like that

EXAMPLES::

<Lots and lots of examples>

AUTHORS:

- Tucker Ervin (2021-10-16): initial version

"""

# ****************************************************************************
#       Copyright (C) 2021 Tucker Ervin tjervin@crimson.ua.edu
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

def mMutation(M, w):
    r'''
    This function takes in a matrix and outputs its image under matrix mutation

    Input:
    - ''M'' -- matrix;

    - ''w'' -- sequence of mutations; The mutations range from 1 to rank(M).

    Output: ''mut'' -- mutated Matrix
    '''
    length = len(w)
    r = M.nrows()
    c = M.ncols()
    mut = copy(M) # Builds a a copy of the original matrix
    if length == 0: # In case no mutation happens
        return mut
    else:
        for i in w:
            if i < 1 or i > r or i > c: # Checks to see if we are mutating at non-vertices
                print("Invalid mutation at vertex: ", i)
                return mut
    if length == 1: # Only one mutation happens
         k = w[0]
         for i in range(r):
             for j in range(c): # Standard matrix mutation here
                 if i == k-1 or j == k-1: # Has the -1's here to account for sage indexing at 0
                     mut[i,j] = -M[i,j]
                 elif M[i,k-1]*M[k-1,j] > 0: # Has the -1's here to account for sage indexing at 0
                     mut[i,j] = M[i,j] + M[i,k-1]*M[k-1,j]*(M[i,k-1]/abs(M[i,k-1])) # Has the -1's here to account for sage indexing at 0
                 else:
                     mut[i,j] = M[i,j]
    else: # Multiple mutations happen
        mut = mMutation(mMutation(M,w[0:length-1]), [w[length-1]]) # Recursively "goes up" the mutation sequence
    return mut

def markovTest(a,b,c):
    r'''
    This function takes in 3 integers and spits out whether the cycle is mutation-cyclic or not
    
    Input:
    - ''a,b,c'' -- positive integers;
    
    Output:
    - ''cyclic'' -- boolean;
    '''
    markov = a^2 + b^2 + c^2 - a*b*c
    cyclic = true

    if min(a,b,c) <= 1:
        print("Is mutation-acyclic.")
        cyclic = false
    elif markov > 4:
        print("Is mutation-acyclic.")
        cyclic = false
    else:
        print("Is mutation-cyclic.")
        
    return cyclic

def successor(M, r):
    r'''
    This function spits out all the vertices succeeding vertex r
    
    Input:
    - ''M'' -- square skew-symmetrizable matrix;
    - ''r'' -- vertex in [1..n];
    
    Output: ''succ'' -- list of vertices succeeding vertex r
    '''
    n = M.nrows()
    succ = []
    
    for i in [0..n-1]:
        if M[r-1,i] > 0:
            succ.append(i)
            
    return succ
            
    
def predecessor(M, r):
    r'''
    This function spits out all the vertices preceding vertex r
    
    Input:
    - ''M'' -- square skew-symmetrizable matrix;
    - ''r'' -- vertex in [1..n];
    
    Output: ''predecessor'' -- list of vertices preceding vertex r
    '''
    n = M.nrows()
    pred = []
    
    for i in [0..n-1]:
        if M[i,r-1] > 0:
            pred.append(i)
            
    return pred

def stopover(j, M, k):
    r'''
    This function spits out the stopover vertices defined in Warkentein's dissertation
    
    Input:
    - ''j'' -- vertex in [1..n];
    - ''M'' -- square skew-symmetrizable matrix;
    - ''k'' -- vertex in [1..n] not equal to j;
    
    Output: ''stover'' -- list of vertices preceding vertex r
    '''
    n = M.nrows()
    succJ = successor(M,j)
    succK = successor(M,k)
    predJ = predecessor(M,j)
    predK = predecessor(M,k)
    
    partOne = [i for i in succJ if i in predK]
    partTwo = [i for i in predJ if i in succK]
    
    stover = list(set(partOne) | set(partTwo))
            
    return stover
    
def sinkVertices(M):
    r'''
    This function spits out the sink vertices of a skew-symmetrizable matrix.
    Assumes positive entry is an arrow from that row to the column
    
    Input:
    - ''M'' -- square skew-symmetrizable matrix;
    
    Output: ''sinks'' -- list of vertices preceding vertex r
    '''
    sinks = []
    n = M.nrows()
    
    for i in [0..n-1]:
        sink = true
        for j in [0..n-1]:
            if M[i,j] < 0:
                sink = false
        
        if sink:
            sinks.append(i)
            
    return sinks

def sourceVertices(M):
    r'''
    This function spits out the source vertices of a skew-symmetrizable matrix.
    Assumes positive entry is an arrow from that row to the column
    
    Input:
    - ''M'' -- square skew-symmetrizable matrix;
    
    Output: ''sources'' -- list of vertices succeeding vertex r
    '''
    sources = []
    n = M.nrows()
    
    for j in [0..n-1]:
        source = true
        for i in [0..n-1]:
            if M[i,j] < 0:
                source = false
        
        if source:
            sources.append(i)
            
    return sources

def hasSinks(M):
    r'''
    This function spits out the number of sink vertices of a skew-symmetrizable matrix.
    Assumes positive entry is an arrow from that row to the column
    
    Input:
    - ''M'' -- square skew-symmetrizable matrix;
    
    Output: number of vertices preceding vertex r
    '''
    return len(sinkVertices(M))
    
def hasSources(M):
    r'''
    This function spits out the number of source vertices of a skew-symmetrizable matrix.
    Assumes positive entry is an arrow from that row to the column
    
    Input:
    - ''M'' -- square skew-symmetrizable matrix;
    
    Output: number of vertices succeeding vertex r
    '''
    return len(sourceVertices(M))

def isAbundantAcyclic(M):
    r'''
    This function tests whether a quiver is abundant acyclic or not.
    
    Input:
    - ''M'' -- square skew-symmetrizable matrix;
    
    Output: value -- false if not abundant acyclic, true if it is
    '''
    value = false
    n = M.nrows()
    for i in [0..n-1]:
        for j in [0..n-1]:
            if i == j:
                continue
            elif M[i,j] < 2:
                value = false
                return value
    
    Q = ClusterQuiver(M)
    if Q.is_acyclic():
        value = true

    return value
    
def isFork(M):
    r'''
    This function takes in a matrix and outputs whether the corresponding quiver is a fork.
    
    Input:
    - ''M'' -- square matrix;
    
    Output: ''por'' -- vertex that is the point of return. 0 if not a fork
    '''
    value = true
    n = M.nrows()
    por = 0
    
    Q = ClusterQuiver(M)
    
    if Q.is_acyclic():
        value = false
        return por
    
    for i in [0..n-1]:
        for j in [0..n-1]:
            if i == j:
                continue
            elif abs(M[i,j]) < 2:
                value = false
                return por
    
    for r in [1..n]:
        suc = successor(M,r)
        pred = predecessor(M, r)
        
        S = M[suc, suc]
        P = M[pred, pred]
        
        Qs = ClusterQuiver(S)
        Qp = ClusterQuiver(P)
        
        if not Qs.is_acyclic() or not Qp.is_acyclic():
            value = false
            continue
        
        value = true
        
        for j in suc:
            for i in pred:
                if M[j,i] <= M[i,r-1] or M[j,i] <= M[r-1,j]:
                    value = false
                    break
            
            if not value:
                break
        
        if value:
            por = r
            break
    
    return por

def isPreFork(M):
    r'''
    This function takes in a matrix and outputs whether the corresponding quiver is a pre-fork.
    
    Input:
    - ''M'' -- square matrix;
    
    Output: ''por'' -- vertex that is the point of return. 0 if not a pre-fork
    '''
    por = 0
    n = M.nrows()
    
    for i in [0..n-1]:
        for j in [0..n-1]:
            if i == j:
                continue
            
            sList = [k for k in [0..n-1] if k != i]
            tList = [k for k in [0..n-1] if k != j]
            
            s = isFork(M[sList, sList])
            t = isFork(M[tList, tList])
            if s >= i + 1:
                s = s + 1
                
            if t >= j + 1:
                t = t + 1
                
            if s == t and s > 0 and stopover(i+1,M,j+1) == []:
                por = s
                return por

    return por

def isDiamondRoot(M):
    r'''
    This function takes in a matrix and outputs whether the corresponding quiver is a Diamond-Root,
    
    Input:
    - ''M'' -- square matrix;
    
    Output: ''kList'' -- list vertices that are r, k, and k' with point of return listed first. [] if not a Diamond-Root
    '''
    por = isPreFork(M)
    numZero = 0
    numOne = 0
    n = M.nrows()
    kList = []
    
    if por == 0:
        return kList
    
    for i in [0..n-1]:
        for j in [0..n-1]:
            if i == j:
                continue
                
            if M[i,j] == 0:
                kList = [por, i+1, j+1]
                return kList
            elif M[i,j] == 1:
                kList = []
                return kList
            
    return kList

def isDiamondWing(M):
    r'''
    This function takes in a matrix and outputs whether the corresponding quiver is a Diamond-Wing,
    
    Input:
    - ''M'' -- square matrix;
    
    Output: ''kList'' -- list vertices that are k and k' with point of return listed first. [] if not a Diamond-Wing
    '''
    por = 0
    n = M.nrows()
    kList = []
    
    for i in [0..n-1]:
        for j in [0..n-1]:
            if i == j:
                continue
            
            sList = [k for k in [0..n-1] if k != i]
            tList = [k for k in [0..n-1] if k != j]
            
            s = isFork(M[sList, sList])
            
            if s > 0 and M[i,j] == 0 and isAbundantAcyclic(M[tList, tList]) and len(stopover(i+1,M,j+1)) == n-2:
                kList = [i+1,j+1]
                return kList

    return kList

def isDiamondTip(M):
    r'''
    This function takes in a matrix and outputs whether the corresponding quiver is a Diamond-Tip,
    
    Input:
    - ''M'' -- square matrix;
    
    Output: ''kList'' -- vertices that are k and k'. [] if not a Diamond-Tip
    '''
    kList = []
    n = M.nrows()
    
    for i in [0..n-1]:
        for j in [0..n-1]:
            if i == j:
                continue
            
            sList = [k for k in [0..n-1] if k != i]
            tList = [k for k in [0..n-1] if k != j]
            
            s = isFork(M[sList, sList])
            t = isFork(M[tList, tList])
            
            if s >= i + 1:
                s = s + 1
                
            if t >= j + 1:
                t = t + 1
                
            if t == i +1 and s == j + 1 and M[i,j] == 0 and stopover(i+1,M,j+1) == []:
                kList = [i+1,j+1]
                return kList

    return kList
            
def forklessMutationGraph(M, depth):
    r'''
    Creates the forkless mutation graph of a quiver 
    
    Current  Bugs:
    If you start with a fork, it will not spit out the correct forkless graph.

    Input:
    - ''M'' -- 'n \\times 2n', initial seed;
    - ''depth'' -- the number of mutations deep we would like to go, if 0 go forever;

    Output: ''G'' -- mutation graph
    '''
    if M.is_mutable():
        N = M
        P = M
    else:
        N = copy(M)
        P = copy(M)
    P.set_immutable()
    N.set_immutable()
    V = [N]
    E = []
    mutations = []
    n = M.nrows()
    
    G = Graph([V,E], loops=true)
    
    por = isFork(M)
    
    if por > 0 and isFork(mMutation(M,[por])) == por:
        return G
    
    for i in [1..depth]:
        mutPrime = []
        
        if i == 1:
            for j in [1..n]:
                N = mMutation(M,[j])
                N.set_immutable()
                if isFork(N) == 0:
                    mutations.append([j])
                    V.append(N)
                    E.append([P,N])
            
            continue
        
        for mut in mutations:
            if len(mut) == i-1:
                mutPrime.append(mut)

        for mut in mutPrime:
            for j in [1..n]:
                if mut[i-2] != j:
                    mutNew = mut.copy()
                    mutNew.append(j)
                    N = mMutation(M,mutNew)
                    N.set_immutable()
                    
                    if isFork(N) == 0:
                        P = mMutation(M,mut)
                        P.set_immutable()
                        
                        if N not in V:
                            mutations.append(mutNew)
                            V.append(N)

                        E.append([P,N])
               
        print(i)

    G = Graph([V,E], loops=true)

    return G

def diamondLessMutationGraph(M,depth):
    r'''
    Creates the forkless mutation graph of a quiver 
    
    Current  Bugs:
    If you start with a fork, it will not spit out the correct forkless graph.

    Input:
    - ''M'' -- 'n \\times 2n', initial seed;
    - ''depth'' -- the number of mutations deep we would like to go, if 0 go forever;

    Output: ''G'' -- mutation graph
    '''
    if M.is_mutable():
        N = M
        P = M
    else:
        N = copy(M)
        P = copy(M)
    P.set_immutable()
    N.set_immutable()
    V = [N]
    E = []
    mutations = []
    n = M.nrows()
    
    G = Graph([V,E], loops=true)
    
    por = isFork(M)
    porList = isDiamondRoot(M)
    
    if por > 0 and isFork(mMutation(M,[por])) == por:
        return G
    elif len(porList) > 0 and isDiamondRoot(mMutation(M, [porList[0]])) == porList[0]:
        return G
    
    for i in [1..depth]:
        mutPrime = []
        
        if i == 1:
            for j in [1..n]:
                N = mMutation(M,[j])
                N.set_immutable()
                if isFork(N) == 0 and len(isDiamondRoot(N)) == 0:
                    mutations.append([j])
                    V.append(N)
                    E.append([P,N])
            
            continue
        
        for mut in mutations:
            if len(mut) == i-1:
                mutPrime.append(mut)

        for mut in mutPrime:
            for j in [1..n]:
                if mut[i-2] != j:
                    mutNew = mut.copy()
                    mutNew.append(j)
                    N = mMutation(M,mutNew)
                    N.set_immutable()
                    
                    if isFork(N) == 0 and len(isDiamondRoot(N)) == 0:
                        P = mMutation(M,mut)
                        P.set_immutable()
                        
                        if N not in V:
                            mutations.append(mutNew)
                            V.append(N)

                        E.append([P,N])
               
        print(i)

    G = Graph([V,E], loops=true)

    return G