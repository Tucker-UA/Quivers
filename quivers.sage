r"""
Code about quivers

EXAMPLES::

<Lots and lots of examples>

AUTHORS:

- Tucker Ervin (2023-03-05): initial version

"""

# ****************************************************************************
#       Copyright (C) 2023 Tucker Ervin tjervin@crimson.ua.edu
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

# Imports
from functools import cached_property

# Code

class Quiver:
    r'''
    The quiver associated to an exchange matrix
    Could be skew-symmetrizable

    The functions usually only work on the mutable subquiver, so keep that in mind

    Note that the properties of the quiver as a rule are immutable.
    If you want to change them, you have to produce another quiver.

    Input:
    - ''M'' -- matrix
    - ''rowsOrCols'' -- string that is cols if a framed matrix would be n x 2n and rows if 2n x n: default = cols
    - ''orientation'' -- string that is positive if i -> j implies a positive entry and negative if not: default = positive


    '''
    
    # Initialization
    def __init__(self, M = matrix(([0])), rowsOrCols = "cols", orientation = "positive"):
        self.numRows = M.nrows()
        self.numCols = M.ncols()
        self.rowsOrCols = rowsOrCols
        self.orientation = orientation
        self.matrix = copy(M)
        self.matrix.set_immutable()
        self.nMut = self.numRows
        self.nFroz = 0
        self._hasFrozen = False
        
        if self.numRows != self.numCols: # Checks to see if the num of rows or cols is off
            if self.numRows * 2 == self.numCols and rowsOrCols == "cols": # We're defining cols here to be the longer one
                self.nMut = self.numRows
                self.nFroz = self.numRows
                self._hasFrozen = True
            elif self.numCols * 2 == self.numRows and rowsOrCols == "rows": # We're defining rows here to be the longer one
                self.nMut = self.numCols
                self.nFroz = self.numCols
                self._hasFrozen = True
            else:
                self.matrix = zero_matrix(self.numRows) # If something is not right, make it the zero quiver on that many vertices
                print("Something went wrong and we got a zero matrix on ", self.nMut, " vertices")
            
        self.n = self.nMut + self.nFroz # Total number of vertices
        assert M[:self.nMut,:self.nMut].is_skew_symmetrizable()
         
    #  Here's the representation functions
    def __repr__(self):
        """Returns a string of an expression that re-creates this object."""
        return f"{self.__class__.__qualname__}({self.matrix})"

    def __str__(self):
        """Prints the matrix itself in a string form"""
        return str(self.matrix)
    
    # Comparison functions
    
    def __eq__(self, other):
        """Replaces the == comparison"""
        '''The equality it tests for is equality of all the values
        May want to eventually replace this with isomorphism or may not'''
        
        equality = False # The default return value is false
        
        if isinstance(other, Quiver):
            if other.matrix == self.matrix and other.rowsOrCols == self.rowsOrCols and other.orientation == self.orientation:
                equality = True
        elif isinstance(other, Matrix):
            P = Quiver(other)
            equality = (self == P)
        
        return equality

    def __ne__(self, other):
        """Replaces the != comparison"""
        
        inequality = not (self == other)
            
        return inequality
   
    # Hash function, because we have to define it if we define the ==
    
    def __hash__(self):
        '''Returns the hash of the matrix as the hash value'''
        
        return hash(self.matrix)

    # The other comparison functions don't make a ton of sense, as we don't usually compare quivers like this
    # Would have to build a partial ordering for this to even work out

    def __lt__(self, other):
        return NotImplemented
            
    def __gt__(self, other):
        return NotImplemented
    
    def __le__(self, other):
        return NotImplemented
    
    def __ge__(self, other):
        return NotImplemented
    
    # Other dunder functions
    
    def __add__(self,other):
        r'''
        Implements the disjoint union of two quivers
        '''
        return self.disjointUnion(other)
    
    # Will build an isomorphism checker function
    
    def isIsomorphicTo(self, P):
        r'''
        Function that checks to see if a given quiver P is isomorphic to our quiver or not
        Returns true if it is, false if it isn't
        '''
        if P.n != self.n: # Vertices don't match
            return False
        elif P.size != self.size: # If the number of arrows does not match up
            return False
        
        
        n = self.nMut
        vertices = self.mutableVertices
        M = self.matrix
        maxNeighbors = 0
        i = 0
        perms = []
        
        for j in vertices:
            neighbor = self.neighborhood(i)
            if len(neighbor) > maxNeighbors:
                maxNeighbors = len(neighbor)
                i = j

        return False
    
    # Building other quivers out of this quiver
    
    def preserveOptions(self, M):
        r'''
        Produces a quiver with the same options as the original quiver, but a different matrix M
        '''
        P = Quiver(M, rowsOrCols = self.rowsOrCols, orientation = self.orientation)
        
        return P
    
    def principalQuiver(self):
        r'''
        Gets the subquiver of mutable vertices
        
        Output:
        - ''P'' -- The desired quiver
        '''
        n = self.nMut # This works because the principal submatrix is always a square matrix
        M = self.matrix[:n,:n] # Easier than converting the tuple back to indices
        P = self.preserveOptions(M)
        
        return P
    
    def framedQuiver(self):
        r'''
        Produces the framed quiver with the current quiver acting as the mutable subquiver
        Adds principal coefficients to matrix if we're thinking of it that way
        
        Output:
        - ''P'' -- The desired quiver
        '''
        n = self.nMut
        B = self.matrix[:n,:n] # Gets the nxn principal submatrix 
        C = identity_matrix(n) # Principal Coefficients
        
        if self.rowsOrCols == "cols": # More Cols
            M = block_matrix([[B,C]])
        elif self.rowsOrCols == "rows": # More Rows
            M = block_matrix([[B],[C]])
      
        P = self.preserveOptions(M)
        
        return P
    
    def oppositeQuiver(self):
        r'''
        Produces the quiver with every arrow flipped
        
        Output:
        - ''P'' -- The desired quiver
        '''
        M = -1 * self.matrix
        Q = self.preserveOptions(M)
        
        return Q
    
    def fullSubquiver(self, vertices = [1]):
        r'''
        Produces the full subquiver induced by the given vertices
        Does not do every full subquiver, as we need a n x n or n x 2n or 2n x n matrix
        
        Input:
        - ''vertices'' -- list of vertices ranging from 1 to n to n+n, where i' = n + i

        Output:
        - ''P'' -- The desired quiver
        '''
        for v in vertices:
            assert v in self.vertices
            
        matrixVertices = [i-1 for i in vertices if i <= self.nMut or i-self.nMut in vertices]
        principalVertices = [i -1 for i in vertices if i <= self.nMut]
        
        if matrixVertices == principalVertices:
            M = copy(self.matrix[matrixVertices,matrixVertices])
        elif self.rowsOrCols == "cols": # More Cols
            M = copy(self.matrix[principalVertices,matrixVertices])
        elif self.rowsOrCols == "rows": # More Rows
            M = copy(self.matrix[matrixVertices, principalVertices])
        else:
            print("Something bad happened")
            
        P = self.preserveOptions(M)
        
        return P
    
    def subquiverRemoveVertex(self, vertex = 1):
        r'''
        Produces the full subquiver P - {vertex}
        
        Input:
        - ''vertex'' -- vertex to remove ranging from 1 to n to n+n, where i' = n + i
        
        Output:
        - ''P'' -- the desired quiver
        '''
        assert vertex in self.vertices
        vertices = list(self.vertices)
        vertices.remove(vertex)
        P = self.fullSubquiver(vertices)
        
        return P
    
    def mutate(self, mutSeq = []):
        r'''
        This function takes in a quiver and outputs its image under mutation
        It does not affect the original quiver

        Input:
        - ''mutSeq'' -- sequence of mutations; The mutations range from 1 to self.n; default is the empty mutation

        Output: ''P'' -- mutated Quiver
        '''
        length = len(mutSeq) # Gets the length of the mutation sequence
        r = self.numRows # Number of rows of matrix
        c = self.numCols # Number of cols of matrix
        mut = copy(self.matrix) # Builds a a copy of the original matrix
        P = self.preserveOptions(mut)
        
        if length == 0: # In case no mutation happens
            return P # If no mutation occurs, then we are just left with the original matrix
        else:
            for i in mutSeq: # Sequentially goes through mutation sequence
                if not (i >= 1 and i <= self.nMut): # Checks to see if we are mutating at non-vertices
                    print("Invalid mutation at vertex: ", i)
                    return P # No mutation happened, so we return the original matrix
        
        if length == 1: # Only one mutation happens
            k = mutSeq[0] - 1 # This is the vertex at which we mutate subtract one as sage indexs at 0
             
            M = copy(self.matrix)
            for i in range(r): # The (i,j) loop over the matrix elements
                for j in range(c): # Standard matrix mutation here
                    if i == k or j == k: 
                        mut[i,j] = -M[i,j] # Same as reversing arrows touching the vertex k
                    elif M[i,k]*M[k,j] > 0:
                        mut[i,j] = M[i,j] + M[i,k]*M[k,j]*(M[i,k]/abs(M[i,k])) # Happens if there is a 2-path passing through k
                    else:
                        mut[i,j] = M[i,j] #If neither of the above occured, we do not modify the number of arrows between i and j 

        else: # Multiple mutations happen
            R = self.mutate(mutSeq[0:length-1])
            S = R.mutate([mutSeq[length-1]])
            mut = S.matrix # Recursively "goes up" the mutation sequence

        P = self.preserveOptions(mut)
        
        return P
    
    def isomorphicQuiver(self, p):
        r'''
        This function returns a quiver that is isomorphic to the given quiver,
        where the isomorphism is produced by the permutation p on the vertices
        
        Input:
        - ''p'' -- list or permutation on n entries
        
        Output: ''P'' -- the isomorphic quiver.
        '''
        p = Permutation(p)
        pM = p.to_matrix() # This makes P a row permutation matrix
        n = self.nMut
        bigPM = block_matrix([[pM, zero_matrix(n)], [zero_matrix(n), identity_matrix(n)]])
        
        if not self.hasFrozen:
            M = pM * self.matrix * pM # Multiplying on left permutes rows, on right permutes columns
        elif self.rowsOrCols == "cols": # More cols than rows
            M = pM * self.matrix * bigPM
        elif self.rowsOrCols == "rows":
            M = bigPM * self.matrix * pM
        else:
            M = self.matrix

        P = self.preserveOptions(M)
        
        return P
    
    def disjointUnion(self, P):
        r'''
        Returns the quiver formed from the disjoint union of the principal subquivers of P and Q.
        If a problem occurs, return just P
        '''
        T = copy(P)
        if P.orientation != self.orientation or P.rowsOrCols != self.rowsOrCols:
            print("Quivers have different orientations or rows or columns")
            return T
        
        n = self.n
        m = P.n
        
        if not P.hasFrozen and not self.hasFrozen: # Want to deal in the case without frozen vertices first
            M = block_matrix([[self.matrix, 0],[0, P.matrix]])
            T = self.preserveOptions(M)
        elif self.hasFrozen and not P.hasFrozen: # Produces a quiver with self.nMut + P.nMut mutable vertices and that many again frozen
            R = self.principalQuiver.disjointUnion(P.principalQuiver)
            n = self.nMut
            if self.rowsOrCols == "cols":
                N = block_matrix([[self.matrix[:,n:],zero_matrix(m)],[zero_matrix(n),zero_matrix(m)]])
                M = block_matrix([[R.matrix,N]])
            elif self.rowsOrCols == "rows":
                N = block_matrix([[self.matrix[n:,:],zero_matrix(m)],[zero_matrix(n),zero_matrix(m)]])
                M = block_matrix([[R.matrix],[N]])
            T = self.preserveOptions(M)
        elif P.hasFrozen and not self.hasFrozen:
            R = self.principalQuiver.disjointUnion(P.principalQuiver)
            m = P.nMut
            if self.rowsOrCols == "cols":
                N = block_matrix([[zero_matrix(n),zero_matrix(m)],[zero_matrix(n),P.matrix[:,m:]]])
                M = block_matrix([[R.matrix,N]])
            elif self.rowsOrCols == "rows":
                N = block_matrix([[zero_matrix(n),zero_matrix(m)],[zero_matrix(n),P.matrix[m:,:]]])
                M = block_matrix([[R.matrix],[N]])
            T = self.preserveOptions(M)
        elif P.hasFrozen and self.hasFrozen:
            R = self.principalQuiver.disjointUnion(P.principalQuiver)
            n = self.nMut
            m = P.nMut
            if self.rowsOrCols == "cols":
                N = block_matrix([[self.matrix[:,n:],zero_matrix(m)],[zero_matrix(n),P.matrix[:,m:]]])
                M = block_matrix([[R.matrix,N]])
            elif self.rowsOrCols == "rows":
                N = block_matrix([[self.matrix[n:,:],zero_matrix(m)],[zero_matrix(n),P.matrix[m:,:]]])
                M = block_matrix([[R.matrix],[N]])
            T = self.preserveOptions(M)
        else:
            T = self.principalQuiver.disjointUnion(P.principalQuiver)
            
        return T
    
    def triangularExtension(self, P, E = {}):
        r'''
        Returns the quiver formed from the triangular extension of principal subquivers of P and Q,
        where E is the set of edges.
        
        Input:
        - ''P'' -- a quiver 
        - ''E'' -- the set of edges to add, where E = {[i,j] : multiplicity} for i in Q.vertices and j in P.vertices
        
        Output:
        - ''T'' -- the desired quiver
        
        Default output will be the disjoint union of P and Q if P and Q are compatible
        
        '''
        T = copy(P)
        if P.orientation != self.orientation or P.rowsOrCols != self.rowsOrCols:
            print("Quivers have different orientations or rows or columns")
            return T
        
        R = self.disjointUnion(P)
        M = copy(R.matrix)
        n = self.nMut
        
        if len(E) == 0:
            return R
        
        signE = 0
        
        for key in E:
            value = E[key]
            if signE == 0 and value != 0:
                signE = sign(value)
            elif value == 0:
                continue
            
            if signE != sign(value):
                print("E does not produce a triangular extension", signE, value)
                return R
            elif key[0] not in self.mutableVertices or key[1] not in P.mutableVertices:
                print("E does not have the right vertices")
                return R
            else:
                i = key[0] - 1
                j = key[1] - 1 + n
                
                M[i,j] = value
                M[j,i] = -1 * value
                
        T = self.preserveOptions(M)
        
        return T
    
    def abundantTriangularExtension(self, P, E = {}):
        r'''
        Returns the triangular extension formed from P and Q if every arrow is >= 2.
        If an edge is not present, adds a default value of \pm 2

        Otherwise works just as the triangularExtension
        '''
        T = copy(P)
        if P.orientation != self.orientation or P.rowsOrCols != self.rowsOrCols:
            print("Quivers have different orientations or rows or columns")
            return T
        
        signE = 0
        
        for key in E:
            value = E[key]
            if value != 0:
                signE = sign(value)
                break
                
        if signE == 0:
            signE = 1
            
        for i in self.mutableVertices:
            for j in P.mutableVertices:
                if (i,j) not in E:
                    E[(i,j)] = signE * 2
                    
        T = self.triangularExtension(P,E)
        
        return T
        
    # Functions that make lists of vertices
    
    @cached_property
    def mutableVertices(self):
        r'''
        Lists the mutable vertices
        '''
        n = self.nMut
        
        return tuple([1..n])
    
    @cached_property
    def vertices(self):
        r'''
        List all vertices
        '''
        n = self.n
        
        return tuple([1..n])
    
    @cached_property
    def frozenVertices(self):
        r'''
        Function that lists all frozen vertices
        '''
        vert = self.vertices
        
        frozen = [v for v in vert if v > self.nMut]
        
        return tuple(frozen)

    def successor(self, r):
        r'''
        This function spits out all the vertices succeeding vertex r, possibly frozen vertices

        Input:
        - ''r'' -- vertex in [1..n,n+1,...,n+n];

        Output: ''succ'' -- list of vertices succeeding vertex r
        '''
        vertices = self.vertices
        assert r in vertices
        succ = [] # Creates an empty list to append to
        
        for i in vertices:
            if self.arrow(r,i) > 0:
                succ.append(i)
                    
        return tuple(succ)

    def predecessor(self, r):
        r'''
        This function spits out all the vertices preceeding vertex r, possibly frozen vertices

        Input:
        - ''r'' -- vertex in [1..n,n+1,...,n+n];

        Output: ''pred'' -- list of vertices preceeding vertex r
        '''
        vertices = self.vertices
        assert r in vertices
        pred = [] # Creates an empty list to append to
        
        for i in vertices:
            if self.arrow(i,r) > 0:
                pred.append(i)
                    
        return tuple(pred)

    def neighborhood(self, v):
        r'''
        This function returns a tuple of all vertices that precede or succede v
        
        Input:
        - ''v'' -- vertex of the quiver
        
        Output: ''neighbor'' -- tuple of vertices touching v
        '''
        return tuple(set(self.successor(v)) | set(self.predecessor(v))) # Union of the both as sets then makes a tuple
    
    def stopover(self,k,l):
        r'''
        This function spits out all the vertices that are pointed to by k and point to l or vice versa
        
        Input:
        - ''k'' -- The starting vertex
        - ''l'' -- The ending vertex
        
        Output:
        - ''vertices'' -- list of stopover vertices
        '''
        succK = set(self.successor(k)) # Gets the successor vertices for k
        succL = set(self.successor(l)) # Gets the successor vertices for l
        predK = set(self.predecessor(k)) # Gets the predecessor vertices for k
        predL = set(self.predecessor(l)) # Gets the predecessor vertices for l
        
        vertices = tuple((succK & predL) | (predK & succL)) # This is the same formula in the documentation
        
        return vertices
   
    @cached_property
    def sinks(self):
        r'''
        Returns a list of all the sink vertices
        Does not return frozen vertices that are sinks
        
        Output:
        - ''sinks'' -- List of all the sinks
        '''
        sinks = []
        vertices = self.mutableVertices

        for i in vertices:
            sink = True
            for j in vertices:
                if i == j:
                    continue
                elif self.arrow(i,j) > 0:
                    sink = False
                    break
            
            if sink:
                sinks.append(i)
        
        return tuple(sinks)
    
    @cached_property
    def sources(self):
        r'''
        This function spits out the source vertices
        Does not return frozen vertices that are sources

        Output: ''sources'' -- list of vertices succeeding vertex r
        '''
        sources = [] # Creates an empty list to append to
        vertices = self.mutableVertices
        
        for i in vertices:
            source = True
            for j in vertices:
                if i == j:
                    continue
                elif self.arrow(j,i) > 0:
                    source = False
                    break
            
            if source:
                sources.append(i)

        return tuple(sources)
    
    @cached_property
    def greens(self):
        r'''
        This function spits out all green vertices
        Only works if the quiver has frozen vertices
        
        Output: ''green'' -- list of vertices that are green
        '''
        green = []
        if not self.hasFrozen:
            return tuple(green)
        
        froz = self.frozenVertices
        mut = self.mutableVertices
        # Next line collects all green vertices v. Throws out ones that would violate sign-coherence (should be impossible)
        green = (v for v in mut if ((set(self.successor(v)) & set(froz) != set()) and (set(self.predecessor(v)) & set(froz) == set())))
        
        return green
    
    @cached_property
    def reds(self):
        r'''
        This function spits out all red vertices
        Only works if the quiver has frozen vertices
        
        Output: ''red'' -- list of vertices that are red
        '''
        red = []
        if not self.hasFrozen:
            return tuple(red)
        
        froz = self.frozenVertices
        mut = self.mutableVertices
        # Next line collects all red vertices v. Throws out ones that would violate sign-coherence (should be impossible)
        red = (v for v in mut if ((set(self.successor(v)) & set(froz) == set()) and (set(self.predecessor(v)) & set(froz) != set())))
        
        return red
    
    # Functions that return information about the quiver
    
    def arrow(self, source, target):
        r'''
        Return the number of arrows pointing from the source vertex to the target vertex

        Input:
        - ''source'' -- a vertex
        - ''target'' -- a vertex

        Output:
        - ''numArrows'' -- number of arrows desired: negative if the arrows point the other way
        '''
        vertices = self.vertices
        mutableVertices = self.mutableVertices
        #print(source,target,vertices)
        assert source in vertices
        assert target in vertices

        i = source - 1
        j = target - 1 # Taking care of indexing at 0
        
        numArrows = 0
        
        if source in mutableVertices and target in mutableVertices:
            numArrows = self.matrix[i,j]
        else:
            if self.rowsOrCols == "cols" and source in mutableVertices:
                numArrows = self.matrix[i,j]
            elif self.rowsOrCols == "cols" and target in mutableVertices:
                numArrows = -1*self.matrix[j,i] # Since we can't access this number of arrows, access the one we can and negate
            elif self.rowsOrCols == "cols":
                numArrows = 0 # No arrows are present between frozen vertices

            if self.rowsOrCols == "rows" and source in mutableVertices:
                numArrows = -1*self.matrix[j,i]
            elif self.rowsOrCols == "rows" and target in mutableVertices:
                numArrows = self.matrix[i,j] # Since we can't access this number of arrows, access the one we can and negate
            elif self.rowsOrCols == "rows":
                numArrows = 0 # No arrows are present between frozen vertices

        if self.orientation == "negative": # We assumed we were working with positive orientation above
            numArrows *= -1
        
        #print(source,target,numArrows) # For testing
        return numArrows
    
    @cached_property
    def sizeMutable(self):
        r'''
        Returns the number of arrows connected to mutable vertices in the quiver
        '''
        vertices = self.mutableVertices
        numArrows = 0
        
        for i in vertices:
            for j in vertices:
                a = arrow(i,j)
                if a > 0:
                    numArrows += a
        
        return numArrows
  
    @cached_property
    def size(self):
        r'''
        Returns the number of arrows in the quiver, counting arrows that are connected to frozen vertices
        '''
        if not self.hasFrozen:
            return self.sizeMutable
        
        numArrows = self.sizeMutable
        mutableVertices = self.mutableVertices
        frozenVertices = self.frozenVertices

        for i in mutableVertices:
            for j in frozenVertices:
                numArrows += abs(arrow(i,j))
        
        return numArrows

    @cached_property
    def _pointOfReturnFork(self):
        r'''
        Returns the point of return for a fork.
        If it doesn't have one, return 0.
        Doesn't worry about checking other vertices for a por,
        as any por is unique.
        '''
        if (not self.isAbundant) or self.isAcyclic:
            return 0
        
        vertices = self.mutableVertices
        R = self.principalQuiver()

        for r in vertices:
            Q = self.subquiverRemoveVertex(r)
            if not Q.isAcyclic:
                continue
                
            P = R.predecessor(r)
            S = R.successor(r)
            
            foundPor = True
            
            for i in P:
                for j in S:
                    if self.arrow(j,i) <= self.arrow(i,r) or self.arrow(j,i) <= self.arrow(r,j):
                        foundPor = False
                        break
                   
                if not foundPor:
                    break
                    
            if foundPor:
                return r
            
        return 0
    
    @cached_property
    def _preForkVertices(self):
        r'''
        Returns the three vertices that make a quiver a pre-fork
        In the format: [r, k, k'], where r is the point of return
        Returns [r,0,0] if the quiver is a fork
        Returns [0,0,0] if the quiver is not a pre-fork
        '''
        if self.isAcyclic:
            return [0,0,0]
        elif self.isFork:
            return [self.pointOfReturn,0,0]
        
        vertices = self.mutableVertices
        
        for i in vertices:
            for j in vertices:
                if i >= j or self.arrow(i,j) not in [-1,0,1]:
                    continue
                    
                Q = self.subquiverRemoveVertex(i).principalQuiver()
                P = self.subquiverRemoveVertex(j).principalQuiver()
                stover = self.stopover(i,j)
                
                if not Q.isFork or not P.isFork or stover != []:
                    return [0,0,0]
                elif Q.pointOfReturn != P.pointOfReturn:
                    return [0,0,0]
                else:
                    return [Q.pointOfReturn, i, j]
                
        return [0,0,0]
                    
    @cached_property
    def _keyVertices(self):
        r'''
        Returns two vertices that make the quiver a key.
        Returns [0,0] if the quiver is not a key
        Returns [1,2] if the quiver is abundant acyclic
        '''
        if not self.isAcyclic:
            return [0,0]
        
        vertices = self.mutableVertices
        
        for i in vertices:
            for j in vertices:
                if i >= j or self.arrow(i,j) not in [-1,0,1]:
                    continue
                
                P = self.subquiverRemoveVertex(i).principalQuiver()
                Q = self.subquiverRemoveVertex(j).principalQuiver()
                
                if P.isAbundantAcyclic and Q.isAbundantAcyclic:
                    return [i,j]
                else:
                    return [0,0]
                
        return [1,2]
    
    @cached_property
    def _wingVertices(self):
        r'''
        Returns the two vertices that make the quiver a wing.
        Returns [0,0] if the quiver is not a wing
        Returns [k,v], where k is the point of return of the wing otherwise
        '''
        if self.isAcyclic or self.isFork or self.isPreFork:
            return [0,0]
        
        vertices = self.mutableVertices
        
        for i in vertices:
            for j in vertices:
                if i >= j or self.arrow(i,j) not in [-1,0,1]:
                    continue
                    
                P = self.subquiverRemoveVertex(i).principalQuiver()
                Q = self.subquiverRemoveVertex(j).principalQuiver()
                stover = self.stopover(i,j)
                remainingVertices = list(vertices)
                remainingVertices.remove([i])
                remainingVertices.remove([j])
                
                if stover != remainingVertices:
                    return [0,0]
                
                k = 0
                v = 0
                
                if P.isFork and Q.isAbundantAcyclic and j == P.pointOfReturn:
                    k = j
                    v = i
                elif Q.isFork and P.isAbundantAcyclic and i == Q.pointOfReturn:
                    k = i
                    v = j
                else:
                    return [0,0]
                
                if self.arrow(k,v) == 0:
                    return [k,v]
                
                for L in vertices:
                    if self.arrow(v,k) > 0 and self.arrow(k,L) > 0:
                        if self.arrow(L,v) < self.arrow(k,L) + 2:
                            return [0,0]
                    elif self.arrow(k,v) > 0 and self.arrow(L,k) > 0:
                        if self.arrow(v,L) < self.arrow(L,k) + 2:
                            return [0,0]
                return [k,v]
            
        return [0,0]
                
    @cached_property
    def _tipVertices(self):
        r'''
        Returns the vertices of a tip
        Returns [0,0] if the quiver is not a tip
        Returns [k,v] if the quiver is a tip with point of return k
        '''
        if self.isAcyclic or self.isFork or self.isPreFork or self.isWing:
            return [0,0]
        
        vertices = self.mutableVertices
        
        for i in vertices:
            for j in vertices:
                if i >= j or self.arrow(i,j) not in [-1,0,1]:
                    continue
                    
                P = self.subquiverRemoveVertex(i).principalQuiver()
                Q = self.subquiverRemoveVertex(j).principalQuiver()
                stover = self.stopover(i,j)
                
                if self.arrow(i,j) == 0 and stover != []:
                    return [0,0]
                elif self.arrow(i,j) == 0:
                    if P.isFork and Q.isFork:
                        if P.pointOfReturn != j or Q.pointOfReturn != i:
                            return [0,0]
                        else:
                            return [i,j]
                    else:
                        return [0,0]
                
                k = 0
                v = 0
                
                if P.isFork and Q.isAbundantAcyclic and j == P.pointOfReturn:
                    k = j
                    v = i
                elif Q.isFork and P.isAbundantAcyclic and i == Q.pointOfReturn:
                    k = i
                    v = j
                else:
                    return [0,0]
                
                if self.arrow(k,v) == 1 and stover != self.predecessor(k):
                    return [0,0]
                elif self.arrow(v,k) == 1 and stover != self.successor(k):
                    return [0,0]
                
                # Need to fix this remaining bit
                
                for L in vertices:
                    if self.arrow(v,k) > 0 and self.arrow(k,L) > 0:
                        if self.arrow(k,L) < self.arrow(L,v) + 2:
                            return [0,0]
                    elif self.arrow(k,v) > 0 and self.arrow(L,k) > 0:
                        if self.arrow(L,k) < self.arrow(v,L) + 2:
                            return [0,0]
                return [k,v]
            
        return [0,0]
    
    @cached_property
    def pointOfReturn(self):
        r'''
        Returns the point of return for the quiver if it has one
        0 otherwise
        '''
        if self.isFork:
            return self._pointOfReturnFork
        elif self.isPreFork:
            return self._preForkVertices[0]
        elif self.isWing:
            return self._wingVertices[0]
        elif self.isTip:
            return self._tipVertices[0]
            
        return 0
    
    # Truth values about the quiver
    # Saved as properties
    
    @cached_property 
    def hasSinks(self):
        r'''
        This property is true if the quiver has mutable sinks
        '''
        return (len(self.sinks) > 0)
    
    @cached_property
    def hasSources(self):
        r'''
        This property is true if the quiver has mutable sinks
        '''
        return (len(self.sources) > 0)
    
    @cached_property
    def hasFrozen(self):
        r'''
        This property is true if the quiver has frozen vertices
        '''
        return self._hasFrozen
    
    @cached_property
    def hasZeroForklessPart(self):
        r'''
        Determines if a quiver has 0 forkless part.
        Property returns true if it does have a 0 forkless part
        '''
        if not self.isFork:
            return False
        
        n = self.nMut # Get number of mutable vertices
        value = False # Default value
        por = self.pointOfReturn # Gets the por

        Q = self.mutate([por]) # Mutate at the point of return

        if not Q.isFork:
            return False

        if Q.pointOfReturn == por: # If the points of return are the same, then we have zero forkless part
            return True

        return Q.hasZeroForklessPart
    
    def hasFiniteForklessPart(self):
        r'''
        Returns true if it does have a finite forkless part
        False otherwise
        Unfinished for now
        '''
        if self.hasZeroForklessPart:
            return self.hasZeroForklessPart
        elif self._hasFiniteForklessPart is not None:
            return self._hasFiniteForklessPart
        
        return False
    
    @cached_property
    def isComplete(self):
        r'''
        This property is true if there is an arrow between every pair of mutable vertices
        '''
        n = self.nMut
        vertices = self.mutableVertices

        for i in vertices:
            for j in vertices:
                if i == j:
                    continue
                elif self.arrow(i,j) == 0:
                    return False
                
        return True
    
    @cached_property
    def isAbundant(self):
        r'''
        This property is true if there is 2 or more arrows between every pair of mutable vertices
        '''
        vertices = self.mutableVertices

        for i in vertices:
            for j in vertices:
                if i == j:
                    continue
                elif abs(self.arrow(i,j)) < 2:
                    return False
                
        return True

    @cached_property
    def isAcyclic(self):
        r'''
        This property is true if the principal quiver is acyclic
        '''
        n = self.nMut
        
        if n <= 2:
            return True
        
        if (not self.hasSinks) or (not self.hasSources):
            return False
        
        s = self.sinks[0]
        vertices = list(self.mutableVertices)
        vertices.remove(s)
        Q = self.fullSubquiver(vertices)
        
        return Q.isAcyclic
   
    @cached_property
    def isAbundantAcyclic(self):
        r'''
        This property is true if the principal quiver is abundant acyclic
        '''
        return (self.isAbundant and self.isAcyclic)
    
    @cached_property
    def isKey(self):
        r'''
        This function tests to see if the principal quiver is a key
        '''
        if not self.isAcyclic:
            return False
        
        vertices = self.mutableVertices
        
        for i in vertices:
            for j in vertices:
                if self.arrow(i,j) not in [-1,0,1] or i == j: # If there are 2 or more arrows between i and j
                    continue
                
                stover = self.stopover(i,j)
                if stover != []:
                    return False
                
                Q = self.subquiverRemoveVertex(i)
                P = self.subquiverRemoveVertex(j)
                
                if not Q.isAbundantAcyclic or not P.isAbundantAcyclic:
                    return False
                else:
                    return True
        
        return True
   
    @cached_property
    def isFork(self):
        '''
        This property is true if the principal quiver is a fork
        '''
        if (not self.isAbundant) or self.isAcyclic:
            return False
        return not (self._pointOfReturnFork == 0)
    
    @cached_property 
    def isPreFork(self):
        r'''
        This property is true if the principal quiver is a fork and is not a pre-fork
        '''
        return not (0 in self._preForkVertices)
        
    @cached_property
    def isWing(self):
        r'''
        This property is true if the principal quiver is a wing
        '''
        return not (0 in self._wingVertices)
    
    @cached_property
    def isTip(self):
        r'''
        This property is true if the principal quiver is a tip
        '''
        return not (o in self._tipVertices)
    
    # Methods that test a given mutation sequence
    
    def isReddeningSequence(self, mutSeq = []):
        r'''
        Takes in a mutation sequence and determines if it is a reddening sequence
        
        Input:
        - ''mutSeq'' -- the mutation sequence we are concerned with
        '''
        if self.hasFrozen:
            Q = self
        else:
            Q = self.framedQuiver()
        
        P = Q.mutate(mutSeq)
        red = set(P.reds)
        vertices = set(Q.mutableVertices)
        
        if vertices == red:
            return True
        
        return False
    
    def isLabeledCycle(self, mutSeq = []):
        r'''
        Takes in a mutation sequence and determines if it is a cycle in the labelled mutation graph
        
        Input:
        - ''mutSeq'' -- the mutation sequence
        '''
        P = self.mutate(mutSeq)
        
        if self.principalQuiver() == P.principalQuiver():
            return True
        
        return False
    
    # Methods that produce a mutation class or part of a mutation class from the quiver
    
    def forklessPart(self, depth):
        r'''
        Creates the forkless part of the labelled mutation graph of a quiver
        If the quiver has frozen vertices, it will only look at the principal quiver

        Input:
        - ''depth'' -- the number of mutations deep we would like to go. When using the function and suspect a quiver with a FFP, shoot for a large number;

        Output: ''G'' -- forkless part
        '''
        done = False # Variable to track if we have finished constructing the graph yet

        G = graphs.EmptyGraph() # This creates the labelled mutation graph using the sage Graph library. The mutation graph can have loops

        if self.hasZeroForklessPart:
            return G
        elif self.isFork: 
            por = self.pointOfReturn
            P = self.mutate([por])
            
            if P.isFork: # If P is a fork and mutating at the point of return produces another fork with different point of return, then we repeat the process
                G = P.forklessPart(depth)
                return G # Returns the forkless mutation graph of P, which will be same as the forkless mutation graph of Q
        elif self.hasFrozen: # Swaps to principal quiver if frozen vertices exist
            P = self.principalQuiver()
            G = P.forklessPart(depth)
            return G
        
        V = [] # This will be our vertex set
        E = {} # No edges yet, but we need a dict
        mutations = [] # We haven't done any mutations yet, but we need a list
        n = self.nMut # Finds our number of mutable vertices
        vertices = self.mutableVertices
        
        for i in [1..depth]:
            mutPrime = [] # Creates an empty list to contain our mutation sequences

            if i == 1: # If we're only looking at quivers one mutation away
                for j in vertices:
                    P = self.mutate([j]) # Sets the quiver P to the quiver formed from mutating at vertex j
                    if not P.isFork:
                        mutations.append([j]) # Adds mutation at j to the list of mutation sequences we have tried
                        V.append(copy(P)) # Appends our mutated quiver to the vertex set

                        if E.get(self, False) == False:
                            E[self] = {P : "mu_"+str(j)} # This adds an edge between vertex P and Q in our mutation graph
                        else:
                            E[self].update({P: "mu_" + str(j)})


                continue # Skips to i = 2 
            
            for mut in mutations: # Checks the mutations we've already done
                if len(mut) == i-1: # If the length of the mutation we picked out of the list is equal to i-1, we will use that as our base
                    mutPrime.append(mut) # Appends the mutation sequence mut to our list mutPrime
            
            newQuiver = False # Checks to see if a new quiver was mutated to

            for mut in mutPrime: # Goes through all the mutation sequences in mutPrime. Each will have length i-1
                for j in vertices:
                    if mut[i-2] != j: # This checks the last element of the mutation sequence. We check that it is not j, as mutating at j twice does nothing.
                        mutNew = mut.copy() # Creates a copy of the mutation sequence mut    
                        mutNew.append(j) # Adds mutation at j to the end of the mutation sequence mut
                        P = self.mutate(mutNew) # Finds the quiver found by mutating Q along the mutation sequence mutNew

                        if not P.isFork:
                            R = self.mutate(mut)

                            if P not in V: # If we mutated to a quiver not present in vertex set, we will add it
                                mutations.append(mutNew) # This adds the mutNew to the set of mutation sequences producing our vertex set
                                V.append(copy(P)) # This adds R to V

                            if E.get(R, False) == False:
                                E[R] = {P : "mu_"+str(j)} # This adds an edge between vertex P and R in our mutation graph
                            else:
                                E[R].update({P: "mu_" + str(j)})

                            newQuiver = True # We created a new quiver, so this becomes true

            if not newQuiver: # If we did not create a new quiver
                done = True # We are done
                print("Quiver has a finite forkless part")
                self._hasFiniteForklessPart = True
                break # Breaks out of the i loop

            print(i, "mutations deep") # This is a visual check that our program is running. It goes slower the further we go, so it is good to be cautious

        G = Graph(E, loops=True,  weighted=True) # This remakes our graph G to be the one with correct vertex set and edge set

        return G # Returns the graph
    
    # Methods involved with displaying the quiver
    
    @cached_property
    def digraph(self):
        r'''
        Returns a digraph from the given information about the quiver.
        Only works if the matrix is skew-symmetric
        
        Output: ''G'' -- digraph
        '''
        mutableVertices = self.mutableVertices
        vertices = self.vertices
        frozenVertices = self.frozenVertices
        
        N = zero_matrix(self.n)
        
        for i in mutableVertices:
            for j in vertices:
                if i == j:
                    continue
                
                a = self.arrow(i,j)
                if j in frozenVertices:
                    if a > 0:
                        N[i-1,j-1] = a
                    else:
                        N[j-1,i-1] = -1*a
                elif a > 0:
                    N[i-1,j-1] = a

        labels = {i-1 : i for i in vertices} # Relabels the graph as sage starts at 0
        G = DiGraph(N, multiedges=True, weighted=True) # Gives us a digraph
        G.relabel(labels,inplace=True) # Relabels the graph

        return G
    
    def plot(self, **kwds):
        r'''
        Builds the graphic object representing the quiver
        
        Output: ''P'' -- graphics object
        '''
        G = self.digraph
        if self.hasFrozen:
            vertices = self.mutableVertices
            greens = self.greens
            reds = self.reds

            vertexLabels = {'blue' : [],
                            'green' : [],
                            'red' : []}
            
            for v in G.vertices(sort=False):
                if v not in vertices:
                    vertexLabels['blue'].append(v)
                elif v in greens:
                    vertexLabels['green'].append(v)
                elif v in reds:
                    vertexLabels['red'].append(v)
                    
            P = G.plot(vertex_colors = vertexLabels, **kwds)
        else:
            P = G.plot(**kwds)
        return P

    def show(self, **kwds):
        r'''
        Displays the quiver, passing the given keywords to Graph.show()

        Output: whatever the return value of the graphics function is.
        '''
        P = self.plot(edge_labels=True, **kwds)
        return P.show()
   
