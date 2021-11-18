
from dataclasses import dataclass
import numpy as np



class ParticleSet(object):
    def __init__(self, N):
        self.N = N
        self.pos  = np.zeros((N,3), dtype=np.float64)
        self.vel  = np.zeros((N,3), dtype=np.float64)
        self.acc  = np.zeros((N,3), dtype=np.float64)
        self.pot  = np.zeros(N,     dtype=np.float64)
        self.mass = np.zeros(N,     dtype=np.float64)
        self.id   = np.zeros(N,     dtype=np.uint32)

    # permute the paticle set according to a permutation array
    # (used for sorting the particles)
    def permute(self, index):
        self.pos[:,:]  = self.pos[index,:]
        self.vel[:,:]  = self.vel[index,:]
        self.acc[:,:]  = self.acc[index,:]
        self.pot[:]  = self.pot[index]
        self.mass[:] = self.mass[index]
        self.id[:]  = self.id[index]

    # find the bounding box for the particles stored in this set
    def boundingCube(self):
        self.boxMin = np.min(self.pos) * np.ones(3) # keep bounding box larger than particles' range
        self.boxMax = np.max(self.pos) * np.ones(3) * (1+3e-10)
        self.center = 0.5*(self.boxMax + self.boxMin)
        self.halfWidth = 0.5*(self.boxMax[0]-self.boxMin[0])

class Morton():
    """
    Morton key class for sorting particles according to Z-order
    """
    def __init__(self, boxMin, boxMax):
        self.boxMin = boxMin
        self.scale = 1/(boxMax - boxMin)
        # some constants
        self.m1  = 0x3ffffffffff
        self.c64 = 0x3ff0000000000000000ffffffff
        self.c32 = 0x3ff00000000ffff00000000ffff
        self.c16 = 0x30000ff0000ff0000ff0000ff0000ff
        self.c8  = 0x300f00f00f00f00f00f00f00f00f00f
        self.c4  = 0x30c30c30c30c30c30c30c30c30c30c3
        self.c2  = 0x9249249249249249249249249249249

        # x-coordinate in Morton key changes most rapidly
        self.mask = np.array([ [-1,-1,-1],
                               [ 1,-1,-1],
                               [-1, 1,-1],
                               [ 1, 1,-1],
                               [-1,-1, 1],
                               [ 1,-1, 1],
                               [-1, 1, 1],
                               [ 1, 1, 1] ])
    def makeKey(self, pos):
        """
        Make a morton key as a binary number using bit-twiddling tricks
        """
        pscaled = (pos-self.boxMin) * self.scale
        # assume that pos is on [0,1)
        p = (pscaled*((1<<42)))
        key = 0
        for d in range(3):
            r = int(p[d])
            r &= self.m1
            r = (r | (r << 64)) & self.c64
            r = (r | (r << 32)) & self.c32
            r = (r | (r << 16)) & self.c16
            r = (r | (r << 8))  & self.c8
            r = (r | (r << 4))  & self.c4
            r = (r | (r << 2))  & self.c2
            key |= (r<<d)
        return key

    # get the octal value for level from key
    def getOct(self, key, level):
        shr = 123-3*(level)
        return (key>>shr)&7

@dataclass
class TreeNode:
    center: np.ndarray
    halfWidth: float
    begin: int
    end:int
    level: int
    parent: int
    firstChild: int
    childCount: int
    whichChild: int

    def __str__(self):
        return \
            f"    center: {self.center}\n" + \
            f" halfWidth: {self.halfWidth}\n" + \
            f"begin, end: [{self.begin},{self.end}]\n" + \
            f"     level: {self.level}\n" + \
            f"    parent: {self.parent}\n" + \
            f"whichChild: {self.whichChild}\n" + \
            f"firstChild: {self.firstChild}\n" + \
            f"childCount: {self.childCount}\n"

class Octree:
    def __init__(self, p, maxLeafSize, check):
        self.p = p
        self.maxLeafSize = maxLeafSize
        self.Morton = Morton(p.boxMin, p.boxMax)
        self.N = p.N

        # make Morton keys
        self.keys = np.empty( self.N, dtype='object')
        for i in range(self.N):
            self.keys[i] = self.Morton.makeKey(p.pos[i])

        # find permutation which sorts the keys
        self.index = np.argsort(self.keys)
        # get the inverse permutation
        self.rindex = np.empty(self.N, dtype=np.uint32)
        self.rindex[self.index] = np.arange(self.N, dtype=np.uint32)

        # permute the particles and the keys
        p.permute(self.index)
        self.keys[:] = self.keys[self.index]

        # perpare for first traversal
        self.nLeaves = 0
        self.nNodes = 1
        self.avgLeafSize = 1
        self.ROOT = 0
        ROOTlevel = 0

        self.levelCount = np.zeros(32, dtype=np.uint32)
        self.levelCount[ROOTlevel] = 1

        self.inOrderTraversal1(ROOTlevel, 0, self.N-1)

        print(f"nNodes: {self.nNodes}")

        for i, v in enumerate(self.levelCount):
            if v == 0:
                print(f"maximum depth of tree: {i}")
                break

        assert self.nNodes == self.levelCount.sum()

        self.avgLeafSize /= self.nLeaves
        print(f"avgLeafSize: {self.avgLeafSize}")
        print(f"nNodes: {self.nNodes}   nLeaves: {self.nLeaves}")

        # prepare for second traversal

        # create a tree array of empty nodes
        self.tree = [ TreeNode(0,0,0,0,0,0,0,0,0) for i in range(self.nNodes)]

        # fill the root node information
        self.tree[self.ROOT].center = p.center
        self.tree[self.ROOT].halfWidth = p.halfWidth
        self.tree[self.ROOT].begin = 0
        self.tree[self.ROOT].end = self.N-1
        self.tree[self.ROOT].parent = np.Inf
        self.tree[self.ROOT].firstChild = 0
        self.tree[self.ROOT].childCount = 0
        self.tree[self.ROOT].level = 0

        # fill in pointers to where first node at each level lies in tree array
        self.levelPointer = np.zeros(32, dtype='int')
        for i in range(1,32):
            self.levelPointer[i] = self.levelPointer[i-1] + self.levelCount[i-1]

        # pointer to location of current node
        self.nodePointer = 1

        self.inOrderTraversal2(self.tree[0].level, 0, self.N-1, 0, self.tree[0].center, 0.5*self.tree[0].halfWidth)
        assert self.nodePointer == self.nNodes

        if check:
            self.checkCount = 0
            if not self.checkTree(self.ROOT, True):
                raise ValueError("tree does not check out")
            elif self.checkCount != self.nNodes:
                raise ValueError("tree does not check out")
            else:
                print(f"octree checks out")

    def inOrderTraversal1(self, parentLevel, begin, end):
        """
        Traverse Morton keys in a depth-first search to count number of nodes at each level
        """
        count = 0

        # get direction of first child from direction of first point in that child
        direction = self.Morton.getOct(self.keys[begin], parentLevel)

        # count the number of points within each child
        while direction <= 7 and (begin<=end):
            count = 0;
            while begin <= end:
                if self.Morton.getOct(self.keys[begin], parentLevel) == direction:
                    begin += 1
                    count += 1
                else:
                    break

            self.levelCount[parentLevel+1] += 1
            self.nNodes += 1

            # if the cell isn't big enough to have children, it is a leaf cell
            if count <= self.maxLeafSize:
                self.nLeaves += 1
                self.avgLeafSize += count
            else:
                # otherwise, keep travwersing
                self.inOrderTraversal1(parentLevel+1, begin-count, begin-1)

            # get the next direction
            if begin <= end:
                direction = self.Morton.getOct(self.keys[begin], parentLevel)

    def inOrderTraversal2(self, parentLevel, begin, end, parent, center, halfWidth):
        """
        Fill all elements of the tree list
        """

        assert parentLevel < 32
        assert begin < end
        assert self.tree[parent].firstChild == 0 # shouldn't be here if parent had any children
        assert self.tree[parent].childCount == 0

        direction = self.Morton.getOct(self.keys[begin], parentLevel)

        while begin <= end:
            assert direction <= 7

            # count number of points in this octant (direction)
            count = 0;
            while begin <= end:
                if self.Morton.getOct(self.keys[begin], parentLevel) == direction:
                    begin += 1
                    count += 1
                else:
                    break
            assert count > 0

            # get child node number in tree
            child = self.levelPointer[parentLevel+1]
            # advance counter to next child node location
            self.levelPointer[parentLevel+1] += 1

            if self.tree[parent].firstChild == 0:
                # first child of parent
                assert self.tree[parent].childCount == 0;
                self.tree[parent].firstChild = child
                self.tree[parent].childCount = 1
            else:
                # subsequent children of parent
                self.tree[parent].childCount += 1
                assert self.tree[parent].childCount <= 8

            self.tree[child].level = parentLevel+1
            self.tree[child].parent = parent
            self.tree[child].begin = begin - count
            self.tree[child].end = begin - 1
            self.tree[child].halfWidth = halfWidth
            self.tree[child].center = center + halfWidth*self.Morton.mask[direction]
            self.tree[child].whichChild = direction
            self.tree[child].firstChild = 0
            self.tree[child].childCount = 0

            self.nodePointer += 1
            assert self.nodePointer <= self.nNodes

            if count > self.maxLeafSize:
                self.inOrderTraversal2(self.tree[child].level, begin-count, begin-1, child, \
                                       self.tree[child].center, 0.5*self.tree[child].halfWidth)


            if begin <= end:
                direction = self.Morton.getOct(self.keys[begin], parentLevel)

    def checkTree(self, node, retval):

        if not retval: return False

        lo = self.tree[node].center - self.tree[node].halfWidth
        hi = self.tree[node].center + self.tree[node].halfWidth

        # need to watch out for floating-point comparisons
        eps = 1e-13*self.tree[node].halfWidth

        b = self.tree[node].begin
        e = self.tree[node].end

        npts = e-b+1

        if npts <= 0:
            print(f"node {node} has no particles: {npts}")

        if self.tree[node].childCount>0 and npts<=self.maxLeafSize:
            print(f"node {node} is a leaf with {npts} > maxLeafSize")

        retval = True
        for d in range(3):
            if not ( np.all(self.p.pos[b:e+1,d]>=lo[d]-eps) and np.all(self.p.pos[b:e+1,d]<=hi[d]+eps) ):
                print(f"node {node} contains an illegal value:")
                print(f"lo: {lo}   hi: {hi}")
                for i in range(b, e+1):
                    if not (self.p.pos[i,d]>=lo[d] and self.p.pos[i,d]<=hi[d]):
                        print(f"[{i},{d}]: {lo[d]} <= {self.p.pos[i,d]} < {hi[d]}")
                        retval = False

        self.checkCount += 1
        for daughter in range(*self.allChildren(node)):
            self.checkTree(daughter, retval)

        return True

    def isLeaf(self, node):
        return (self.tree[node].childCount == 0)

    def nodeFirst(self, node):
        return self.tree[node].begin

    def nodeLast(self, node):
        return self.tree[node].end

    def nodeSize(self, node):
        return self.tree[node].end - self.tree[node].begin + 1

    def allPoints(self, node):
        """
        Return the range of particle indices in node, suitable for the range function
        (i.e. the end is one past the last particle)
        """
        return self.tree[node].begin, self.tree[node].end + 1

    def allChildren(self, node):
        """
        Return the range of child node indices, suitable for the range function
        (i.e. the end is one past the last node)
        """
        return self.tree[node].firstChild, self.tree[node].firstChild + self.tree[node].childCount


class BHtree:

    def __init__(self, p, maxLeafSize, epsSmooth):
        """
        Build the octree
        Compute the corresponding monopole information (masses, COM's, and Bmax's of each node)
        Define the traversals of the tree to compute acceleration and potential
        """

        self.p = p
        self.N = p.N
        self.eps2 = epsSmooth**2

        # start by building the octree
        self.octree = Octree(p, maxLeafSize, False)

        # get space for the data
        self.com = np.zeros((self.octree.nNodes,3), dtype=np.float64)
        self.mass = np.zeros(self.octree.nNodes, dtype=np.float64)
        self.Bmax = np.zeros(self.octree.nNodes, dtype=np.float64)

        self.propagateCOM(self.octree.ROOT)
        self.propagateBmax(self.octree.ROOT)

    def propagateCOM(self, node):

        if self.octree.isLeaf(node):
            for i in range(*self.octree.allPoints(node)):
                self.com[node] += self.p.pos[i] * self.p.mass[i, np.newaxis]
                self.mass[node] += self.p.mass[i]

        else:
            for d in range(*self.octree.allChildren(node)):
                self.propagateCOM(d)
                self.com[node] += self.com[d] * self.mass[d]
                self.mass[node] += self.mass[d]

        self.com[node] /= self.mass[node]

    def propagateBmax(self, node):
        beg, end = self.octree.allPoints(node)
        self.Bmax[node] = np.max( np.linalg.norm(self.p.pos[beg:end] - self.com[node]) )
        if not self.octree.isLeaf(node):
            for d in range(*self.octree.allChildren(node)):
                self.propagateBmax(d)

    def bhSingle(self, sink, node, particleList, nodeList):
        if self.octree.isLeaf(node):
            particleList += list(range(*self.octree.allPoints(node)))
        else:
            dr = np.linalg.norm(self.p.pos[sink] - self.com[node])
            if self.Bmax[node] < self.theta * dr:
                nodeList.append(node)
            else:
                for d in range(*self.octree.allChildren(node)):
                    self.bhSingle(sink, d, particleList, nodeList)

    def direct(self, sinkpos, srcpos, srcmass, index, eps2):
        drvec = sinkpos - srcpos[index]                               # drvec = r_sink - rsrc_i
        dr2 = np.sum(drvec*drvec, axis=1)                             # |dr|^2 = drvec dot drvec
        ir = 1/np.sqrt(dr2 + eps2)                                    # (|dr|^2 + eps^2)^(-3/2)
        return -np.sum( ir*ir*ir * drvec.T * srcmass[index], axis=1), -np.sum(srcmass[index]*ir)

    def acc1(self, sink):
        particleList = []
        nodeList = []
        self.bhSingle(sink, self.octree.ROOT, particleList, nodeList)
        acc,pot = self.direct(self.p.pos[sink], self.p.pos, self.p.mass, np.asarray(particleList, dtype=np.int32), self.eps2)
        dacc,dpot = self.direct(self.p.pos[sink], self.com, self.mass, np.asarray(nodeList, dtype=np.int32), self.eps2)
        self.p.acc[sink] = acc + dacc
        self.p.pot[sink] = pot + dpot

    def accAll(self, theta):
        self.theta = theta
        for i in range(self.N):
            self.acc1(i)
