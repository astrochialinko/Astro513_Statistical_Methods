import numpy as np
import ctypes

"""
Wrapper classes for ParticleSet and BHTree
"""

# This is a Python function called by the C++ code
# in order to print in a notebook
@ctypes.CFUNCTYPE(None, ctypes.c_char_p)
def print_callback(a):
    print(a.decode(), end='')

# some abbreviations for types
c_double_array = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags=['C_CONTIGUOUS','ALIGNED'])
c_double2_array = np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags=['C_CONTIGUOUS','ALIGNED'])
c_int_array = np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags=['C_CONTIGUOUS','ALIGNED'])
c_fpointer = ctypes.POINTER(ctypes.c_char)
c_int = ctypes.c_int
c_double = ctypes.c_double
c_bool = ctypes.c_bool
c_double_pointer = ctypes.POINTER(c_double)
c_char_p = ctypes.c_char_p


def createDtype(descString, _verbose):
    """
    Take the string returned from ParticleSet::getTypeString and
    create a dictionary of variable names, types, and widths from
    which the corresponding Numpy arrays can be made.
    """

    # C++ types likely to be found in output, mapped to numpy type and length of each element
    # code will throw an error if an unknown type is found.
    _cpptypes = { "SmallVec<double, 3>":['float64',3], "SmallVec<double, 2>":['float64',2], "SmallVec<double, 1>":['float64',1],
                  "SmallVec<float, 3>":['float32',3], "SmallVec<float, 2>":['float32',2], "SmallVec<float, 1>":['float32',1],
                  "SmallVec<int, 3>":['int32',3], "SmallVec<int, 2>":['int32',2], "SmallVec<int, 1>":['int32',1],
                  "double":['float64',1], "float":['float32',1], "int":['int32',1], "unsigned int":['uint32',1], "long int":['uint64',1], "long unsigned int":['uint64',1] }

    _desc = descString.split(";")

    # dictionary of type descriptors
    _data = {}
    for _d in _desc:
        _td, _variable = _d.split(":")
        try:
            _t, _n = _cpptypes[_td]                    # type and multiplicity of each variable
        except KeyError:
            print(f"\nParticleSet.py: Error reading {_fname}: unknown type '{_td}'\n")
            raise

        _data[_variable] = [_t, _n]                    # numpy type and multiplicity of each variable

    _nvars = len(_data)
    if _verbose: print(f"    nvars = {_nvars}")
    if _verbose: print(f"     dict = {_data}")

    return _data


class ParticleSet:
    """
    Wrapper for the ParticleSet class.

    Numpy arrays are created in the reserve(N) function and their pointers are then passed
    to the C++ code.

    NB: One must be very careful not to redefine these Numpy ndarrays:
        use, e.g.: ps.r[:,:] = something
        not      : ps.r = something, which redefines ps.r, which will no longer point to the C++ array
                   ( and which may be garbage-collected, causing the C++ code to segfault)

    """
    def __init__(self):

        # load the shared library; windows users will need to modify this:
        try:
            self.dll = ctypes.cdll.LoadLibrary('./treelib.so')
        except OSError:
            self.dll = ctypes.cdll.LoadLibrary('./treelib.dylib')

        # the constructor; returns a pointer to the class instance
        self.dll.newParticleSetDouble.restype = c_fpointer
        self.dll.newParticleSetDouble.argtypes = None;

        # call the constructor
        # need to pass length for plain C arrays
        self.me = self.dll.newParticleSetDouble()

        self.dll.freeParticleSet.restype = None
        self.dll.freeParticleSet.argtypes = [ c_fpointer ]

        self.dll.setPSDStorage.restype = None
        self.dll.setPSDStorage.argtypes = [ c_fpointer,
                                            c_double2_array, c_double2_array, c_double2_array,
                                            c_double_array, c_double_array, c_int_array,
                                            c_int ]

        self.dll.computeBoundingBox.restype = None
        self.dll.computeBoundingBox.argtypes = [ c_fpointer ]
        self.dll.getBoundingBox.restype = None
        self.dll.getBoundingBox.argtypes = [ c_fpointer, c_double_array, ctypes.POINTER(c_double) ]

        self.dll.getTypeString.restype = c_char_p
        self.dll.getTypeString.argtypes = [ c_fpointer ]

        self.dll.read.restype = None
        self.dll.read.argtypes = [ c_fpointer, c_char_p ]

        self.dll.write.restype = None
        self.dll.write.argtypes = [ c_fpointer, c_char_p ]

    def free(self):
        self.dll.freeParticleSet(self.me)

    def getTypeString(self):
        tmp = self.dll.getTypeString(self.me)
        return tmp.decode("utf-8")

    def reserve(self, N, printVars=False):
        """
        Create the variables in the particle set and send their pointers to the C++ library
        """
        self.N = N

        dtype = self.getTypeString()
        _data = createDtype(dtype, False)

        for _v, _d in _data.items():
            if _d[1] == 1:
                cmd = f"self.{_v} = np.zeros(N, dtype='{_d[0]}')"
            else:
                cmd = f"self.{_v} = np.zeros((N,{_d[1]}), dtype='{_d[0]}')"
            if printVars:
                print(cmd)
            exec(cmd)

        self.dll.setPSDStorage(self.me, self.r, self.v, self.a, self.pot, self.mass, self.id, self.N)

    def computeBoundingBox(self):
        self.dll.computeBoundingBox(self.me)

    def getBoundingBox(self):
        halfWidth = c_double()
        halfWidth_p = ctypes.pointer(halfWidth)
        centre = np.zeros(3, dtype='float64')
        self.dll.getBoundingBox(self.me, centre, halfWidth_p)
        return centre, halfWidth.value

    def read(self, fname):
        b_fname = fname.encode('utf-8')
        self.dll.read(self.me, b_fname)

    def write(self, fname):
        b_fname = fname.encode('utf-8')
        self.dll.write(self.me, b_fname)


class BHTree:
    """
    Python wrapper for the BHTree class.

    NB: While the C++ code will automatically destroy a BHTree instance when it goes out of scope,
        there is no way to make this happen in the wrapped code.
        -> You must explicitly call BHTree.free() to do this when you are finished with the class
           instance. Otherwise, you will create a memory leak!

    """
    def __init__(self, pset, maxLeafSize, epsSmooth, maxSources):

        # load the shared library; windows users will need to modify this:
        try:
            self.dll = ctypes.cdll.LoadLibrary('./treelib.so')
        except OSError:
            self.dll = ctypes.cdll.LoadLibrary('./treelib.dylib')

        # the constructor; returns a pointer to the class instance
        self.dll.newBHTreeDouble.restype = c_fpointer
        self.dll.newBHTreeDouble.argtypes = [c_fpointer, c_int, c_double, c_int ]

        # call the constructor
        # need to pass length for plain C arrays
        self.me = self.dll.newBHTreeDouble(pset.me, maxLeafSize, epsSmooth, maxSources)

        self.dll.BHfree.restype = None
        self.dll.BHfree.argtypes = [ c_fpointer ]

        self.dll.makeTree.restype = None
        self.dll.makeTree.argtypes = [ c_fpointer, c_bool, c_bool, c_bool ]

        self.dll.BHsubsets.restype = None
        self.dll.BHsubsets.argtypes = [ c_fpointer, c_double, c_int ]

        self.dll.rspAcc.restype = None
        self.dll.rspAcc.argtypes = [ c_fpointer, c_double, c_int ]

        self.dll.accAll.restype = None
        self.dll.accAll.argtypes = [ c_fpointer, c_double ]

        # send the address of the callback print function to the C++ code
        self.dll.registerPrint(self.me, print_callback)

        self.dll.putParticlesInOriginalOrder.restype = None
        self.dll.putParticlesInOriginalOrder.argtypes = [ c_fpointer ]


    def free(self):
        self.dll.BHfree(self.me)

    def makeTree(self, verbose=False, check=False, getStats=False):
        self.dll.makeTree(self.me, verbose, check, getStats)

    def rspAcc(self, theta, maxInteractionLength):
        self.dll.rspAcc(self.me, theta, maxInteractionLength)

    def BHsubsets(self, theta, maxInteractionLength):
        self.dll.BHsubsets(self.me, theta, maxInteractionLength)

    def accAll(self, theta):
        self.dll.accAll(self.me, theta);

    def putParticlesInOriginalOrder(self):
        self.dll.putParticlesInOriginalOrder(self.me)

class Conservation:
    def __init__(self, ps):
        self.ps = ps
        self.K0, self.P0, self.E0 = self.getEnergy()

    def getEnergy(self):
        K = 0.5*np.sum(self.ps.mass * np.sum(self.ps.v**2, axis=1))
        P = 0.5*np.sum(self.ps.mass * self.ps.pot)
        return K, P, K+P

class Hamiltonian:
    def __init__(self, theta, epsSmooth, maxLeafSize, maxSrcLen):
        self.theta = theta
        self.epsSmooth = epsSmooth
        self.maxLeafSize = maxLeafSize
        self.maxSrcLen = maxSrcLen

        self.verbose = False
        self.check = False
        self.getStats = False

    def getAcceleration(self, ps):

        ps.computeBoundingBox()
        center, halfWidth = ps.getBoundingBox()

        BH = BHTree(ps, self.maxLeafSize, self.epsSmooth, self.maxSrcLen)

        BH.makeTree(self.verbose, self.check, self.getStats)

        BH.BHsubsets(self.theta, ps.N)     # use leaves as sink sets
        #BH.accAll(self.theta)             # use single particles as sinks

        BH.putParticlesInOriginalOrder()  # so we know which particles are which galaxy

        BH.free()                         # important to free memory on the C side...

    def positionEquation(self, ps, newvel):
        return ps.v

    def momentumEquation(self, ps, newacc):
        if newacc: self.getAcceleration(ps)
        return ps.a


class State:
    def __init__(self, ps, hamilton ):
        self.time = 0.0
        self.step = 0
        self.ps = ps

    def kick(self, h, hamilton, recalculate):
        # Update velocity
        self.ps.v += h * hamilton.momentumEquation(self.ps, recalculate)
        return self

    def drift(self, h, hamilton, recalculate):
        # Update positions
        self.ps.r += h * hamilton.positionEquation(self.ps, recalculate)
        return self

def KDK(dt, hamilton, s):
    s = s.kick(dt/2, hamilton, False).drift(dt, hamilton, False).kick(dt/2, hamilton, True)
    s.time += dt
    s.step += 1
    return s

if __name__ == "__main__":

    import matplotlib.pyplot as plt
    import numpy as np
    import numpy.random as rng
    from plummer import PlummerModel

    rng.seed(2394820)

    N = 2**20
    PS = ParticleSet()
    PS.reserve(N)

    PS.r[:,:], PS.v[:,:], PS.mass[:] = PlummerModel(N)

    PS.computeBoundingBox()
    centre, halfWidth = PS.getBoundingBox()
    print(f"centre: {centre}   halfWidth: {halfWidth}")

    maxLeafSize = 32
    epsSmooth = 0.98*N**(-0.26)
    maxSources = N

    BH = BHTree(PS, maxLeafSize, epsSmooth, maxSources)

    BH.makeTree()

    BH.accAll(0.6)

    BH.free()

    fig, ax = plt.subplots()
    r = np.linalg.norm(PS.r, axis=1)
    v = np.linalg.norm(PS.v, axis=1)
    a = np.linalg.norm(PS.a, axis=1)
    ax.loglog(r, a, ',')
    plt.show()
