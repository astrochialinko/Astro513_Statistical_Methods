import numpy as np
import ctypes

"""
Wrapper class for MCRT3
"""


# This is a Python function called by the C++ code
# in order to print in a notebook
@ctypes.CFUNCTYPE(None, ctypes.c_char_p)
def print_callback(a):
    print(a.decode(), end='')


class RT:

    def __init__(self, r, v, rho):

        # some abbreviations for types
        c_double_array = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS')
        c_fpointer = ctypes.POINTER(ctypes.c_char)
        c_int = ctypes.c_int
        c_double = ctypes.c_double
        c_double_pointer = ctypes.POINTER(c_double)

        # load the shared library
        self.dll = ctypes.cdll.LoadLibrary('./libmcrt.so')

        # the constructor; returns a pointer to the class instance
        self.dll.newRT.restype = c_fpointer
        self.dll.newRT.argtypes = [c_double_array, c_double_array, c_double_array, c_int]

        # the remaining wrapper functions
        self.dll.setSpectrumGrid.restype = None
        self.dll.setSpectrumGrid.argtypes = [c_fpointer, c_double, c_double, c_int]

        self.dll.doSimulation.restype = None
        self.dll.doSimulation.argtypes = [c_fpointer, c_int, ctypes.c_ulonglong]

        self.dll.getSpectrum.restype = None
        self.dll.getSpectrum.argtypes = [c_fpointer, c_double_array,c_double_array, c_double_array]

        self.dll.getDeposition.restype = None
        self.dll.getDeposition.argtypes = [c_fpointer, c_double_array, c_double_array]

        self.dll.getPdV.restype = None
        self.dll.getPdV.argtypes = [c_fpointer, c_double_array, c_double_array]

        self.dll.getEmitSpec.restype = None
        self.dll.getEmitSpec.argtypes = [c_fpointer, c_double_array, c_double_array, c_double_array]

        self.dll.getInner.restype = None
        self.dll.getInner.argtypes = [c_fpointer, c_double_pointer, c_double_pointer]

        # call the constructor
        # need to pass length for plain C arrays
        self.nz = r.shape[0]
        self.me = self.dll.newRT(r, v, rho, self.nz)

        # send the address of the callback print function to the C++ code
        self.dll.registerPrint(self.me, print_callback)

    #--------------------------------------------------------------------------------
    #
    # Python functions corresponding to the wrapped C functions
    #
    #--------------------------------------------------------------------------------

    def setSpectrumGrid(self, eMin, eMax, ne):
        self.ne = ne
        self.dll.setSpectrumGrid(self.me, eMin, eMax, ne);

    def doSimulation(self, N, seed):
        self.dll.doSimulation(self.me, N, seed)

    def getSpectrum(self):
        # allocate space to recieve data from C++ code
        egrid = np.zeros(self.ne)
        spectrum = np.zeros(self.ne)
        spectrumsig = np.zeros(self.ne)
        self.dll.getSpectrum(self.me, egrid, spectrum, spectrumsig)
        return egrid, spectrum, spectrumsig

    def getDeposition(self):
        deposition = np.zeros(self.nz)
        depositionsig = np.zeros(self.nz)
        self.dll.getDeposition(self.me, deposition, depositionsig)
        return deposition, depositionsig

    def getPdV(self):
        PdV = np.zeros(self.nz)
        PdVsig = np.zeros(self.nz)
        self.dll.getPdV(self.me, PdV, PdVsig)
        return PdV, PdVsig

    def getEmitSpec(self):
        egrid = np.zeros(self.ne)
        emitSpec = np.zeros(self.ne)
        emitSpecsig = np.zeros(self.ne)
        self.dll.getEmitSpec(self.me, egrid, emitSpec, emitSpecsig)
        return egrid, emitSpec, emitSpecsig

    def getInner(self):
        inner = 0
        sig = 0
        self.dll.getInner(self.me, inner, sig)


if __name__ == "__main__":

    import matplotlib.pyplot as plt

    clight = 2.99792458e8

    nz = 11
    rmin, rmax = 0, 10
    vmin, vmax = 0, 0#.001*clight
    rhomin, rhomax = 1,1
    r = np.linspace(rmin, rmax, nz)
    v = np.linspace(vmin, vmax, nz)
    rho = np.linspace(rhomin, rhomin, nz)

    ne = 101
    eMin, eMax = 0.8, 1.2

    seed = 2873642343;
    N0 = 2**20

    sim = RT(r, v, rho)
    sim.setSpectrumGrid(eMin, eMax, ne)
    sim.doSimulation(N0,seed)

    egrid, spec, sig = sim.getSpectrum()
    egrid = egrid[:-1]
    spec = spec[:-1]
    sig = sig[:-1]

    fig, ax = plt.subplots()
    offset = (egrid[1]-egrid[0])/2
    ax.step(egrid+offset, spec, 'b')
    plt.show()
