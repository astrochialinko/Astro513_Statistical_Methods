"""
Class to read a .pset file.

See __main__ below for an example of use.

"""



import numpy as np
import struct


def createDtype(descString, _verbose):

    # C++ types likely to be found in output, mapped to numpy type and length of each element
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
            _t, _n = _cpptypes[_td]                    # flag and multiplicity of each variable
        except KeyError:
            print(f"\nreadParticleSet: Error reading {_fname}: unknown type '{_td}'\n")
            raise

        _data[_variable] = [_t, _n]                    # numpy type and multiplicity of each variable

    _nvars = len(_data)
    if _verbose: print(f"    nvars = {_nvars}")
    if _verbose: print(f"     dict = {_data}")

    _npdtype = []                                          # numpy structured array dtype
    for _variable in _data:
        _s = _data[_variable][1]                           # length of each element
        _dtype = _data[_variable][0]                       # dtype of element
        if _s > 1:
            _npdtype.append( ( _variable, _dtype, (_s) ) )
        else:
            _npdtype.append( ( _variable, _dtype ) )

    return _npdtype, _variable, _data


def readParticleSet(_fname, _verbose=False):
    """
    Read the data from the particleset in the file "fname"
    The number of particles N and the centre and halfWidth are returned, as is
    a numpy structured array whose entries are named the same as the variables in a Ptype.

    The function call is then:
       N, centre, halfWidth, data = readParticleSet(filename)

    Calling readParticleSet with verbose=True as an additional parameter will print the numpy calls used
    to retrieve the data.
    """

    # C++ types likely to be found in output, mapped to numpy type and length of each element
    _cpptypes = { "SmallVec<double, 3>":['float64',3], "SmallVec<double, 2>":['float64',2], "SmallVec<double, 1>":['float64',1],
                  "SmallVec<float, 3>":['float32',3], "SmallVec<float, 2>":['float32',2], "SmallVec<float, 1>":['float32',1],
                  "SmallVec<int, 3>":['int32',3], "SmallVec<int, 2>":['int32',2], "SmallVec<int, 1>":['int32',1],
                  "double":['float64',1], "float":['float32',1], "int":['int32',1], "unsigned int":['uint32',1], "long int":['uint64',1], "long unsigned int":['uint64',1] }

    if _verbose: print(f'\nreadPset reading from file "{_fname}":')

    bytes = 0
    with open(_fname, "rb") as f:
        _strlen, = struct.unpack('i',f.read(4))     # length of arrays descriptor string from file
        bytes += 4
        _desc = f.read(_strlen+1)[:-1].decode()       # list of arrays descriptor entries
        bytes += _strlen+1
        if _verbose: print(f"\nDescriptor from file:\n     desc = {_desc}");

        # number of particles
        N, = struct.unpack('i',f.read(4))
        bytes+=4

        # center of volume
        centre = np.array(struct.unpack('ddd',f.read(3*8)), dtype=np.float64)
        bytes+=3*8

        # half-width
        halfWidth, = struct.unpack('d', f.read(8))
        bytes+=8
        print(f"header length: {bytes}")

        _npdtype, _variable, _data = createDtype(_desc, _verbose)

        if _verbose: print(f"\nNumpy dtype for data array:\n    dtype = {_npdtype}\n")
        bigdata = np.empty(N, dtype=_npdtype)                  # create the array


        if _verbose: print("Data read from file with lines:")
        for _variable in _data:
            _s = _data[_variable][1]
            _dtype = _data[_variable][0]
            if _s > 1:
                bigdata[_variable] = np.fromfile(f, dtype=_dtype, count=_s*N).reshape((N,_s))
                bytes += _s*N*8
                if _verbose: print(f"   data['{_variable}'] = np.fromfile(f, dtype='{_dtype}', count={_s}*N).reshape((N,{_s}))")
            else:
                bigdata[_variable] = np.fromfile(f, dtype=_dtype, count=N)
                if(_variable != 'id'):
                    bytes += N*8
                else:
                    bytes += N*4
                if _verbose: print(f"   data['{_variable}'] = np.fromfile(f, dtype='{_dtype}', count=N)")

        print(f"read {bytes} bytes in all")
        if _verbose: print("\n")


    return N, centre, halfWidth, bigdata

if __name__ == '__main__':
    import sys

    def pnp(a, form):
        return np.array2string(a, formatter={'float': lambda x: f'{x:{form}}'}, separator=',')

    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <pset filename>")
        raise ValueError("bad input")

    N, centre, halfWidth, pset = readParticleSet(sys.argv[1], True)

    print(f"N: {N}")
    print(f"centre: {centre}")
    print(f"halfWidth: {halfWidth}")
    for i in range(10):
        print(f"{i: 4d}  {pnp(pset['r'][i],' .4e')}  {pnp(pset['v'][i],' .4e')}  {pnp(pset['a'][i],' .4e')} {pset['mass'][i]:8.4f}  {pset['id'][i]: 4d}")
