import numpy as np


def transformToCon(filename, data, double=False):
    shapeofdata = np.shape(data)
    reslst = []
    lx = shapeofdata[4]
    ly = shapeofdata[3]
    lz = shapeofdata[2]
    lt = shapeofdata[1]
    for x in range(lx):
        for y in range(ly):
            for z in range(lz):
                for t in range(lt):
                    for l in range(4):
                        for mx in range(3):
                            for my in range(3):
                                reslst.append(np.real(data[l, t, z, y, x, mx, my]))
                                reslst.append(np.imag(data[l, t, z, y, x, mx, my]))
    reslst = np.array(reslst)
    data_float32 = reslst.astype(np.float32 if not double else np.float64)
    data_float32.tofile(filename)

def transformFromCon(filename, lx, ly, lz, lt, double=False):
    data = np.fromfile(filename, dtype=np.float32 if not double else np.float64)
    print(np.shape(data))
    arr = np.zeros((4, lt, lz, ly, lx, 3, 3), dtype=np.complex64)
    for x in range(lx):
        for y in range(ly):
            for z in range(lz):
                for t in range(lt):
                    for l in range(4):
                        idxnow = x * (ly * lz * lt * 4 * 9) + y * (lz * lt * 4 * 9) + z * (lt * 4 * 9) + t * 36 + l * 9
                        for mx in range(3):
                            for my in range(3):
                                arr[l, t, z, y, x, mx, my] = data[(idxnow + mx * 3 + my) * 2] + data[(idxnow + mx * 3 + my) * 2 + 1] * 1j
    return arr

def transformStaggeredFermionFromCon(filename, lx, ly, lz, lt, double=False):
    data = np.fromfile(filename, dtype=np.float32 if not double else np.float64)
    print(np.shape(data))
    arr = np.zeros((lt, lz, ly, lx, 3), dtype=np.complex64 if not double else np.complex128)
    for x in range(lx):
        for y in range(ly):
            for z in range(lz):
                for t in range(lt):
                    idxnow = x * (ly * lz * lt * 3) + y * (lz * lt * 3) + z * (lt * 3) + t * 3
                    for c in range(3):
                        arr[t, z, y, x, c] = data[(idxnow + c) * 2] + data[(idxnow + c) * 2 + 1] * 1j
    return arr

def transformStaggeredFermionToCon(filename, data, double=False):
    lt, lz, ly, lx, nc = data.shape
    assert nc == 3
    reslst = []
    for x in range(lx):
        for y in range(ly):
            for z in range(lz):
                for t in range(lt):
                    for c in range(3):
                        reslst.append(np.real(data[t, z, y, x, c]))
                        reslst.append(np.imag(data[t, z, y, x, c]))
    reslst = np.array(reslst)
    out = reslst.astype(np.float32 if not double else np.float64)
    out.tofile(filename)