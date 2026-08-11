"""
* 0: 1
* 1: (-1)^nx
* 2: (-1)^ny
* 3: (-1)^nz
* 4: (-1)^nt
* 5: (-1)^(nx+ny)
* 6: (-1)^(nx+nz)
* 7: (-1)^(nx+nt)
* 8: (-1)^(ny+nz)
* 9: (-1)^(ny+nt)
* 10: (-1)^(nz+nt)
* 11: (-1)^(ny+nz+nt)
* 12: (-1)^(nx+nz+nt)
* 13: (-1)^(nx+ny+nt)
* 14: (-1)^(nx+ny+nz)
* 15: (-1)^(nx+ny+nz+nt)

"""


def all_flavours():
    return [["uu"], ["dd"], ["du", "ud"]]

def all_case():
    return [[0], [11], [2, 3], [4], [8], [9, 10]]

def channel_pm():
    return [1, 2, 1, 1, 2, 2]