import numpy as np

rawdata = np.load('data/cfg_100.npy')
# onelinktest = np.load('data/cfg_100_onelink_test.npy')
fat3data = np.load('data/cfg_100_fat3.npy')

# link, t, z, y, x, nc, nc
print(np.shape(rawdata))
# print(np.shape(onelinktest))

xdir = 0
ydir = 1
zdir = 2
tdir = 3

site_x, site_y, site_z, site_t = 1, 2, 3, 4

# print(rawdata[xdir, site_t, site_z, site_y, site_x])
# print(onelinktest[xdir, site_t, site_z, site_y, site_x])

print(fat3data[xdir, site_t, site_z, site_y, site_x])
#
def forward(m1, m2, m3):
    l = m1
    l = np.dot(l, m2)
    l = np.dot(l, m3.conj().transpose())
    return l

def backward(m1, m2, m3):
    l = m1.conj().transpose()
    l = np.dot(l, m2)
    l = np.dot(l, m3)
    return l

# yx-y
m11 = rawdata[ydir, site_t, site_z, site_y, site_x]
m12 = rawdata[xdir, site_t, site_z, site_y + 1, site_x]
m13 = rawdata[ydir, site_t, site_z, site_y, site_x + 1]
res = forward(m11, m12, m13)

m21 = rawdata[ydir, site_t, site_z, site_y - 1, site_x]
m22 = rawdata[xdir, site_t, site_z, site_y - 1, site_x]
m23 = rawdata[ydir, site_t, site_z, site_y - 1, site_x + 1]
res = res + backward(m21, m22, m23)

# zx-z
m31 = rawdata[zdir, site_t, site_z, site_y, site_x]
m32 = rawdata[xdir, site_t, site_z + 1, site_y, site_x]
m33 = rawdata[zdir, site_t, site_z, site_y, site_x + 1]
res = res + forward(m31, m32, m33)

m41 = rawdata[zdir, site_t, site_z - 1, site_y, site_x]
m42 = rawdata[xdir, site_t, site_z - 1, site_y, site_x]
m43 = rawdata[zdir, site_t, site_z - 1, site_y, site_x + 1]
res = res + backward(m41, m42, m43)

# tx-t
m51 = rawdata[tdir, site_t, site_z, site_y, site_x]
m52 = rawdata[xdir, site_t + 1, site_z, site_y, site_x]
m53 = rawdata[tdir, site_t, site_z, site_y, site_x + 1]
res = res + forward(m51, m52, m53)

m61 = rawdata[tdir, site_t - 1, site_z, site_y, site_x]
m62 = rawdata[xdir, site_t - 1, site_z, site_y, site_x]
m63 = rawdata[tdir, site_t - 1, site_z, site_y, site_x + 1]
res = res + backward(m61, m62, m63)

print(res)



