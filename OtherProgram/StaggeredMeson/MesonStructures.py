import numpy as np


def eta_i(i):
    return np.array([1 if (j + 1) < i else 0 for j in range(4)])

def xi_i(i):
    return np.array([1 if (j + 1) > i else 0 for j in range(4)])

def eta_i3(i):
    return np.array([1 if (j + 1) < i else 0 for j in range(3)])

def xi_i3(i):
    return np.array([1 if (j + 1) > i else 0 for j in range(3)])

def eps3():
    return np.array([1, 1, 1])

def spatial_sign(arr):
    return np.array([arr[0], arr[1], arr[2]])

def all_signs():
    all_sign_funcs0 = []
    all_sign_funcs1 = []
    all_sign_funcs2 = []
    all_sign_funcs3 = []
    all_sign_funcs4 = []
    all_sign_funcs5 = []
    all_sign_funcs6 = []
    all_sign_funcs7 = []
    all_sign_funcs8 = []
    all_sign_funcs9 = []
    all_sign_funcs10 = []
    all_sign_funcs11 = []
    all_sign_funcs12 = []
    all_sign_funcs13 = []
    all_sign_funcs14 = []
    all_sign_funcs15 = []
    all_sign_funcs16 = []
    all_sign_funcs17 = []
    all_sign_funcs18 = []
    all_sign_funcs19 = []

    all_sign_funcs0.append(np.array([0, 0, 0]))

    all_sign_funcs1.append(spatial_sign((eta_i(4) + xi_i(4))%2))

    all_sign_funcs2.append(spatial_sign((eta_i(1) + xi_i(1) + eta_i(5))%2))
    all_sign_funcs2.append(spatial_sign((eta_i(2) + xi_i(2) + eta_i(5))%2))
    all_sign_funcs2.append(spatial_sign((eta_i(3) + xi_i(3) + eta_i(5))%2))

    all_sign_funcs3.append(spatial_sign((eta_i(4) + xi_i(4) + eta_i(1) + xi_i(1) + eta_i(5))%2))
    all_sign_funcs3.append(spatial_sign((eta_i(4) + xi_i(4) + eta_i(2) + xi_i(2) + eta_i(5))%2))
    all_sign_funcs3.append(spatial_sign((eta_i(4) + xi_i(4) + eta_i(3) + xi_i(3) + eta_i(5))%2))

    all_sign_funcs4.append(spatial_sign((eta_i(1))%2))
    all_sign_funcs4.append(spatial_sign((eta_i(2))%2))
    all_sign_funcs4.append(spatial_sign((eta_i(3))%2))

    all_sign_funcs5.append(spatial_sign((eta_i(4) + xi_i(4) + eta_i(1))%2))
    all_sign_funcs5.append(spatial_sign((eta_i(4) + xi_i(4) + eta_i(2))%2))
    all_sign_funcs5.append(spatial_sign((eta_i(4) + xi_i(4) + eta_i(3))%2))

    all_sign_funcs6.append(spatial_sign((xi_i(1) + eta_i(5))%2))
    all_sign_funcs6.append(spatial_sign((xi_i(2) + eta_i(5))%2))
    all_sign_funcs6.append(spatial_sign((xi_i(3) + eta_i(5))%2))

    all_sign_funcs7.append(spatial_sign((eta_i(4) + xi_i(4) + xi_i(1) + eta_i(5))%2))
    all_sign_funcs7.append(spatial_sign((eta_i(4) + xi_i(4) + xi_i(2) + eta_i(5))%2))
    all_sign_funcs7.append(spatial_sign((eta_i(4) + xi_i(4) + xi_i(3) + eta_i(5))%2))

    all_sign_funcs8.append(spatial_sign((eta_i(1) + xi_i(1) + eta_i(2) + eta_i(5))%2))
    all_sign_funcs8.append(spatial_sign((eta_i(2) + xi_i(2) + eta_i(3) + eta_i(5))%2))
    all_sign_funcs8.append(spatial_sign((eta_i(3) + xi_i(3) + eta_i(1) + eta_i(5))%2))

    all_sign_funcs9.append(spatial_sign((eta_i(4) + xi_i(4) + eta_i(1) + xi_i(1) + eta_i(2) + eta_i(5))%2))
    all_sign_funcs9.append(spatial_sign((eta_i(4) + xi_i(4) + eta_i(2) + xi_i(2) + eta_i(3) + eta_i(5))%2))
    all_sign_funcs9.append(spatial_sign((eta_i(4) + xi_i(4) + eta_i(3) + xi_i(3) + eta_i(1) + eta_i(5))%2))

    all_sign_funcs10.append(spatial_sign((eta_i(1) + eta_i(2))%2))
    all_sign_funcs10.append(spatial_sign((eta_i(2) + eta_i(3))%2))
    all_sign_funcs10.append(spatial_sign((eta_i(3) + eta_i(1))%2))

    all_sign_funcs11.append(spatial_sign((eta_i(4) + xi_i(4) + eta_i(1) + eta_i(2))%2))
    all_sign_funcs11.append(spatial_sign((eta_i(4) + xi_i(4) + eta_i(2) + eta_i(3))%2))
    all_sign_funcs11.append(spatial_sign((eta_i(4) + xi_i(4) + eta_i(3) + eta_i(1))%2))

    all_sign_funcs12.append(spatial_sign((xi_i(1) + xi_i(2))%2))
    all_sign_funcs12.append(spatial_sign((xi_i(2) + xi_i(3))%2))
    all_sign_funcs12.append(spatial_sign((xi_i(3) + xi_i(1))%2))

    all_sign_funcs13.append(spatial_sign((eta_i(4) + xi_i(4) + xi_i(1) + xi_i(2))%2))
    all_sign_funcs13.append(spatial_sign((eta_i(4) + xi_i(4) + xi_i(2) + xi_i(3))%2))
    all_sign_funcs13.append(spatial_sign((eta_i(4) + xi_i(4) + xi_i(3) + xi_i(1))%2))

    all_sign_funcs14.append(spatial_sign((eta_i(3) + xi_i(3) + eta_i(1) + xi_i(2))%2))
    all_sign_funcs14.append(spatial_sign((eta_i(1) + xi_i(1) + eta_i(2) + xi_i(3))%2))
    all_sign_funcs14.append(spatial_sign((eta_i(2) + xi_i(2) + eta_i(3) + xi_i(1))%2))

    all_sign_funcs15.append(spatial_sign((eta_i(4) + xi_i(4) + eta_i(3) + xi_i(3) + eta_i(1) + xi_i(2))%2))
    all_sign_funcs15.append(spatial_sign((eta_i(4) + xi_i(4) + eta_i(1) + xi_i(1) + eta_i(2) + xi_i(3))%2))
    all_sign_funcs15.append(spatial_sign((eta_i(4) + xi_i(4) + eta_i(2) + xi_i(2) + eta_i(3) + xi_i(1))%2))

    all_sign_funcs16.append(spatial_sign((eta_i(1) + eta_i(2) + eta_i(3))%2))

    all_sign_funcs17.append(spatial_sign((eta_i(4) + xi_i(4) + eta_i(1) + eta_i(2) + eta_i(3))%2))

    all_sign_funcs18.append(spatial_sign((eta_i(1) + xi_i(1) + eta_i(5) + eta_i(1) + eta_i(2) + eta_i(3))%2))
    all_sign_funcs18.append(spatial_sign((eta_i(2) + xi_i(2) + eta_i(5) + eta_i(1) + eta_i(2) + eta_i(3))%2))
    all_sign_funcs18.append(spatial_sign((eta_i(3) + xi_i(3) + eta_i(5) + eta_i(1) + eta_i(2) + eta_i(3))%2))

    all_sign_funcs19.append(spatial_sign((eta_i(4) + xi_i(4) + eta_i(1) + xi_i(1) + eta_i(5) + eta_i(1) + eta_i(2) + eta_i(3))%2))
    all_sign_funcs19.append(spatial_sign((eta_i(4) + xi_i(4) + eta_i(2) + xi_i(2) + eta_i(5) + eta_i(1) + eta_i(2) + eta_i(3))%2))
    all_sign_funcs19.append(spatial_sign((eta_i(4) + xi_i(4) + eta_i(3) + xi_i(3) + eta_i(5) + eta_i(1) + eta_i(2) + eta_i(3))%2))

    all_signs = []
    all_signs.append(all_sign_funcs0)
    all_signs.append(all_sign_funcs1)
    all_signs.append(all_sign_funcs2)
    all_signs.append(all_sign_funcs3)
    all_signs.append(all_sign_funcs4)
    all_signs.append(all_sign_funcs5)
    all_signs.append(all_sign_funcs6)
    all_signs.append(all_sign_funcs7)
    all_signs.append(all_sign_funcs8)
    all_signs.append(all_sign_funcs9)
    all_signs.append(all_sign_funcs10)
    all_signs.append(all_sign_funcs11)
    all_signs.append(all_sign_funcs12)
    all_signs.append(all_sign_funcs13)
    all_signs.append(all_sign_funcs14)
    all_signs.append(all_sign_funcs15)
    all_signs.append(all_sign_funcs16)
    all_signs.append(all_sign_funcs17)
    all_signs.append(all_sign_funcs18)
    all_signs.append(all_sign_funcs19)
    return all_signs

def all_signs3():
    all_sign_funcs0 = []
    all_sign_funcs1 = []
    all_sign_funcs2 = []
    all_sign_funcs3 = []
    all_sign_funcs4 = []
    all_sign_funcs5 = []
    all_sign_funcs6 = []
    all_sign_funcs7 = []
    all_sign_funcs8 = []
    all_sign_funcs9 = []
    all_sign_funcs10 = []
    all_sign_funcs11 = []
    all_sign_funcs12 = []
    all_sign_funcs13 = []
    all_sign_funcs14 = []
    all_sign_funcs15 = []
    all_sign_funcs16 = []
    all_sign_funcs17 = []
    all_sign_funcs18 = []
    all_sign_funcs19 = []

    all_sign_funcs0.append(np.array([0, 0, 0]))

    all_sign_funcs1.append(spatial_sign(eps3()))

    all_sign_funcs2.append(spatial_sign((eta_i3(1) + xi_i3(1) + eps3())%2))
    all_sign_funcs2.append(spatial_sign((eta_i3(2) + xi_i3(2) + eps3())%2))
    all_sign_funcs2.append(spatial_sign((eta_i3(3) + xi_i3(3) + eps3())%2))

    all_sign_funcs3.append(spatial_sign((eta_i3(1) + xi_i3(1))%2))
    all_sign_funcs3.append(spatial_sign((eta_i3(2) + xi_i3(2))%2))
    all_sign_funcs3.append(spatial_sign((eta_i3(3) + xi_i3(3))%2))

    all_sign_funcs4.append(spatial_sign((eta_i3(1))%2))
    all_sign_funcs4.append(spatial_sign((eta_i3(2))%2))
    all_sign_funcs4.append(spatial_sign((eta_i3(3))%2))

    all_sign_funcs5.append(spatial_sign((eta_i3(1) + eps3())%2))
    all_sign_funcs5.append(spatial_sign((eta_i3(2) + eps3())%2))
    all_sign_funcs5.append(spatial_sign((eta_i3(3) + eps3())%2))

    all_sign_funcs6.append(spatial_sign((xi_i3(1) + eps3())%2))
    all_sign_funcs6.append(spatial_sign((xi_i3(2) + eps3())%2))
    all_sign_funcs6.append(spatial_sign((xi_i3(3) + eps3())%2))

    all_sign_funcs7.append(spatial_sign((xi_i3(1))%2))
    all_sign_funcs7.append(spatial_sign((xi_i3(2))%2))
    all_sign_funcs7.append(spatial_sign((xi_i3(3))%2))

    all_sign_funcs8.append(spatial_sign((eta_i3(1) + xi_i3(1) + eta_i3(2) + eps3())%2))
    all_sign_funcs8.append(spatial_sign((eta_i3(2) + xi_i3(2) + eta_i3(3) + eps3())%2))
    all_sign_funcs8.append(spatial_sign((eta_i3(3) + xi_i3(3) + eta_i3(1) + eps3())%2))

    all_sign_funcs9.append(spatial_sign((eta_i3(1) + xi_i3(1) + eta_i3(2))%2))
    all_sign_funcs9.append(spatial_sign((eta_i3(2) + xi_i3(2) + eta_i3(3))%2))
    all_sign_funcs9.append(spatial_sign((eta_i3(3) + xi_i3(3) + eta_i3(1))%2))

    all_sign_funcs10.append(spatial_sign((eta_i3(1) + eta_i3(2))%2))
    all_sign_funcs10.append(spatial_sign((eta_i3(2) + eta_i3(3))%2))
    all_sign_funcs10.append(spatial_sign((eta_i3(3) + eta_i3(1))%2))

    all_sign_funcs11.append(spatial_sign((eta_i3(1) + eta_i3(2) + eps3())%2))
    all_sign_funcs11.append(spatial_sign((eta_i3(2) + eta_i3(3) + eps3())%2))
    all_sign_funcs11.append(spatial_sign((eta_i3(3) + eta_i3(1) + eps3())%2))

    all_sign_funcs12.append(spatial_sign((xi_i3(1) + xi_i3(2))%2))
    all_sign_funcs12.append(spatial_sign((xi_i3(2) + xi_i3(3))%2))
    all_sign_funcs12.append(spatial_sign((xi_i3(3) + xi_i3(1))%2))

    all_sign_funcs13.append(spatial_sign((xi_i3(1) + xi_i3(2) + eps3())%2))
    all_sign_funcs13.append(spatial_sign((xi_i3(2) + xi_i3(3) + eps3())%2))
    all_sign_funcs13.append(spatial_sign((xi_i3(3) + xi_i3(1) + eps3())%2))

    all_sign_funcs14.append(spatial_sign((eta_i3(3) + xi_i3(3) + eta_i3(1) + xi_i3(2))%2))
    all_sign_funcs14.append(spatial_sign((eta_i3(1) + xi_i3(1) + eta_i3(2) + xi_i3(3))%2))
    all_sign_funcs14.append(spatial_sign((eta_i3(2) + xi_i3(2) + eta_i3(3) + xi_i3(1))%2))

    all_sign_funcs15.append(spatial_sign((eta_i3(3) + xi_i3(3) + eta_i3(1) + xi_i3(2) + eps3())%2))
    all_sign_funcs15.append(spatial_sign((eta_i3(1) + xi_i3(1) + eta_i3(2) + xi_i3(3) + eps3())%2))
    all_sign_funcs15.append(spatial_sign((eta_i3(2) + xi_i3(2) + eta_i3(3) + xi_i3(1) + eps3())%2))

    all_sign_funcs16.append(spatial_sign((eta_i3(1) + eta_i3(2) + eta_i3(3))%2))

    all_sign_funcs17.append(spatial_sign((eta_i3(1) + eta_i3(2) + eta_i3(3) + eps3())%2))

    all_sign_funcs18.append(spatial_sign((eta_i3(1) + xi_i3(1) + eta_i3(1) + eta_i3(2) + eta_i3(3) + eps3())%2))
    all_sign_funcs18.append(spatial_sign((eta_i3(2) + xi_i3(2) + eta_i3(1) + eta_i3(2) + eta_i3(3) + eps3())%2))
    all_sign_funcs18.append(spatial_sign((eta_i3(3) + xi_i3(3) + eta_i3(1) + eta_i3(2) + eta_i3(3) + eps3())%2))

    all_sign_funcs19.append(spatial_sign((eta_i3(1) + xi_i3(1) + eta_i3(1) + eta_i3(2) + eta_i3(3))%2))
    all_sign_funcs19.append(spatial_sign((eta_i3(2) + xi_i3(2) + eta_i3(1) + eta_i3(2) + eta_i3(3))%2))
    all_sign_funcs19.append(spatial_sign((eta_i3(3) + xi_i3(3) + eta_i3(1) + eta_i3(2) + eta_i3(3))%2))

    all_signs = []
    all_signs.append(all_sign_funcs0)
    all_signs.append(all_sign_funcs1)
    all_signs.append(all_sign_funcs2)
    all_signs.append(all_sign_funcs3)
    all_signs.append(all_sign_funcs4)
    all_signs.append(all_sign_funcs5)
    all_signs.append(all_sign_funcs6)
    all_signs.append(all_sign_funcs7)
    all_signs.append(all_sign_funcs8)
    all_signs.append(all_sign_funcs9)
    all_signs.append(all_sign_funcs10)
    all_signs.append(all_sign_funcs11)
    all_signs.append(all_sign_funcs12)
    all_signs.append(all_sign_funcs13)
    all_signs.append(all_sign_funcs14)
    all_signs.append(all_sign_funcs15)
    all_signs.append(all_sign_funcs16)
    all_signs.append(all_sign_funcs17)
    all_signs.append(all_sign_funcs18)
    all_signs.append(all_sign_funcs19)
    return all_signs

def third_sign():
    return np.array([ 1, 1, 1, 1, 1,
                      1, 1, 1, 1, 1,
                     -1,-1,-1,-1,-1,
                     -1, 1, 1, 1, 1])

def deltaidx_to_delta(delta: int):
    return np.array([(delta&1), ((delta>>1)&1), ((delta>>2)&1)])

def delta_to_deltaidx(delta):
    return (delta[2]<<2)|(delta[1]<<1)|delta[0]

def all_deltas():
    all_delta0 = []
    all_delta1 = []
    all_delta2 = []
    all_delta3 = []
    all_delta4 = []
    all_delta5 = []
    all_delta6 = []
    all_delta7 = []
    all_delta8 = []
    all_delta9 = []
    all_delta10 = []
    all_delta11 = []
    all_delta12 = []
    all_delta13 = []
    all_delta14 = []
    all_delta15 = []
    all_delta16 = []
    all_delta17 = []
    all_delta18 = []
    all_delta19 = []

    all_delta0.append(np.array([0, 0, 0]))

    all_delta1.append(np.array([0, 0, 0]))

    all_delta2.append(np.array([0, 0, 0]))
    all_delta2.append(np.array([0, 0, 0]))
    all_delta2.append(np.array([0, 0, 0]))

    all_delta3.append(np.array([0, 0, 0]))
    all_delta3.append(np.array([0, 0, 0]))
    all_delta3.append(np.array([0, 0, 0]))

    all_delta4.append(np.array([1, 0, 0]))
    all_delta4.append(np.array([0, 1, 0]))
    all_delta4.append(np.array([0, 0, 1]))

    all_delta5.append(np.array([1, 0, 0]))
    all_delta5.append(np.array([0, 1, 0]))
    all_delta5.append(np.array([0, 0, 1]))

    all_delta6.append(np.array([1, 0, 0]))
    all_delta6.append(np.array([0, 1, 0]))
    all_delta6.append(np.array([0, 0, 1]))

    all_delta7.append(np.array([1, 0, 0]))
    all_delta7.append(np.array([0, 1, 0]))
    all_delta7.append(np.array([0, 0, 1]))

    all_delta8.append(np.array([0, 1, 0]))
    all_delta8.append(np.array([0, 0, 1]))
    all_delta8.append(np.array([1, 0, 0]))

    all_delta9.append(np.array([0, 1, 0]))
    all_delta9.append(np.array([0, 0, 1]))
    all_delta9.append(np.array([1, 0, 0]))

    all_delta10.append(np.array([1, 1, 0]))
    all_delta10.append(np.array([0, 1, 1]))
    all_delta10.append(np.array([1, 0, 1]))

    all_delta11.append(np.array([1, 1, 0]))
    all_delta11.append(np.array([0, 1, 1]))
    all_delta11.append(np.array([1, 0, 1]))

    all_delta12.append(np.array([1, 1, 0]))
    all_delta12.append(np.array([0, 1, 1]))
    all_delta12.append(np.array([1, 0, 1]))

    all_delta13.append(np.array([1, 1, 0]))
    all_delta13.append(np.array([0, 1, 1]))
    all_delta13.append(np.array([1, 0, 1]))

    all_delta14.append(np.array([1, 1, 0]))
    all_delta14.append(np.array([0, 1, 1]))
    all_delta14.append(np.array([1, 0, 1]))

    all_delta15.append(np.array([1, 1, 0]))
    all_delta15.append(np.array([0, 1, 1]))
    all_delta15.append(np.array([1, 0, 1]))

    all_delta16.append(np.array([1, 1, 1]))

    all_delta17.append(np.array([1, 1, 1]))

    all_delta18.append(np.array([1, 1, 1]))
    all_delta18.append(np.array([1, 1, 1]))
    all_delta18.append(np.array([1, 1, 1]))

    all_delta19.append(np.array([1, 1, 1]))
    all_delta19.append(np.array([1, 1, 1]))
    all_delta19.append(np.array([1, 1, 1]))

    all_deltas = []
    all_deltas.append(all_delta0)
    all_deltas.append(all_delta1)
    all_deltas.append(all_delta2)
    all_deltas.append(all_delta3)
    all_deltas.append(all_delta4)
    all_deltas.append(all_delta5)
    all_deltas.append(all_delta6)
    all_deltas.append(all_delta7)
    all_deltas.append(all_delta8)
    all_deltas.append(all_delta9)
    all_deltas.append(all_delta10)
    all_deltas.append(all_delta11)
    all_deltas.append(all_delta12)
    all_deltas.append(all_delta13)
    all_deltas.append(all_delta14)
    all_deltas.append(all_delta15)
    all_deltas.append(all_delta16)
    all_deltas.append(all_delta17)
    all_deltas.append(all_delta18)
    all_deltas.append(all_delta19)
    return all_deltas

