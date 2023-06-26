import numpy as np
from qecsim import paulitools as pt
from qecsim.models.generic import BitFlipErrorModel
from qecsim.models.toric import ToricCode, ToricMWPMDecoder
import math

# set number of iterations and torus size
n = 10**4
L = 2;

# initialise models
my_code = ToricCode(L,L)
my_error_model = BitFlipErrorModel()
my_decoder = ToricMWPMDecoder()

# set physical error probability to some percentage% (note threshold for toric code is around 11%)
error_probability = 0.1

# init variables, we hold counts in a dictionary
error_counts = dict()
syndrome_counts = dict()
succ_counts = dict()
syndrome_succ_counts = dict()
succ_counts["success"] = 0
succ_counts["fail"] = 0

for i in range(n):

    # seed random number generator for repeatability if necessary
    rng = np.random.default_rng()

    # error: random error based on error probability
    error = my_error_model.generate(my_code, error_probability, rng)

    # syndrome: stabilizers that do not commute with the error
    syndrome = pt.bsp(error, my_code.stabilizers.T)

    recovery = my_decoder.decode(my_code, syndrome)
    is_correct = sum(pt.bsp(recovery ^ error, my_code.logicals.T))

    error = str(error)
    syndrome = str(syndrome)

    if error in error_counts:
        error_counts[error] = error_counts[error]+1
    else:
        error_counts[error] = 1

    if syndrome in syndrome_counts:
        syndrome_counts[syndrome] = syndrome_counts[syndrome]+1
    else:
        syndrome_counts[syndrome] = 1

    if is_correct == 0:
        succ_counts["success"] = succ_counts["success"] + 1
        syndrome_succ = syndrome + "1"
    else:
       succ_counts["fail"] = succ_counts["fail"] + 1
       syndrome_succ = syndrome + "0"

    if syndrome_succ in syndrome_succ_counts:
        syndrome_succ_counts[syndrome_succ] = syndrome_succ_counts[syndrome_succ] + 1
    else:
        syndrome_succ_counts[syndrome_succ] = 0

else:
    He = 0
    for error in error_counts:
        p = error_counts[error]/n
        He = He - p*math.log2(p)

    Hs = 0
    for syndrome in syndrome_counts:
        p = syndrome_counts[syndrome]/n
        Hs = Hs - p*math.log2(p)
    
    Hec = 0
    for succ in succ_counts:
        p = succ_counts[succ]/n
        Hec = Hec - p*math.log2(p)

    Hsec = 0
    for syndrome_succ in syndrome_succ_counts:
        p = syndrome_succ_counts[syndrome_succ]/n
        Hsec = Hsec - p*math.log2(p)

    print("L = " + str(L))
    print("n = " + str(n))
    print("H(E) = " + str(He))
    print("H(S) = " + str(Hs))
    print("H(EC) = " + str(Hec))
    print("H(S,EC) = " + str(Hsec))
    print("I(S;EC) = " + str(Hs + Hec - Hsec))