import numpy as np
from qecsim import paulitools as pt
from qecsim.models.generic import BitFlipErrorModel
from qecsim.models.toric import ToricCode, ToricMWPMDecoder
import math
import itertools
import more_itertools as mit
import csv

# set size of lattice
L = 9
num_qubits = 2*(L**2)
num_iter = 10**7

# initialise models
my_code = ToricCode(L,L)
my_error_model = BitFlipErrorModel()
my_decoder = ToricMWPMDecoder()

# set physical error probability to 2%
error_probability = 0.02

# init variables
syndrome_count = dict()
error_correction_count = dict()
syndrome_error_correction_count = dict()
error_correction_count[0] = 0
error_correction_count[1] = 0

# seed random number generator for repeatability (leave blank for random)
rng = np.random.default_rng()

for values in range(num_iter):
    error = my_error_model.generate(my_code, error_probability, rng)

    # syndrome: stabilizers that do not commute with the error
    syndrome = pt.bsp(error, my_code.stabilizers.T)

    recovery = my_decoder.decode(my_code, syndrome)
    is_correct = int(sum(pt.bsp(recovery ^ error, my_code.logicals.T)) == 0)

    # convert the syndrome binary-valued array into equivalent integer
    syndrome = "".join(str(x) for x in syndrome)

    # update syndrome count
    if syndrome in syndrome_count:
        syndrome_count[syndrome] = syndrome_count[syndrome] + 1
    else:
        syndrome_count[syndrome] = 1

    # update EC count
    if is_correct == 0:
        error_correction_count[0] = error_correction_count[0] + 1
        syndrome_ec = syndrome + "0"
    else:
        error_correction_count[1] = error_correction_count[1] + 1
        syndrome_ec = syndrome + "1"

    # update joint S, EC count
    if syndrome_ec in syndrome_error_correction_count:
        syndrome_error_correction_count[syndrome_ec] = syndrome_error_correction_count[syndrome_ec] + 1
    else:
        syndrome_error_correction_count[syndrome_ec] = 1
else:
    # now write to csv
    header = ["Syndrome", "S-Count"]
    with open('MI-L9-S.csv', 'w', encoding='UTF8', newline='') as f:
        writer = csv.writer(f)

        # write the header
        writer.writerow(header)
        # write the data
        writer.writerows(zip(list(syndrome_count.keys()), list(syndrome_count.values())))
    
    header = ["EC", "EC-Count"]
    with open('MI-L9-EC.csv', 'w', encoding='UTF8', newline='') as f:
        writer = csv.writer(f)

        # write the header
        writer.writerow(header)
        # write the data
        writer.writerows(zip(list(error_correction_count.keys()), list(error_correction_count.values())))

    header = ["S_EC", "S_EC-Count"]
    with open('MI-L9-S_EC.csv', 'w', encoding='UTF8', newline='') as f:
        writer = csv.writer(f)

        # write the header
        writer.writerow(header)
        # write the data
        writer.writerows(zip(list(syndrome_error_correction_count.keys()), list(syndrome_error_correction_count.values())))    
    