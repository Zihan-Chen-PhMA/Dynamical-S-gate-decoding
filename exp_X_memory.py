from dem_parsor import *
from hypergraph import *
import pymatching



# X-memory circuit experiments



# The following function is a code snippet directly from 
#           Stim/doc/getting_started.ipynb
# and is used to decode X memory circuits. 

def count_logical_errors(circuit: stim.Circuit, num_shots: int) -> int:
    # Sample the circuit.
    sampler = circuit.compile_detector_sampler()
    detection_events, observable_flips = sampler.sample(num_shots, separate_observables=True)

    # Configure a decoder using the circuit.
    detector_error_model = circuit.detector_error_model(decompose_errors=True)
    matcher = pymatching.Matching.from_detector_error_model(detector_error_model)

    # Run the decoder.
    predictions = matcher.decode_batch(detection_events)

    # Count the mistakes.
    num_errors = 0
    for shot in range(num_shots):
        actual_for_shot = observable_flips[shot]
        predicted_for_shot = predictions[shot]
        if not np.array_equal(actual_for_shot, predicted_for_shot):
            num_errors += 1
    return num_errors





distance_list = [3,5,7,9]
num_batch_list = [10,20*3,20*10,60*15]

noise = 0.001


for distance, num_batch in zip(distance_list,
                               num_batch_list):

    num_pad = (distance//2)+1
    rotated = Rotated_surface_code(distance)

    memory_circ = FT_Circuit_dynamical_phase(rotated,
                                            'X_memory.stim',noise)

    memory_circ.initialize_in_x()

    for i in range(num_pad):
        memory_circ.stab_measurement_functional()




    for i in range((distance+4)):
        memory_circ.stab_measurement_functional()


    for i in range(num_pad):
        memory_circ.stab_measurement_functional()


    memory_circ.logic_x_measurement()



    stim_circuit = memory_circ.get_stim_circuit()


    num_shots = 50_000
    err_rate_list = []
    err_acc = 0
    err_max = 800

    for i in range(num_batch):
        num_logical_errors = count_logical_errors(stim_circuit, num_shots)
        err_rate_list.append(num_logical_errors/num_shots)
        with open('X_memory.txt','a') as f:
            f.write('surface_X: ')
            f.write('distance: ')
            f.write(str(distance))
            f.write(' ')
            f.write('noise_prob: ')
            f.write(str(noise))
            f.write(' ')
            f.write('num_shots: ')
            f.write(str(num_shots))
            f.write(' ')
            f.write('num_hits: ')
            f.write(str(num_logical_errors))
            f.write('\n')
        print('progress: ',i+1,'/',num_batch)
        err_acc += num_logical_errors
        if err_acc>err_max:
            break
    print("error rate is", np.mean(np.array(err_rate_list)))





