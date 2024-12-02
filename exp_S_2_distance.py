from dem_decompose_corr import *


# logic error rate - distance scaling 
# S-2 circuits with n_pad = d//2 + 1 
# n_m is currently chosen to be d+1
# plain decoding is used (can be switched to FR decoding for 
# better logical performance) 




noise_list = [0.001]
distance_list = [3,5,7,9]
num_batch_list = [10,20*2,20*5,60*10]
noise = 0.001



# noise_list = [0.001]
# distance_list = [3,5,7,9,11,13]
# num_batch_list = [10,20*2,20*5,60*10,3000,200*60]
# noise = 0.001



for distance, num_batch in zip(distance_list,
                               num_batch_list):

    num_pad = (distance//2)+1
    rotated = Rotated_surface_code(distance)

    ft_circuit = FT_Circuit_dynamical_phase(rotated,
                                            'S_2_mid_full.stim',noise)

    # X-basis initialization followed by an I-SE round
    ft_circuit.initialize_in_x()

    for i in range(num_pad):
        ft_circuit.stab_measurement_functional()


    # contains a single S-SE round and a single I-SE round
    ft_circuit.phase_gate_dynamical()


    # (number-1) of rounds between two S-SE rounds 
    for i in range((distance)):
        ft_circuit.stab_measurement_functional()


    ft_circuit.phase_gate_dynamical()

    for i in range(num_pad):
        ft_circuit.stab_measurement_functional()


    ft_circuit.logic_x_measurement()



    stim_circuit = ft_circuit.get_stim_circuit()
    dem = DEM(stim_circuit)

    decomposer = Decomposer(ft_circuit, dem)
    decomposer._hpg_initialize()
    decomposer._fund_det_x_errs()
    decomposer._fund_det_z_errs()
    decomposer._det_z_errs_weight_cate()
    decomposer._check_z_foot_print()
    decomposer._attach_hyperedges()
    decomposer.z_matcher_compilation()
    decomposer.x_matcher_compilation()
    decomposer.reweight_helper()

    num_shots = 50_000
    if distance >= 15:
        num_shots = 30_000
    err_rate_list = []
    err_acc = 0
    max_err = 800   # maximal number of hits

    for i in range(num_batch):
        start = timer()
        res = decomposer.decoding_batch_plain(num_shots)
        end = timer()
        print('sparse: ',res[0], '     time: ',end-start)
        num_logical_errors = res[1]
        err_acc += num_logical_errors
        err_rate_list.append(num_logical_errors/num_shots)
        with open('S_2_mid_full.txt','a') as f:
            f.write('surface_S: ')
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

        if err_acc > max_err:
            break