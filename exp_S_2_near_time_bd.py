from dem_decompose_corr import *


# S-SE rounds near time boundaries
# S-2 circuit with n_pad = 2 and n_m = d+1
# default decoding: VTB (can be changed to other VTB-based 
# decoding methods)





noise_list = [0.001]
distance_list = [3,5,7,9]
num_batch_list = [10,20*2,20*5,60*10]

noise = 0.001

for distance, num_batch in zip(distance_list,
                               num_batch_list):

    num_pad = 1
    rotated = Rotated_surface_code(distance)

    ft_circuit = FT_Circuit_dynamical_phase(rotated,
                                            'S_2_timebd.stim',noise)

    ft_circuit.initialize_in_x()

    for i in range(num_pad):
        ft_circuit.stab_measurement_functional()


    ft_circuit.phase_gate_dynamical()


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

    # VTB-based methods require special z matcher compilation step
    # because we need to add virtual vertices to the decoding graph.
    decomposer.z_matcher_compilation_VTB()  
    decomposer.x_matcher_compilation()
    decomposer.reweight_helper()

        

    num_shots = 50_000
    if distance >= 15:
        num_shots = 30_000
    err_rate_list = []
    err_acc = 0
    max_err = 800

    for i in range(num_batch):
        start = timer()
        # VTB decoding:
        res = decomposer.decoding_batch_VTB(num_shots)
        # VTB-PR decoding:
        # res = decomposer.decoding_batch_VTB_PR(num_shots)  
        # VTB-FR decoding:
        # res = decomposer.decoding_batch_VTB_FR(num_shots)
        end = timer()
        print('sparse: ',res[0], '     time: ',end-start)
        num_logical_errors = res[1]
        err_acc += num_logical_errors
        err_rate_list.append(num_logical_errors/num_shots)
        with open('S_2_timebd_VTB.txt','a') as f:
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