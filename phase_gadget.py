from surface_codes import *
from typing import Callable
import copy
import pymatching

      
# Build S-2 or X-memory circuits from I-SE or S-SE rounds. 

class FT_Circuit_dynamical_phase():

    def __init__(self, code: Rotated_surface_code, 
                 circuit_filename: str,
                 noise_probability: float = 0, 
                 control_functional: Callable = None) -> None:
        self.code = code
        self.circuit_filename = circuit_filename
        self.flag_displacement = self.code.distance*5
        self.flag_switch = True
        self.circ = Circuit_helper(self.code.data_qubits_collection, 
                            self.code.x_check_collection,
                            self.code.z_check_collection,
                            circuit_filename)
        self.circ.initialize()
        self.z_check_ctrl_dict: dict[int,list[Qubit]] = \
                                        self._z_check_ctrl_dict()
        self.x_check_targ_dict: dict[int,list[Qubit]] = \
                                        self._x_check_targ_dict()
        self.outer_clock = -1
        self.noise_tag = ['RESET', 'MEASUREMENT', 'PREROUND_DEPOLAR']
        self.noise_probability = noise_probability
        if control_functional == None:
            self.control_functional = self._control_functional()
        else:
            self.control_functional = control_functional
        pass


    def _control_functional(self) -> Callable:

        def control() -> list[
                              dict[int,Union[list[Qubit],None]],
                              dict[int,Union[list[Qubit],None]],
                              dict[int,Union[list[Qubit],None]],
                              Union[None,list[Qubit]]]:
            ctrl_qubit_dict: dict[int,list[Qubit]] = {}
            targ_qubit_dict: dict[int,list[Qubit]] = {}
            pre_transform_dict: dict[int,list[Qubit]] = {}
            final_transform_list: list[Qubit] = []
            for i in [1,2,3,4]:
                ctrl_qubit_dict[i] = self.z_check_ctrl_dict[i] \
                                        + self.code.x_check_collection
                targ_qubit_dict[i] = self.code.z_check_collection \
                                        + self.x_check_targ_dict[i] 
                    
            pre_transform_dict[1] = self.code.x_check_collection 
            pre_transform_dict[2] = None
            pre_transform_dict[3] = None
            pre_transform_dict[4] = None
            final_transform_list = self.code.x_check_collection


            return [ctrl_qubit_dict, targ_qubit_dict, pre_transform_dict,
                    final_transform_list]
        return control



                

        

    def _z_check_ctrl_dict(self) -> dict[str,list[Qubit]]:
        z_check_ctrl_dict = {}
        z_check_ctrl_dict[1] = [qubit.neighbor['nw'] 
                                  for qubit in self.code.z_check_collection]
        z_check_ctrl_dict[2] = [qubit.neighbor['ne']
                                  for qubit in self.code.z_check_collection]
        z_check_ctrl_dict[3] = [qubit.neighbor['sw'] 
                                  for qubit in self.code.z_check_collection]
        z_check_ctrl_dict[4] = [qubit.neighbor['se'] 
                                  for qubit in self.code.z_check_collection]
        return z_check_ctrl_dict
    

    def _x_check_targ_dict(self) -> dict[str,list[Qubit]]:
        x_check_targ_dict = {}
        x_check_targ_dict[1] = [qubit.neighbor['nw'] 
                                  for qubit in self.code.x_check_collection]
        x_check_targ_dict[2] = [qubit.neighbor['sw'] 
                                for qubit in self.code.x_check_collection]
        x_check_targ_dict[3] = [qubit.neighbor['ne'] 
                                for qubit in self.code.x_check_collection]
        x_check_targ_dict[4] = [qubit.neighbor['se'] 
                                for qubit in self.code.x_check_collection]
        return x_check_targ_dict
        
    
    def initialize_in_z(self) -> None:
        # Initialization in logical |0> and first round of stab measurement.
        # Z-Detector is set using the first round of Z-stab measurements.
        if 'RESET' in self.noise_tag:
            reset_error = ['X_ERROR', self.noise_probability]
        else: 
            reset_error = None
        if 'MEASUREMENT' in self.noise_tag:
            measurement_error = ['X_ERROR', self.noise_probability]
        else:
            measurement_error = None
        self.circ.reset(self.code.data_qubits_collection
                        +self.code.x_check_collection
                        +self.code.z_check_collection, 
                        noise = reset_error)
        self.circ.single_qubit_gate('H', self.code.x_check_collection,
                                    noise=['DEPOLARIZE1',self.noise_probability])
        self.circ.tick()
        self.circ.tick()
        for i in [1,2,3,4]:
            self.circ.two_qubit_gate('CX', self.z_check_ctrl_dict[i]
                                           + self.code.x_check_collection,
                                           self.code.z_check_collection
                                           + self.x_check_targ_dict[i],
                                        noise=['DEPOLARIZE2',self.noise_probability])
            active_qubit_list = self.z_check_ctrl_dict[i] \
                                + self.code.x_check_collection \
                                + self.code.z_check_collection \
                                + self.x_check_targ_dict[i]
            
            idle_qubit_list = []
            for qubit in self.circ.qubit_collection:
                if qubit not in active_qubit_list:
                    idle_qubit_list.append(qubit)

            self.circ.single_qubit_noise(idle_qubit_list,
                                         noise=['DEPOLARIZE1', self.noise_probability])

            self.circ.tick()
        
        self.circ.single_qubit_gate('H', self.code.x_check_collection,
                                    noise=['DEPOLARIZE1',self.noise_probability])
        self.circ.tick()
        self.circ.measure_reset(self.code.x_check_collection 
                                + self.code.z_check_collection, 
                                noise = measurement_error)
        
        for qubit in self.code.z_check_collection:
            self.circ.detector([qubit], self.circ.time_temp, [], None, 
                        (qubit.pos[0], qubit.pos[1],self.circ.time_temp))


        self.circ.clock_plus_one()
        self.circ.tick()
        pass
        pass

    def initialize_in_x(self) -> None:
        # Initialization in logical |+> and first round of stab measurement.
        # X-Detector is set using the first round of X-stab measurements.
        if 'RESET' in self.noise_tag:
            reset_error = ['X_ERROR', self.noise_probability]
        else: 
            reset_error = None
        if 'MEASUREMENT' in self.noise_tag:
            measurement_error = ['X_ERROR', self.noise_probability]
        else:
            measurement_error = None
        self.circ.reset(self.code.data_qubits_collection
                        +self.code.x_check_collection
                        +self.code.z_check_collection, 
                        noise = reset_error)
        self.circ.single_qubit_gate('H', self.code.x_check_collection
                                         + self.code.data_qubits_collection,
                                         noise=['DEPOLARIZE1',self.noise_probability])
        self.circ.tick()
        for i in [1,2,3,4]:
            self.circ.two_qubit_gate('CX', self.z_check_ctrl_dict[i]
                                           + self.code.x_check_collection,
                                           self.code.z_check_collection
                                           + self.x_check_targ_dict[i],
                                           noise=['DEPOLARIZE2',self.noise_probability])
            active_qubit_list = self.z_check_ctrl_dict[i] \
                                + self.code.x_check_collection \
                                + self.code.z_check_collection \
                                + self.x_check_targ_dict[i]
            idle_qubit_list = []
            for qubit in self.circ.qubit_collection:
                if qubit not in active_qubit_list:
                    idle_qubit_list.append(qubit)

            self.circ.single_qubit_noise(idle_qubit_list,
                                         noise=['DEPOLARIZE1', self.noise_probability])
            
            self.circ.tick()
        
        self.circ.single_qubit_gate('H', self.code.x_check_collection,
                                    noise=['DEPOLARIZE1',self.noise_probability])
        self.circ.tick()
        self.circ.measure_reset(self.code.x_check_collection 
                                + self.code.z_check_collection, 
                                noise = measurement_error)
        for qubit in self.code.x_check_collection:
            self.circ.detector([qubit], self.circ.time_temp, [], None, 
                        (qubit.pos[0], qubit.pos[1],self.circ.time_temp))


        self.circ.clock_plus_one()
        self.circ.tick()
        pass





    def stab_measurement_functional(self) -> None:
        if 'PREROUND_DEPOLAR' in self.noise_tag:
            preround_error = ['DEPOLARIZE1', self.noise_probability]
        else:
            preround_error = None
        if 'MEASUREMENT' in self.noise_tag:
            measurement_error = ['X_ERROR', self.noise_probability]
        else:
            measurement_error = None
        control_signals = self.control_functional()
        ctrl_qubit_dict: dict[int,Union[list[Qubit],None]] = control_signals[0]
        targ_qubit_dict: dict[int,Union[list[Qubit],None]] = control_signals[1]
        pre_transform_dict: dict[int,Union[list[Qubit],None]] = control_signals[2]
        final_transform_list: Union[None,list[Qubit]] = control_signals[3]

        self.circ.single_qubit_gate('I', self.code.data_qubits_collection,
                                    noise = preround_error)
        for i in [1,2,3,4]:
            self.circ.single_qubit_gate('H', pre_transform_dict[i],
                                        ['DEPOLARIZE1',self.noise_probability])
            self.circ.tick()
            self.circ.two_qubit_gate('CX', ctrl_qubit_dict[i],
                                           targ_qubit_dict[i],
                                           ['DEPOLARIZE2', self.noise_probability])
            active_qubit_list = ctrl_qubit_dict[i] + targ_qubit_dict[i]

            idle_qubit_list = []
            for qubit in self.circ.qubit_collection:
                if qubit not in active_qubit_list:
                    idle_qubit_list.append(qubit)

            self.circ.single_qubit_noise(idle_qubit_list,
                                         noise=['DEPOLARIZE1', self.noise_probability])

            self.circ.tick()
        
        self.circ.single_qubit_gate('H',final_transform_list,
                                    noise=['DEPOLARIZE1',self.noise_probability])
        self.circ.tick()
        self.circ.measure_reset(self.code.x_check_collection
                                + self.code.z_check_collection,
                                noise = measurement_error)
        for check in self.code.x_check_collection + self.code.z_check_collection:
            self.circ.detector([check], self.circ.time_temp - 1,
                               [check], self.circ.time_temp,
                        (check.pos[0],check.pos[1], self.circ.time_temp))

        self.circ.clock_plus_one()
        self.circ.tick()
        pass


    

    def logic_z_measurement(self) -> None:
        self.circ.measure(self.code.data_qubits_collection,
                          noise = ['X_ERROR', self.noise_probability])
        for qubit in self.code.z_check_collection:
            self.circ.detector([qubit], self.circ.time_temp - 1,
                        [data_qubit for data_qubit in qubit.neighbor.values()],
                        self.circ.time_temp,
                        (qubit.pos[0],qubit.pos[1], self.circ.time_temp))
        self.circ.observable(self.code.logic_z_collection,self.circ.time_temp)

        pass

    def logic_x_measurement(self) -> None:
        self.circ.single_qubit_gate('H',self.code.data_qubits_collection,
                                    noise= ['DEPOLARIZE1',self.noise_probability])
        self.circ.measure(self.code.data_qubits_collection, 
                          noise = ['X_ERROR', self.noise_probability])
        for qubit in self.code.x_check_collection:
            self.circ.detector([qubit], self.circ.time_temp - 1,
                        [data_qubit for data_qubit in qubit.neighbor.values()],
                        self.circ.time_temp,
                        (qubit.pos[0],qubit.pos[1], self.circ.time_temp))
        self.circ.observable(self.code.logic_x_collection,self.circ.time_temp)

        pass







    def phase_gate_dynamical(self) -> None:
        control_signals = self.control_functional()
        ctrl_qubit_dict: dict[int,Union[list[Qubit],None]] = control_signals[0]
        targ_qubit_dict: dict[int,Union[list[Qubit],None]] = control_signals[1]
        pre_transform_dict: dict[int,Union[list[Qubit],None]] = control_signals[2]
        final_transform_list: Union[None,list[Qubit]] = control_signals[3]
        # ctrl_qubit_dict, \
        # targ_qubit_dict, \
        # pre_transform_dict, \
        # final_transform_list = self.control_functional()
        self.circ.single_qubit_gate('I', self.code.data_qubits_collection,
                                    noise = ['DEPOLARIZE1', self.noise_probability])
        
        for i in [1,2]:
            self.circ.single_qubit_gate('H', pre_transform_dict[i],
                                        noise=['DEPOLARIZE1',self.noise_probability])
            self.circ.tick()
            self.circ.two_qubit_gate('CX', ctrl_qubit_dict[i],
                                           targ_qubit_dict[i],
                                    ['DEPOLARIZE2', self.noise_probability])
            active_qubit_list = ctrl_qubit_dict[i] + targ_qubit_dict[i]

            idle_qubit_list = []
            for qubit in self.circ.qubit_collection:
                if qubit not in active_qubit_list:
                    idle_qubit_list.append(qubit)

            self.circ.single_qubit_noise(idle_qubit_list,
                                         noise=['DEPOLARIZE1', self.noise_probability])

            self.circ.tick()

        mid_line = []
        mid_line_even = []
        mid_line_odd = []
        upper_triangle = []
        lower_triangle = []

        for qubit in self.circ.qubit_collection:
            qubit_pos = qubit.pos
            # no boundary qubit check
            if qubit_pos[0] < 0 or qubit_pos[1] < 0:
                continue
            if qubit_pos[0] > (self.code.distance - 1) or \
               qubit_pos[1] > (self.code.distance - 1):
                continue
            # check for mid line qubits:
            if qubit_pos[0] == self.code.distance - 1 - qubit_pos[1]:
                if qubit_pos[0] - (qubit_pos[0] // 1) < 0.0001:
                    mid_line_even.append(qubit)
                else:
                    mid_line_odd.append(qubit)
                mid_line.append(qubit)
            # check for upper triangle
            elif qubit_pos[0] < self.code.distance - 1 - qubit_pos[1]:
                upper_triangle.append(qubit)
                lower_pos = (self.code.distance-1-qubit_pos[1],
                             self.code.distance-1-qubit_pos[0])
                lower_qubit = self.circ.pos_qubit_dict[lower_pos]
                lower_triangle.append(lower_qubit)

        self.circ.two_qubit_gate('CZ', upper_triangle,
                                 lower_triangle,
                                 ['DEPOLARIZE2', self.noise_probability])
        self.circ.single_qubit_gate('S',mid_line_even,
                                    ['DEPOLARIZE1', self.noise_probability])
        self.circ.single_qubit_gate('S_DAG',mid_line_odd,
                                    ['DEPOLARIZE1', self.noise_probability])
        
        active_qubit_list = upper_triangle + lower_triangle \
                            + mid_line_even + mid_line_odd
        
        idle_qubit_list = []
        for qubit in self.circ.qubit_collection:
            if qubit not in active_qubit_list:
                idle_qubit_list.append(qubit)

        self.circ.single_qubit_noise(idle_qubit_list,
                                        noise=['DEPOLARIZE1', self.noise_probability])

        
        for i in [3,4]:
            self.circ.single_qubit_gate('H', pre_transform_dict[i],
                                        noise=['DEPOLARIZE1',self.noise_probability])
            self.circ.tick()
            self.circ.two_qubit_gate('CX', ctrl_qubit_dict[i],
                                           targ_qubit_dict[i],
                                           ['DEPOLARIZE2', self.noise_probability])

            active_qubit_list = ctrl_qubit_dict[i] + targ_qubit_dict[i]

            idle_qubit_list = []
            for qubit in self.circ.qubit_collection:
                if qubit not in active_qubit_list:
                    idle_qubit_list.append(qubit)

            self.circ.single_qubit_noise(idle_qubit_list,
                                         noise=['DEPOLARIZE1', self.noise_probability])

            self.circ.tick()
        
        self.circ.single_qubit_gate('H',final_transform_list,
                                    noise=['DEPOLARIZE1',self.noise_probability])
        self.circ.measure_reset(self.code.x_check_collection
                                + self.code.z_check_collection,
                                ['X_ERROR', self.noise_probability])

        for z_check in self.code.z_check_collection:
            self.circ.detector([z_check], self.circ.time_temp - 1,
                               [z_check], self.circ.time_temp,
                               (z_check.pos[0], z_check.pos[1], 
                                self.circ.time_temp))
        
        for x_check in self.code.x_check_collection:
            # right boundary
            if x_check.pos[0] == self.code.distance - 0.5:
                self.circ.detector([x_check], self.circ.time_temp-1,
                                   [x_check], self.circ.time_temp,
                                   (x_check.pos[0], x_check.pos[1],
                                    self.circ.time_temp))
            # not inner right boudary
            elif x_check.pos[0] != self.code.distance - 1.5:
                dual_pos = (self.code.distance-1-x_check.pos[1],
                            self.code.distance-2-x_check.pos[0])
                dual_check = self.circ.pos_qubit_dict[dual_pos]
                self.circ.detector([x_check], self.circ.time_temp-1,
                                   [x_check, dual_check], 
                                    self.circ.time_temp,
                                   (x_check.pos[0], x_check.pos[1],
                                    self.circ.time_temp))

        # self.circ.detector(stab_prior, self.circ.time_temp - 1,
        #                        stab_current, self.circ.time_temp,
        #                 (qubit.pos[0],qubit.pos[1], self.circ.time_temp))


        
        self.circ.clock_plus_one()
        self.circ.tick()

        # Next round. 

        self.circ.single_qubit_gate('I', self.code.data_qubits_collection,
                                    noise = ['DEPOLARIZE1', self.noise_probability])
        
        for i in [1,2,3,4]:
            self.circ.single_qubit_gate('H', pre_transform_dict[i],
                                        noise=['DEPOLARIZE1',self.noise_probability])
            self.circ.tick()
            self.circ.two_qubit_gate('CX', ctrl_qubit_dict[i],
                                           targ_qubit_dict[i],
                                           ['DEPOLARIZE2', self.noise_probability])
            active_qubit_list = ctrl_qubit_dict[i] + targ_qubit_dict[i]

            idle_qubit_list = []
            for qubit in self.circ.qubit_collection:
                if qubit not in active_qubit_list:
                    idle_qubit_list.append(qubit)

            self.circ.single_qubit_noise(idle_qubit_list,
                                         noise=['DEPOLARIZE1', self.noise_probability])
            
            self.circ.tick()

        self.circ.single_qubit_gate('H',final_transform_list,
                                    noise=['DEPOLARIZE1',self.noise_probability])
        self.circ.measure_reset(self.code.x_check_collection
                                + self.code.z_check_collection,
                                ['X_ERROR', self.noise_probability])

        for z_check in self.code.z_check_collection:
            self.circ.detector([z_check], self.circ.time_temp - 1,
                               [z_check], self.circ.time_temp,
                               (z_check.pos[0], z_check.pos[1], 
                                self.circ.time_temp))
        
        for x_check in self.code.x_check_collection:
            # left boundary
            if x_check.pos[0] == - 0.5:
                self.circ.detector([x_check], self.circ.time_temp-1,
                                   [x_check], self.circ.time_temp,
                                   (x_check.pos[0], x_check.pos[1],
                                    self.circ.time_temp))
            # inner right boudary 
            elif x_check.pos[0] == self.code.distance - 1.5:
                # (prev round)
                dual_pos = (self.code.distance-1-x_check.pos[1],
                            self.code.distance-2-x_check.pos[0])
                dual_check = self.circ.pos_qubit_dict[dual_pos]
                self.circ.detector_qt_list(
                                    [x_check, x_check, dual_check], 
                                    [self.circ.time_temp-2,
                                     self.circ.time_temp-1,
                                     self.circ.time_temp-1],
                                   (x_check.pos[0], x_check.pos[1],
                                    self.circ.time_temp-1))
                # (This round)
                self.circ.detector([x_check], self.circ.time_temp-1,
                                   [x_check], self.circ.time_temp,
                                   (x_check.pos[0], x_check.pos[1],
                                    self.circ.time_temp))
            # inner left boundary
            elif x_check.pos[0] == 0.5:
                dual_pos = (self.code.distance-1-x_check.pos[1],
                            self.code.distance-x_check.pos[0])
                dual_check = self.circ.pos_qubit_dict[dual_pos]
                self.circ.detector_qt_list(
                                    [x_check, x_check, dual_check],
                                    [self.circ.time_temp-1,
                                     self.circ.time_temp,
                                     self.circ.time_temp],
                                    (x_check.pos[0], x_check.pos[1],
                                     self.circ.time_temp)
                                    )
            else:
                self.circ.detector([x_check], self.circ.time_temp-1,
                                   [x_check], self.circ.time_temp,
                                   (x_check.pos[0], x_check.pos[1],
                                    self.circ.time_temp))

        
        self.circ.clock_plus_one()
        self.circ.tick()



    def sample_result_parity(self, qubit_list: list[Qubit],
                             time_list: list[int], num_shots: int = 5000) -> float:
        if len(qubit_list) != len(time_list):
            raise ValueError('Parameter length mismatch.')
        circuit = stim.Circuit.from_file(self.circuit_filename)
        sampler = circuit.compile_sampler(skip_reference_sample=False)
        result_samples = sampler.sample(shots=num_shots)
        parity_result_samples = []
        result_neg_index_list = []
        for qubit, time in zip(qubit_list, time_list):
            index = -1
            for i, time_meas_id_tuple in enumerate(qubit.time_meas_id_seq):
                if time == time_meas_id_tuple[0]:
                    index = i 
            if index == -1:
                raise ValueError('undefined time step')
            meas_id = qubit.time_meas_id_seq[index][-1]
            neg_pos = (self.circ.meas_id_list.index(meas_id) 
                        -len(self.circ.meas_id_list))
            result_neg_index_list.append(neg_pos)
        for result_line in result_samples:
            subset = [result_line[i] for i in result_neg_index_list]
            parity_result_samples.append(sum(subset) % 2)
        parity_avg = sum(parity_result_samples)/num_shots
        print('parity_avg: ',parity_avg)
        return parity_avg

    def check_circuit_distance(self) -> int:
        circuit = stim.Circuit.from_file(self.circuit_filename)
        explained_error = circuit.search_for_undetectable_logical_errors(
                            dont_explore_detection_event_sets_with_size_above=10,
                            dont_explore_edges_with_degree_above=5,
                            dont_explore_edges_increasing_symptom_degree=False)
        print(len(explained_error))
        for item in explained_error:
            print(item)
        return len(explained_error)
    

        
    def get_stim_circuit(self) -> stim.Circuit:
        circuit = stim.Circuit.from_file(self.circuit_filename)
        return circuit
    







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





