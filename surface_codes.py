from __future__ import annotations

import numpy as np
from numpy import ndarray
from typing import Any, Union, TypeVar, Tuple
from warnings import warn
import stim


oddint = TypeVar('oddint')

class Operation():

    def __init__(self, operation: str, target: list[Qubit],
                 operation_id: int,
                 tick_id: int,
                 meas_id: Union[int, None],
                 clocktime: int,
                 arg: Union[float, None] = None
                 ) -> None:
        self.operation = operation
        self.target_id_list = [qubit.circuit_id for qubit in target]
        self.operation_id = operation_id
        self.tick_id = tick_id
        self.meas_id = meas_id
        self.clocktime = clocktime
        self.arg = None

class Measurement():

    def __init__(self, target: Qubit,
                 meas_id: Union[int, None],
                 clocktime: int
                 ) -> None:
        self.target = target
        self.meas_id = meas_id
        self.clocktime = clocktime

class Detector():

    def __init__(self, spacetime_coord: tuple[float], 
                 detector_id: int, 
                 list_meas_id: list[int] = None
                ) -> None:
        self.spacetime_coord = spacetime_coord
        self.list_meas_id = []
        self.detector_id = detector_id
        self.flag = False
        self.cz_flag = False
        pass

    def _append_meas_id(self, meas_id: int) -> None:
        self.list_meas_id.append(meas_id)
        pass

     

        
class Circuit_helper():
    # helper class for the stim package. 
    def __init__(self, data_qubit_collection: list[Qubit],
                 X_qubit_collection: list[Qubit],
                 Z_qubit_collection: list[Qubit],
                 stim_filename: str
                 ) -> None:
        self.data_qubit_collection = data_qubit_collection
        self.X_qubit_collection = X_qubit_collection
        self.Z_qubit_collection = Z_qubit_collection
        self.stim_filename = stim_filename
        open(self.stim_filename, "w").close()

        self.qubit_coll_dict: dict[str,list[Qubit]] = {}
        self.qubit_coll_dict['Data'] = self.data_qubit_collection
        self.qubit_coll_dict['X'] = self.X_qubit_collection
        self.qubit_coll_dict['Z'] = self.Z_qubit_collection

        # self.qubit_collection, self.id_qubit_dict = self._stim_qubit_initialize()
        # self.num_qubits = len(self.qubit_collection)
        self.time_temp: int = 0
        self.clock: list[int] = []
        self.time_meas_dict: dict[int,list[Measurement]] = {}
        self.tick_temp = 0
        self.tick_seq = [0]
        self.operation_id_temp = 0
        self.operation_id_seq = []
        self.operation_id_dict: dict[int,Operation] = {}
        self.meas_num = 0
        self.meas_id_list = []
        self.meas_oper_dict: dict[int,Operation] = {}
        self.detector_id_temp = 0
        self.detector_id_list = []
        self.detector_dict: dict[int, Detector] = {}
        self.spacetime_detector_dict: dict[Tuple,Detector] = {}

        self.qubit_collection: Union[list[Qubit],None] = None
        self.id_qubit_dict: dict[int,Qubit] = {}
        self.num_qubits = None
        self.pos_qubit_dict: dict[Tuple,Qubit] = {}
        pass


    def initialize(self) -> None:

        num_qubits_temp = 0
        qubit_collection = []
        for qubit_coll in self.qubit_coll_dict.values():
            num_qubits_temp += len(qubit_coll)
            qubit_collection = qubit_collection + qubit_coll
        
        self.qubit_collection = qubit_collection
        self.num_qubits =  num_qubits_temp
        circuit_id_list = [i for i in range(self.num_qubits)]
        for circuit_id, qubit in zip(circuit_id_list,
                                     self.qubit_collection):
            qubit.set_circuit_id(circuit_id)
            self.id_qubit_dict[circuit_id] = qubit
            self.pos_qubit_dict[qubit.pos] = qubit
        
        with open(self.stim_filename,'a') as f:
            for qubit in self.qubit_collection:
                f.write('QUBIT_COORDS'+qubit.pos_str+' '+str(qubit.circuit_id))
                f.write('\n')
        


        pass


    def add_qubit_coll(self, qubit_coll: list[Qubit],
                        qtype: str) -> None:
        self.qubit_coll_dict[qtype] = qubit_coll



    def reset(self, qubit_collection: list[Qubit],
                noise: Union[list[str,float],None] = None):
        with open(self.stim_filename,'a') as f:
            f.write('R')
            for qubit in qubit_collection:
                oper = Operation('R', [qubit], self.operation_id_temp,
                                 self.tick_temp, None, self.time_temp)
                self.operation_id_dict[self.operation_id_temp] = oper
                qubit.append_operation(self.tick_temp, self.operation_id_temp)
                self.operation_id_seq.append(self.operation_id_temp)
                f.write(' ')
                f.write(str(qubit.circuit_id))
                self.operation_id_temp = self.operation_id_temp + 1
            f.write('\n')
        if noise != None:
            # No matter the noise type, do the X_ERROR. 
            # However, do give a warning if the noise type is not 
            # the X_ERROR.
            noise_str = noise[0]
            noise_probability = noise[1]
            if noise_str != 'X_ERROR':
                warn('Not an X_ERROR after reset.')
            with open(self.stim_filename,'a') as f:
                f.write('X_ERROR')
                f.write('(' + str(noise_probability) + ')')
                for qubit in qubit_collection:
                    f.write(' ')
                    f.write(str(qubit.circuit_id))
                f.write('\n')
    
    def single_qubit_gate(self, gate_name: str, qubit_collection: list[Qubit],
                          noise: Union[list[str,float],None] = None):
        if qubit_collection == None:
            return None
        if len(qubit_collection) == 0:
            return None
        with open(self.stim_filename,'a') as f:
            f.write(gate_name)
            for qubit in qubit_collection:
                oper = Operation(gate_name, [qubit], self.operation_id_temp,
                                 self.tick_temp, None, self.time_temp)
                self.operation_id_dict[self.operation_id_temp] = oper
                self.operation_id_seq.append(self.operation_id_temp)
                qubit.append_operation(self.tick_temp, self.operation_id_temp)
                f.write(' ')
                f.write(str(qubit.circuit_id))
                self.operation_id_temp = self.operation_id_temp + 1
            f.write('\n')
        if noise != None:
            noise_str = noise[0]
            noise_probability = noise[1]
            with open(self.stim_filename,'a') as f:
                f.write(noise_str)
                f.write('(' + str(noise_probability) + ')')
                for qubit in qubit_collection:
                    if qubit != None:
                        f.write(' ')
                        f.write(str(qubit.circuit_id))
                f.write('\n')
    
    def single_qubit_noise(self, qubit_collection: list[Qubit],
                           noise: Union[list[str,float],None] = None):
        if qubit_collection == None:
            return None
        if len(qubit_collection) == 0:
            return None
        if noise != None:
            noise_str = noise[0]
            noise_probability = noise[1]
            with open(self.stim_filename,'a') as f:
                f.write(noise_str)
                f.write('(' + str(noise_probability) + ')')
                for qubit in qubit_collection:
                    if qubit != None:
                        f.write(' ')
                        f.write(str(qubit.circuit_id))
                f.write('\n')
        
    
    def two_qubit_gate(self, gate_name: str, 
                       control_qubit_collection: list[Qubit],
                       target_qubit_collection: list[Qubit],
                       noise: Union[list[str,float],None] = None):
        if control_qubit_collection == None:
            return None
        if len(control_qubit_collection) == 0:
            return None
        if len(control_qubit_collection) != len(target_qubit_collection):
            raise ValueError('control-target mismatch')
        with open(self.stim_filename,'a') as f:
            f.write(gate_name)
            for ctrl_qubit, targ_qubit in zip(control_qubit_collection,
                                              target_qubit_collection):
                if ctrl_qubit != None and targ_qubit != None:
                    oper = Operation(gate_name, [ctrl_qubit, targ_qubit], 
                                     self.operation_id_temp,
                                     self.tick_temp, None, self.time_temp)
                    self.operation_id_dict[self.operation_id_temp] = oper
                    self.operation_id_seq.append(self.operation_id_temp)
                    ctrl_qubit.append_operation(self.tick_temp, 
                                                self.operation_id_temp)
                    targ_qubit.append_operation(self.tick_temp,
                                                self.operation_id_temp)
                    f.write(' ')
                    f.write(str(ctrl_qubit.circuit_id))
                    f.write(' ')
                    f.write(str(targ_qubit.circuit_id))
                    self.operation_id_temp = self.operation_id_temp + 1
            f.write('\n')
        if noise != None:
            noise_str = noise[0]
            noise_probability = noise[1]
            with open(self.stim_filename,'a') as f:
                f.write(noise_str)
                f.write('(' + str(noise_probability) + ')')
                for ctrl_qubit, targ_qubit in zip(control_qubit_collection,
                                                  target_qubit_collection):
                    if ctrl_qubit != None and targ_qubit != None:         
                        f.write(' ')
                        f.write(str(ctrl_qubit.circuit_id))
                        f.write(' ')
                        f.write(str(targ_qubit.circuit_id))
                f.write('\n')

        pass

    def two_qubit_noise(self, 
                       control_qubit_collection: list[Qubit],
                       target_qubit_collection: list[Qubit],
                       noise: Union[list[str,float],None] = None):
        if control_qubit_collection == None:
            return None
        if len(control_qubit_collection) == 0:
            return None
        if len(control_qubit_collection) != len(target_qubit_collection):
            raise ValueError('control-target mismatch')
        if noise != None:
            noise_str = noise[0]
            noise_probability = noise[1]
            with open(self.stim_filename,'a') as f:
                f.write(noise_str)
                f.write('(' + str(noise_probability) + ')')
                for ctrl_qubit, targ_qubit in zip(control_qubit_collection,
                                                  target_qubit_collection):
                    if ctrl_qubit != None and targ_qubit != None:         
                        f.write(' ')
                        f.write(str(ctrl_qubit.circuit_id))
                        f.write(' ')
                        f.write(str(targ_qubit.circuit_id))
                f.write('\n')

        pass

    def measure_reset(self, qubit_collection: list[Qubit],
                      noise: Union[list[str,float],None] = None) -> None:
        self.clock.append(self.time_temp)
        meas_oper_list = []
        if noise != None:
            # No matter the noise type, do the X_ERROR. 
            # However, do give a warning if the noise type is not 
            # the X_ERROR.
            noise_str = noise[0]
            noise_probability = noise[1]
            if noise_str != 'X_ERROR':
                warn('Not an X_ERROR after reset.')
            with open(self.stim_filename,'a') as f:
                f.write('X_ERROR')
                f.write('(' + str(noise_probability) + ')')
                for qubit in qubit_collection:
                    f.write(' ')
                    f.write(str(qubit.circuit_id))
                f.write('\n')
        with open(self.stim_filename,'a') as f:
            f.write('MR')
            for qubit in qubit_collection:
                self.meas_num = self.meas_num + 1
                meas_id = self.meas_num
                self.meas_id_list.append(meas_id)
                qubit.meas_id_seq.append(meas_id)
                qubit.time_meas_id_seq.append((self.time_temp, meas_id))
                qubit.time_meas_id_dict[self.time_temp] = meas_id
                meas_oper = Operation('MR', [qubit], self.operation_id_temp,
                                      self.tick_temp, meas_id, self.time_temp)
                self.operation_id_dict[self.operation_id_temp] = meas_oper
                self.operation_id_seq.append(self.operation_id_temp)
                qubit.append_operation(self.tick_temp, self.operation_id_temp)
                meas_oper_list.append(meas_oper)
                self.meas_oper_dict[meas_id] = meas_oper
                f.write(' ')
                f.write(str(qubit.circuit_id))
                self.operation_id_temp = self.operation_id_temp + 1
            f.write('\n')
        if noise != None:
            # No matter the noise type, do the X_ERROR. 
            # However, do give a warning if the noise type is not 
            # the X_ERROR.
            noise_str = noise[0]
            noise_probability = noise[1]
            if noise_str != 'X_ERROR':
                warn('Not an X_ERROR after reset.')
            with open(self.stim_filename,'a') as f:
                f.write('X_ERROR')
                f.write('(' + str(noise_probability) + ')')
                for qubit in qubit_collection:
                    f.write(' ')
                    f.write(str(qubit.circuit_id))
                f.write('\n')
        self.time_meas_dict[self.time_temp] = meas_oper_list
        # self.time_temp = self.time_temp + 1
        pass

    def measure(self, qubit_collection: list[Qubit],
                noise: Union[list[str,float],None] = None) -> None:
        self.clock.append(self.time_temp)
        meas_oper_list = []
        if noise != None:
            # No matter the noise type, do the X_ERROR. 
            # However, do give a warning if the noise type is not 
            # the X_ERROR.
            noise_str = noise[0]
            noise_probability = noise[1]
            if noise_str != 'X_ERROR':
                warn('Not an X_ERROR after reset.')
            with open(self.stim_filename,'a') as f:
                f.write('X_ERROR')
                f.write('(' + str(noise_probability) + ')')
                for qubit in qubit_collection:
                    f.write(' ')
                    f.write(str(qubit.circuit_id))
                f.write('\n')
        with open(self.stim_filename,'a') as f:
            f.write('M')
            for qubit in qubit_collection:
                self.meas_num = self.meas_num + 1
                meas_id = self.meas_num
                self.meas_id_list.append(meas_id)
                qubit.meas_id_seq.append(meas_id)
                qubit.time_meas_id_seq.append((self.time_temp, meas_id))
                qubit.time_meas_id_dict[self.time_temp] = meas_id
                meas_oper = Operation('MR', [qubit], self.operation_id_temp,
                                      self.tick_temp, meas_id, self.time_temp)
                self.operation_id_dict[self.operation_id_temp] = meas_oper
                self.operation_id_seq.append(self.operation_id_temp)
                qubit.append_operation(self.tick_temp, self.operation_id_temp)
                meas_oper_list.append(meas_oper)
                self.meas_oper_dict[meas_id] = meas_oper
                f.write(' ')
                f.write(str(qubit.circuit_id))
                self.operation_id_temp = self.operation_id_temp + 1
            f.write('\n')
        self.time_meas_dict[self.time_temp] = meas_oper_list
        # self.time_temp = self.time_temp + 1
        pass

    def clock_plus_one(self) -> None:
        self.time_temp = self.time_temp + 1

    def detector(self, qubit_front_list: list[Qubit], time_front: int,
                 qubit_back_list: list[Qubit], time_back: int, 
                 spacetime_coord: tuple[float],
                 flag: bool = False,
                 cz_flag: bool = False) -> None:
        if len(qubit_front_list) + len(qubit_back_list) == 0:
            raise ValueError('operation not supported on anything')
        ret_detector = Detector(spacetime_coord, 
                                self.detector_id_temp)
        if flag == True:
            ret_detector.flag = True
        if cz_flag == True:
            ret_detector.cz_flag = True
        self.detector_id_list.append(self.detector_id_temp)
        spacetime_coord_str = pos_str_conversion_general([coord for
                                                    coord in spacetime_coord])
        with open(self.stim_filename,'a') as f:
            f.write('DETECTOR')
            f.write(spacetime_coord_str)
            for qubit in qubit_front_list:
                if qubit == None:
                    continue
                f.write(' rec')
                # search for meas_id at time_front:
                index = -1
                for i, time_meas_id_tuple in enumerate(qubit.time_meas_id_seq):
                    if time_front == time_meas_id_tuple[0]:
                        index = i 
                if index == -1:
                    raise ValueError('undefined time step')
                meas_id = qubit.time_meas_id_seq[index][-1]
                ret_detector._append_meas_id(meas_id)
                neg_pos = (self.meas_id_list.index(meas_id) 
                           -len(self.meas_id_list))
                f.write(str([neg_pos]))
            for qubit in qubit_back_list:
                if qubit == None:
                    continue
                f.write(' rec')
                # search for meas_id at time_front:
                index = -1
                for i, time_meas_id_tuple in enumerate(qubit.time_meas_id_seq):
                    if time_back == time_meas_id_tuple[0]:
                        index = i 
                if index == -1:
                    raise ValueError('undefined time step')
                meas_id = qubit.time_meas_id_seq[index][-1]
                ret_detector._append_meas_id(meas_id)
                neg_pos = (self.meas_id_list.index(meas_id) 
                           -len(self.meas_id_list))
                f.write(str([neg_pos]))
            f.write('\n')
        self.detector_dict[self.detector_id_temp] = ret_detector
        self.spacetime_detector_dict[ret_detector.spacetime_coord] = \
                                                    ret_detector
        self.detector_id_temp = self.detector_id_temp + 1
        pass
    
    def detector_qt_list(self, 
                 qubit_list: list[Qubit], 
                 time_list: list[int],
                 spacetime_coord: tuple[float],
                 flag: bool = False,
                 cz_flag: bool = False) -> None:
        if len(qubit_list) == 0:
            raise ValueError('detector on null set')
        if len(qubit_list) != len(time_list):
            raise ValueError('length mismatch')
        ret_detector = Detector(spacetime_coord, 
                                self.detector_id_temp)
        if flag == True:
            ret_detector.flag = True
        if cz_flag == True:
            ret_detector.cz_flag = True
        self.detector_id_list.append(self.detector_id_temp)
        spacetime_coord_str = pos_str_conversion_general([coord for
                                                    coord in spacetime_coord])
        with open(self.stim_filename,'a') as f:
            f.write('DETECTOR')
            f.write(spacetime_coord_str)
            for qubit, time_q in zip(qubit_list, time_list):
                if qubit == None:
                    continue
                f.write(' rec')
                # search for meas_id at time_q:
                index = -1
                for i, time_meas_id_tuple in enumerate(qubit.time_meas_id_seq):
                    if time_q == time_meas_id_tuple[0]:
                        index = i 
                if index == -1:
                    raise ValueError('no meas result at this time step')
                meas_id = qubit.time_meas_id_seq[index][-1]
                ret_detector._append_meas_id(meas_id)
                neg_pos = (self.meas_id_list.index(meas_id) 
                           -len(self.meas_id_list))
                f.write(str([neg_pos]))

            f.write('\n')
        self.detector_dict[self.detector_id_temp] = ret_detector
        self.spacetime_detector_dict[ret_detector.spacetime_coord] = \
                                                    ret_detector
        self.detector_id_temp = self.detector_id_temp + 1
        pass


    def observable(self, qubit_list: list[Qubit], time: int) -> None:
        with open(self.stim_filename,'a') as f:
            f.write('OBSERVABLE_INCLUDE')
            f.write('(0)')
            for qubit in qubit_list:
                f.write(' rec')
                # search for meas_id at time_front:
                index = -1
                for i, time_meas_id_tuple in enumerate(qubit.time_meas_id_seq):
                    if time == time_meas_id_tuple[0]:
                        index = i 
                if index == -1:
                    raise ValueError('undefined time step')
                meas_id = qubit.time_meas_id_seq[index][-1]
                neg_pos = (self.meas_id_list.index(meas_id) 
                           -len(self.meas_id_list))
                f.write(str([neg_pos]))
        pass


    def tick(self) -> None:
        with open(self.stim_filename,'a') as f:
            f.write('TICK\n')
        self.tick_temp = self.tick_temp + 1
        pass
    
    
    


class Flag():

    def __init__(self, dict_qubit_coll: dict[str,list[Qubit]],
                 flag_displacement: int,
                 flag_type: str) -> None:
        self.dict_qubit_coll = dict_qubit_coll
        self.num_qubits = self._sanity_check()
        self.flag_list: list[Qubit] = []
        for i in range(self.num_qubits):
            list_temp = list(self.dict_qubit_coll.values())[0]
            qubit = list_temp[i]
            pos = qubit.pos
            flag_pos = (pos[0] + flag_displacement,
                        pos[1] + flag_displacement)
            flag = Qubit(id = i, pos=flag_pos, tribe=flag_type)
            self.flag_list.append(flag)

        for qtype in self.dict_qubit_coll.keys():
            list_temp = self.dict_qubit_coll[qtype]
            for qubit, flag in zip(list_temp, 
                                   self.flag_list):
                flag.append_neighbor(qubit,qtype)
                qubit.append_neighbor(flag,flag_type)
        

        pass

    def _sanity_check(self) -> int:
        num_qubits_temp = len(list(self.dict_qubit_coll.values())[0])
        for qubit_coll in self.dict_qubit_coll.values():
            if len(qubit_coll) != num_qubits_temp:
                raise ValueError('size mismatch')
                return False
        return num_qubits_temp


            

class Qubit():
    
    def __init__(self, id: int, pos: tuple[float, float], tribe: str,
                 tag: str = None
                 ) -> None:
        self.id = id
        self.pos = pos
        self.pos_str = pos_str_conversion(self.pos)
        self.tribe = tribe     # 'D', 'X', 'Z', 'FLAG'
        self.neighbor: dict[str,Union[Qubit,None]] = {}
        self.check = {}
        self.css = None
        self.tag = tag
        # self.oper_seq = []
        self.meas_id_seq = []
        self.time_meas_id_seq = []
        self.time_meas_id_dict: dict[int,int] = {}
        self.circuit_id = None
        self.tick_oper_id_dict: dict[int, int] = {}
    
    def append_neighbor(self, qubit: Union[Qubit,None], rel_pos_str: str):
        self.neighbor[rel_pos_str] = qubit

    def append_check(self, qubit: Qubit, single_check: str):
        self.check[pos_str_conversion(qubit.pos)] = [qubit, single_check]
        pass

    def set_circuit_id(self, circuit_id: int) -> None:
        self.circuit_id = circuit_id

    def append_operation(self, tick: int, oper_id: int) -> None:
        self.tick_oper_id_dict[tick] = oper_id
        pass


class Rotated_surface_code():


    def __init__(self, distance: oddint) -> None:
        if (distance%2) != 1:
            raise ValueError('distance is not an odd number') 
        self.distance = distance
        self.num_qubits = distance**2
        self.id_collection = [i for i in range(self.num_qubits)]
        self.data_qubits_collection = self.initialize_data_qubits()
        self.data_qubits_network, self.data_qubits_coord = grid_network(
            self.distance, self.distance, 1, 1, 0, 0,
            self.data_qubits_collection)
        self.z_check_collection = self.initialize_z_checks()
        self.z_check_network, self.z_check_coord \
                = self.generate_z_check_network()
        self.num_z_checks = len(self.z_check_collection)
        self.x_check_collection = self.initilize_x_checks()
        self.num_x_checks = len(self.x_check_collection)
        self.x_check_network, self.x_check_coord \
                = self.generate_x_check_network()
        self.Hx = self.compose_Hx()
        self.Hz = self.compose_Hz()
        self.logic_z_collection = self._logic_z_assemble()
        self.logic_x_collection = self._logic_x_assemble()
        self.pos_y_data_qubits_dict: dict[str,list[Qubit]] = \
                                        self._pos_y_data_qubits_dict_gen()
        self.pos_y_ancilla_qubits_dict: dict[str,list[Qubit]] = \
                                        self._pos_y_ancilla_qubits_dict_gen()
        pass

    def _pos_y_data_qubits_dict_gen(self) -> dict[str,list[Qubit]]:
        pos_y_data_qubits_dict = {}
        for y_pos in range(self.distance):
            row_data_qubits = []
            for x_pos in range(self.distance):
                pos_str = pos_str_conversion_general([x_pos,y_pos])
                row_data_qubits.append(self.data_qubits_network[pos_str])
            pos_y_data_qubits_dict[str(y_pos)] = row_data_qubits
        return pos_y_data_qubits_dict
    
    def _pos_y_ancilla_qubits_dict_gen(self) -> dict[str,list[Qubit]]:
        pos_y_data_qubits_dict = {}
        y_pos_list = [i-0.5 for i in range(self.distance)]
        y_pos_list = y_pos_list[1:]
        x_pos_list = [i-0.5 for i in range(self.distance)]
        x_check_pos_list = [pos_str_conversion(qubit.pos)
                                 for qubit in self.x_check_collection]
        z_check_pos_list = [pos_str_conversion(qubit.pos)
                                 for qubit in self.z_check_collection]
        for y_pos in y_pos_list:
            row_data_qubits = []
            for x_pos in x_pos_list:
                pos_str = pos_str_conversion_general([x_pos,y_pos])
                if pos_str in x_check_pos_list:
                    row_data_qubits.append(self.x_check_network[pos_str])
                elif pos_str in z_check_pos_list:
                    row_data_qubits.append(self.z_check_network[pos_str])
            pos_y_data_qubits_dict[
                pos_str_conversion_general([y_pos])] = row_data_qubits
        return pos_y_data_qubits_dict

    def _logic_z_assemble(self) -> list[Qubit]:
        x_pos = (self.distance-1)/2
        ret_z_list = []
        for qubit in self.data_qubits_collection:
            if x_pos == qubit.pos[0]:
                ret_z_list.append(qubit)
        return ret_z_list

    def _logic_x_assemble(self) -> list[Qubit]:
        y_pos = (self.distance-1)/2
        ret_x_list = []
        for qubit in self.data_qubits_collection:
            if y_pos == qubit.pos[1]:
                ret_x_list.append(qubit)
        return ret_x_list

    def compose_Hx(self) -> ndarray[int]:
        Hx = np.zeros((self.num_x_checks, self.num_qubits), dtype=np.int64)
        for check in self.x_check_collection:
            row_num = check.id
            for pos_str in check.check:
                qubit_supp = check.check[pos_str][0].id
                Hx[row_num][qubit_supp] = 1
        return Hx

    def compose_Hz(self) -> ndarray[int]:
        Hz = np.zeros((self.num_z_checks, self.num_qubits), dtype=np.int64)
        for check in self.z_check_collection:
            row_num = check.id
            for pos_str in check.check:
                qubit_supp = check.check[pos_str][0].id
                Hz[row_num][qubit_supp] = 1
        return Hz

    def initialize_data_qubits(self) -> list[Qubit]:
        data_qubits_collection = []
        for i in range(self.distance):
            for j in range(self.distance):
                order = i*self.distance + j
                pos = (i,j)
                id = self.id_collection[order]
                data_qubits_collection.append(Qubit(id, pos, 'D'))
        return data_qubits_collection
    
    def initilize_x_checks(self) -> list[Qubit]: # Unfinished. 
        x_check_collection = []
        id_sub_A_collection = []
        id_sub_B_collection = []
        id_left_bd_collection = []
        id_right_bd_collection = []
        id_collection = []
        id_index = -1 
        for i in range(self.distance//2):
            for j in range(self.distance//2):
                id_index = id_index + 1
                id_sub_A_collection.append(id_index)
                id_collection.append(id_index)
                pos_x = i*2 + 0.5
                pos_y = j*2 + 0.5
                x_check_temp = Qubit(id_index, (pos_x, pos_y), 'X', 'sub_A')
                x_check_collection.append(x_check_temp)
        for i in range(self.distance//2):
            for j in range(self.distance//2):
                id_index = id_index + 1
                id_sub_B_collection.append(id_index)
                id_collection.append(id_index)
                pos_x = i*2 + 1.5
                pos_y = j*2 + 1.5
                x_check_temp = Qubit(id_index, (pos_x, pos_y), 'X', 'sub_B')
                x_check_collection.append(x_check_temp)
        for j in range(self.distance//2):
            id_index = id_index + 1
            id_left_bd_collection.append(id_index)
            id_collection.append(id_index)
            pos_x = -0.5
            pos_y = 2*j + 1.5
            x_check_temp = Qubit(id_index, (pos_x, pos_y), 'X', 'left_bd')
            x_check_collection.append(x_check_temp)
        # print(id_collection)
        for j in range(self.distance//2):
            id_index = id_index + 1
            id_right_bd_collection.append(id_index)
            id_collection.append(id_index)
            pos_x = self.distance - 0.5
            pos_y = 2*j + 0.5
            x_check_temp = Qubit(id_index,(pos_x, pos_y), 'X', 'right_bd')
            x_check_collection.append(x_check_temp)
        return x_check_collection

    
    def initialize_z_checks(self) -> list[Qubit]:
        z_check_collection = []
        id_upper_bd_collection = []
        id_sub_A_collection = []
        id_sub_B_collection = []
        id_lower_bd_collection = []
        id_collection = []
        id_index = -1 
        for i in range(self.distance//2):
            id_index = id_index + 1
            id_upper_bd_collection.append(id_index)
            id_collection.append(id_index)
            pos_x = i*2 + 0.5
            pos_y = -0.5
            z_check_temp = Qubit(id_index, (pos_x, pos_y), 'Z', 'upper_bd')
            z_check_collection.append(z_check_temp)
        for i in range(self.distance//2):
            for j in range(self.distance//2):
                id_index = id_index + 1
                id_sub_A_collection.append(id_index)
                id_collection.append(id_index)
                pos_x = i*2 + 0.5
                pos_y = j*2 + 1.5
                z_check_temp = Qubit(id_index, (pos_x, pos_y), 'Z', 'sub_A')
                z_check_collection.append(z_check_temp)
        for i in range(self.distance//2):
            for j in range(self.distance//2):
                id_index = id_index + 1
                id_sub_B_collection.append(id_index)
                id_collection.append(id_index)
                pos_x = i*2 + 1.5
                pos_y = j*2 + 0.5
                z_check_temp = Qubit(id_index, (pos_x, pos_y), 'Z', 'sub_B')
                z_check_collection.append(z_check_temp)
        for i in range(self.distance//2):
            id_index = id_index + 1
            id_lower_bd_collection.append(id_index)
            id_collection.append(id_index)
            pos_x = i*2 + 1.5
            pos_y = self.distance - 0.5
            z_check_temp = Qubit(id_index, (pos_x, pos_y), 'Z', 'lower_bd')
            z_check_collection.append(z_check_temp)
        # print(id_collection)
        return z_check_collection
    


    def generate_data_qubit_network(self): 
        
        pass

    def generate_x_check_network(self):
        x_check_network = {}
        x_check_coord = {}
        for check in self.x_check_collection:
            x_check_network[pos_str_conversion(check.pos)] = check
            x_check_coord[pos_str_conversion(check.pos)] = check.pos
        for check in self.x_check_collection:
            if check.tag=='sub_A' or check.tag=='sub_B':
                sw_pos = (check.pos[0]-0.5, check.pos[1]+0.5)
                sw_pos_str = pos_str_conversion(sw_pos)
                se_pos = (check.pos[0]+0.5, check.pos[1]+0.5)
                se_pos_str = pos_str_conversion(se_pos)
                nw_pos = (check.pos[0]-0.5, check.pos[1]-0.5)
                nw_pos_str = pos_str_conversion(nw_pos)
                ne_pos = (check.pos[0]+0.5, check.pos[1]-0.5)
                ne_pos_str = pos_str_conversion(ne_pos)
                check.append_neighbor(self.data_qubits_network[sw_pos_str],'sw')
                check.append_neighbor(self.data_qubits_network[se_pos_str],'se')
                check.append_neighbor(self.data_qubits_network[nw_pos_str],'nw')
                check.append_neighbor(self.data_qubits_network[ne_pos_str],'ne')
                check.append_check(self.data_qubits_network[sw_pos_str],'X')
                check.append_check(self.data_qubits_network[se_pos_str],'X')
                check.append_check(self.data_qubits_network[nw_pos_str],'X')
                check.append_check(self.data_qubits_network[ne_pos_str],'X')
            elif check.tag == 'left_bd':
                se_pos = (check.pos[0]+0.5, check.pos[1]+0.5)
                se_pos_str = pos_str_conversion(se_pos)
                ne_pos = (check.pos[0]+0.5, check.pos[1]-0.5)
                ne_pos_str = pos_str_conversion(ne_pos)
                check.append_neighbor(None,'sw')
                check.append_neighbor(self.data_qubits_network[se_pos_str],'se')
                check.append_neighbor(None,'nw')
                check.append_neighbor(self.data_qubits_network[ne_pos_str],'ne')
                check.append_check(self.data_qubits_network[se_pos_str],'X')
                check.append_check(self.data_qubits_network[ne_pos_str],'X')
            elif check.tag == 'right_bd':
                sw_pos = (check.pos[0]-0.5, check.pos[1]+0.5)
                sw_pos_str = pos_str_conversion(sw_pos)
                nw_pos = (check.pos[0]-0.5, check.pos[1]-0.5)
                nw_pos_str = pos_str_conversion(nw_pos)
                check.append_neighbor(self.data_qubits_network[sw_pos_str],'sw')
                check.append_neighbor(None,'se')
                check.append_neighbor(self.data_qubits_network[nw_pos_str],'nw')
                check.append_neighbor(None,'ne')
                check.append_check(self.data_qubits_network[sw_pos_str],'X')
                check.append_check(self.data_qubits_network[nw_pos_str],'X')
                pass
        return x_check_network, x_check_coord


    def generate_z_check_network(self):
        z_check_network = {}
        z_check_coord = {}
        for check in self.z_check_collection:
            z_check_network[pos_str_conversion(check.pos)] = check
            z_check_coord[pos_str_conversion(check.pos)] = check.pos
        for check in self.z_check_collection:
            pos_str = pos_str_conversion(check.pos)
            if check.tag == 'upper_bd':
                sw_pos = (check.pos[0]-0.5, 0)
                sw_pos_str = pos_str_conversion(sw_pos)
                se_pos = (check.pos[0]+0.5, 0)
                se_pos_str = pos_str_conversion(se_pos)
                check.append_neighbor(self.data_qubits_network[sw_pos_str],'sw')
                check.append_neighbor(self.data_qubits_network[se_pos_str],'se')
                check.append_neighbor(None,'nw')
                check.append_neighbor(None,'ne')
                check.append_check(self.data_qubits_network[sw_pos_str],'Z')
                check.append_check(self.data_qubits_network[se_pos_str],'Z')
            elif check.tag == 'lower_bd':
                nw_pos = (check.pos[0]-0.5, check.pos[1]-0.5)
                nw_pos_str = pos_str_conversion(nw_pos)
                ne_pos = (check.pos[0]+0.5, check.pos[1]-0.5)
                ne_pos_str = pos_str_conversion(ne_pos)
                check.append_neighbor(self.data_qubits_network[nw_pos_str],'nw')
                check.append_neighbor(self.data_qubits_network[ne_pos_str],'ne')
                check.append_neighbor(None,'sw')
                check.append_neighbor(None,'se')
                check.append_check(self.data_qubits_network[nw_pos_str],'Z')
                check.append_check(self.data_qubits_network[ne_pos_str],'Z')
            elif check.tag=='sub_A' or check.tag=='sub_B':
                sw_pos = (check.pos[0]-0.5, check.pos[1]+0.5)
                sw_pos_str = pos_str_conversion(sw_pos)
                se_pos = (check.pos[0]+0.5, check.pos[1]+0.5)
                se_pos_str = pos_str_conversion(se_pos)
                nw_pos = (check.pos[0]-0.5, check.pos[1]-0.5)
                nw_pos_str = pos_str_conversion(nw_pos)
                ne_pos = (check.pos[0]+0.5, check.pos[1]-0.5)
                ne_pos_str = pos_str_conversion(ne_pos)
                check.append_neighbor(self.data_qubits_network[sw_pos_str],'sw')
                check.append_neighbor(self.data_qubits_network[se_pos_str],'se')
                check.append_neighbor(self.data_qubits_network[nw_pos_str],'nw')
                check.append_neighbor(self.data_qubits_network[ne_pos_str],'ne')
                check.append_check(self.data_qubits_network[sw_pos_str],'Z')
                check.append_check(self.data_qubits_network[se_pos_str],'Z')
                check.append_check(self.data_qubits_network[nw_pos_str],'Z')
                check.append_check(self.data_qubits_network[ne_pos_str],'Z')
        return z_check_network, z_check_coord







def pos_str_conversion(pos: tuple[Union[int,oddint], 
                                  Union[int,oddint]]) -> str:
    pos_x = pos[0]
    pos_y = pos[1]
    if abs(pos_x-round(pos_x)) < 10**(-8):
        pos_x_str = "{:.0f}".format(pos_x)
    elif abs(abs(pos_x-round(pos_x))-0.5) < 10**(-8):
        pos_x_str = "{:.1f}".format(pos_x)
    else:
        return ValueError('pos_x has a strange value')
    if abs(pos_y-round(pos_y)) < 10**(-8):
        pos_y_str = "{:.0f}".format(pos_y)
    elif abs(abs(pos_y-round(pos_y))-0.5) < 10**(-8):
        pos_y_str = "{:.1f}".format(pos_y)
    else:
        return ValueError('pos_y has a strange value')
    pos_str = '(' + pos_x_str + ',' + pos_y_str + ')'
    return pos_str


def pos_str_conversion_general(pos: list[Union[int, oddint]]) -> str:
    ret_str = '('
    for coord, index in zip(pos,range(len(pos))):
        if abs(coord-round(coord)) < 10**(-8):
            coord_str = "{:.0f}".format(coord)
        elif abs(abs(coord-round(coord))-0.5) < 10**(-8):
            coord_str = "{:.1f}".format(coord)
        else:
            return ValueError('coord has a strange value')
        ret_str = ret_str + coord_str
        if index != len(pos)-1:
            ret_str = ret_str + ','
    ret_str = ret_str + ')'
    return ret_str



def grid_network(lx: float, ly: float,
                 x_gap: float, y_gap: float,
                 x_offset: float, y_offset: float, 
                 nodes: list[Qubit]) -> Tuple[dict[str, Qubit],
                                           dict[str, tuple[float, float]]]: 
    if len(nodes) != lx*ly:
        raise ValueError('length of the list of the nodes is incompatible')
    network = {}
    network_coord = {}
    for i in range(lx):
        for j in range(ly):
            pos_x = i*x_gap + x_offset
            pos_y = j*y_gap + y_offset
            pos_str = pos_str_conversion((pos_x, pos_y))
            network[pos_str] = nodes[i*ly + j]
            network_coord[pos_str] = (pos_x,pos_y)
    return [network, network_coord]






