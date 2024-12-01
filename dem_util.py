import pymatching.matching
from dem_parsor import *
from hypergraph import *
import pymatching
from joblib import Parallel, delayed
from scipy.sparse import csr_array, csc_array, csr_matrix, vstack, hstack
from timeit import default_timer as timer


# Decomposition of errors in dem into edges in decoding (hyper-)graph. 



class Qubit_ops():

    def __init__(self, qubit_id: int):
        self.qubit_id = qubit_id
        self.oper_dict = {}
        self.sq_fault_dict: dict[tuple,int] = {}
        self.oper_line_list = []
        pass

    def _previous_line(self, line_num: int) -> Union[int,None]:
        if self.oper_line_list.index(line_num) == 0:
            return None
        ret_value = self.oper_line_list[self.oper_line_list.index(line_num)-1]
        return ret_value
        
    def _resort(self) -> None:
        self.oper_line_list = sorted(self.oper_line_list)
        pass

    def _cp_fault_ids(self, line_num: int) -> list[int]:
        fault_id_list = []
        for err_type in ['X','Z']:
            try:
                fault_id_list.append(self.sq_fault_dict[(line_num,err_type)])
            except:
                continue
        if len(fault_id_list) == 0:
            raise ValueError('current line does not have any single qubit err')
        oper_index = self.oper_line_list.index(line_num)
        fault_id_curr = copy.deepcopy(fault_id_list)
        for i in range(oper_index):
            try:
                fault_id_x = self.sq_fault_dict[(self.oper_line_list[oper_index-1-i],
                                                'X')]
            except:
                continue
            if fault_id_x not in fault_id_curr:
                fault_id_list.append(fault_id_x)
                break
        
        for i in range(oper_index):
            try:
                fault_id_z = self.sq_fault_dict[(self.oper_line_list[oper_index-1-i],
                                                'Z')]
            except:
                continue
            if fault_id_z not in fault_id_curr:
                fault_id_list.append(fault_id_z)
                break
        
        return fault_id_list




class Decompose_DEM():


    def __init__(self, dem: DEM, num_dets: int):
        self.dem = dem
        self.num_dets = num_dets
        self.H: Union[csr_array,None] = None
        self.error_channel: Union[ndarray,None] = None
        self._construct_H() 
        self.qubit_ops_dict: dict[int,Qubit_ops] = {}
        self.decomp_dict:dict[int,list[int]] = {}
        self.sq_fault_id_list = []
        self.comp_fault_id_list = []
        self.explicit_fault_id = []
        self.nonexplicit_fault_id = []
        self.fault_det_dict: dict[int,frozenset] = {}
        self.basis_dict: dict[int,list] = {}
        self.decomp_csr: Union[csr_array,None] = None
        self._single_qubit_faults()
        self._fault_det_init()
        self._explicit_decomp()
        self._nonexplicit_decomp()
        



    def _construct_H(self):

        # error location -> classical bits 
        # detectors -> checks

        error_channel = []

        data_list = []
        row_list = []
        col_list = []
        num_row = self.num_dets
        num_col = len(self.dem.fault_id_list)

        for fault_id in range(num_col):
            try:
                fault_unit = self.dem.fault_id_dict[fault_id]
            except:
                raise ValueError('fault id not found.')
            det_ids = fault_unit.flipped_detector_id_list
            error_channel.append(fault_unit.err_probability)
            for det_id in det_ids:
                data_list.append(1)
                row_list.append(det_id)
                col_list.append(fault_id)

        self.H = csr_matrix((data_list,(row_list,col_list)),
                           shape=(num_row,num_col),
                           dtype=np.uint8)
        
        self.error_channel = np.array(error_channel)


    def _fault_det_init(self) -> None:

        for fault_unit in self.dem.fault_unit_list:
            self.fault_det_dict[fault_unit.fault_id] \
                = frozenset(fault_unit.flipped_detector_id_list)


    def _single_qubit_faults(self) -> None:
        for fault_unit in self.dem.fault_unit_list:
            sq_check = False
            for err_unit in fault_unit.err_unit_list:
                if len(err_unit.error_type) == 1 and \
                    'Y' not in err_unit.error_type:
                    qubit_id = err_unit.error_supp[0]
                    err_type = err_unit.error_type[0]
                    oper_line = err_unit.instruction_line
                    try:
                        Q_ops = self.qubit_ops_dict[qubit_id]
                        Q_ops.sq_fault_dict[(oper_line,err_type)] \
                                = fault_unit.fault_id
                        if oper_line not in Q_ops.oper_line_list:
                            Q_ops.oper_line_list.append(oper_line)
                    except:
                        self.qubit_ops_dict[qubit_id] = Qubit_ops(qubit_id)
                        self.qubit_ops_dict[qubit_id].sq_fault_dict\
                            [(oper_line,err_type)] = fault_unit.fault_id
                        Q_ops = self.qubit_ops_dict[qubit_id]
                        if oper_line not in Q_ops.oper_line_list:
                            Q_ops.oper_line_list.append(oper_line)
                    sq_check = True
            if sq_check == True:
                self.sq_fault_id_list.append(fault_unit.fault_id)
                self.decomp_dict[fault_unit.fault_id] = [fault_unit.fault_id]
            else:
                self.comp_fault_id_list.append(fault_unit.fault_id)
    
    def _explicit_decomp(self) -> None:
        for fault_id in self.comp_fault_id_list:
            fault_unit = self.dem.fault_id_dict[fault_id]
            explicitness = False
            for err_unit in fault_unit.err_unit_list:
                if len(err_unit.error_type) == 2 and \
                    'Y' not in err_unit.error_type:
                    fault_id_list = []
                    for err_type, qubit_id in zip(err_unit.error_type,
                                                  err_unit.error_supp):
                        fault_id_list.append(self.qubit_ops_dict[qubit_id]\
                                .sq_fault_dict[(err_unit.instruction_line,
                                                err_type)])
                    self.decomp_dict[fault_unit.fault_id] = fault_id_list
                    explicitness = True
                    break
                if len(err_unit.error_type) == 1 and \
                    'Y' in err_unit.error_type:
                    fault_id_list = []
                    qubit_id = err_unit.error_supp[0]
                    for err_type in ['X','Z']:
                        fault_id_list.append(self.qubit_ops_dict[qubit_id]\
                                .sq_fault_dict[(err_unit.instruction_line,
                                                err_type)])
                    self.decomp_dict[fault_unit.fault_id] = fault_id_list
                    explicitness = True
                    break
            if explicitness == False:
                self.nonexplicit_fault_id.append(fault_unit.fault_id)

    def _nonexplicit_decomp(self) -> None:
        for qubit_ops in self.qubit_ops_dict.values():
            qubit_ops._resort()
        for fault_id in self.nonexplicit_fault_id:
            fault_unit = self.dem.fault_id_dict[fault_id]
            dets = self.fault_det_dict[fault_id]
            decompse_check = False
            for err_unit in fault_unit.err_unit_list:
                fault_ids_1 = []
                fault_ids_2 = []
                oper_line = err_unit.instruction_line
                qubit_id_1 = err_unit.error_supp[0]
                qubit_id_2 = err_unit.error_supp[1]
                fault_ids_1 = self.qubit_ops_dict[qubit_id_1]._cp_fault_ids(oper_line)
                fault_ids_2 = self.qubit_ops_dict[qubit_id_2]._cp_fault_ids(oper_line)
                for f1 in fault_ids_1:
                    for f2 in fault_ids_2:
                        dets_1 = self.fault_det_dict[f1]
                        dets_2 = self.fault_det_dict[f2]
                        if dets_1.symmetric_difference(dets_2) == dets:
                            self.decomp_dict[fault_id] = [f1,f2]
                            decompse_check = True
                            break
                    if decompse_check == True:
                        break
                if decompse_check == True:
                    break
            if decompse_check == False:
                # print(fault_id)
                raise ValueError('decomposition failed')
            

    def futher_decompose(self, attached_faults: list[int],
                         unattached_faults: list[int],
                         unattached_decomp:dict[int,list[int]]):
        attached_fault_id_set = frozenset(attached_faults)
        unattached_fault_id_set = frozenset(unattached_faults)
        for fault_id in self.decomp_dict:
            decomp_ids = self.decomp_dict[fault_id]
            decomp_further_set = frozenset([])
            for id in decomp_ids:
                if id in unattached_fault_id_set:
                    decomp_further_set = decomp_further_set.symmetric_difference(
                                frozenset(unattached_decomp[id]))
                else:
                    decomp_further_set = decomp_further_set.symmetric_difference(
                                frozenset([id]))
            for basis_id in decomp_further_set:
                try: 
                    self.basis_dict[basis_id].append(fault_id)
                except:
                    self.basis_dict[basis_id] = [fault_id]
            self.decomp_dict[fault_id] = list(decomp_further_set)

        for fault_id in attached_fault_id_set:
            decomp_ids = self.decomp_dict[fault_id]
            if frozenset(decomp_ids) != frozenset([fault_id]):
                self.decomp_dict[fault_id] = frozenset([fault_id])
                try:
                    self.basis_dict[fault_id].append(fault_id)
                except:
                    self.basis_dict[fault_id] = [fault_id]
        
        # check:
        for decomped_ids in self.decomp_dict.values():
            for id in decomped_ids:
                if id not in attached_fault_id_set:
                    raise ValueError('conversion failure')
        pass
    

    def gen_decompose_csr(self, attached_faults: list[int]):

        data_list = []
        row_list = []
        col_list = []
        num_rows = len(attached_faults)
        num_cols = len(self.decomp_dict.values())

        for i in range(num_cols):
            decomp_ids = self.decomp_dict[i]
            for id in decomp_ids:
                data_list.append(1)
                row_list.append(attached_faults.index(id))
                col_list.append(i)
        
        ret_csr = csr_array((data_list,(row_list,col_list)),
                            shape=(num_rows,num_cols),dtype=int)
        
        return ret_csr


            
                            
        
            




