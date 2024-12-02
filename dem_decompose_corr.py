import pymatching.matching
from dem_parsor import *
from hypergraph import *
import pymatching
from joblib import Parallel, delayed
from scipy.sparse import csr_array, csc_array, vstack, hstack
from timeit import default_timer as timer
from dem_util import *




# A meta class for constructing decoding (hyper-)graphs and 
# carrying out various decoding tasks. (Admittedly, the name of the class
# is slightly confusing.)

# Currently supported decoders:
# 1. plain decoder,
# 2. VTB decoder,
# 3. VTB-PR decoder,
# 4. VTB-FR decoder,
# 5. FR decoder.
# No. 1,3,4 decoders are used in this work. 

# Decoding is run in parrallel using the joblib package
# The number of processes is set to be 5 by default, and 
# should be reset according to your own computing machine 
# via the num_batch parameter in each decoding method. 

# with current setup, it would also be straight forward 
# to implement a BP-assisted decoder. 


class Decomposer():

    def __init__(self,
                 ft_circuit: FT_Circuit_dynamical_phase,
                 dem: DEM) -> None:
        self.ft_circuit = ft_circuit
        self.dem = dem
        self.hypergraph = HyperGraph()
        self.x_det_id_list = []
        self.z_det_id_list = []
        self.fund_det_z_fault_ids = []
        self.fund_det_x_fault_ids = []
        self.composite_fault_ids = []
        self.composite_weight_dict: dict[int,list[int]] = {}
        self.attached_z_fault_ids = []
        self.unattached_z_fault_ids = []
        self.unattached_z_decom_dict = {}
        self.ill_decomposable_fault_ids = []
        self.soft_hpe_core: list[int] = []
        self.sampler = self.ft_circuit.get_stim_circuit().compile_detector_sampler()
        pass


    def _hpg_initialize(self) -> None:
        '''
        add all detectors as vertices to the hypergraph
        '''
        id_temp = 0
        for det_id in self.ft_circuit.circ.detector_dict.keys():
            if det_id != id_temp:
                raise ValueError('det_id order is strange')
            detector =  self.ft_circuit.circ.detector_dict[id_temp]
            det_space_coord = (detector.spacetime_coord[0],
                               detector.spacetime_coord[1])
            # figure out the type of detector
            try:
                self.ft_circuit.code.x_check_network[pos_str_conversion(
                                                    det_space_coord)]
                det_type = 'X'
                self.hypergraph.append_vertex_set([detector.detector_id],
                                                  [det_type],
                                                  [detector.spacetime_coord])
                self.x_det_id_list.append(det_id)
            except:
                try:
                    self.ft_circuit.code.z_check_network[pos_str_conversion(
                                                        det_space_coord)]
                    det_type = 'Z'
                    self.hypergraph.append_vertex_set([detector.detector_id],
                                                      [det_type],
                                                      [detector.spacetime_coord])
                    self.z_det_id_list.append(det_id)
                except:
                    raise ValueError('det space coord not on stab ancilla')

            id_temp += 1

        pass


    def _fund_det_z_errs(self) -> None:
        '''
        All Z-type or X-type fault units (error locations) that triggers Z detectors only
        and triggers less or equal than two Z detectors.
        '''
        for fault_unit in self.dem.x_faults + self.dem.z_faults:
            if fault_unit.fault_id in self.hypergraph.id_hyperedge_dict.keys():
                continue
            x_det_check = 0
            for det_id in fault_unit.flipped_detector_id_list:
                if det_id in self.x_det_id_list:
                    x_det_check = 1
            if x_det_check == 0 \
                and len(fault_unit.flipped_detector_id_list) <= 2:
                self.fund_det_z_fault_ids.append(fault_unit.fault_id)
            

    def _fund_det_x_errs(self) -> None:
        '''
        All Z-type or X-type that triggers X detectors only
        and triggers less or equal than two X detectors

        '''
        for fault_unit in self.dem.x_faults + self.dem.z_faults:
            if fault_unit.fault_id in self.hypergraph.id_hyperedge_dict.keys():
                continue
            x_det_check = 0
            for det_id in fault_unit.flipped_detector_id_list:
                if det_id in self.z_det_id_list:
                    x_det_check = 1
            if x_det_check == 0 \
                and len(fault_unit.flipped_detector_id_list) <= 2:
                self.hypergraph.append_hyperedge_set([fault_unit.fault_id],
                                [fault_unit.err_probability],
                                [fault_unit.logical_flag])
                conn_det_ids = fault_unit.flipped_detector_id_list
                hpe_ids = [fault_unit.fault_id for i in conn_det_ids]
                self.hypergraph.attach(conn_det_ids,
                                       hpe_ids)
                self.hypergraph.hyperedge_csc_gen(fault_unit.fault_id)
                self.fund_det_x_fault_ids.append(fault_unit.fault_id)


    def _det_z_errs_weight_cate(self) -> None:
        '''
        Sort out all fault units that flips <= 2 Z detectors.
        And sort them by weight (the number of flipped detectors).
        '''
        for fault_unit in self.dem.x_faults + self.dem.z_faults:
            if fault_unit.fault_id in self.fund_det_x_fault_ids \
                or fault_unit.fault_id in self.fund_det_z_fault_ids \
                or fault_unit.fault_id in self.composite_fault_ids:
                continue
            z_det_count = 0
            for det_id in fault_unit.flipped_detector_id_list:
                if det_id in self.z_det_id_list:
                    z_det_count += 1
            if z_det_count == 0:
                continue
            # All registered composite has at least one z det 
            # and one x det triggered.
            # Also, the number of triggered z dets should not be larger
            # than 2. (So that matching can still be used.)
            if z_det_count <= 2:
                self.composite_fault_ids.append(fault_unit.fault_id)
                det_weight = len(fault_unit.flipped_detector_id_list)
                if det_weight in self.composite_weight_dict.keys():
                    self.composite_weight_dict[det_weight].append(
                                                    fault_unit.fault_id)
                else:
                    self.composite_weight_dict[det_weight] = [
                                                    fault_unit.fault_id]
            
    

    def _fault_id_detector_type_spacetime(self, fault_id: int) -> tuple[
                                    list[int],str,list[tuple[float]]]:
        fault_unit = self.dem.fault_id_dict[fault_id]
        det_id_list = fault_unit.flipped_detector_id_list
        det_type_list = []
        det_spacetime_list = []
        for det_id in det_id_list:
            if det_id in self.x_det_id_list:
                det_type_list.append('X')
            elif det_id in self.z_det_id_list:
                det_type_list.append('Z')
            
            det_vet = self.hypergraph.id_vertex_dict[det_id]
            det_spacetime_list.append(det_vet.spacetime_coords)
        
        return (det_id_list, det_type_list, det_spacetime_list)


    

    def _check_z_foot_print(self) -> bool:
        """
        set up a dict where keys are partial_{Z} E
        (the set of z dets triggered by E) and the corresponding
        value is the list of errors {E'} with 
        partial_{Z} {E'} = partial_{Z} E. 
        """

        z_foot_print: list[frozenset] = []
        self.z_fp_fault_id_dict: dict[frozenset,list[int]] = {}
        for fault_id in self.fund_det_z_fault_ids:
            fault_unit = self.dem.fault_id_dict[fault_id]
            det_ids_frozen = frozenset(fault_unit.flipped_detector_id_list)
            z_foot_print.append(det_ids_frozen)
            if det_ids_frozen in self.z_fp_fault_id_dict.keys():
                self.z_fp_fault_id_dict[det_ids_frozen].append(fault_id)
            else:
                self.z_fp_fault_id_dict[det_ids_frozen] = [fault_id]
        
        for fault_id in self.composite_fault_ids:
            fault_unit = self.dem.fault_id_dict[fault_id]
            det_ids = fault_unit.flipped_detector_id_list
            det_z_ids = []
            for det_id in det_ids:
                if det_id in self.z_det_id_list:
                    det_z_ids.append(det_id)
            z_foot_print.append(frozenset(det_z_ids))
            frozen_det_z = frozenset(det_z_ids)
            if frozen_det_z in self.z_fp_fault_id_dict.keys():
                self.z_fp_fault_id_dict[frozen_det_z].append(fault_id)
            else:
                self.z_fp_fault_id_dict[frozen_det_z] = [fault_id]


        return True


    
    def _attach_hyperedges(self) -> None:
        for det_id_frozen in self.z_fp_fault_id_dict.keys():
            fault_ids = self.z_fp_fault_id_dict[det_id_frozen]
            if len(fault_ids) == 1:
                if fault_ids[0] in self.hypergraph.id_hyperedge_dict.keys():
                    continue
                else:
                    fault_unit = self.dem.fault_id_dict[fault_ids[0]]
                    det_id_list = fault_unit.flipped_detector_id_list
                    self.hypergraph.append_hyperedge_set([fault_ids[0]],
                                        [fault_unit.err_probability],
                                        [fault_unit.logical_flag])
                    hpe_ids = [fault_ids[0] for i in det_id_list]
                    self.hypergraph.attach(det_id_list,
                                                hpe_ids)
                    self.attached_z_fault_ids.append(fault_ids[0])
                    self.hypergraph.hyperedge_csc_gen(fault_unit.fault_id)
            else:
                weight_list = []
                for fault_id in fault_ids:
                    if fault_id in self.hypergraph.id_hyperedge_dict.keys():
                        raise ValueError('what?')
                    fault_unit = self.dem.fault_id_dict[fault_id]
                    weight_list.append(self._probability_to_weight(
                                        fault_unit.err_probability))
                min_index = weight_list.index(min(weight_list))

                fault_id_min = fault_ids[min_index]
                self.soft_hpe_core.append(fault_id_min)
                fault_unit = self.dem.fault_id_dict[fault_id_min]
                det_id_list = fault_unit.flipped_detector_id_list
                self.hypergraph.append_hyperedge_set([fault_id_min],
                                    [fault_unit.err_probability],
                                    [fault_unit.logical_flag])
                hpe_ids = [fault_id_min for i in det_id_list]
                self.hypergraph.attach(det_id_list,
                                            hpe_ids)
                self.attached_z_fault_ids.append(fault_id_min)
                self.hypergraph.hyperedge_csc_gen(fault_unit.fault_id)

                for fault_id in fault_ids:
                    if fault_id == fault_id_min:
                        continue
                    fault_unit = self.dem.fault_id_dict[fault_id]
                    det_id_list_2 = fault_unit.flipped_detector_id_list
                    det_id_diff_list = []
                    for det_id in det_id_list:
                        if det_id not in det_id_list_2:
                            det_id_diff_list.append(det_id)
                    for det_id in det_id_list_2:
                        if det_id not in det_id_list:
                            det_id_diff_list.append(det_id)
                    is_there_hpe, hpe_id_diff = self.hypergraph.find_hpe(
                                                det_id_diff_list)
                    if is_there_hpe == False:
                        print(fault_id_min, fault_id)
                        raise ValueError('cannot be decomposed to a ' 
                                         + 'fundamental error')
                    else:
                        hpe = self.hypergraph.id_hyperedge_dict[fault_id_min]
                        hpe.soft_hpe_id.append(hpe_id_diff)
                        hpe.soft_fault_id.append(fault_id)
                        hpe_diff = self.hypergraph.id_hyperedge_dict[hpe_id_diff]
                        hpe_diff.soft_fault_id_core.append(fault_id_min)
                        self.unattached_z_fault_ids.append(fault_id)
                        self.unattached_z_decom_dict[fault_id]=[fault_id_min,hpe_id_diff]


    def decoding_batch_plain(self, num_shots: int, num_batch=5) -> None:

        detection_events_full, observable_flips_full = self.sampler\
                                    .sample(num_shots, separate_observables=True)


        detection_events_list = np.array_split(detection_events_full,num_batch)
        observable_flips_list = np.array_split(observable_flips_full,num_batch)

        edge_prob = np.transpose(self.edge_decompose_csr \
                                @ np.reshape(self.decomposer_dem.error_channel,(-1,1)))
        edge_prob = np.where(edge_prob>1e-12,edge_prob,1e-12)
        z_prob = np.transpose(self.z_reweight_csr \
                                @ np.transpose(edge_prob))
        z_prob = np.where(z_prob>1e-12,z_prob,1e-12)
        z_prob = np.where(z_prob>0.5,0.5,z_prob)
        z_weights = np.log(1-z_prob) - np.log(z_prob)
        z_weights = np.where(z_weights>0,z_weights,0)
        z_weights = np.reshape(z_weights,(-1))
        x_prob = np.transpose(self.x_reweight_csr \
                                @ np.transpose(edge_prob))
        x_prob = np.where(x_prob>1e-12,x_prob,1e-12)
        x_prob = np.where(x_prob>0.5,0.5,x_prob)
        x_weights = np.log(1-x_prob) - np.log(x_prob)
        x_weights = np.where(x_weights>0,x_weights,0)
        x_weights = np.reshape(x_weights,(-1))
        
        res_list = Parallel(n_jobs=num_batch)(delayed(
                                    plain_decoder)(
                                        detection_events,
                                        observable_flips,
                                        self.z_ext_csc,
                                        self.logic_flip_z_mat,
                                        self.Hz_total,
                                        self.x_ext_csc,
                                        self.logic_flip_x_mat,
                                        self.Hz, z_weights,
                                        self.Hx, x_weights) 
                    for detection_events, observable_flips in 
                                zip(detection_events_list,
                                    observable_flips_list))  
        
        tot_err_shot = 0
        for res in res_list:
            tot_err_shot += res[1]

        return (tot_err_shot/num_shots, tot_err_shot, num_shots)
    


    def decoding_batch_VTB(self, num_shots: int, num_batch=5) -> None:

        detection_events_full, observable_flips_full = self.sampler\
                                    .sample(num_shots, separate_observables=True)


        detection_events_list = np.array_split(detection_events_full,num_batch)
        observable_flips_list = np.array_split(observable_flips_full,num_batch)

        edge_prob = np.transpose(self.edge_decompose_csr \
                                @ np.reshape(self.decomposer_dem.error_channel,(-1,1)))
        edge_prob = np.where(edge_prob>1e-12,edge_prob,1e-12)
        z_prob = np.transpose(self.z_reweight_csr \
                                @ np.transpose(edge_prob))
        z_prob = np.where(z_prob>1e-12,z_prob,1e-12)
        z_prob = np.where(z_prob>0.5,0.5,z_prob)
        z_weights = np.log(1-z_prob) - np.log(z_prob)
        z_weights = np.where(z_weights>0,z_weights,0)
        z_weights = np.reshape(z_weights,(-1))
        x_prob = np.transpose(self.x_reweight_csr \
                                @ np.transpose(edge_prob))
        x_prob = np.where(x_prob>1e-12,x_prob,1e-12)
        x_prob = np.where(x_prob>0.5,0.5,x_prob)
        x_weights = np.log(1-x_prob) - np.log(x_prob)
        x_weights = np.where(x_weights>0,x_weights,0)
        x_weights = np.reshape(x_weights,(-1))


        
        res_list = Parallel(n_jobs=num_batch)(delayed(
                                    VTB_decoder)(
                                        detection_events,
                                        observable_flips,
                                        self.z_ext_csc,
                                        self.logic_flip_z_mat,
                                        self.Hz_total,
                                        self.x_ext_csc,
                                        self.logic_flip_x_mat,
                                        self.Hz, z_weights,
                                        self.Hx, x_weights) 
                    for detection_events, observable_flips in 
                                zip(detection_events_list,
                                    observable_flips_list))  
        
        tot_err_shot = 0
        for res in res_list:
            tot_err_shot += res[1]

        return (tot_err_shot/num_shots, tot_err_shot, num_shots)
    


    def decoding_batch_VTB_PR(self, num_shots: int, num_batch=5) -> None:

        detection_events_full, observable_flips_full = self.sampler\
                                    .sample(num_shots, separate_observables=True)


        detection_events_list = np.array_split(detection_events_full,num_batch)
        observable_flips_list = np.array_split(observable_flips_full,num_batch)


        edge_prob = np.transpose(self.edge_decompose_csr \
                                @ np.reshape(self.decomposer_dem.error_channel,(-1,1)))
        edge_prob = np.where(edge_prob>1e-12,edge_prob,1e-12)
        z_prob = np.transpose(self.z_reweight_csr \
                                @ np.transpose(edge_prob))
        z_prob = np.where(z_prob>1e-12,z_prob,1e-12)
        z_prob = np.where(z_prob>0.5,0.5,z_prob)
        z_weights = np.log(1-z_prob) - np.log(z_prob)
        z_weights = np.where(z_weights>0,z_weights,0)
        z_weights = np.reshape(z_weights,(-1))
        x_prob = np.transpose(self.x_reweight_csr \
                                @ np.transpose(edge_prob))
        x_prob = np.where(x_prob>1e-12,x_prob,1e-12)
        x_prob = np.where(x_prob>0.5,0.5,x_prob)
        x_weights = np.log(1-x_prob) - np.log(x_prob)
        x_weights = np.where(x_weights>0,x_weights,0)
        x_weights = np.reshape(x_weights,(-1))
        
        res_list = Parallel(n_jobs=num_batch)(delayed(
                                    VTB_PR_decoder)(
                                        detection_events,
                                        observable_flips,
                                        self.z_ext_csc,
                                        self.logic_flip_z_mat,
                                        self.Hz_total,
                                        self.x_ext_csc,
                                        self.logic_flip_x_mat,
                                        self.Hz, z_weights,
                                        self.Hx, x_weights,
                                        self.full_err_reweight_from_z_csr,
                                        self.edge_decompose_csr,
                                        self.z_reweight_csr,
                                        self.x_reweight_csr,
                                        self.decomposer_dem.error_channel) 
                    for detection_events, observable_flips in 
                                zip(detection_events_list,
                                    observable_flips_list))  
        
        tot_err_shot = 0
        for res in res_list:
            tot_err_shot += res[1]

        return (tot_err_shot/num_shots, tot_err_shot, num_shots)




    def decoding_batch_VTB_FR(self, num_shots: int, num_batch = 5):

        detection_events_full, observable_flips_full = self.sampler\
                                    .sample(num_shots, separate_observables=True)


        detection_events_list = np.array_split(detection_events_full,num_batch)
        observable_flips_list = np.array_split(observable_flips_full,num_batch)


        edge_prob = np.transpose(self.edge_decompose_csr \
                                @ np.reshape(self.decomposer_dem.error_channel,(-1,1)))
        edge_prob = np.where(edge_prob>1e-12,edge_prob,1e-12)
        z_prob = np.transpose(self.z_reweight_csr \
                                @ np.transpose(edge_prob))
        z_prob = np.where(z_prob>1e-12,z_prob,1e-12)
        z_prob = np.where(z_prob>0.5,0.5,z_prob)
        z_weights = np.log(1-z_prob) - np.log(z_prob)
        z_weights = np.where(z_weights>0,z_weights,0)
        z_weights = np.reshape(z_weights,(-1))
        x_prob = np.transpose(self.x_reweight_csr \
                                @ np.transpose(edge_prob))
        x_prob = np.where(x_prob>1e-12,x_prob,1e-12)
        x_prob = np.where(x_prob>0.5,0.5,x_prob)
        x_weights = np.log(1-x_prob) - np.log(x_prob)
        x_weights = np.where(x_weights>0,x_weights,0)
        x_weights = np.reshape(x_weights,(-1))
        
        res_list = Parallel(n_jobs=num_batch)(delayed(
                                    VTB_FR_decoder)(
                                        detection_events,
                                        observable_flips,
                                        self.z_ext_csc,
                                        self.logic_flip_z_mat,
                                        self.Hz_total,
                                        self.x_ext_csc,
                                        self.logic_flip_x_mat,
                                        self.Hz, z_weights,
                                        self.Hx, x_weights,
                                        self.full_err_reweight_from_z_csr,
                                        self.edge_decompose_csr,
                                        self.z_reweight_csr,
                                        self.x_reweight_csr,
                                        self.decomposer_dem.error_channel) 
                    for detection_events, observable_flips in 
                                zip(detection_events_list,
                                    observable_flips_list))  
        
        tot_err_shot = 0
        for res in res_list:
            tot_err_shot += res[1]

        return (tot_err_shot/num_shots, tot_err_shot, num_shots)
    

    def decoding_batch_FR(self, num_shots: int, num_batch = 5):

        detection_events_full, observable_flips_full = self.sampler\
                                    .sample(num_shots, separate_observables=True)


        detection_events_list = np.array_split(detection_events_full,num_batch)
        observable_flips_list = np.array_split(observable_flips_full,num_batch)


        edge_prob = np.transpose(self.edge_decompose_csr \
                                @ np.reshape(self.decomposer_dem.error_channel,(-1,1)))
        edge_prob = np.where(edge_prob>1e-12,edge_prob,1e-12)
        z_prob = np.transpose(self.z_reweight_csr \
                                @ np.transpose(edge_prob))
        z_prob = np.where(z_prob>1e-12,z_prob,1e-12)
        z_prob = np.where(z_prob>0.5,0.5,z_prob)
        z_weights = np.log(1-z_prob) - np.log(z_prob)
        z_weights = np.where(z_weights>0,z_weights,0)
        z_weights = np.reshape(z_weights,(-1))
        x_prob = np.transpose(self.x_reweight_csr \
                                @ np.transpose(edge_prob))
        x_prob = np.where(x_prob>1e-12,x_prob,1e-12)
        x_prob = np.where(x_prob>0.5,0.5,x_prob)
        x_weights = np.log(1-x_prob) - np.log(x_prob)
        x_weights = np.where(x_weights>0,x_weights,0)
        x_weights = np.reshape(x_weights,(-1))
        
        res_list = Parallel(n_jobs=num_batch)(delayed(
                                    FR_decoder)(
                                        detection_events,
                                        observable_flips,
                                        self.z_ext_csc,
                                        self.logic_flip_z_mat,
                                        self.Hz_total,
                                        self.x_ext_csc,
                                        self.logic_flip_x_mat,
                                        self.Hz, z_weights,
                                        self.Hx, x_weights,
                                        self.full_err_reweight_from_z_csr,
                                        self.edge_decompose_csr,
                                        self.z_reweight_csr,
                                        self.x_reweight_csr,
                                        self.decomposer_dem.error_channel) 
                    for detection_events, observable_flips in 
                                zip(detection_events_list,
                                    observable_flips_list))  
        
        tot_err_shot = 0
        for res in res_list:
            tot_err_shot += res[1]

        return (tot_err_shot/num_shots, tot_err_shot, num_shots)





    def z_matcher_compilation(self) -> None:
        '''
        Construction of the Gz decoding graph for the plain decoder
        '''
        # first extract the z-detection events
        data_list = []
        row_list = []
        col_list = []
        z_det_id_list_sorted = sorted(self.z_det_id_list)
        num_rows = len(z_det_id_list_sorted)
        num_cols = len(self.hypergraph.vertices)
        for i in range(num_rows):
            data_list.append(1)
            row_list.append(i)
            col_list.append(z_det_id_list_sorted[i])
        self.z_ext_csc = csr_array((data_list, (row_list,col_list)),
                               shape=(num_rows,num_cols),
                               dtype=np.int64)
        
        # extract list of fault_ids with >= 1 z-dets
        fault_id_list = sorted(list(self.hypergraph.id_hyperedge_dict.keys()))
        fault_z_id_list = []
        
        for fault_id in fault_id_list:
            fault_unit = self.dem.fault_id_dict[fault_id]
            det_ids = fault_unit.flipped_detector_id_list
            z_det_counter = 0
            for det_id in det_ids:
                if det_id in self.z_det_id_list:
                    z_det_counter += 1
            if z_det_counter > 0:
                fault_z_id_list.append(fault_id)
            pass


        # construct z-syndromes triggered by z-edges
        # construct weight for z-edges (This is not used later,
        # as a better way to calculate the weights requires accounting
        # for correlated errors.) 
        # construct logical flips by z-edges

        data_list = []
        row_list = []
        col_list = []
        num_rows = len(z_det_id_list_sorted)
        num_cols = len(fault_z_id_list)
        self.fault_z_id_list = copy.deepcopy(fault_z_id_list)
        self.z_weight_list = []
        logic_list = []

        for fault_id, col in zip(fault_z_id_list,
                                 range(num_cols)):
            fault_unit = self.dem.fault_id_dict[fault_id]
            self.z_weight_list.append(self._probability_to_weight(
                                fault_unit.err_probability))
            if fault_unit.logical_flag == True:
                logic_list.append(1)
            elif fault_unit.logical_flag == False:
                logic_list.append(0)
            else:
                raise ValueError('not even wrong')
            
            det_ids = fault_unit.flipped_detector_id_list
            z_det_ids = []
            for det_id in det_ids:
                if det_id in self.z_det_id_list:
                    z_det_ids.append(det_id)
            
            for det_id in z_det_ids:
                row = z_det_id_list_sorted.index(det_id)
                data_list.append(1)
                row_list.append(row)
                col_list.append(col)
        
        self.Hz = csr_array((data_list, (row_list,col_list)),
                        shape=(num_rows,num_cols),
                        dtype=np.int64)

        
        self.logic_flip_z_mat = np.reshape(np.array(logic_list,
                                    dtype=np.int64),(1,-1))
        
        # create z-matcher

        self.z_matcher = pymatching.Matching(self.Hz, 
                            weights=np.array(self.z_weight_list))
        


        # construct syndromes triggered by z-edges 

        data_list = []
        row_list = []
        col_list = []
        num_rows = len(self.hypergraph.vertices)
        num_cols = len(fault_z_id_list)

        for fault_id, col in zip(fault_z_id_list,
                                 range(num_cols)):
            fault_unit = self.dem.fault_id_dict[fault_id]
            det_ids = fault_unit.flipped_detector_id_list
            for det_id in det_ids:
                row = det_id
                data_list.append(1)
                row_list.append(row)
                col_list.append(col)
        
        self.Hz_total = csr_array((data_list, (row_list,col_list)),
                        shape=(num_rows,num_cols),
                        dtype=np.int64)
        


        pass








    def z_matcher_compilation_VTB(self) -> None:
        '''
        Construction of the Gz decoding graph for 
        the VTB-assisted decoders
        '''
        # max z det time coord:
        time_coord_list = []
        for det_id in self.z_det_id_list:
            det_vert = self.hypergraph.id_vertex_dict[det_id]
            time_coord_list.append(det_vert.spacetime_coords[-1])

        max_time = max(time_coord_list)
        min_time = min(time_coord_list)


        # first extract the z-detection events
        data_list = []
        row_list = []
        col_list = []
        z_det_id_list_sorted = sorted(self.z_det_id_list)
        num_rows = len(z_det_id_list_sorted)
        num_cols = len(self.hypergraph.vertices)
        for i in range(num_rows):
            data_list.append(1)
            row_list.append(i)
            col_list.append(z_det_id_list_sorted[i])
        self.z_ext_csc = csr_array((data_list, (row_list,col_list)),
                               shape=(num_rows,num_cols),
                               dtype=np.int64)
        
        # extract list of fault_ids with >= 1 z-dets
        # upper boundary connected to a virtual det
        # lower boundary connected to another virtual det

        fault_id_list = sorted(list(self.hypergraph.id_hyperedge_dict.keys()))
        fault_z_id_list = []
        fault_z_top_bd = []
        fault_z_bottom_bd = []
        
        for fault_id in fault_id_list:
            fault_unit = self.dem.fault_id_dict[fault_id]
            det_ids = fault_unit.flipped_detector_id_list
            z_det_counter = 0
            z_det_spacetime_coords = []
            top_bd_check = False
            lower_bd_check = False
            for det_id in det_ids:
                if det_id in self.z_det_id_list:
                    z_det_counter += 1
                    z_det_vert = self.hypergraph.id_vertex_dict[det_id]
                    z_det_spacetime_coords.append(z_det_vert.spacetime_coords)
            if z_det_counter > 0:
                fault_z_id_list.append(fault_id)
            if z_det_counter == 1:
                time_coord = z_det_spacetime_coords[0][-1]
                x_coord = z_det_spacetime_coords[0][0]
                if time_coord == min_time:
                    top_bd_check = True
                    fault_z_top_bd.append(True)
                elif time_coord == max_time:
                    lower_bd_check = True
                    fault_z_bottom_bd.append(True)


            if z_det_counter > 0 and top_bd_check == False:
                fault_z_top_bd.append(False)
            if z_det_counter > 0 and lower_bd_check == False:
                fault_z_bottom_bd.append(False)
        
                    
            pass





        # construct z-matching H for det-z-faults
        # construct weight
        # construct logical flips

        data_list = []
        row_list = []
        col_list = []
        num_rows = len(z_det_id_list_sorted) + 2
        num_cols = len(fault_z_id_list)
        self.fault_z_id_list = copy.deepcopy(fault_z_id_list)
        self.z_weight_list = []
        logic_list = []

        for fault_id, col, top_bd, bottom_bd in zip(fault_z_id_list,
                                 range(num_cols),
                                 fault_z_top_bd,
                                 fault_z_bottom_bd):
            fault_unit = self.dem.fault_id_dict[fault_id]
            self.z_weight_list.append(self._probability_to_weight(
                                fault_unit.err_probability))
            if fault_unit.logical_flag == True:
                logic_list.append(1)
            elif fault_unit.logical_flag == False:
                logic_list.append(0)
            else:
                raise ValueError('not even wrong')
            
            det_ids = fault_unit.flipped_detector_id_list
            z_det_ids = []
            for det_id in det_ids:
                if det_id in self.z_det_id_list:
                    z_det_ids.append(det_id)
            
            for det_id in z_det_ids:
                row = z_det_id_list_sorted.index(det_id)
                data_list.append(1)
                row_list.append(row)
                col_list.append(col)
            
            if top_bd == True:
                row = num_rows-2
                data_list.append(1)
                row_list.append(row)
                col_list.append(col)
            
            if bottom_bd == True:
                row = num_rows-1
                data_list.append(1)
                row_list.append(row)
                col_list.append(col)
            

        
        self.Hz = csr_array((data_list, (row_list,col_list)),
                        shape=(num_rows,num_cols),
                        dtype=np.int64)
        
        self.logic_flip_z_mat = np.reshape(np.array(logic_list,
                                    dtype=np.int64),(1,-1))
        
        # create z-matcher

        self.z_matcher = pymatching.Matching(self.Hz, 
                            weights=np.array(self.z_weight_list))

        

        # construct total H for det-z-faults 

        data_list = []
        row_list = []
        col_list = []
        num_rows = len(self.hypergraph.vertices)
        num_cols = len(fault_z_id_list)

        for fault_id, col in zip(fault_z_id_list,
                                 range(num_cols)):
            fault_unit = self.dem.fault_id_dict[fault_id]
            det_ids = fault_unit.flipped_detector_id_list
            for det_id in det_ids:
                row = det_id
                data_list.append(1)
                row_list.append(row)
                col_list.append(col)
        
        self.Hz_total = csr_array((data_list, (row_list,col_list)),
                        shape=(num_rows,num_cols),
                        dtype=np.int64)
        


        pass


    

    def reweight_helper(self):
        

        stim_circuit = self.ft_circuit.get_stim_circuit()
        dem = DEM(stim_circuit)
        decomposer_dem = Decompose_DEM(dem,len(self.hypergraph.vertices))
        decomposer_dem.futher_decompose(self.attached_z_fault_ids \
                                    + self.fund_det_x_fault_ids,
                                    self.unattached_z_fault_ids,
                                    self.unattached_z_decom_dict)
        
        attached_edges = self.fault_z_id_list + self.fault_x_id_list
        self.decomposer_dem = decomposer_dem
        self.edge_decompose_csr = copy.deepcopy(self.decomposer_dem.gen_decompose_csr(attached_edges))


        data_list = []
        row_list = []
        col_list = []
        num_rows = len(self.fault_z_id_list)
        num_cols = len(attached_edges)

        for z_fault_id, i in zip(self.fault_z_id_list,range(num_rows)):
            data_list.append(1)
            row_list.append(i)
            col_list.append(attached_edges.index(z_fault_id))
        
        self.z_reweight_csr = csr_array((data_list,(row_list,col_list)),
                                        shape=(num_rows,num_cols),dtype=int)
        

        # use z-(hyper)edges to infer general faults:

        data_list = []
        row_list = []
        col_list = []
        num_rows = len(self.decomposer_dem.decomp_dict.keys())  
        num_cols = len(self.fault_z_id_list)

        for z_fault_id, col in zip(self.fault_z_id_list,
                                   range(num_cols)):
            conn_fault_ids = self.decomposer_dem.basis_dict[z_fault_id]
            reweight_prob = 1/len(conn_fault_ids)
            for i in conn_fault_ids:
                data_list.append(reweight_prob)
                row_list.append(i)
                col_list.append(col)
        
        self.full_err_reweight_from_z_csr = csr_array(
                            (data_list,(row_list,col_list)),
                            shape=(num_rows,num_cols))
        

        



        data_list = []
        row_list = []
        col_list = []
        num_rows = len(self.fault_x_id_list)
        num_cols = len(attached_edges)

        for x_fault_id, i in zip(self.fault_x_id_list,range(num_rows)):
            data_list.append(1)
            row_list.append(i)
            col_list.append(attached_edges.index(x_fault_id))
        
        self.x_reweight_csr = csr_array((data_list,(row_list,col_list)),
                                        shape=(num_rows,num_cols),dtype=int)
        



     

        pass




        







    def x_matcher_compilation(self) -> None:

        # first extract the x-detection events

        data_list = []
        row_list = []
        col_list = []
        x_det_id_list_sorted = sorted(self.x_det_id_list)
        num_rows = len(x_det_id_list_sorted)
        num_cols = len(self.hypergraph.vertices)
        for i in range(num_rows):
            data_list.append(1)
            row_list.append(i)
            col_list.append(x_det_id_list_sorted[i])
        self.x_ext_csc = csr_array((data_list, (row_list,col_list)),
                               shape=(num_rows,num_cols),
                               dtype=np.int64)
        
        # create Hx, weight list and logic flips
        # here only consider the fundamental faults that
        # triggers X detectors

        data_list = []
        row_list = []
        col_list = []
        x_det_id_list_sorted = sorted(self.x_det_id_list)
        num_rows = len(x_det_id_list_sorted)
        num_cols = len(self.fund_det_x_fault_ids)
        self.fault_x_id_list = copy.deepcopy(self.fund_det_x_fault_ids)

        self.x_weight_list = []
        logic_flip_list = []

        for fault_id, col in zip(self.fund_det_x_fault_ids,
                                 range(num_cols)):
            fault_unit = self.dem.fault_id_dict[fault_id]
            self.x_weight_list.append(self._probability_to_weight(
                                fault_unit.err_probability))
            if fault_unit.logical_flag == True:
                logic_flip_list.append(1)
            elif fault_unit.logical_flag == False:
                logic_flip_list.append(0)
            else:
                raise ValueError('not even wrong')
            
            det_ids = fault_unit.flipped_detector_id_list
            for det_id in det_ids:
                row = x_det_id_list_sorted.index(det_id)
                data_list.append(1)
                row_list.append(row)
                col_list.append(col)
            
        self.Hx = csr_array((data_list, (row_list,col_list)),
                             shape=(num_rows,num_cols),
                             dtype=np.int64)
        
        self.logic_flip_x_mat = np.reshape(np.array(logic_flip_list,
                                            dtype=np.int64),
                                    (1,-1))


        # create X-matcher

        self.x_matcher = pymatching.Matching(self.Hx,
                                weights=np.array(self.x_weight_list))

        pass


    def _probability_to_weight(self, probability: float) -> float:
        if probability < 0:
            raise ValueError('negative probability')
        ret_weight = np.log((1-probability)/probability)
        if ret_weight < 0:
            ret_weight = 0
        
        return ret_weight
    


    def _update_matchers(self, dem: DEM, stim_circuit: stim.Circuit):
        fault_id_list = sorted(list(self.hypergraph.id_hyperedge_dict.keys()))
        fault_z_id_list = []
        
        for fault_id in fault_id_list:
            fault_unit = self.dem.fault_id_dict[fault_id]
            det_ids = fault_unit.flipped_detector_id_list
            z_det_counter = 0
            for det_id in det_ids:
                if det_id in self.z_det_id_list:
                    z_det_counter += 1
            if z_det_counter > 0:
                fault_z_id_list.append(fault_id)
            pass

        self.z_weight_list = []

        for fault_id in fault_z_id_list:
            fault_unit = dem.fault_id_dict[fault_id]
            self.z_weight_list.append(self._probability_to_weight(
                                fault_unit.err_probability))

        self.x_weight_list = []
        for fault_id in self.fund_det_x_fault_ids:
            fault_unit = dem.fault_id_dict[fault_id]
            self.x_weight_list.append(self._probability_to_weight(
                                fault_unit.err_probability))
            

        decomposer_dem = Decompose_DEM(dem,len(self.hypergraph.vertices))
        decomposer_dem.futher_decompose(self.attached_z_fault_ids \
                                    + self.fund_det_x_fault_ids,
                                    self.unattached_z_fault_ids,
                                    self.unattached_z_decom_dict)
        
        attached_edges = self.fault_z_id_list + self.fault_x_id_list
        self.decomposer_dem = decomposer_dem

        self.z_matcher = pymatching.Matching(self.Hz, 
                            weights=np.array(self.z_weight_list))
        self.x_matcher = pymatching.Matching(self.Hx,
                                weights=np.array(self.x_weight_list))
        

        self.sampler = stim_circuit.compile_detector_sampler()
        
        pass

        







def dem_generator(distance, noise, num_pad,
                  filename = 'dem_gen.stim') -> tuple[DEM,
                                                      stim.Circuit]:

    rotated = Rotated_surface_code(distance)

    ft_circuit = FT_Circuit_dynamical_phase(rotated,
                                            filename,noise)

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

    return (dem, stim_circuit)






def plain_decoder(detection_events,
                    observable_flips,
                    z_ext_csr,logic_flip_z_mat,Hz_total,
                    x_ext_csr, logic_flip_x_mat,
                    Hz, z_weight_list,
                    Hx, x_weight_list
                    ) -> tuple[float, int, int]:
    obs_logic_arr = np.reshape(observable_flips,(-1))*1
    num_shots = len(obs_logic_arr)

    detection_events_csr = csr_array(detection_events)

    z_matcher = pymatching.Matching(Hz, 
                        weights=np.array(z_weight_list))

    z_detection_events = (z_ext_csr \
                            @ detection_events_csr.transpose()).transpose()
    det_z_decoded_faults = z_matcher.decode_batch(z_detection_events.toarray())
    det_z_dec_logic_flips = np.reshape(logic_flip_z_mat \
                                @ np.transpose(det_z_decoded_faults),
                                (-1)) % 2

    det_z_decoded_faults_csr = csr_array(det_z_decoded_faults)

    feedback_det_events_csr = (Hz_total @ 
                det_z_decoded_faults_csr.transpose()).transpose()
    
    res_det_events_csr = detection_events_csr + feedback_det_events_csr
    res_det_events_csr.data = res_det_events_csr.data % 2
    
    x_matcher = pymatching.Matching(Hx,
                            weights=np.array(x_weight_list))

    x_detection_events = (x_ext_csr \
                            @ res_det_events_csr.transpose()).transpose()
    det_x_decoded_faults = x_matcher.decode_batch(x_detection_events.toarray())
    det_x_dec_logic_flips = np.reshape(logic_flip_x_mat \
                                @ np.transpose(det_x_decoded_faults),
                                (-1)) % 2
    
    combined_logic_flips = (det_z_dec_logic_flips \
                                + det_x_dec_logic_flips) % 2
    err_shots = np.sum((combined_logic_flips-obs_logic_arr)%2)
    err_rate = err_shots/num_shots

    return (err_rate, err_shots, num_shots)








def VTB_decoder(detection_events,
                    observable_flips,
                    z_ext_csr,logic_flip_z_mat,Hz_total,
                    x_ext_csr, logic_flip_x_mat,
                    Hz, z_weight_arr,
                    Hx, x_weight_arr
                    ) -> tuple[float, int, int,
                               np.ndarray[int],
                               csc_array,csr_array]:
    obs_logic_arr = np.reshape(observable_flips,(-1))*1
    num_shots = len(obs_logic_arr)


    detection_events_csr = csr_array(detection_events)

    x_matcher = pymatching.Matching(Hx,
                            weights=np.array(x_weight_arr))

    z_detection_events = (z_ext_csr \
                            @ detection_events_csr.transpose()).transpose()

    z_matcher = pymatching.Matching(Hz, 
                    weights=np.reshape(z_weight_arr,(-1)))


    
    fw_list = []
    logic_list = []
    for i in [0,1]:
        if i == 0:
            z_detection_events_temp = hstack((z_detection_events,
                                    np.zeros((num_shots,
                                            1),dtype=np.int64)),
                                            dtype=np.int64)
        else:
            z_detection_events_temp = hstack((z_detection_events,
                                    np.ones((num_shots,
                                            1),dtype=np.int64)),
                                            dtype=np.int64)
        for j in [0,1]:
            if j == 0:
                z_detection_events_bd = hstack((z_detection_events_temp,
                                    np.zeros((num_shots,
                                            1),dtype=np.int64)),
                                            dtype=np.int64)
            else:
                z_detection_events_bd = hstack((z_detection_events_temp,
                                    np.ones((num_shots,
                                            1),dtype=np.int64)),
                                            dtype=np.int64)
            
            det_z_decoded_faults = z_matcher.decode_batch(z_detection_events_bd.toarray())
            det_z_dec_logic_flips = np.reshape(logic_flip_z_mat \
                                        @ np.transpose(det_z_decoded_faults),
                                        (-1)) % 2
            det_z_decoded_faults_csr = csr_array(det_z_decoded_faults)

            feedback_det_events_csr = (Hz_total @ 
                        det_z_decoded_faults_csr.transpose()).transpose()
            
            res_det_events_csr = detection_events_csr + feedback_det_events_csr
            res_det_events_csr.data = res_det_events_csr.data % 2
            
            

            x_detection_events = (x_ext_csr \
                                    @ res_det_events_csr.transpose()).transpose()
            det_x_decoded_faults = x_matcher.decode_batch(x_detection_events.toarray())
            det_x_dec_logic_flips = np.reshape(logic_flip_x_mat \
                                        @ np.transpose(det_x_decoded_faults),
                                        (-1)) % 2
            
            z_faults_weight = np.reshape(np.reshape(z_weight_arr,(1,-1)) \
                                @ np.transpose(det_z_decoded_faults),(-1))
            x_faults_weight = np.reshape(np.reshape(x_weight_arr,(1,-1)) \
                                @ np.transpose(det_x_decoded_faults),(-1))
            
        

            combined_logic_flips = (det_z_dec_logic_flips \
                                        + det_x_dec_logic_flips) % 2
            
            
            fault_weight = z_faults_weight + x_faults_weight
            fw_list.append(fault_weight)
            logic_list.append(combined_logic_flips)


    fw_arr = np.transpose(np.vstack(fw_list))
    logic_arr = np.transpose(np.vstack(logic_list))
    final_logic_flips = []
    for fw_row, logic_row in zip(fw_arr,logic_arr):
        
        fw_row_list = list(fw_row)
        fw_min_index = fw_row_list.index(min(fw_row_list))
        final_logic_flips.append(logic_row[fw_min_index])


    final_logic_flips_arr = np.array(final_logic_flips,
                                        dtype = np.int64)


    err_shots = np.sum((final_logic_flips_arr-obs_logic_arr)%2)
    err_rate = err_shots/num_shots

    return (err_rate, err_shots, num_shots)



class reweight_calc():

    def __init__(self,
                 full_err_reweight_csr: csr_array,
                 bp_error_channel: ndarray,
                 edge_decomposer_csr: csr_array,
                 x_reweight_csr: csr_array
                 ):
        self.full_err_reweight_csr = full_err_reweight_csr
        self.bp_error_channel = bp_error_channel
        self.edge_decomposer_csr = edge_decomposer_csr
        self.x_reweight_csr = x_reweight_csr

    def x_weight_cal(self, z_decoded_faults,
                     x_decoded_faults) -> float:
        full_err_reweight = np.transpose(self.full_err_reweight_csr @ 
                            np.reshape(z_decoded_faults,(-1,1)))
        full_err_signal = np.where(full_err_reweight!=0,0,1)
        full_err = np.reshape(self.bp_error_channel,(1,-1)) * full_err_signal
        full_err = full_err + full_err_reweight

        edge_prob = np.transpose(self.edge_decomposer_csr \
                                    @ np.transpose(full_err))
        edge_prob = np.where(edge_prob>1e-12,edge_prob,1e-12)
        x_prob = np.transpose(self.x_reweight_csr \
                                @ np.transpose(edge_prob))
        x_prob = np.where(x_prob>1e-12,x_prob,1e-12)
        x_prob = np.where(x_prob>0.5,0.5,x_prob)
        x_reweights = np.log(1-x_prob) - np.log(x_prob)
        x_reweights = np.reshape(np.where(x_reweights>0,x_reweights,0),(-1))

        x_faults_weight = np.sum(x_reweights \
                                * np.reshape(x_decoded_faults,(-1)))
        return x_faults_weight
    
        




class reweight_calc_full():

    def __init__(self,
                 full_err_reweight_csr: csr_array,
                 bp_error_channel: ndarray,
                 edge_decomposer_csr: csr_array,
                 x_reweight_csr: csr_array,
                 Hx: csr_array,
                 logic_flip_x_mat
                 ):
        self.full_err_reweight_csr = full_err_reweight_csr
        self.bp_error_channel = bp_error_channel
        self.edge_decomposer_csr = edge_decomposer_csr
        self.x_reweight_csr = x_reweight_csr
        self.logic_flip_x_mat = logic_flip_x_mat
        self.Hx = Hx

    def x_weight_cal(self, z_decoded_faults,
                     x_decoded_faults) -> float:
        full_err_reweight = np.transpose(self.full_err_reweight_csr @ 
                            np.reshape(z_decoded_faults,(-1,1)))
        full_err_signal = np.where(full_err_reweight!=0,0,1)
        full_err = np.reshape(self.bp_error_channel,(1,-1)) * full_err_signal
        full_err = full_err + full_err_reweight

        edge_prob = np.transpose(self.edge_decomposer_csr \
                                    @ np.transpose(full_err))
        edge_prob = np.where(edge_prob>1e-12,edge_prob,1e-12)
        x_prob = np.transpose(self.x_reweight_csr \
                                @ np.transpose(edge_prob))
        x_prob = np.where(x_prob>1e-12,x_prob,1e-12)
        x_prob = np.where(x_prob>0.5,0.5,x_prob)
        x_reweights = np.log(1-x_prob) - np.log(x_prob)
        x_reweights = np.reshape(np.where(x_reweights>0,x_reweights,0),(-1))

        x_faults_weight = np.sum(x_reweights \
                                * np.reshape(x_decoded_faults,(-1)))
        return x_faults_weight
    
    def x_decoding_res(self, z_decoded_faults,
                       x_detection_events):
        full_err_reweight = np.transpose(self.full_err_reweight_csr @ 
                            np.reshape(z_decoded_faults,(-1,1)))
        full_err_signal = np.where(full_err_reweight!=0,0,1)
        full_err = np.reshape(self.bp_error_channel,(1,-1)) * full_err_signal
        full_err = full_err + full_err_reweight

        edge_prob = np.transpose(self.edge_decomposer_csr \
                                    @ np.transpose(full_err))
        edge_prob = np.where(edge_prob>1e-12,edge_prob,1e-12)
        x_prob = np.transpose(self.x_reweight_csr \
                                @ np.transpose(edge_prob))
        x_prob = np.where(x_prob>1e-12,x_prob,1e-12)
        x_prob = np.where(x_prob>0.5,0.5,x_prob)
        x_reweights = np.log(1-x_prob) - np.log(x_prob)
        x_reweights = np.reshape(np.where(x_reweights>0,x_reweights,0),(-1))
        x_matcher = pymatching.Matching(self.Hx,weights=x_reweights)
        x_faults = x_matcher.decode_batch(x_detection_events.toarray())
        x_faults_weight = np.sum(x_reweights \
                                 * np.reshape(x_faults,(-1)))
        logic_flip = np.reshape(self.logic_flip_x_mat \
                                        @ np.transpose(x_faults),
                                        (-1)) % 2
        return (logic_flip[0],x_faults_weight)



def VTB_PR_decoder(detection_events,
                    observable_flips,
                    z_ext_csr,logic_flip_z_mat,Hz_total,
                    x_ext_csr, logic_flip_x_mat,
                    Hz, z_weight_list,
                    Hx, x_weight_list,
                    full_err_reweight_csr,
                    edge_decomposer_csr,
                    z_reweight_csr,
                    x_reweight_csr,
                    bp_error_channel
                    ) -> tuple[float, int, int,
                               np.ndarray[int],
                               csc_array,csr_array]:
    obs_logic_arr = np.reshape(observable_flips,(-1))*1
    num_shots = len(obs_logic_arr)
    detection_events_csr = csr_array(detection_events)

    x_matcher = pymatching.Matching(Hx,
                            weights=np.array(x_weight_list))

    z_detection_events = (z_ext_csr \
                            @ detection_events_csr.transpose()).transpose()
    
    rwc = reweight_calc(full_err_reweight_csr,
                        bp_error_channel,
                        edge_decomposer_csr,
                        x_reweight_csr)

    z_matcher = pymatching.Matching(Hz, 
                    weights=np.reshape(z_weight_list,(-1)))
    
    fw_list = []
    logic_list = []
    for i in [0,1]:
        if i == 0:
            z_detection_events_temp = hstack((z_detection_events,
                                    np.zeros((num_shots,
                                            1),dtype=np.int64)),
                                            dtype=np.int64)
        else:
            z_detection_events_temp = hstack((z_detection_events,
                                    np.ones((num_shots,
                                            1),dtype=np.int64)),
                                            dtype=np.int64)
        for j in [0,1]:
            if j == 0:
                z_detection_events_bd = hstack((z_detection_events_temp,
                                    np.zeros((num_shots,
                                            1),dtype=np.int64)),
                                            dtype=np.int64)
            else:
                z_detection_events_bd = hstack((z_detection_events_temp,
                                    np.ones((num_shots,
                                            1),dtype=np.int64)),
                                            dtype=np.int64)
            
            det_z_decoded_faults = z_matcher.decode_batch(z_detection_events_bd.toarray())
            det_z_dec_logic_flips = np.reshape(logic_flip_z_mat \
                                        @ np.transpose(det_z_decoded_faults),
                                        (-1)) % 2
            det_z_decoded_faults_csr = csr_array(det_z_decoded_faults)

            feedback_det_events_csr = (Hz_total @ 
                        det_z_decoded_faults_csr.transpose()).transpose()
            
            res_det_events_csr = detection_events_csr + feedback_det_events_csr
            res_det_events_csr.data = res_det_events_csr.data % 2
            
            

            x_detection_events = (x_ext_csr \
                                    @ res_det_events_csr.transpose()).transpose()
            det_x_decoded_faults = x_matcher.decode_batch(x_detection_events.toarray())
            det_x_dec_logic_flips = np.reshape(logic_flip_x_mat \
                                        @ np.transpose(det_x_decoded_faults),
                                        (-1)) % 2
            
            z_faults_weight = np.reshape(np.reshape(z_weight_list,(1,-1)) \
                                @ np.transpose(det_z_decoded_faults),(-1))
            
        

            combined_logic_flips = (det_z_dec_logic_flips \
                                        + det_x_dec_logic_flips) % 2
            
            x_faults_weight = [rwc.x_weight_cal(z_decoded_faults,
                                                x_decoded_faults)
                    for z_decoded_faults, x_decoded_faults in 
                    zip(det_z_decoded_faults,det_x_decoded_faults)]
            
            fault_weight = z_faults_weight + np.array(x_faults_weight)
            fw_list.append(fault_weight)
            logic_list.append(combined_logic_flips)


    fw_arr = np.transpose(np.vstack(fw_list))
    logic_arr = np.transpose(np.vstack(logic_list))
            
            

    


    final_logic_flips = []


    for fw_row, logic_row in zip(fw_arr,logic_arr):
        
        fw_row_list = list(fw_row)
        fw_min_index = fw_row_list.index(min(fw_row_list))
        final_logic_flips.append(logic_row[fw_min_index])


    final_logic_flips_arr = np.array(final_logic_flips,
                                        dtype = np.int64)


    err_shots = np.sum((final_logic_flips_arr-obs_logic_arr)%2)
    err_rate = err_shots/num_shots

    return (err_rate, err_shots, num_shots)




def VTB_FR_decoder(detection_events,
                    observable_flips,
                    z_ext_csr,logic_flip_z_mat,Hz_total,
                    x_ext_csr, logic_flip_x_mat,
                    Hz, z_weight_list,
                    Hx, x_weight_list,
                    full_err_reweight_csr,
                    edge_decomposer_csr,
                    z_reweight_csr,
                    x_reweight_csr,
                    bp_error_channel
                    ) -> tuple[float, int, int,
                               np.ndarray[int],
                               csc_array,csr_array]:
    obs_logic_arr = np.reshape(observable_flips,(-1))*1
    num_shots = len(obs_logic_arr)

    detection_events_csr = csr_array(detection_events)

    x_matcher = pymatching.Matching(Hx,
                            weights=np.array(x_weight_list))

    z_detection_events = (z_ext_csr \
                            @ detection_events_csr.transpose()).transpose()
    
    rwc = reweight_calc_full(full_err_reweight_csr,
                        bp_error_channel,
                        edge_decomposer_csr,
                        x_reweight_csr,Hx,logic_flip_x_mat)
    
    z_matcher = pymatching.Matching(Hz, 
                    weights=np.reshape(z_weight_list,(-1)))


    
    fw_list = []
    logic_list = []
    for i in [0,1]:
        if i == 0:
            z_detection_events_temp = hstack((z_detection_events,
                                    np.zeros((num_shots,
                                            1),dtype=np.int64)),
                                            dtype=np.int64)
        else:
            z_detection_events_temp = hstack((z_detection_events,
                                    np.ones((num_shots,
                                            1),dtype=np.int64)),
                                            dtype=np.int64)
        for j in [0,1]:
            if j == 0:
                z_detection_events_bd = hstack((z_detection_events_temp,
                                    np.zeros((num_shots,
                                            1),dtype=np.int64)),
                                            dtype=np.int64)
            else:
                z_detection_events_bd = hstack((z_detection_events_temp,
                                    np.ones((num_shots,
                                            1),dtype=np.int64)),
                                            dtype=np.int64)
            
            det_z_decoded_faults = z_matcher.decode_batch(z_detection_events_bd.toarray())
            det_z_dec_logic_flips = np.reshape(logic_flip_z_mat \
                                        @ np.transpose(det_z_decoded_faults),
                                        (-1)) % 2
            det_z_decoded_faults_csr = csr_array(det_z_decoded_faults)

            feedback_det_events_csr = (Hz_total @ 
                        det_z_decoded_faults_csr.transpose()).transpose()
            
            res_det_events_csr = detection_events_csr + feedback_det_events_csr
            res_det_events_csr.data = res_det_events_csr.data % 2
            
            

            x_detection_events = (x_ext_csr \
                                    @ res_det_events_csr.transpose()).transpose()
            x_res_list = [rwc.x_decoding_res(z_decoded_faults, x_detection_events)
                          for z_decoded_faults, x_detection_events in 
                          zip(det_z_decoded_faults, x_detection_events)]
            
            det_x_dec_logic_flips = [item[0] for item in x_res_list]
            x_faults_weight = [item[1] for item in x_res_list]
            det_x_dec_logic_flips = np.array(det_x_dec_logic_flips)
            x_faults_weight = np.array(x_faults_weight)

            
            z_faults_weight = np.reshape(np.reshape(z_weight_list,(1,-1)) \
                                @ np.transpose(det_z_decoded_faults),(-1))
            


            combined_logic_flips = (det_z_dec_logic_flips \
                                        + det_x_dec_logic_flips) % 2
            
            
            fault_weight = z_faults_weight + np.array(x_faults_weight)
            fw_list.append(fault_weight)
            logic_list.append(combined_logic_flips)

    fw_arr = np.transpose(np.vstack(fw_list))
    logic_arr = np.transpose(np.vstack(logic_list))
    final_logic_flips = []

    for fw_row, logic_row in zip(fw_arr,logic_arr):
        
        fw_row_list = list(fw_row)
        fw_min_index = fw_row_list.index(min(fw_row_list))
        final_logic_flips.append(logic_row[fw_min_index])

    final_logic_flips_arr = np.array(final_logic_flips,
                                        dtype = np.int64)
    err_shots = np.sum((final_logic_flips_arr-obs_logic_arr)%2)
    err_rate = err_shots/num_shots

    return (err_rate, err_shots, num_shots)






def FR_decoder(detection_events,
                    observable_flips,
                    z_ext_csr,logic_flip_z_mat,Hz_total,
                    x_ext_csr, logic_flip_x_mat,
                    Hz, z_weight_list,
                    Hx, x_weight_list,
                    full_err_reweight_csr,
                    edge_decomposer_csr,
                    z_reweight_csr,
                    x_reweight_csr,
                    bp_error_channel
                    ) -> tuple[float, int, int,
                               np.ndarray[int],
                               csc_array,csr_array]:
    obs_logic_arr = np.reshape(observable_flips,(-1))*1
    num_shots = len(obs_logic_arr)

    detection_events_csr = csr_array(detection_events)
    z_detection_events = (z_ext_csr \
                            @ detection_events_csr.transpose()).transpose()
    
    rwc = reweight_calc_full(full_err_reweight_csr,
                        bp_error_channel,
                        edge_decomposer_csr,
                        x_reweight_csr,Hx,logic_flip_x_mat)


    z_matcher = pymatching.Matching(Hz, 
                    weights=np.reshape(z_weight_list,(-1)))
  
    det_z_decoded_faults = z_matcher.decode_batch(z_detection_events.toarray())
    det_z_dec_logic_flips = np.reshape(logic_flip_z_mat \
                                @ np.transpose(det_z_decoded_faults),
                                (-1)) % 2
    det_z_decoded_faults_csr = csr_array(det_z_decoded_faults)

    feedback_det_events_csr = (Hz_total @ 
                det_z_decoded_faults_csr.transpose()).transpose()
    
    res_det_events_csr = detection_events_csr + feedback_det_events_csr
    res_det_events_csr.data = res_det_events_csr.data % 2
    
    x_detection_events = (x_ext_csr \
                            @ res_det_events_csr.transpose()).transpose()
    x_res_list = [rwc.x_decoding_res(z_decoded_faults, x_detection_events)
                    for z_decoded_faults, x_detection_events in 
                    zip(det_z_decoded_faults, x_detection_events)]
    
    det_x_dec_logic_flips = [item[0] for item in x_res_list]
    det_x_dec_logic_flips = np.array(det_x_dec_logic_flips)
    combined_logic_flips = (det_z_dec_logic_flips \
                                + det_x_dec_logic_flips) % 2
            
    err_shots = np.sum((combined_logic_flips-obs_logic_arr)%2)
    err_rate = err_shots/num_shots

    return (err_rate, err_shots, num_shots)








