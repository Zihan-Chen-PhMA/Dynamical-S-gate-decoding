from __future__ import annotations
from typing import Union
from scipy.sparse import csc_matrix
from abc import ABC
import numpy as np

class Vertex():

    def __init__(self, id: int, vtype: str,
                 spacetime_coords: tuple[float]) -> None:
        self.id = id    # Detector id 
        self.vtype: str = vtype    # Detector type
        self.spacetime_coords = spacetime_coords
        self.connected_hyperedges_id = []
        self.subgraph_id: Union[int,None] = None
        pass

class HyperEdge():

    def __init__(self, id: int, probability: Union[float,None],
                 logical_flag: bool) -> None:
        self.id = id    # Fault id 
        self.probability = probability
        self.logical_flag = logical_flag
        self.connected_vertices_id = []
        self.err_type: Union[str,None] = None
        self.hpe_id_composition: Union[list[int],None] = None
        self.connected_id_csc_row: Union[csc_matrix, None] = None
        self.subgraph_connected_id_csc_row: Union[csc_matrix, None] = None
        self.soft_hpe_id: list[int] = []
        self.soft_fault_id: list[int] = []
        self.soft_fault_id_core: list[int] = []
        pass

class HyperGraph(ABC):

    def __init__(self) -> None:
        self.vertices: list[Vertex] = []
        self.id_vertex_dict: dict[int, Vertex] = {}
        self.hyperedges: list[HyperEdge] = []
        self.id_hyperedge_dict: dict[int, HyperEdge] = {}
        pass

    def append_vertex_set(self, vertex_id_list: list[int], 
                          vertex_type_list: list[str],
                          spacetime_coords_list: list[tuple[float]]):
        if len(vertex_id_list) != len(vertex_type_list):
            raise ValueError('length mismatch')
        for vertex_id, vertex_type, spacetime_coords\
                                 in zip(vertex_id_list, 
                                        vertex_type_list,
                                        spacetime_coords_list):
            self.id_vertex_dict[vertex_id] = Vertex(vertex_id, vertex_type,
                                                    spacetime_coords)
            self.vertices.append(self.id_vertex_dict[vertex_id])
        pass

    def append_hyperedge_set(self, hyperedge_id_list: list[int],
            hyperedge_probability_list: Union[list[float], None], 
            logical_flag_list: Union[list[bool],None]):
        for hyperedge_id, probability, logical_flag in zip(
                    hyperedge_id_list, 
                    hyperedge_probability_list,
                    logical_flag_list):
            self.id_hyperedge_dict[hyperedge_id] = HyperEdge(
                hyperedge_id, probability, logical_flag)
            self.hyperedges.append(self.id_hyperedge_dict[hyperedge_id])
        pass

    def set_hyperedge_probabiltiy(self, hyperedge_id_list: list[int],
            hyperedge_probability_list: list[float]):
        pass

    def attach(self, vertex_id_list: list[int], 
               hyperedge_id_list: list[int]) -> None:
        if len(vertex_id_list) != len(hyperedge_id_list):
            raise ValueError('length mismatch')
        for vertex_id, hyperedge_id in zip(vertex_id_list, 
                                           hyperedge_id_list):
            if vertex_id not in self.id_vertex_dict.keys():
                raise ValueError('Vertex id not found')
            if hyperedge_id not in self.id_hyperedge_dict.keys():
                raise ValueError('Hyperedge id not found')
            self.id_vertex_dict[vertex_id].connected_hyperedges_id\
                                                .append(hyperedge_id)
            self.id_hyperedge_dict[hyperedge_id].connected_vertices_id\
                                                .append(vertex_id)
        pass


    def hyperedge_csc_gen(self, hyperedge_id: int) -> None:
        hyperedge = self.id_hyperedge_dict[hyperedge_id]
        col_list = hyperedge.connected_vertices_id
        row_list = [0 for i in col_list]
        data_list = [1 for i in col_list]
        hyperedge.connected_id_csc_row = csc_matrix((data_list,
                                        (row_list,col_list)),
                            shape=(1,len(self.vertices)), 
                            dtype= np.int64)
    

    def subgraph_gen(self, vertex_list: list[int]):
        pass

    def decompose_into_hpe_wrong(self, vertex_ids_1: list[int],
                           vertex_ids_2: list[int]) -> tuple[bool,
                                                Union[list[int],None]]:
        '''
        We assume only two-body decomposition
        '''
        conn_hpes_1: list[HyperEdge] = []
        conn_hpe_1_ids = []
        conn_hpes_2: list[HyperEdge] = []
        conn_hpe_2_ids = []
        for vertex_id in vertex_ids_1:
            vertex = self.id_vertex_dict[vertex_id]
            for hpe_id in vertex.connected_hyperedges_id:
                if hpe_id not in conn_hpe_1_ids:
                    conn_hpes_1.append(self.id_hyperedge_dict[hpe_id])
                    conn_hpe_1_ids.append(hpe_id)
        for vertex_id in vertex_ids_2:
            vertex = self.id_vertex_dict[vertex_id]
            for hpe_id in vertex.connected_hyperedges_id:
                if hpe_id not in conn_hpe_2_ids:
                    conn_hpes_2.append(self.id_hyperedge_dict[hpe_id])
                    conn_hpe_2_ids.append(hpe_id)
        
        for hpe_1 in conn_hpes_1:
            for hpe_2 in conn_hpes_2:
                if self.is_zero_mod_2(hpe_1.connected_id_csc_row + 
                                      hpe_2.connected_id_csc_row) is True:
                    return (True, [hpe_1.id, hpe_2.id])
        
        
        return (False, None)
    
    def decompose_into_hpe(self, vertex_ids: list[int]) -> tuple[bool,
                                                Union[list[int],None]]:
        conn_hpes: list[HyperEdge] = []
        conn_hpes_ids = []
        ver_csc_row = self.vset_to_csc(vertex_ids)
        for vertex_id in vertex_ids:
            vertex = self.id_vertex_dict[vertex_id]
            for hpe_id in vertex.connected_hyperedges_id:
                if hpe_id not in conn_hpes_ids:
                    conn_hpes_ids.append(hpe_id)
                    conn_hpes.append(self.id_hyperedge_dict[hpe_id])
        
        for i in range(len(conn_hpes_ids)):
            for j in range(len(conn_hpes_ids)):
                if j > i:
                    hpe_1 = conn_hpes[i]
                    hpe_2 = conn_hpes[j]
                    hpg_matching = hpe_1.connected_id_csc_row \
                                   + hpe_2.connected_id_csc_row \
                                   + ver_csc_row
                    if self.is_zero_mod_2(hpg_matching) is True:
                        return (True, [hpe_1.id, hpe_2.id])
        
        return (False, None)
    
    def find_hpe(self, vertex_ids: list[int]) -> tuple[bool,
                                            Union[int,None]]:
        conn_hpes: list[HyperEdge] = []
        conn_hpes_ids = []
        ver_csc_row = self.vset_to_csc(vertex_ids)
        vertex = self.id_vertex_dict[vertex_ids[0]]
        for hpe_id in vertex.connected_hyperedges_id:
            if hpe_id not in conn_hpes_ids:
                conn_hpes_ids.append(hpe_id)
                conn_hpes.append(self.id_hyperedge_dict[hpe_id])
        
        for i in range(len(conn_hpes_ids)):
            hpe_1 = conn_hpes[i]
            hpg_matching = hpe_1.connected_id_csc_row \
                            + ver_csc_row
            if self.is_zero_mod_2(hpg_matching) is True:
                return (True, hpe_1.id)
        
        return (False, None)


    def is_zero_mod_2(self, csc_mat: csc_matrix) -> bool:
        arr = csc_mat.toarray() % 2
        sum_all = np.sum(arr) 
        if sum_all == 0:
            return True
        else:
            return False
        
    def vset_to_csc(self, vertex_ids: list[int]) -> csc_matrix:
        data_list = [1 for i in vertex_ids]
        row_list = [0 for i in vertex_ids]
        ret_csc = csc_matrix((data_list,(row_list,vertex_ids)),
                            shape=(1,len(self.vertices)), 
                            dtype= np.int64)
        return ret_csc

