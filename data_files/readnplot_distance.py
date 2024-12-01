from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import sinter

class noise_single_result:

    def __init__(self, noise_prob) -> None:
        self.noise_prob = noise_prob
        self.num_shots = 0
        self.num_hits = 0
        self.max_likelihood = 1000

    def append_result(self, num_shots, num_hits):
        self.num_shots = self.num_shots + num_shots
        self.num_hits = self.num_hits + num_hits
    
    def err_rate(self) -> float:
        sinter_res = sinter.fit_binomial(num_shots=self.num_shots,
                            num_hits=self.num_hits,
                            max_likelihood_factor=self.max_likelihood)
        return sinter_res.best
    
    def err_lower(self) -> float:
        sinter_res = sinter.fit_binomial(num_shots=self.num_shots,
                            num_hits=self.num_hits,
                            max_likelihood_factor=self.max_likelihood)
        return sinter_res.low
    
    def err_higher(self) -> float:
        sinter_res = sinter.fit_binomial(num_shots=self.num_shots,
                            num_hits=self.num_hits,
                            max_likelihood_factor=self.max_likelihood)
        return sinter_res.high
        
    
    


class dict_noise_result:

    def __init__(self, distance) -> None:
        self.distance = distance
        self.noise_list = []
        self.noise_dict: dict[float,noise_single_result] = {}
    
    def add_result(self, noise_prob, num_shots, num_hits):
        if noise_prob not in self.noise_list:
            self.noise_list.append(noise_prob)
            self.noise_dict[noise_prob] = noise_single_result(noise_prob)
        self.noise_dict[noise_prob].append_result(num_shots,num_hits)

noise_list = []




X_memory_dict: dict[int,dict_noise_result] = {}

with open("X_memory.txt") as file:
    file_line = [line.rstrip() for line in file]
    for line_str in file_line:
        line_parse = line_str.split(' ')
        distance = int(line_parse[2])
        noise_prob = float(line_parse[4])
        num_shots = int(line_parse[6])
        num_hits = int(line_parse[8])
        if distance not in X_memory_dict.keys():
            X_memory_dict[distance] = dict_noise_result(distance)
        X_memory_dict[distance].add_result(noise_prob,
                                           num_shots, num_hits)





S_2_mid_full_dict: dict[int,dict_noise_result] = {}

with open("S_2_mid_full.txt") as file:
    file_line = [line.rstrip() for line in file]
    for line_str in file_line:
        line_parse = line_str.split(' ')
        distance = int(line_parse[2])
        noise_prob = float(line_parse[4])
        num_shots = int(line_parse[6])
        num_hits = int(line_parse[8])
        if distance not in S_2_mid_full_dict.keys():
            S_2_mid_full_dict[distance] = dict_noise_result(distance)
        S_2_mid_full_dict[distance].add_result(noise_prob,
                                           num_shots, num_hits)




S_2_mid_2_dict: dict[int,dict_noise_result] = {}

with open("S_2_mid_2.txt") as file:
    file_line = [line.rstrip() for line in file]
    for line_str in file_line:
        line_parse = line_str.split(' ')
        distance = int(line_parse[2])
        noise_prob = float(line_parse[4])
        num_shots = int(line_parse[6])
        num_hits = int(line_parse[8])
        if distance not in S_2_mid_2_dict.keys():
            S_2_mid_2_dict[distance] = dict_noise_result(distance)
        S_2_mid_2_dict[distance].add_result(noise_prob,
                                           num_shots, num_hits)






color_dict = {
    'Blue': '#3366CC',
    'SkyBlue': '#0099C6',
    'Teal': '#22AA99',
    'Red': '#DC3912',
    'FireBrick': '#B82E2E',
    'Pink': '#DD4477',
    'Orange': '#FF9900',
    'DeepOrange': '#E67300',
    'Green': '#109618',
    'LightGreen': '#66AA00',
    'Purple': '#990099'
}


blue_gradient = ['#99c4ff',
                 '#5ca1ff',
                 '#1f7eff',
                 '#005fe0',
                 '#0045a3',
                 '#002c66']

orange_gradient = ['#ffd199',
                   '#ffb55c',
                   '#ff9a1f',
                   '#e07c00',
                   '#a35900',
                   '#663800']




plt.figure(figsize=(8, 6))


single_noise_distance_list = [3,5,7,9,11,13]

plt.yscale('log')
# plt.xscale('log')
plt.grid(True,which='major',axis='y',color='gray')
plt.grid(True,which='minor',axis='y',color='#cccccc')
plt.xticks([3,5,7,9,11,13])
plt.xlabel(r'Code distance $d$')
plt.ylabel('Logical error rate (per shot)')







distance_list = sorted(X_memory_dict.keys())


noise = 0.001
logic_err_list = []
logic_err_floor_list = []
logic_err_ceil_list = []
distance_list_temp = []


for distance in distance_list:
    
    try:
        single_noise_result = X_memory_dict[distance].noise_dict[noise]
        logic_err_list.append(single_noise_result.err_rate())
        logic_err_floor_list.append(single_noise_result.err_lower())
        logic_err_ceil_list.append(single_noise_result.err_higher())
        distance_list_temp.append(distance)
        err_lower = np.array(logic_err_list) - np.array(logic_err_floor_list)
        err_upper = np.array(logic_err_ceil_list) - np.array(logic_err_list)
        plt_err = np.array([err_lower,
                            err_upper])
    except:
        continue



plt.plot(distance_list_temp,logic_err_list,marker='s',color=blue_gradient[2]
                ,markersize=7)


plt.fill_between(distance_list_temp, logic_err_floor_list,
                 logic_err_ceil_list,
                 color=blue_gradient[2],
                 alpha=0.4)






distance_list = sorted(S_2_mid_full_dict.keys())


noise = 0.001
logic_err_list = []
logic_err_floor_list = []
logic_err_ceil_list = []
distance_list_temp = []


for distance in distance_list:
    
    try:
        single_noise_result = S_2_mid_full_dict[distance].noise_dict[noise]
        logic_err_list.append(single_noise_result.err_rate())
        logic_err_floor_list.append(single_noise_result.err_lower())
        logic_err_ceil_list.append(single_noise_result.err_higher())
        distance_list_temp.append(distance)
        err_lower = np.array(logic_err_list) - np.array(logic_err_floor_list)
        err_upper = np.array(logic_err_ceil_list) - np.array(logic_err_list)
        plt_err = np.array([err_lower,
                            err_upper])
    except:
        continue



plt.plot(distance_list_temp,logic_err_list,marker='s',color=orange_gradient[2]
                ,markersize=7)


plt.fill_between(distance_list_temp, logic_err_floor_list,
                 logic_err_ceil_list,
                 color=orange_gradient[2],
                 alpha=0.4)






distance_list = sorted(S_2_mid_2_dict.keys())


noise = 0.001
logic_err_list = []
logic_err_floor_list = []
logic_err_ceil_list = []
distance_list_temp = []

for distance in distance_list:
    

    try:
        single_noise_result = S_2_mid_2_dict[distance].noise_dict[noise]
        logic_err_list.append(single_noise_result.err_rate())
        logic_err_floor_list.append(single_noise_result.err_lower())
        logic_err_ceil_list.append(single_noise_result.err_higher())
        distance_list_temp.append(distance)
        err_lower = np.array(logic_err_list) - np.array(logic_err_floor_list)
        err_upper = np.array(logic_err_ceil_list) - np.array(logic_err_list)
        plt_err = np.array([err_lower,
                            err_upper])
    except:
        continue




plt.plot(distance_list_temp,logic_err_list,marker='s',color='#ce1fff'
                ,markersize=7)

plt.fill_between(distance_list_temp, logic_err_floor_list,
                 logic_err_ceil_list,
                 color='#ce1fff',
                 alpha=0.4)







distance_list = sorted(S_2_mid_full_dict.keys())


noise = 0.003
logic_err_list = []
logic_err_floor_list = []
logic_err_ceil_list = []
distance_list_temp = []


for distance in distance_list:
    
    try:
        single_noise_result = S_2_mid_full_dict[distance].noise_dict[noise]
        logic_err_list.append(single_noise_result.err_rate())
        logic_err_floor_list.append(single_noise_result.err_lower())
        logic_err_ceil_list.append(single_noise_result.err_higher())
        distance_list_temp.append(distance)
        err_lower = np.array(logic_err_list) - np.array(logic_err_floor_list)
        err_upper = np.array(logic_err_ceil_list) - np.array(logic_err_list)
        plt_err = np.array([err_lower,
                            err_upper])
    except:
        continue



plt.plot(distance_list_temp,logic_err_list,marker='o',color=orange_gradient[2]
                ,markersize=7)


plt.fill_between(distance_list_temp, logic_err_floor_list,
                 logic_err_ceil_list,
                 color=orange_gradient[2],
                 alpha=0.4)










distance_list = sorted(X_memory_dict.keys())


noise = 0.003
logic_err_list = []
logic_err_floor_list = []
logic_err_ceil_list = []
distance_list_temp = []


for distance in distance_list:
    
    try:
        single_noise_result = X_memory_dict[distance].noise_dict[noise]
        logic_err_list.append(single_noise_result.err_rate())
        logic_err_floor_list.append(single_noise_result.err_lower())
        logic_err_ceil_list.append(single_noise_result.err_higher())
        distance_list_temp.append(distance)
        err_lower = np.array(logic_err_list) - np.array(logic_err_floor_list)
        err_upper = np.array(logic_err_ceil_list) - np.array(logic_err_list)
        plt_err = np.array([err_lower,
                            err_upper])
    except:
        continue


plt.plot(distance_list_temp,logic_err_list,marker='o',color=blue_gradient[2]
                ,markersize=7)

plt.fill_between(distance_list_temp, logic_err_floor_list,
                 logic_err_ceil_list,
                 color=blue_gradient[2],
                 alpha=0.4)














distance_list = sorted(S_2_mid_2_dict.keys())


noise = 0.003
logic_err_list = []
logic_err_floor_list = []
logic_err_ceil_list = []
distance_list_temp = []

for distance in distance_list:
    

    try:
        single_noise_result = S_2_mid_2_dict[distance].noise_dict[noise]
        logic_err_list.append(single_noise_result.err_rate())
        logic_err_floor_list.append(single_noise_result.err_lower())
        logic_err_ceil_list.append(single_noise_result.err_higher())
        distance_list_temp.append(distance)
        err_lower = np.array(logic_err_list) - np.array(logic_err_floor_list)
        err_upper = np.array(logic_err_ceil_list) - np.array(logic_err_list)
        plt_err = np.array([err_lower,
                            err_upper])
    except:
        continue




plt.plot(distance_list_temp,logic_err_list,marker='o',color='#ce1fff'
                ,markersize=7)

plt.fill_between(distance_list_temp, logic_err_floor_list,
                 logic_err_ceil_list,
                 color='#ce1fff',
                 alpha=0.4)


    





plt.rcParams.update({'font.size': 16})
plt.savefig('logic_err_distance_scaling.pdf',bbox_inches='tight')