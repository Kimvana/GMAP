import numpy as np
import matplotlib.pyplot as plt

# Only works for Vesna's system!

class graphBuilder():
    def __init__(self):

        # matrix_data = graphBuilder.read_files() # read
        # # print([row[0] for row in matrix_data])
        # # print(len(matrix_data[0]) / 121)

        # # Each Nx_list contains a list of lists. The sublists contain the potentials of each xth N atom in each residue. The lists 
        # # correspond to spheresizes, ranging from 10 to 200 a    

        # N1_list = graphBuilder.extract_N(matrix_data, 2)
        # N2_list = graphBuilder.extract_N(matrix_data, 10)
        # N3_list = graphBuilder.extract_N(matrix_data, 16)
        # N4_list = graphBuilder.extract_N(matrix_data, 21)
        # emptyList = [[0] * len(N1_list[0])] * len(N1_list) 

        # # print(emptyList)
        # # print(N1_list)

        # # print(N4_list[-1])

        # to_plot = graphBuilder.potential_plotter(N1_list, N3_list)

        # print(len(to_plot))

        type_list = ["aa_short/subbox_aa_linear_", "ma_short/subbox_ma_linear_", "perres_short/perres_mm_linear_"] # This list contains all the names of folders that contain data
        distance_list = np.arange(0, 10, 1)
        # distance_list = [0, 1, 2, 5, 9]

        # print(distance_list)
        radii = [spheresize for spheresize in range(10, 50)]

        N_pairs = [[2,10],[2,16],[2,21],[10,16],[10,21],[16,21]]

        # Calculates per smoothing distance

        # for atom_of_interest in range(1): 
        #     for folder in type_list:
        #         match folder:
        #             case "aa_short/subbox_aa_linear_":
        #                 method = "aa"
        #             case "ma_short/subbox_ma_linear_":
        #                 method = "ma"
        #             case "perres_short/perres_mm_linear_":
        #                 method = "perres"
                
        #         for pair in N_pairs:
        #             for distance in distance_list:
        #                 matrix_data = graphBuilder.read_files(folder, distance)
        #                 N1_list = graphBuilder.extract_N(matrix_data, pair[0])
        #                 N2_list = graphBuilder.extract_N(matrix_data, pair[1])
        #                 to_plot = graphBuilder.potential_plotter(N1_list, N2_list, True, atom_of_interest)
        #                 plt.plot(radii, to_plot, label=f"{distance}")
        #             plt.legend()
        #             plt.xlabel("Distance (A)")
        #             plt.ylabel("Potential(V)")
        #             plt.savefig(f"molecule_{atom_of_interest}_pair_{pair}_method_{method}")
        #             plt.clf()
        

        # Calculates per method

        for atom_of_interest in range(1):
            for pair in N_pairs:
                for distance in distance_list:
                    for folder in type_list:
                        matrix_data = graphBuilder.read_files(folder, distance)
                        N1_list = graphBuilder.extract_N(matrix_data, pair[0])
                        N2_list = graphBuilder.extract_N(matrix_data, pair[1])
                        to_plot = graphBuilder.potential_plotter(N1_list, N2_list, True, atom_of_interest)
                        plt.plot(radii, to_plot, label=f"{folder}{distance}")
                    plt.legend()
                    plt.xlabel("Distance (A)")
                    plt.ylabel("Potential(V)")
                    plt.savefig(f"atom_{atom_of_interest}_dist_{distance}_pairing_{pair[0]}_{pair[1]}")
                    plt.clf()

        # Graphs the potential vs. smoothing distance

        # for atom_of_interest in range(1):
        #     for pair in N_pairs:
        #         for folder in type_list:
        #             match folder:
        #                 case "aa_short/subbox_aa_linear_":
        #                     method = "aa"
        #                 case "ma_short/subbox_ma_linear_":
        #                     method = "ma"
        #                 case "perres_short/perres_mm_linear_":
        #                     method = "perres"
        #             for spheresize_index in [0, 2, 5, 10, 15, 20, 30]: # Each of these numbers refer to an index in radii
        #                 potential_list = []
        #                 for distance in distance_list:
        #                     matrix_data = graphBuilder.read_files(folder, distance)
        #                     N1_list = graphBuilder.extract_N(matrix_data, pair[0])
        #                     N2_list = graphBuilder.extract_N(matrix_data, pair[1])
        #                     to_plot = graphBuilder.potential_plotter(N1_list, N2_list, True, atom_of_interest)
        #                     potential_list.append(to_plot[spheresize_index])
        #                 plt.plot(distance_list, potential_list, label=f"spheresize {spheresize_index + 10}")
        #             plt.legend()
        #             plt.xlabel("Smoothing Distance (A)")
        #             plt.ylabel("Potential(V)")
        #             plt.savefig(f"molecule_{atom_of_interest}_pairing_{pair[0]}_{pair[1]}_method_{method}")
        #             plt.clf()


        # plt.plot([spheresize for spheresize in range(10, 200)], to_plot, label=input_file_name)
        # to_plot = graphBuilder.potential_plotter(N2_list, N4_list)
        # plt.plot([spheresize for spheresize in range(10, 200)], to_plot, label=f"{input_file_name}!!!" )
        # plt.legend()
        # plt.savefig("plot.png")

        


    def read_files(folder, distance): # Reads the file and turns it into a matrix which is returned
        return np.array(list(map(lambda x: np.array(x.split(" ")[1:-1]), np.array(open(f"{folder}{distance}.txt", "r").readlines()))), dtype=float)

    def extract_N(matrix_data, atom): # turns the text file into a list of lists
        N_list = []
        list_ = 0
        while list_ < len(matrix_data): 
            sublist = 0
            N_list.append([])
            while sublist < 11:
                N_list[list_].append(matrix_data[list_][sublist * 121 + atom]) 
                sublist += 1
            list_ += 1
        return N_list
    
    def potential_plotter(first_list, second_list, average, atom_of_interest = 0):
        # Takes two output lists, averages the differences of the items. Each item on the output list is the average potential difference.

        difference_list = [[first_pot - second_pot for first_pot, second_pot in zip(first_sublist, second_sublist)] for first_sublist, second_sublist in zip(first_list, second_list)]
        if average:
            return list(map(lambda pot_diff_list: sum(pot_diff_list) / len(pot_diff_list) , difference_list)) # This returns the averages the differences
        else:
            return list(map(lambda pot_diff_list: pot_diff_list[atom_of_interest] , difference_list)) #This returns one of the differences
    

if __name__ == "__main__":
    graphBuilder()