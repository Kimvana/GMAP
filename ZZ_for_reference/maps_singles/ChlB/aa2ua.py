import pandas as pd

file = pd.read_csv("chlb.csv", sep = ",", keep_default_na= False)  # Read in csv file with difference charges

columns_to_sum = file.columns.to_list()[2:]  # From which column are we summing the values

atom_names = file["atom_name"].to_list()  # Grab the name of the atoms as a list

indices_to_delete = []  # Create blank list where indices of 'absorbed atoms' will be saved into

for index, name in enumerate(atom_names):  # enumerate() function adds a counter to each item in an iterable
    correct_name = None 
    if name[0] == "H":
        cont = 1
        while True:
            if atom_names[index - cont][0] == 'C' or atom_names[index - cont][0] == 'O':
                correct_name = atom_names[index - cont]
                break
            else:
                cont += 1
        
        if correct_name:
            for col in columns_to_sum:
                value = file.iloc[index][col]
                file.loc[file["atom_name"] == correct_name,col] += value
            indices_to_delete.append(index)

file.drop(labels=indices_to_delete,inplace=True, axis=0)

file.to_csv("chlb_result.csv", index=False)