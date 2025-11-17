import os
import csv

# Set the path to the parent directory containing directories 1 to 34
parent_dir = "/home/obiwan/Projects/Stuff/Mofrad/5-VirtualScreening/Separate-sdf"

# Initialize a list to store information for CSV
csv_data = []

# List to store coconut IDs for the first 50,000 ligands
coconut_ids_list = []

# Traverse directories 1 to 34
for dir_num in range(1, 35):
    dir_path = os.path.join(parent_dir, str(dir_num))
    
    # Check if directory exists
    if not os.path.isdir(dir_path):
        print(f"Directory {dir_path} does not exist.")
        continue
    
    # Traverse all files in the current directory
    for file_name in os.listdir(dir_path):
        if file_name.endswith("-out.pdbqt"):  # Check for Vina output file
            coconut_id = file_name.split("-")[0]  # Extract coconut ID
            
            file_path = os.path.join(dir_path, file_name)
            lowest_energy = None
            best_model = None
            
            try:
                with open(file_path, 'r') as file:
                    model_number = None
                    for line in file:
                        if line.startswith("MODEL"):
                            model_number = int(line.split()[1])  # Get model number
                        elif line.startswith("REMARK VINA RESULT:"):
                            energy = float(line.split()[3])  # Get affinity energy
                            
                            # Check if this model has the lowest energy so far
                            if lowest_energy is None or energy < lowest_energy:
                                lowest_energy = energy
                                best_model = model_number
            except FileNotFoundError:
                print(f"File {file_path} not found.")
                continue
            
            # Add the results to csv_data
            if lowest_energy is not None and best_model is not None:
                csv_data.append([coconut_id, best_model, lowest_energy])

# Sort the csv_data list by the lowest_energy (third element)
csv_data.sort(key=lambda x: x[2])  # Sort by lowest_energy (third element)

# Extract the sorted coconut IDs for the first 50,000 entries
sorted_coconut_ids = [row[0] for row in csv_data]

# Write the results to a CSV file
csv_output_path = os.path.join(parent_dir, "docking_results.csv")
with open(csv_output_path, 'w', newline='') as csvfile:
    csv_writer = csv.writer(csvfile)
    csv_writer.writerow(["Coconut_ID", "Best_Model", "Lowest_Energy"])
    csv_writer.writerows(csv_data)

# Write the first 30,000 sorted coconut IDs to a text file
ligand_output_path = os.path.join(parent_dir, "first_30000_ligands.txt")
with open(ligand_output_path, 'w') as ligand_file:
    for coconut_id in sorted_coconut_ids[:min(30000, len(sorted_coconut_ids))]:
        ligand_file.write(coconut_id + "\n")

print("CSV and sorted ligand text files generated successfully.")

