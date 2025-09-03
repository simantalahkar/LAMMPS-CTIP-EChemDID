import os
import glob
import re
import time
from lammpskit.ecellmodel.filament_layer_analysis import track_filament_evolution

step = 1000
Nchunks = 7
timestep=0.0005  # in ps, for 0.5fs timestep
analysis_time = 250     # in ps

data_dir = os.path.join(os.path.dirname(__file__))
file_list = glob.glob(os.path.join(data_dir, "run.*.lammpstrj"))


r_cmo_cutoff_list = (2.3, 2.4, 2.5, 2.6)
filament_cmo_cutoff: float = 2.3

## Forming

max_timestep_forming = 7100000
min_timestep_forming = max_timestep_forming - 250/(timestep)
print("min_timestep_forming = ", min_timestep_forming)
print("max_timestep_forming = ", max_timestep_forming)

filtered_sorted_filelist_forming = sorted(
    [f for f in file_list if min_timestep_forming <= int(re.search(r"run\.(\d+)\.lammpstrj", os.path.basename(f)).group(1)) <= max_timestep_forming],
    key=lambda f: int(re.search(r"run\.(\d+)\.lammpstrj", os.path.basename(f)).group(1))
)

(print("filtered_sorted_filelist_forming = \n", filtered_sorted_filelist_forming[:5],  filtered_sorted_filelist_forming[-5:]))

time.sleep(5)


filament_using_locpot = False
filament_using_all_cmo = True

output_dir = os.path.join(data_dir,"forming", "with_all_cmo", "output_using_coordination")
print('data directory path to file_list',data_dir)
for r_cmo_cutoff in r_cmo_cutoff_list:
    analysis_name = f'{r_cmo_cutoff}_'
    # Call the function
    if file_list:
        track_filament_evolution(filtered_sorted_filelist_forming, analysis_name, timestep, step, 
                                 output_dir=output_dir, filament_using_locpot=filament_using_locpot, 
                                 r_cmo_cutoff=r_cmo_cutoff,
                                 filament_using_all_cmo=filament_using_all_cmo,
                                 filament_cmo_cutoff=filament_cmo_cutoff)
        print("Filament evolution tracking completed.")
    else:
        print("No trajectory files found for filament evolution analysis.")


filament_using_locpot = False
filament_using_all_cmo = False

output_dir = os.path.join(data_dir,"forming", "with_undercoord_cmo", "output_using_coordination")
print('data directory path to file_list',data_dir)
for r_cmo_cutoff in r_cmo_cutoff_list:
    analysis_name = f'{r_cmo_cutoff}_'
    # Call the function
    if file_list:
        track_filament_evolution(filtered_sorted_filelist_forming, analysis_name, timestep, step, 
                                 output_dir=output_dir, filament_using_locpot=filament_using_locpot, 
                                 r_cmo_cutoff=r_cmo_cutoff,
                                 filament_using_all_cmo=filament_using_all_cmo,
                                 filament_cmo_cutoff=filament_cmo_cutoff)
        print("Filament evolution tracking completed.")
    else:
        print("No trajectory files found for filament evolution analysis.")