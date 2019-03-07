Before running a simulation with new maps, it is advised to preprocess them with the command line tool 'loadmap'.
An example is available in this directory at processmaps.sh

All the available options are listed below:

| shortcut | identifier | description | type |
| -------- |:----------:|:---------------------------------------------------------------------------------------------------------------------:|-----:|
| e        | efield     | Electric field map file | string|
|          | efield_sc  | Electric field map skip columns | int|
|          | efield_sr  | Electric field map skip rows | int|
|          | efield_scale  | Electric field scaling factor | int|
|          | efield_scalespace  | Electric field coordinates scaling factor | double|
|          | efield_output      | Electric field output file | string |
|  E       | emob               | Electron mobility map file | string |
|          | emob_sc | Electron mobility map skip columns | int |
|          | emob_sr | Electron mobility map skip rows | int |
|          | emob_scale | Electron mobility map scaling factor | double |
|          | emob_scalespace | Electron mobility map coordinates scaling factor | double |
|          | emob_out | Electron mobility map output file | string |
|  H       | hmob               |Hole mobility map file | string |
|          | hmob_sc | Hole mobility map skip columns | int |
|          | hmob_sr | Hole mobility map skip rows | int |
|          | hmob_scale | Hole mobility map scaling factor | double |
|          | hmob_scalespace | Hole mobility map coordinates scaling factor | double |
|          | hmob_out | Hole mobility map output file | string |
| w        | wfield     | Weighting field map file | string|
|          | wfield_sc  | Weighting field map skip columns | int|
|          | wfield_sr  | Weighting field map skip rows | int|
|          | wfield_scale | Weighting field scaling factor | int|
|          | wfield_scalespace  | Weighting field coordinates scaling factor | double|
|          | wfield_output      | Weighting field output file | string |
| o        | OutputDirectory    | Output directory | string |


    
If a single map with all quantitites is given (the order must be x,y,z, Ex, Ey, Ez, emob, hmob, Wx, Wy, Wz) the options are:

| shortcut | identifier | description | type |
| -------- |:----------:|:-----------:|-----:|
| p | physmap | Full map of fields | string |
|   | physmap_sc | Full map skip columns | int |
|   | physmap_sr | Full map skip rows | int |
|   | physmap_scale | Full map scaling factor | double |
|   | physmap_scalespace | Full map space scaling factor | double |
|   | physmap_out | Full map output file | string |
| o | OutputDirectory    | Output directory | string |
