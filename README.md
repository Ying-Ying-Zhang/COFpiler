# COFpiler

Build the statistical structure for layered materials by stacking layers following the Maxwell-Boltzmann energy distribution of their stacking modes.

## usage:

```
python COFpiler.py [-h] [--path PATH] [--data DATA] [--i I] [--o O] [--T T] [--s S] [--mirror MIRROR] [--mplane MPLANE] [--L L] [--M M]
```

## Optional arguments: 

```
-h, --help           Show this help message and exit  
-data DATA           The input file include data needed, form as: form as: the 1st column is the stacking_type, the 2nd is the 
                   Erel (relative energies), the 3rd, 4th and 5th columns are the x, y, z (vector of the shift_vectors)  
-i --instr           The input structure  
-T --tem             The experimental synthesize temperature  
-path PATH           The path where the infile and instr are, the output structure will also save there  
-o --outstr_format   The output structure format  
-s --symmetry        The symmetry for the structure, C3, C6 or C4  
-m --mirror          Enable the consider of mirror shift or not  
-mp --mplane         The mirror plane for s2 shift  
-L L                 The layer number  
-M M                 The model number
```
  
  ## License
  MIT © Richard McRichface
