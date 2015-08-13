`.` is the project directory. All paths must be defined locally (except toolbox path).


# Preprocessing

A first step reads raw data in the following files:

```
./data/ (cond1/cond2/) plate_id1/plate_id1_well_id1.fcs
                                 plate_id1_well_id2.fcs
                                 ...
                                 plate_id1_cond.od.txt  (e.g. plate_id1_2hgrowth.od.txt, or plate_id1_pbs_10.od.txt)
```

and writes control plots to the directory of raw fcs files (here `data`).
Cache files (and gated data in .tab or .fcs format if any) are written to the directory of your choice (e.g. `preproc`):

```
./data/ (cond1/cond2/) plate_id1/plate_id1.pdf
./preproc/ (cond1/cond2/) plate_id1/plate_id1.Rdata
                                    plate_id1_well_id1.txt
                                    plate_id1_well_id2.txt
                                    ...
```

NB: The directory structured of the raw data directory will be replicated in the preproc directory.

NB: data2preproc is a function that must be defined; it converts raw data path to preproc path. If you want cache files to be stored along with raw data, simply use `identity`

preproc_facs_plate() returns a list containing the following 3 dataframes:
- gates (coordinates of the fsc/ssc gates)
- preproc (all gated events: one line per event)
- stats (stats per well: one line per well, including od values. Extra columns from the index file can be subsequently appended).

## FACS related parameters

### overriding FACS parameters for one plate

## Index files

Index files can be defined in 2 ways, per plate or per well (or set of wells).

### Plates index

Table with one row per plate. It must have a column `dir` with the path to the plate raw data, and a column `discard` with a comma separated list of wells to be discarded (typically edited after producing control plots).
IMPORTANT: All paths must be relative to the project directory and have no trailing slash.

### Wells index
Table with one row per well or set of wells. It must have a column `dir` with the path to the plate raw data, and a column `well` with a comma separated list of wells corresponding to one condition/strain. Only wells described in this index file will be kept for further analysis.
IMPORTANT: All paths must be relative to the project directory and have no trailing slash.

## Data reloading 

When analysis is resumed on already preprocessed data, `gates`, `preproc`, and `stats` are taken from cache files if they exists. Info from the index file still need to be appended.

# Data structures

Control plots can be checked manually to specify which well must be excluded in the index file.

A second step subsets preprocessed data following the index file in order to produce the following data structures:
- f_preproc (all gated events: one line per event)
- f_od (od values of all wells: one line per well)
- f_stats (stats per well: one line per well with extra columns from the index file)


