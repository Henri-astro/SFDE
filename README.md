# Star Formation Duration Estimator

## Description

SFDE is designed to compute the duration for which star formation (SF) must have lasted during a given globular cluster's (GC) formation.
To this end it takes the current physical and orbit parameters, to determine the size of the GC at birth and then uses the iron spread within the GC members to determine, how much iron must have been expelled by super novae (SNe) befor star formation stopped.
Using a given function of failed SNe and SNe ejecta depending on the stellar mass SFDE computes which stars exploded before SF ended and from these stars life expectancies, how long SF must have lasted.

## Usage

Currently SFDE uses 4 input files: one for the properties of the GCs, one to describe which stars explode in a SN and which don't, one file showing the amount of iron ejected by a star depending on its initial stellar mass and one file showing the remnant masses depending on the initial masses.
It is called taking the paths to these 3 files and a path to the output folder as input parameters:

```
python main.py <GC_property_file> <SN_file> <Ejecta_file> <Remnant_file> <output_folder>
```

The output folder may not exist at the time of calling the script.

All parameters in the input files are organised on columns.
The order of these columns is not relevant, however, each column should have the correct column header.
The column headers and expected input units are listed in the tables below.
The columns are seperated by spaces.

### The GC property file

The GC property file contains the following GC parameters:

| Quantity                     | Unit | Name in File |
| :--------------------------- | :--- | :----------- |
| Name                         | -    | Name         |
| present-day cluster Mass     | Msun | Mass         |
| Apocentre                    | kpc  | R_a          |
| Pericentre                   | kpc  | R_p          |
| Star formation efficiency    | -    | SFE          |
| [Fe/H]                       | dex  | Fe-H         |
| iron spread                  | dex  | FeSpread     |
| cluster age                  | Gyr  | Age          |

### The SN file

The SN file contains the following parameters

| Quantity     | Unit   | Name in File |
| :----------- | :----- | :----------- |
| stellar mass | Msun   | mass[Msun]   |
| SN (yes/no)  | 1 or 0 | SN           |

The program interpolates between the stellar masses.
That means that if a stellar mass is closest to a mass given in the file for which an explosion occurs, a star of that mass is seen as exploding and vice versa.

### The ejecta file

A file containing the amount of iron ejected by a SN depending on the initial stellar mass.

| Quantity                | Unit | Name in File |
| :---------------------- | :--- | :----------- |
| stellar mass            | Msun | mass[Msun]   |
| amount of iron ejected  | Msun | Fe[Msun]     |

The program interpolates between the stellar masses.

### The remnant file

A file containing the life times and remnant masses of stars for different initial masses and metallicities.

| Quantity                | Unit          | Name in File | Comment                      |
| :---------------------- | :------------ | :----------- | :--------------------------- |
| stellar mass            | log10( Msun ) | Mstar        |                              |
| life time of the star   | log10( yr )   | t_x          | x can be any value for [Z/H] |
| remant mass             | log10( Msun ) | Mfin_x       | x can be any value for [Z/H] |

Any number of columns for the life time and remnant masses can exist in this file for different metallicities.
To find the best values for a given star the program will inter- and extrapolate using the given metallicities and initial stellar masses.
