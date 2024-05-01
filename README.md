# Star Formation Duration Estimator

## Description

SFDE is designed to compute the duration for which star formation (SF) must have lasted during a given globular cluster's (GC) formation.
To this end it takes the current physical and orbit parameters, to determine the size of the GC at birth and then uses the iron spread within the GC members to determine, how much iron must have been expelled by super novae (SNe) befor star formation stopped.
Using a given function of failed SNe and SNe ejecta depending on the stellar mass SFDE computes which stars exploded before SF ended and from these stars life expectancies, how long SF must have lasted.

## Requirements

The SFDE was tested using the following software packages. We recommend running pytest before trying to run the code with different versions of these packages:

| Package         | Version |
| :-------------- | ------: |
| python          |  3.8.10 |
| cycler          |  0.10.0 |
| kiwisolver      |   1.3.1 |
| matplotlib      |   3.3.3 |
| numpy           |  1.23.5 |
| pandas          |   1.1.5 |
| pillow          |   7.0.0 |
| pyparsing       |   2.4.6 |
| python-dateutil |   2.8.1 |
| pytz            |  2019.3 |
| six             |  1.14.0 |

For unittesting:

| Package         | Version |
| :-------------- | ------: |
| pytest          |   7.2.1 |
| attrs           |  19.3.0 |
| exceptiongroup  |   1.1.0 |
| iniconfig       |   2.0.0 |
| packaging       |    20.3 |
| pluggy          |   1.0.0 |
| tomli           |   2.0.1 |

## Usage

Currently SFDE uses 4 input files: one for the properties of the GCs, one to describe which stars explode in a SN and which don't, one file showing the amount of iron ejected by a star depending on its initial stellar mass and one file showing the remnant masses depending on the initial masses.
It is called taking the paths to these 3 files and a path to the output folder as input parameters:

```
python main.py <GC_property_file> <SN_file> <Ejecta_file> <Remnant_file> <output_folder>
```

Alternatively the SN, Ejecta and Remnant file can be combined into one file. All columns for all these files must be present in the combined file:

```
python main.py <GC_property_file> <Combined_file> <output_folder>
```

The output folder may not exist at the time of calling the script.

All parameters in the input files are organised on columns.
The order of these columns is not relevant, however, each column should have the correct column header.
The column headers and expected input units are listed in the tables below.
The columns are seperated by spaces.
'#' is used to indicate comments.
Everything after this symbol to the end of the line will not be parsed.

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
| stellar mass            | Msun          | mass[Msun]   |                              |
| life time of the star   | log10( yr )   | t_x          | x can be any value for [Z/H] |
| remant mass             | Msun          | Mfin_x       | x can be any value for [Z/H] |

Any number of columns for the life time and remnant masses can exist in this file for different metallicities.
To find the best values for a given star the program will inter- and extrapolate using the given metallicities and initial stellar masses.

## Output

All output is stored in the specified output folder.
This folder is expected to be non-existend and is created by the program itself.
The program outputs a file with the GC properties.
This contains all the GC properties that were read in and additionally the following:

| Quantity                                        | Unit | Name in File |
| :---------------------------------------------- | :--- | :----------- |
| the total mass of iron produced before SF ends  | Msun | ProducedIron |
| number of SN exploding before SF ends           |      | NSN          |
| number of SN exploding in total                 |      | NSNPos       |
| mass of the last star to explode before SF ends | Msun | mlast        |
| star formation duration                         | Gyr  | SFD          |

Additionally the program creates a barcode plot for each cluster that produces enough SNe to generate its observed iron spread.
These plots are generated in the subfolder 'Barcodes'.
The files are named after the cluster names.
