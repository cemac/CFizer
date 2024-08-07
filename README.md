# CFizer
CFizer is a tool to make NetCDF output files from MONC (Met Office NERC cloud model) CF-1.10 compliant (Climate and Forecast metadata convention). It is also able to merge and split datasets according to the number of spatial dimensions, and apply compression to a user-specified level. Please note that although CFizer generates CF-1.10 compliant MONC output, the metadata is labelled as CF-1.8 so that - if required - it can be run through CEDA's cf-checker which, as of July 2024, checks up to CF-1.8 only.

CEDA provides an [overview of the CF metadata convention](https://help.ceda.ac.uk/article/4507-the-cf-metadata-convention).

[Full details of the convention are available here.](http://cfconventions.org/Data/cf-conventions/cf-conventions-1.10/cf-conventions.html)

# Installation

Clone this repository and cd into it.

Create an empty Anaconda environment in which to install CFizer, and activate it:
```
conda create -n cfizer
conda activate cfizer
```
Alternatively, a pip venv may be used:
```
python -m venv /path/to/cfizer_env
source /path/to/cfizer_env/bin/activate
```
**If using on JASMIN**, load the `jaspy` environment: `module load jaspy/3.10/r20230718`. This contains all the required dependencies.

**If not using on JASMIN**, load the following packages in your cfizer conda environment (this has been tested on both ARC4 and ARCHER2):

- netcdf4=1.5.7
- xarray = 2022.3.0
- pyyaml = 6.0
- cfunits = 3.3.4
- numpy = 1.26.3
- dask = 2023.11.0
- openpyxl = 3.0.10

**Important**: *Before installing CFizer*, update the configuration file, `cfizer/config.yml` and, if necessary, the vocabulary (`cfizer/vocabulary.yml`) - see [Setup](#setup) below. Installation creates a copy of these in the virtual environment, so any subsequent changes either need to be made to the version in that environment, or require re-installation.

Install CFizer:
From the top level of the repository type `pip install .`

# Setup
The most important setup for users is to check/complete the `cfizer/vocabulary.yml` and `cfizer/config.yml` files. They define, respectively, how the MONC variables are to be modified for CF compliance, and parameters that should be uniform across a set of files, including the source of the original data.

At the beginning of `config.yml` are the global attributes that [CF Conventions recommend](http://cfconventions.org/Data/cf-conventions/cf-conventions-1.10/cf-conventions.html#description-of-file-contents) be present in all datasets:
- `title`:          Set automatically by CFizer, based on file names and any merge/split operations.
- `institution`:    This should be provided in `config.yml`. "Specifies where the original data was produced."
- `source`:         Details of how the data were generated, e.g. MONC version ... To be provided in `config.yml`.
- `history`:        Any modifications to the data after generation. Set automatically by CFizer.
- `references`:     To be provided in `config.yml`. "Published or web-based references that describe the data or methods used to produce it."
- `comment`:        To be provided in `config.yml`. "Miscellaneous information about the data or methods used to produce it."
- `conventions`:    Set automatically by CFizer.

# Running
To run the CFizer on a directory of MONC outputs:

`cfize [options] <source_directory>`

where `<source_directory>` is the path to the directory containing the NC files to be processed.

## Options
Option | Required (Y/N)? | Argument | Function | Example
:---|:---|:---|:---|:---
`--target_dir`, `-t`|N|Target directory|Specify directory for processed NetCDF files. Default is to create a sibling directory to the source directory, appending `+processed` to the same name.|`-t /work/project/diagnostic_outputs/230228`
`--reference_time`, `-r` |Y| Datetime in ISO format|Date or date-time of origin for time units, in ISO format and as a string (`"yyyy-mm-dd[{T/ }hh:mm[:ss][{+/-}hh:mm]]"`). If not specified, CFizer tries to find it in any input file it finds, and requests user input if it fails. If no input file is present and no reference time is specified or is found in existing time units, the software will exit with an error message.|`-r "2020-01-25 00:00+00:00"`
`--calendar`, `-c`|N|Calendar|Calendar to use for time units. See [CF Conventions 4.4.1](https://cfconventions.org/Data/cf-conventions/cf-conventions-1.7/build/ch04s04.html) for valid options. Default: `proleptic_gregorian`|`--calendar proleptic_gregorian`
`--cpus`, `-p`|N|Integer|Total number of cores to use for parallel processing (including controller process). Default: 1 (serial)|`-p 16`
`--keep-interim`,`-i`|N| |Set if multiple dimension sets are being merged, but you want to keep the single-dimension-set files as well. By default, these redundant files are deleted.|
`--verbose`, `-v`|N| |Set to report all progress and function timings to stdout.|
`--quiet`, `-q`|N| |Set to suppress any warnings normally printed to stdout upon completion.|

## Example Slurm Submission Script (JASMIN)
```bash
#!/bin/bash
#SBATCH --partition=test
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -N 1
#SBATCH --job-name=cfizer_example
#SBATCH -o <run_directory>/%j.out
#SBATCH -e <run_directory>/%j.err
#SBATCH --time=60
#SBATCH --mem-per-cpu=16384

# module add jaspy/3.10/r20230718  # If not already running
# conda activate cfizer  # If not already active

srun cfize -p 4 -v -i -t <output_directory> <directory_to_process>
```

## Creating Vocabulary File from Excel Spreadsheet

The `vocab_from_xlsx` tool creates the required `vocabulary.yml` file from a Microsoft Excel spreadsheet, providing it is in the expected format:

- Columns:
    - `updated_name`
    - `units` (CF-compliant)
    - `axis` (CF-compliant)
    - `standard_name` (CF-compliant)
    - `long_name`
    - `dimension_changes` (in form `current:new`)
    - `perturbation_to_absolute` (value either `True` or `False`)
    - `reference_variable` (for `perturbation_to_absolute` only)

Any field can be left blank if not needed (e.g. if field is already CF-compliant or, in the case of variable names, if no change is needed).

Any additional columns will be ignored.

Because the vocabulary is organised by the number of spatial dimensions for a given variable (and output file), the spreadsheet should likewise be divided by dimensions. Before each set of variables should be a row containing only a dimensionality indicator in the first cell, e.g. `0D` for files/variables with only a time dimension and `3D` for file/variables with 3 spatial dimensions as well as time.

### Usage

`xlvocab <path_to_spreadsheet>`

The new vocabulary YAML file will be created in the CFizer source code directory.

# Version History

1.0.0  Initial release.

1.0.1  Allows backwards compatibility with Python 3.9, by removing `root_dir` option from `iglob` function calls.

1.0.2  Code reformatting using *Black*, for consistency, and added maintainer to pyproject.toml.

1.0.3  Includes the following fixes: 1) Adds checking for empty DsGroup filepath collections; 2) Adds check for availability of reference variables for converting perturbation values; 3) Corrects Boolean values in vocabulary in Excel spreadsheet to vocabulary tool; 4) Makes some minor corrections and adds some debugging clarifications. README.md updated to include instructions for installation on HPC systems other than JASMIN. Added maintainer to pyproject.toml.

1.0.4 Although CFizer generates CF-1.10 compliant MONC output, `CF_VERSION` in startup.py has been relabelled as 1.8 instead of 1.10 (this label becomes part of the CFizer'd output metadata) so that output can be run through CEDA's cf-checker, which currently checks up to CF-1.8 - see issue #58 for further details.

1.0.5 Included fix to cfize_ds.py in order to preserve time dimension when splitting datasets.
