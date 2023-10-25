# CFizer
CFizer is a tool to make NetCDF output files from MONC (Met Office NERC cloud model) CF-compliant (Climate and Forecast metadata convention). It is also able to merge and split datasets according to the number of spatial dimensions, and apply compression to a user-specified level.

CEDA provides an [overview of the CF metadata convention](https://help.ceda.ac.uk/article/4507-the-cf-metadata-convention).

[Full details of the convention are available here.](http://cfconventions.org/Data/cf-conventions/cf-conventions-1.10/cf-conventions.html)

# Installation

Clone this repository.

**If using on JASMIN**, load the `jaspy` environment: `module load jaspy/3.10/r20230718`. This contains all the required dependencies.

Create an empty Anaconda environment in which to install CFizer, and activate it:
```
conda env create -n cfizer
conda activate cfizer
```
Alternatively, a pip venv may be used:
```
python -m venv /path/to/cfizer_env
source /path/to/cfizer_env/bin/activate
```

**Important**: *Before installing CFizer*, update the configuration file, `cfizer/config.yml` and, if necessary, the vocabulary (`cfizer/vocabulary.yml`) - see [[Setup]] below. Installation creates a copy of these in the virtual environment, so any subsequent changes either need to be made to the version in that environment, or require re-installation.

Install CFizer:
`pip install .`

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
To run this version, from this directory:

`cfize [options] <source_directory>`

where `<source_directory>` is the path for the directory containing the MONC NC files to be processed.

## Options
Option | Argument | Function | Example
:---|:---|:---|:---
`--target_dir`, `-t`|Target directory|Specify directory for processed NetCDF files. Default is to create a sibling directory to the source directory, appending `+processed` to the same name.|`-t /work/project/diagnostic_outputs/230228`
`--reference_time`, `-r` | Datetime in ISO format|Date or date-time of origin for time units, in ISO format (`yyyy-mm-dd[{T/ }hh:mm[:ss][{+/-}hh:mm]]`). If not specified, CFizer tries to find it in any input file it finds, and requests user input if it fails. If no input file is present and no reference time is specified or is found in existing time units, the software will exit with an error message.|`-r 2020-01-25 00:00+00:00`
`--calendar`, `-c`|Calendar|Calendar to use for time units. See [CF Conventions 4.4.1](https://cfconventions.org/Data/cf-conventions/cf-conventions-1.7/build/ch04s04.html) for valid options. Default: `proleptic_gregorian`|`--calendar proleptic_gregorian`
`--cpus`, `-p`|Integer|Total number of cores to use for parallel processing (including controller process). Default: 1 (serial)|`-p 16`
`--keep-interim`,`-i`||Set if multiple dimension sets are being merged, but you want to keep the single-dimension-set files as well. By default, these redundant files are deleted.|
`--verbose`, `-v`||Set to report all progress and function timings to stdout.|
`--quiet`, `-q`||Set to suppress any warnings normally printed to stdout upon completion.|

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

## Usage

From the root directory:
`xlvocab <path_to_spreadsheet>`

The new vocabulary YAML file will be created in the CFizer root directory.
