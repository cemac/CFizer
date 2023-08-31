# CFizer
Tools to make NetCDF files CF-compliant (Climate and Forecast metadata convention), initially working with MONC (Met Office NERC cloud model) output.

CEDA provides an [overview of the CF metadata convention](https://help.ceda.ac.uk/article/4507-the-cf-metadata-convention).

[Full details of the convention are available here.](http://cfconventions.org/Data/cf-conventions/cf-conventions-1.10/cf-conventions.html)

# Installation
Clone this repository.

If running on JASMIN, load the `jaspy` environment: `module load jaspy/3.10/r20230718`. This contains all the required dependencies.

# Running
To run this development version, from this directory:
`./dev/algorithm.py [options] <source_directory>`
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
