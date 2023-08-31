#!/usr/bin/env python

import openpyxl as xl
import yaml
import re
import os
import os.path as op
import argparse
from cfunits import Units
import sys
# os.chdir(op.join(os.getcwd(), 'dev'))
# print(os.getcwd())
# Workaround to find package directory (https://stackoverflow.com/a/61571300)
if os.getcwd() not in sys.path:
    sys.path.append(os.getcwd())  # This shouldn't be necessary once packaged
# print(sys.path)
from dev import VOCAB_FIELDS, app_dir


def vocab_from_xls(filepath: str) -> dict:
    '''
    This will for Excel file to convert to vocabulary dictionary, if it 
    contains the right columns (fields).
    It will attempt to fill any missing data based on CF conventions, including:
    - inferring cell_method from presence of terms such as "mean" and "max" in variable name.
    - looking up standard_name, if supplied.
    It will check any units given using `cfunits`, and check against CF standard names definitions.

    The intention is to compile a full lookup table, such that the main 
    application does not need to infer anything or look anywhere else for the 
    required data.
    '''
    # TODO: Load or download standard_names table.

    vocabulary = {}
    try:
        wb = xl.load_workbook(filepath)
        ws = wb['Vocabulary']
    except:
        exit("Workbook can't be loaded, or does not contain 'Vocabulary' worksheet.")

    header_row = {c.value: i for i, c in enumerate(ws['1']) 
                  if c.value in VOCAB_FIELDS.keys()}
    
    subset = 'all'  # As backup, in case headings not found, or not found at beginning

    for i, row in enumerate(ws.rows):
        if i==0:  # header row
            continue
        elif not row[0].value:  # Blank row
            continue
        elif row[0].font.b or (
            len(row[0].value)==2 and 
            row[0].value[-1].upper() == 'D' and 
            row[0].value[0].isdigit()
        ):  # heading indicated by bold type or a 2-character dimension code
            if (len(row[0].value)==2 and 
                row[0].value[-1].upper() == 'D' and 
                row[0].value[0].isdigit()):  # Dimensionality
                subset=int(row[0].value[0])
            else:  # Non-dimension heading
                subset = row[0].value
            # print(subset)
            vocabulary.update({subset: {}})
        else:  # variable
            if subset == 'all':
                vocabulary.update()
            else:
                vocabulary[subset].update({
                    row[0].value.split('(')[0].strip(): 
                    {field: row[index].value 
                    for field, index in header_row.items() 
                    if row[index].value is not None}
                })
                # print(row[0].value.split('(')[0].strip(), ':', 
                #     vocabulary[subset][row[0].value.split('(')[0].strip()])

    empty = []
    for subset, var_set in vocabulary.items():
        if not var_set:  # empty subset
            empty.append(subset)
            continue
        for var, attrs in var_set.items():
            if 'dimension_changes' in attrs:
                if not isinstance(vocabulary[subset][var]['dimension_changes'], dict):
                    
                    vocabulary[subset][var]['dimension_changes'] = {c.split(':')[0].strip(): c.split(':')[1].strip() 
                            for c in re.split(
                                '[,;]', 
                                vocabulary[subset][var]['dimension_changes']
                                )
                    }
                    # print(vocabulary[subset][var]['dimension_changes'])
            if 'long_name' in attrs:
                vocabulary[subset][var]['long_name'] = vocabulary[subset][var]['long_name'].title()
    [vocabulary.pop(e) for e in empty]
    return vocabulary


def vocab_to_yaml(vocab: dict) -> str:
    filepath = op.join(app_dir, 'vocabulary.yml')
    if op.exists(filepath):
        while True:
            overwrite = input(
                f"Vocabulary file {filepath} already exists. Overwrite?"
            )
            if overwrite[0].lower() == 'y': break
            elif overwrite[0].lower() == 'n':
                print("Exiting without updating vocabulary.yml")
                sys.exit(0)
    with open(filepath, "w") as f:
        yaml.safe_dump(vocab, f)
    return filepath


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'filepath',
        help='Path to Microsoft Excel spreadsheet (.xlsx) file containing the vocabulary in a sheet of that name.'
    )
    args = parser.parse_args()
    excel_file = op.abspath(args.filepath)
    if not op.exists(excel_file):
        raise OSError(f"Filepath {excel_file} not found.")
    if not 'xls' in excel_file[-4:]:  # Finds 'xls' or 'xlsx':
        raise TypeError("Filepath must be to an Excel spreadsheet file.")
    vocabulary = vocab_from_xls(filepath=excel_file)
    yaml_file = vocab_to_yaml(vocab=vocabulary)
    print(f"New vocabulary file created: {yaml_file}")
    sys.exit(0)
    