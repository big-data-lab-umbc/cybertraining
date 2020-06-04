# Standard library imports.
import datetime
import re
from os import walk
from os.path import join

# Third party imports.
import pandas as pd

# Constant definitions.
TIMESTAMP_FORMAT = "%Y-%m-%dT%H-%M-%SZ"
DATA_FILE_REGEX  = (r"CAL_LID_L2_(?P<resolution>\d\d)kmCLay-"
                    r"(?P<production_strategy>.*)-(?P<version>V\d-\d\d)."
                    r"(?P<timestamp>\d\d\d\d-\d\d-\d\dT\d\d-\d\d-\d\dZ)"
                    r"(N|D).*")


def _get_files(directories, directory_key):
    
    # TODO: Account for no files found.
    
    path = directories[directory_key]
    
    files = pd.DataFrame()
    
    filenames, paths = \
        zip(*[(file, join(root, file))
              for root, _, files in walk(path)
              for file in files])
        
    files["filename"] = filenames
    files["path"]     = paths
    
    regex_matches = \
        files.filename.map(lambda filename: re.match(DATA_FILE_REGEX,
                                                     filename))
        
    files["timestamp"] = regex_matches.map(_regex_match_to_datetime)
    
    # Filter out those filenames which do not match the regular expression.
    mask = regex_matches.map(lambda regex_match: bool(regex_match))
    
    files = files[mask]
    
    files.drop(columns = "filename", inplace = True)
    
    files.set_index("timestamp", inplace = True)
    
    column_names = \
        {column : f"{directory_key}" for column in files.columns}
        
    files.rename(columns = column_names, inplace = True)
    
    return files


def _regex_match_to_datetime(m):
    
    if m is not None:
    
        timestamp_string = m.groupdict()["timestamp"]
        timestamp        = datetime.datetime.strptime(timestamp_string,
                                                      TIMESTAMP_FORMAT)
        
        return timestamp


def get_filenames(directories, start_date = None, end_date = None):
    
    filenames = []
    
    for key in directories.keys():
        
        filename_array = _get_files(directories, key)
        
        # Drops duplicate indices.
        filename_array = filename_array.loc[~filename_array.index.duplicated()]
        
        filenames.append(filename_array)
    
    filenames = pd.concat(filenames, axis = 1)
    
    filenames.dropna(inplace = True)
            
    if start_date is not None:
        
        filenames = filenames[start_date <= filenames.index]
        
    if end_date is not None:
        
        filenames = filenames[filenames.index < end_date + datetime.timedelta(days = 1)]
        
    return filenames
