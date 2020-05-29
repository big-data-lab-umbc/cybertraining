# Standard library imports.
import datetime

# Third party imports.
import numpy as np

# Constant definitions.
REFERENCE_DATE = datetime.datetime(1993, 1, 1)


def get_timestamp(seconds : int) -> datetime.datetime:
    
    '''
    
    Accepts as input the number of seconds since the reference date and returns
    the corresponding datetime.
    
    '''
    
    seconds      = np.array(seconds).astype(object)
    time_elapsed = np.vectorize(datetime.timedelta)(seconds = seconds)
    timestamp    = REFERENCE_DATE + time_elapsed
    
    return timestamp
