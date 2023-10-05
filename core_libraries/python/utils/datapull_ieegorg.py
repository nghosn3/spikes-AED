import os
import argparse

def file_methods(filetype):
    
    if filetype in ['pickle','pkl']:
        method = 

def access_type(filepath):
    """
    Determine the type of the data we are accessing by rules and user-input.

    Parameters
    ----------
    filename : STR
        Filepath.

    Returns
    -------
    Bound method for reading in the data.

    """
    
    # Get the filename and filetype
    filename = os.path.basename(filepath)
    filetype = filename.split('.')[-1]
    
    # Get the most likely read in statements
    

def main():
    
    # Command line options needed to obtain data.
    parser   = argparse.ArgumentParser()
    parser.add_argument('-f', '--local_path', type=str, help='Read in user provided file directly.')
    parser.add_argument('-u', '--user', help='username')
    parser.add_argument('-p', '--password', help='password (will be prompted if omitted)')
    parser.add_argument('--dataset', help='dataset name')
    parser.add_argument('--start', type=int, help='start offset in usec')
    parser.add_argument('--duration', type=int, help='number of usec to request')
    parser.add_argument('--datapointer_path', default=None, type=str, help='Path to an alternate data pointer file.')
    parser.add_argument('--silent', dest='verbose', default=True, action='store_false', help='Silent Verbose Output. Default=False.')
    parser.add_argument('--nchan', type=int, help='Number of channels')
    parser.add_argument('--szfile', help='Filename containing seizure times created using Dr. Conrads xlsx.')
    args = parser.parse_args()
    return args