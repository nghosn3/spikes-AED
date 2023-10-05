# Utility Scripts

Below we provide a brief description of each utility script.

## stream_ssh.py
This enables streaming data between a remote host and a python instance. You can either read data from a remote source into memory using the data access method of your choice (i.e. Pandas, numpy, pyedflib, etc.) with the requisite options or save the data to the 
remote system.

This should be considered when code is close to release and data needs to be stored on systems such as Lief or Borel.

### Examples

```
import pandas as PD
from functools import partial
import stream_ssh as STRMSSH
DF = STRMSSH.read('path/to/file/on/remote/system.csv','borel.seas.upenn.edu','USERNAME','USERPASSWORD',partial(PD.read_csv,names=['a','b'],usecols=[1,2])
```

will stream the data on borel into the partial function instantiation of a pandas dataframe and return your data with the associated key:value pairs applied.
