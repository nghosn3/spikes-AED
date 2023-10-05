# Purpose
This example brew is meant to demonstrate how to locally install a python environment in the case that you do not have access to the normal install directory

## Optional First Steps
1. Create a folder for your repository/workspace
2. Clone the workspace to the current directory using:
    - `git clone <git-url>`

## Process
1. Create a virtual environment in your desired directory
    - `virtualenv /path/to/your/env`
2. Activate the environment
    - `source /path/to/your/env/bin/activate`
3. Pip install from yaml file. This will now store the results in your envionment folder, outside of the normal install directory
    - `pip install -r requirements.yml`
