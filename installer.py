import subprocess

# List of libraries to install
libraries = ['numpy', 'pandas', 'matplotlib',
             'requests', 'selenium', 'openpyxl',
             'biopython', 'scikit-learn']

# Possible pip commands to try
pip_commands = ['pip',
                'pip3',
                'python -m pip',
                'python3 -m pip',
                'py -m pip']

# Function to check if a pip command is valid
def is_valid_pip_command(pip_command):
    try:
        subprocess.check_call(f'{pip_command} --version', stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    except subprocess.CalledProcessError:
        return False

# Find the valid pip command
valid_pip_command = None
for pip_command in pip_commands:
    if is_valid_pip_command(pip_command):
        valid_pip_command = pip_command
        break

if valid_pip_command:
    # Iterate over the libraries and install each one
    for library in libraries:
        #try:
        subprocess.check_call(f'{valid_pip_command} install {library}')
        print(f"Successfully installed {library}")
        #except subprocess.CalledProcessError:
        #    print(f"Failed to install {library}")
else:
    print("No valid pip command found. Please ensure that pip is installed and accessible.")
