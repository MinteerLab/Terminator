import os, requests, json, shutil, inspect, sys, traceback, time
import pandas as pd
import xml.etree.ElementTree as ET
from openpyxl import load_workbook
from openpyxl.formatting.rule import ColorScaleRule

def str_to_intFloat(value):
    try:
        return float(value) if '.' in value else int(value)
    except ValueError:
        return value
##############################################
def is_online():
    try:
        # Send a GET request to a reliable website (e.g., Google)
        response = requests.get("https://www.google.com")
        statusCode = response.status_code
        print(statusCode)
        # Check if the response status code indicates success (2xx or 3xx)
        if statusCode // 100 == 2 or statusCode // 100 == 3:
            return True
        else:
            print('connection failure')
            time.sleep(2)   #to avoid DDOSing google
            return False
    except requests.ConnectionError:
        time.sleep(2)
        return False
def issue(text, projDir=None):
    def get_script_function_names():
        # Get the frame of the caller's caller
        caller_frame = inspect.stack()[2]
        # Extract the module object from the frame
        module = inspect.getmodule(caller_frame[0])
        # Get the script name from the module's __file__ attribute, using basename for just the file name
        script_name = os.path.basename(module.__file__) if module and hasattr(module, '__file__') else None
        # Get the name of the caller's caller function
        function_name = caller_frame[3]

        return script_name, function_name

    script_name, function_name = get_script_function_names()
    issueLine = f"in {script_name}, {function_name}(), {text}"
    if projDir:
        printWrite(issueLine, projDir)
    else:
        print(issueLine)
def check_subdir(projDir, subdirectory):
    """Check if a subdirectory exists within a given directory, and create it if not."""
    subdirectory_path = os.path.join(projDir, subdirectory)
    if not os.path.exists(subdirectory_path):
        os.makedirs(subdirectory_path)
    return subdirectory_path
def get_traceback():
    """to print full tracebacks in case of error"""
    exc_type, exc_obj, exc_tb = sys.exc_info()
    full_traceback = "".join(traceback.format_exception(exc_type, exc_obj, exc_tb))
    return full_traceback
def printWrite(text, projDir, firstWrite = False):
    if firstWrite:
        uniprotID = projDir.split('\\')[-1].split(' ')[0]
        timeNow   = time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())
        txtName   = f'{uniprotID}_{timeNow}_outputFile.txt'
        fileName  = name_from_projDir(txtName, projDir, exactName=True)
    else:
        fileName  = name_from_projDir('outputFile.txt', projDir)
    with open(fileName, 'a') as outTxt:
        outTxt.write(f'\n{text}')
    print(text)

##############################################

def make_tempDir():
    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    tempDir = os.path.join(parent_dir, "temporary")
    # Check if the "temporary" directory already exists, and delete it if so
    if os.path.exists(tempDir):
        shutil.rmtree(tempDir)
    os.makedirs(tempDir)
    return tempDir
def fetchUrl(url, attribute=''):
    response    = requests.get(url)
    attributes  = {
        '' : response,
        'text': response.text,
        'status_code': response.status_code}
    if attribute not in attributes:
        issue(f"{attribute} is not a valid attribute, \n available attributes are {[key for key in attributes.keys()]}")
    return attributes.get(attribute, False)
def projDir_if_projDir(uniprotID):
    """Check if a project directory exists for the specified UniProt ID. returns it."""
    parentDir    = os.pardir
    subdirs      = os.listdir(parentDir)
    dirs         = [subdir for subdir in subdirs if os.path.isdir(os.path.join(parentDir, subdir)) and subdir.startswith(uniprotID)]
    if dirs:
        dirName     = dirs[0]
        rel_path   = os.path.join(parentDir, dirName)
        full_path   = os.path.abspath(rel_path)
        return full_path
    return False
def if_not_in_projDir(aFunction, outputFileName, targetDir, exactName = False, save = True, passProjDir = False, **kwargs):
    """if outputFile isn't in projDir, runs aFunction to generate outputFile, passing kwargs to the function, then saves it to projDir\n passProjDir True when aFunction needs to receive projDir as a kwarg as well"""
    def projDirWrapper(aFunction, projDir, **kwargs):
        return aFunction(**kwargs, projDir = projDir)
    print(f'\nchecking for {outputFileName} in project directory')
    data = data_from_projDir(outputFileName, targetDir, exactName=exactName)
    if not data:
        print(f"running {aFunction.__name__} to generate {outputFileName}")
        if passProjDir:
            functionOutput = projDirWrapper(aFunction, targetDir, **kwargs)
        else:
            functionOutput = aFunction(**kwargs)
        if save:
            save_to_dir(functionOutput, outputFileName, targetDir, exactName = exactName)
        return functionOutput
    else:
        return data
def data_from_projDir(fileName, projDir, exactName = False):
    """Get a file from the project directory as an object."""
    fileName    = name_from_projDir(fileName, projDir, exactName)
    shortName   = fileName.split('\\')[-1]

    if is_in_Dir(fileName, projDir, exactName = False):
        print(f"extracting data from {shortName} in project directory")
        fileType = '.'+fileName.split('.')[-1]
        if fileType == '.xml':            # Try to get_consfDicts as XML
            try:
                tree = ET.parse(fileName)
                return tree.getroot()           # Return the root element of the XML tree
            except ET.ParseError:
                issue(f"Error: {shortName} is not a valid XML.", projDir)
                return False
        elif fileType == '.json':
            with open(fileName, 'r') as file:
                return json.load(file)
        elif fileType in ['.txt', '.pdb']:
            with open(fileName, 'r') as file:
                return file.read()
        elif fileType == '.xlsx':
            data = list(pd.read_excel(fileName))
            return data
        else:
            issue(f"data_from_projDir, unsupported filetype {fileType} in {shortName}", projDir)
            return False
    else:
        print(f"{shortName} not found in project directory")
        return False
def is_in_Dir(fileName, projDir, exactName = False):
    """Check if a file exists in the specified directory."""
    file_path = name_from_projDir(fileName, projDir, exactName=exactName)
    return os.path.exists(file_path)
def name_from_projDir(fileName, projDir, exactName = False):
    """Find the full path to a file in the project directory that ends with fileName, returns last match.
    If exactName is True, returns the full path to the file with the exact name.
     Else returns what the path would be if the file existed."""

    if exactName:
        return os.path.join(projDir, fileName)

    fileNames = []
    for file in os.listdir(projDir):
        if file.endswith(fileName):
            fileNames.append(file)
    if fileNames:
        return os.path.join(projDir, fileNames[-1])
    else:
        return os.path.join(projDir, fileName)
def save_to_dir(data, fileName, projDir, force = False, formatExcel = False, exactName = False):
    """download a file to the project directory, or overwrite if force is True."""

    def format_xlsx(fileName):
        """freezes top row and conditional formats an excel file"""
        wb = load_workbook(fileName)
        ws = wb.active  # select sheet
        ws.freeze_panes = 'B2'  # Freeze the first row and column
        color_scale_rule = ColorScaleRule(start_type='min', start_color='367d5f',
                                          mid_type='percentile', mid_value=50, mid_color='ffffe0',
                                          end_type='max', end_color='8c0035')
        ws.conditional_formatting.add('B2:AA1000', color_scale_rule)

        # Save the workbook with changes
        wb.save(fileName)

    if is_in_Dir(fileName, projDir, exactName = exactName) and not force:
        print(f"tried to save, but found an existing {fileName} in project directory")
    else:
        if data:
            print(f"saving {fileName} to project directory")
            fileName = name_from_projDir(fileName, projDir)
            if fileName.endswith('.xml'):
                # Save as XML
                if isinstance(data, ET.Element):
                    tree = ET.ElementTree(data)
                    tree.write(fileName)
                    return True
                else:
                    shortName = fileName.split('\\')[-1]
                    issue(f"Error: data sent to {shortName} is not an XML element, but a {type(data)}", projDir)
            elif fileName.endswith('.txt'):
                # Save as text file
                try:
                    with open(fileName, 'w') as outFile:
                        outFile.write(str(data))
                    return True
                except:
                    issue(f"Error: data sent to {fileName} is not a string.", projDir)
                    return False
            elif fileName.endswith('.xlsx'):
                # Save as Excel file
                data = pd.DataFrame(data)
                data.to_excel(fileName, index=False)
                if formatExcel:
                    format_xlsx(fileName)
                return True
            elif fileName.endswith('.json'):
                # Save as JSON file
                with open(fileName, 'w') as outFile:
                    json.dump(data, outFile, indent=4)
                return True
            elif fileName.endswith('.pdb'):
                # Save as PDB file
                with open(fileName, 'wb') as outFile:
                    outFile.write(data)
                return True
            else:
                issue(f"Unsupported file type in {fileName}", projDir)
                return False
        else:
            issue(f'tried to save {fileName} but object received was not savable data.', projDir)
def delete_from_projDir(fileName, projDir):
    """Deletes a file from the project directory if it exists.

    Args:
        fileName (str): The name of the file to be deleted.
        projDir (str): The path to the project directory.

    Returns:
        bool: True if the file was deleted successfully, False otherwise."""
    file_path = name_from_projDir(fileName, projDir)
    if os.path.exists(file_path):
        try:
            os.remove(file_path)
            print(f"{fileName} has been deleted from the project directory.")
            return True
        except Exception as e:
            issue(f"Error deleting {fileName}: {e}", projDir)
            return False
    else:
        print(f"{fileName} does not exist in the project directory.")
        return False
##############################################
