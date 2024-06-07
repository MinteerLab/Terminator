import fetch_uniprot

from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import TimeoutException, NoSuchElementException
from QOL_functions import *

#DB request
def req_consurfDB(pdbID, targetSeq, projDir, debug = False, consurfPDB = False, consurfChain = False):
    """inputs requestInput into consurfDB, downloads consurf summary to projDir"""
    chromeOptions  = Options()
    consurfURL      = 'https://consurfdb.tau.ac.il/index.php'
    chromeOptions  .add_experimental_option("prefs", {
        "download.default_directory": projDir+'\\',
        "download.prompt_for_download": False,  # To disable the download prompt
        "download.directory_upgrade": True,     # To ensure file downloads work in headless mode
        "safebrowsing.enabled": True})
    if not debug:
        chromeOptions.add_argument("--headless=new")  # don't visibly open window while making consurf

    def consurf_input_by_sequence(targetSeq, chromeOptions, consurfURL, projDir):
        """easier to input by sequence"""
        try:
            driver = webdriver.Chrome(options=chromeOptions)
            driver.get(consurfURL)  # Navigate to consurfDB
        except:
            issue("Failed to access ConsurfDB to download conservation information.", projDir)
            return 'failed', None
        wait = WebDriverWait(driver, 10)
        #input target sequence
        textarea = wait.until(EC.visibility_of_element_located((By.ID, "protein_seq")))
        textarea .clear()
        textarea.send_keys(targetSeq)
        #make request
        applyButton = wait.until(EC.element_to_be_clickable((By.XPATH, "//input[@value='Apply']")))
        applyButton .click()
        #await results page
        try:
            success = "https://consurfdb.tau.ac.il/main_output.php?pdb_ID="
            failure = "https://consurfdb.tau.ac.il/scripts/error.php?reason=Sequence%20not%20found%20in%20ConSurf%20database"

            urlWait = wait.until(
                lambda d: d.current_url.startswith(success) or d.current_url.startswith(failure))
            url = driver.current_url
            if url.startswith(failure):
                driver.quit()
                return False, None
            else:
                return urlWait, driver
        except TimeoutException:
            driver.quit()
            return False, None

    def consurf_input_by_pdbID(pdbID, chromeOptions, consurfURL, consurfPDB, consurfChain, projDir):
        """When a truncated structure is in consurfDB, a sequence search will fail, must search by PDBID"""
        driver      = webdriver.Chrome(options=chromeOptions)
        driver      .get(consurfURL)  # Navigate to consurfDB
        wait        = WebDriverWait(driver, 10)  # an object containing page info

        #input pdbID
        inputField = wait.until(EC.visibility_of_element_located((By.ID, "pdb_ID_field")))
        inputField.clear()
        if not consurfPDB:
            inputField.send_keys(pdbID)
        else:
            inputField.send_keys(consurfPDB)
        # set chain
        try:
            dropdown = WebDriverWait(driver, 4).until(EC.element_to_be_clickable((By.CSS_SELECTOR, ".dropdown-select.wide.select_box")))
        except:
            issue(f'no chains offered for {pdbID}', projDir)
            return False, None
        dropdown.click()
        wait.until(EC.visibility_of_element_located(
            (By.CSS_SELECTOR, ".list ul")))  # Wait for the dropdown options to become visible
        if not consurfChain:
            option = wait.until(EC.element_to_be_clickable((By.XPATH,
                            "//li[contains(@class, 'option') and @data-value != ''][1]")))
        else:
            option = wait.until(EC.element_to_be_clickable(
                ((By.XPATH, f"//li[starts-with(@data-value, '{consurfChain} ')]"))))
        option.click()
        #make request
        applyButton = wait.until(EC.element_to_be_clickable((By.XPATH, "//input[@value='Apply']")))
        applyButton .click()
        #await results page
        try:
            success = "https://consurfdb.tau.ac.il/main_output.php?pdb_ID="
            failure = "https://consurfdb.tau.ac.il/scripts/error.php?reason=Sequence%20not%20found%20in%20ConSurf%20database"

            urlWait = wait.until(
                lambda d: d.current_url.startswith(success) or d.current_url.startswith(failure))
            url = driver.current_url
            if url.startswith(failure):
                return False, None
            else:
                return urlWait, driver
        except TimeoutException:
            driver.quit()
            return False, None

    try:
        driverResult, outDriver = consurf_input_by_sequence(targetSeq, chromeOptions, consurfURL, projDir)
        if driverResult == 'failed': #failed to even access the database
            issue("Failed to access ConsurfDB to download conservation information. check if https://consurfdb.tau.ac.il/ is online.", projDir)
            return None
        if not driverResult:         #failed to retrieve by sequence
            if pdbID:
                driverResult, outDriver = consurf_input_by_pdbID(pdbID, chromeOptions, consurfURL, consurfPDB, consurfChain)
            if not driverResult:
                issue("\nOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO\n"
                      "Either failed to find a consurf result for this protein by both sequence and pdbID or no pdbID was available\n"
                              "Repeat the analysis with consurfJob = True\n"
                      "OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO", projDir)
                return None

        driver = outDriver
            #await results page

            #download consurf summary
        downloadLink = WebDriverWait(outDriver, 10).until(EC.element_to_be_clickable((By.XPATH, "//a[contains(@href, 'consurf_summary.txt') and contains(@class, 'download')]")))
        downloadLink.click()

        #wait for file to download
        file = ''
        for i in range(30):
            file = [i for i in os.listdir(projDir) if i.endswith('consurf_summary.txt')]
            if file:
                break
            time.sleep(1)
        if file:
            return data_from_projDir('consurf_summary.txt', projDir)
        else:
            issue("consurf download failed in 30 seconds", projDir)

    except NoSuchElementException as e:
        issue(f"Failed while navigating ConsurfDB.\n failed element: {e}", projDir)
        return None
    except TimeoutException as e:
        issue(f"Failed while navigating ConsurfDB.\n failed element: {e}", projDir)
        return None
    driver.quit()
def main_consurfDB(projDir, pdbID = None, targetSeq = None, debug = False, consurfPDB = False, consurfChain = False):
    """runs if_not_in_projDir on req_consurfDB(), then retrieves consurfData from projDir"""
    #because req_consurfDB can't actually return
    consurfData = if_not_in_projDir(req_consurfDB, 'consurf_summary.txt', targetDir= projDir, save=False,
                                    passProjDir=True, pdbID = pdbID, targetSeq=targetSeq, debug = debug,
                                    consurfPDB = consurfPDB, consurfChain = consurfChain)
    return consurfData
#job request
def req_consurfJob(targetSeq, uniprotID, email, projDir, debug = False):
    consurfURL = 'https://consurf.tau.ac.il/consurf_index.php'
    chrome_options = Options()
    if not debug:
        chrome_options.add_argument("--headless=new")  # Run in headless mode

    driver = webdriver.Chrome(options=chrome_options)
    driver.get(consurfURL)

    # Check if the page has loaded correctly
    try:
        wait = WebDriverWait(driver, 10)
        wait.until(EC.presence_of_element_located((By.ID, "FASTA_text")))
    except TimeoutException:
        issue('Failed to load Consurf job page, check if it is online at:\n'
              'https://consurf.tau.ac.il/consurf_index.php', projDir)
        driver.quit()
        return None

    # Select the checkbox
    checkbox_span = wait.until(EC.presence_of_element_located((By.XPATH, "//label[@for='no_structure']/span")))
    checkbox_span.click()

    # Input the target sequence
    textbox = wait.until(EC.presence_of_element_located((By.ID, "FASTA_text")))
    textbox.send_keys(targetSeq)

    # Input the job name (uniprotID)
    job_name_input = wait.until(EC.presence_of_element_located((By.ID, "JOB_TITLE")))
    job_name_input.send_keys(uniprotID)

    # Input the email
    email_input = wait.until(EC.presence_of_element_located((By.ID, "user_email")))
    email_input.send_keys(email)

    # Click the submit button
    submit_button = wait.until(EC.presence_of_element_located((By.CSS_SELECTOR, ".btn_sub.big_btn")))
    submit_button.click()

    # Wait for the page to navigate to the progress page
    printWrite(f'waiting for consurf progress URL page to load.\n'
          f' DO NOT turn off the program until you see:\n'
          f'"running download_job_results..." message below:', projDir)
    wait.until(EC.url_contains("progress.php"))
    progress_url = driver.current_url

    driver.quit()
    return progress_url
def download_job_results(progress_url, projDir, debug = False):
    chrome_options = Options()
    chrome_options.add_experimental_option("prefs", {
        "download.default_directory": projDir+'\\',
        "download.prompt_for_download": False,  # To disable the download prompt
        "download.directory_upgrade": True,     # To ensure file downloads work in headless mode
        "safebrowsing.enabled": True})
    if not debug:
        chrome_options.add_argument("--headless=new")  # Run in headless mode

    driver = webdriver.Chrome(options=chrome_options)
    driver.get(progress_url)

    while True:
        try:
            finish_element = WebDriverWait(driver, 5).until(
                EC.presence_of_element_located((By.XPATH, "//span[contains(text(), 'FINISHED')]"))
            )
            if finish_element:
                break
        except TimeoutException:
            time.sleep(25)  # Wait for 30 seconds before checking again
        printWrite("Consurf jobs often take at least 20 minutes to complete.\n"
              "if you close this program and submit the same uniprotID later, the program will re-query this same consurf job.\n"
              "It may be complete then.", projDir)
    # Click the "Go to Results" button
    results_button = WebDriverWait(driver, 10).until(
        EC.presence_of_element_located((By.CSS_SELECTOR, ".btn"))
    )
    results_button.click()

    # Wait for the results page to load
    WebDriverWait(driver, 10).until(EC.url_contains("final_output.php"))

    # Download the sonsurf_summary.txt file
    download_link = WebDriverWait(driver, 10).until(EC.element_to_be_clickable(
        (By.XPATH, "//a[contains(@href, 'consurf_grades.txt') and contains(@class, 'download')]")))
    download_link.click()
    # Wait for the file to download
    file = ''
    for i in range(30):
        file = [f for f in os.listdir(projDir) if f.endswith('consurf_grades.txt')]
        if file:
            break
        time.sleep(1)
    driver.quit()
    if file:
        return data_from_projDir('consurf_grades.txt', projDir)
    else:
        issue("Consurf job download failed after 30 seconds.", projDir)
def main_ConsurfJob(targetSeq, uniprotID, email, projDir, debug = False):
    """queries consurf for a full job on the uniprotID rather than querying the database"""
    progress_url = if_not_in_projDir(req_consurfJob, 'consurf_progress_url.txt', projDir, save = True, passProjDir = True,
                                     targetSeq = targetSeq, uniprotID = uniprotID, email = email, debug = debug)
    jobData      = if_not_in_projDir(download_job_results, 'consurf_grades.txt', projDir, save = False, passProjDir = True,
                                     progress_url = progress_url, debug = debug)
    return jobData
#pick job or DB
def process_query(consurfJob, allPdbIDs, targetSeq, consurfPDB, consurfChain, uniprotID, email, projDir, debug):

    # else only req consurfDB
    consfData, consfID = '',''
    if not consurfJob:

        #if a job was performed previously
        progressURL = data_from_projDir('consurf_progress_url.txt',projDir)
        if progressURL:
            consfID = uniprotID
            consfData = data_from_projDir('consurf_grades.txt', projDir)
            if not consfData:
                consfData = download_job_results(progressURL, projDir, debug)

        #if no previous job
        else:
            #request DB
            firstPdbID = allPdbIDs[0] if allPdbIDs else None
            consfData  = main_consurfDB(projDir, pdbID=firstPdbID,
                                                     targetSeq=targetSeq, debug=debug,
                                                     consurfPDB=consurfPDB, consurfChain=consurfChain)  # query consurfDB

            consfPath = name_from_projDir('consurf_summary.txt', projDir)
            consfID = consfPath.split('\\')[-1].split('_')[0][:-1]  # pdbID used by consurf

            #if DB request fails, request a job
            if not consfData:
                printWrite('failed to fetch from consurfDB with both sequence and structures. Requesting consurfJob...', projDir)
                consurfJob = True

    # if a consurfJob is requested outside of function OR in previous block
    if consurfJob:
        consfData = main_ConsurfJob(targetSeq, uniprotID, email, projDir,
                                                  debug=debug)  # request a consurf Job
        consfID = uniprotID
    return consfData, consfID
#process
def get_consfDicts(consfData, targetSeq, projDir):
    """ parses consurf_summary.txt into:
            consfScores = {index: conservation score}
            colorScores = {index: color score}
            residues    = {index: residue}
            indexes are remapped against original target seq since sometimes consurf aligns to structure sequence and not genetic sequence

            returns: consfScores, colorScores, residues"""

    def find_ishift(targetSeq, consfSeq, consfIndexes):
        """find the index shift between the consurf sequence and the target sequence"""
        def reindex_Nones(indexMap):
            """assigns pre-1 indexes to neg values"""
            firstNonNoneIndex = next((i for i, v in indexMap.items() if v is not None), None)
            if firstNonNoneIndex is not None:
                decrementValue = 0
                for i in range(firstNonNoneIndex - 1, 0, -1):
                    indexMap[i] = decrementValue - 1
                    decrementValue -= 1
            return indexMap
        consfSeqString                              = ''.join(consfSeq)
        consfIndexMap,_,percMismatch,percGaps       = fetch_uniprot. get_index_consf_maps_percMismatch(targetSeq, consfSeqString, consfIndexes)
        finalConsfIndexMap                          = reindex_Nones(consfIndexMap)
        keyiShift                                   = min([key for key in consfIndexMap.keys() if consfIndexMap[key] != None])-1
        valiShift                                   = min([val for val in consfIndexMap.values() if val != None][:-1])-1
        iShift = valiShift - keyiShift -1
        return percMismatch, percGaps, finalConsfIndexMap, iShift

    lines  = consfData.split("\n")

    # Dynamically identify the header line based on specific keywords

    headerKeywords   = ["POS", "SEQ", "SCORE", "COLOR", "CONFIDENCE"]
    headerLineIndex = next((i for i, line in enumerate(lines) if
                              all(keyword in line for keyword in headerKeywords)), None)

    if headerLineIndex is None:
        issue(f"issue reading consurf_summary.txt\n did you modify the file?\n"
                      f" if the file is open in another program, close it and press ENTER. \n",projDir)
        input()
        headerLineIndex = next((i for i, line in enumerate(lines) if all(keyword in line for keyword in headerKeywords)), None)
        if headerLineIndex is None:
            raise ValueError("Header line not found in consurf_summary.txt.\n"
                             "first, check if https://consurfdb.tau.ac.il/ is online\n"
                             "if so, try deleting it and rerunning. If consurfDB doesn't load (it is often unavailable),\n"
                             "if you have another consurf output file, try repairing the file manually in excel.")
    # Adjusting to handle variable number of columns
    headerLine = lines[headerLineIndex].split("\t")
    headerLine = [i.strip() for i in headerLine]
    dataLines  = lines[headerLineIndex + 1:]
    # Split each line by tab and only take as many elements as there are headers
    data        = [line.split("\t")[:len(headerLine)] for line in dataLines if line]
    # Convert the adjusted data into a DataFrame
    df          = pd.DataFrame(data, columns=headerLine)

    oriConsurfDict          = df.to_dict("list")
    colsConsurfDict         = {key:value for key,value in oriConsurfDict.items() if
                               key in ['POS','SEQ','SCORE','COLOR']}                                                #grab relevant columns
    consurfDict             = {key: [str_to_intFloat(i.strip()) for i in value[1:] if i] for
                               key, value in colsConsurfDict.items()}                                                   #convert to types, skip empty row

    unshiftedConsfSeq                           = consurfDict['SEQ'][:len(consurfDict['SCORE'])]                     #omit comments at bottom SEQ line
    unshiftedConsfPositions                     = consurfDict['POS'][:len(consurfDict['SCORE'])]                     #omit comments at bottom POS line
    percMismatch, percGap, index_map, iShift    = find_ishift(targetSeq, unshiftedConsfSeq, unshiftedConsfPositions)
    consurfDict['POS']                          = [index_map.get(i,None) for i in unshiftedConsfPositions]

    errorDict               = {'%mismatch':percMismatch,'%gap':percGap}
    consurfDict['COLOR']    = [999 if type(color) == str else color for color in consurfDict['COLOR']]
    consfScores             = {consurfDict['POS'][i]: consurfDict['SCORE'][i] for i in range(len(consurfDict['POS']))}
    colorScores             = {consurfDict['POS'][i]: consurfDict['COLOR'][i] for i in range(len(consurfDict['POS']))}
    residues                = {consurfDict['POS'][i]: consurfDict['SEQ'][i] for i in range(len(consurfDict['POS']))}

    consfScores[None]       = 999                                                                                                  #missing residues are "unconserved"
    colorScores[None]       = 999
    return [consfScores, colorScores, residues, errorDict]
