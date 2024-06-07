import fetch_uniprot

from QOL_functions import *


### QOL functions
def request_job_id(projDir, email, targetSeq):
    """submit a blast request, returns the jobID"""
    requestData = {
        'email': (None, email),
        'program': (None, 'blastp'),
        'matrix': (None, 'BLOSUM62'),
        'alignments': (None, '50'),
        'scores': (None, '50'),
        'exp': (None, '10'),
        'filter': (None, 'F'),
        'gapalign': (None, 'false'),
        'compstats': (None, 'F'),
        'align': (None, '0'),
        'stype': (None, 'protein'),
        'sequence': (None, targetSeq),
        'database': (None, 'uniprotkb'),
    }
    output = requests.post('https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/run', files=requestData)
    jobID = output.text
    printWrite(f"BLAST job submitted. Job ID: {jobID}", projDir)
    return jobID
def check_blast_loop(jobID, projDir):
    """if if_not_in_projDir(request_job_id()) returned a job_id, check status until finished, save to blastResults.xml,
     \nand return the xml, else (if the function returned an xml) return the xml"""
    def check_blast_status(jobID):
        return fetchUrl(f"https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/status/{jobID}", 'text')

    def fetch_blast_results(jobID):
        strBlast = fetchUrl(f"https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/{jobID}/xml", 'text')
        xmlBlast = ET.fromstring(strBlast)
        return xmlBlast

    checkInterval = 30  # Interval in seconds to check the status
    if type(jobID) == str: # If jobID is a string, it's the job ID
        totalTime = 0
        while True:
            status = check_blast_status(jobID)
            if status == "FINISHED":
                printWrite("BLAST job finished. Fetching results...", projDir)
                results = fetch_blast_results(jobID)
                save_to_dir(results, 'blastResults.xml', projDir)
                return results
            elif status in ["QUEUED", "RUNNING"]:
                printWrite(f"Job is {status}. Checking again in 30 seconds...", projDir)
                time.sleep(checkInterval)
                totalTime += checkInterval
                if totalTime > 119:
                    printWrite(f"BLAST job has not completed in over 2 minutes, if this is a uniprotID based search, and you need to end this run,\n"
                          f"restarting later will query this same BLAST request, and it may be complete then.", projDir)
            elif status == "NOT_FOUND":
                printWrite(f"Job not found. Old job likely expired. Requesting a new job id.", projDir)
                return False
            else:
                printWrite(f"in request_BLAST, Error or unknown BLAST status: {status}", projDir)
                return False
    else:                # If jobID is not a string, request_job_id returned the actual blast xml from projDir
        return jobID
def search__get_BLAST(projDir, email, targetSeq):
    """get BLAST results from project directory, or request them if not found"""
    if is_in_Dir('blastResults.xml', projDir):                                      #if blast results in the project directory
        blastResults = data_from_projDir('blastResults.xml', projDir)                       #get the blast results from the project directory
    else:                                                                                       #else no blast results in the project directory
        job_id          = if_not_in_projDir(request_job_id, 'blastJobID.txt',  #check for a job id
                                            targetDir=projDir, email=email, targetSeq=targetSeq, passProjDir=True)
        blastResults    = check_blast_loop(job_id, projDir)                                          #check for blast results with that job id
        if blastResults is False:                                                                       #if unable to get blast results
            job_id          = request_job_id(email, targetSeq)                                              #request a new job id
            save_to_dir(job_id, 'blastJobID.txt', projDir, force = True)                             #save the job id to the project directory
            blastResults    = check_blast_loop(job_id, projDir)                                             #check for blast results with that job id
            if blastResults is False:                                                                       #if unable to get blast results again
                raise Exception("Error: unable to get BLAST results on second attempt.")                        #raise an exception

    return blastResults
