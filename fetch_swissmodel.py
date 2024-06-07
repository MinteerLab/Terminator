from QOL_functions import *

def search__query_swissmodel(uniprotID, projDir):
    """query SwissModel for a given UniProt ID if swissmodel file isn't in the project directory"""
    def query_swissmodel(uniprotID):
        """
        Query SwissModel for a given UniProt ID.

        Parameters:
        uniprot_id (str): The UniProt ID to query.

        Returns:
        dict: The JSON response from the SwissModel API.
        """
        base_url = "https://swissmodel.expasy.org/repository/uniprot/"
        response = requests.get(f"{base_url}{uniprotID}.json")
        if response.status_code == 200:
            return response.json()
        else:
            return {"error": "Could not retrieve data from SwissModel."}
    fileName    = f"{uniprotID}_swissmodel.json"
    swmData  = if_not_in_projDir(query_swissmodel, outputFileName=fileName, targetDir=projDir, uniprotID=uniprotID)
    return swmData
def search__pdb_from_swissmodel(swmData, pdbDir, projDir):

    def pdb_from_swissmodel(structure, gmqe, qsqe, _projDir):
        """
        Download a PDB file from a URL contained in a JSON object and save it to a specified directory.

        Parameters:
        json_data (dict): JSON object containing the URL under a "coordinates" key.
        save_directory (str): The directory where the PDB file should be saved.
        """
        oligo       = structure["oligo-state"]
        url         = structure.get("coordinates")

        if url:
            try:
                response        = requests.get(url)
                response        . raise_for_status()  # Raises an HTTPError if the response status code is 4XX/5XX
                swmData         = response.content.decode('utf-8')
                quatLine        = f'REMARK 350 SOFTWARE DETERMINED QUATERNARY STRUCTURE: {oligo}'
                gmqeLine        = f'gmqe: {gmqe}'
                qsqeLine        = f'qsqe: {qsqe}'
                pdbData         = f'{quatLine}\n{gmqeLine}\n{qsqeLine}\n{swmData}'
                pdbData_bytes   = pdbData.encode('utf-8')
                return pdbData_bytes
            except requests.exceptions.RequestException as e:
                issue(f"Error downloading file from swissmodel: {e}", _projDir)
        else:
            issue("swissmodel result does not contain a url for the PDB file.", _projDir)

    targetStructure = None
    bestGMQE        = 0
    bestQSQE        = 0
    for aStructure in swmData['result']['structures']:
        if 'template_qsqe' in aStructure and aStructure['template_qsqe'] > bestQSQE:
            bestQSQE = aStructure['template_qsqe']
        if 'gmqe' in aStructure and aStructure['gmqe'] > bestGMQE:
            bestGMQE = aStructure['gmqe']
            targetStructure = aStructure

    if not targetStructure:
        issue("No structures found in SwissModel data.", projDir)
        return [None, None]
    pdbName     = targetStructure.get("template").split('.')[0].upper()
    gmqe        = targetStructure["gmqe"]
    qsqe        = targetStructure["template_qsqe"]
    filename    = f"SWISSMODEL_{pdbName}.pdb"
    if_not_in_projDir(pdb_from_swissmodel, outputFileName=filename, targetDir=pdbDir,
                      structure=targetStructure, gmqe = gmqe, qsqe = qsqe, _projDir=projDir)
    return [gmqe, pdbName]
def main_swissmodel(allUnipIDs, projDir, pdbDir, resetSwm):
    """querys swissmodel with all sibling uniprotIDs until a structure is returned"""
    for unipID in allUnipIDs:
        if resetSwm:
            delete_from_projDir('_swissmodel.json', projDir)
        swmData        = search__query_swissmodel(unipID, projDir)
        gmqe, swmPdbID = search__pdb_from_swissmodel(swmData, pdbDir, projDir)
        if swmPdbID:
            return gmqe, swmPdbID
        else:
            issue(f'{unipID} failed to retrieve a SWISSMODEL.', projDir)
            resetSwm = True
    return None, None