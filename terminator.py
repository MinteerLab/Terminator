import fetch_uniprot, request_BLAST, parse_BLAST, fetch_consurf, parse_PDB, fetch_swissmodel, process_result
import traceback, sys, pandas, os, copy
from sortedcontainers import SortedDict
from QOL_functions import *

#required inputs
uniprotIDs     = ['G0SGU4', 'P22256', 'Q96EY8', 'A0A892IHP6', 'Q47PU3', 'P53704', 'D0VWQ0', 'G8LK72', 'B1VK30', 'H7C697', 'G9F1Y9', 'P02787', 'P76342', 'O46414']#['P94559', 'O31590', 'P71018', 'Q06795', 'P69816', 'P0C0L2', 'P55343', 'Q9X0G9', 'P75467', 'P75223', 'P75232', 'Q9WZB9']
#OR
sequence      = ''

#required
email         = 'r.gerulskis@utah.edu'

#optional inputs
rushes        = ''#['SWISSMODEL_2GSD', '3FN4', '2GSD']#['1W3N','SWISSMODEL_1W3I']      #PDBIDs to rush process. If a PDBID is not found in uniprot entries, it will not be processed.
forcePDB      = ''      #to force a PDBID to be run, even if it's not in the uniprot entries
consurfPDB    = ''      #PDBID to search consurfDB with
consurfChain  = ''      #chainID to search consurfDB with
resetSwm      = False   #delete old SWISSMODEL data to retry if a previous attempt failed
resetConsurf  = False   #for inputing a new consurfPDB and consurfChain, if true, deletes previous consurfDB results and queries a new DB request.
consurfJob    = False    #see below
cavitySearch  = True
debug         = False
swmOnly       = True
targetChainName = ''   #chain to process of target PDB
resetCavities   = False
getColorCmds    = False #prints series of color commands to recolor pdb file by consurf colors
#if false, requests from consurfDB, if true, requests a full consurf job. full jobs take much longer and use consurf computing power.
#a full job should only be requested if a completed termini_scores.xlsx shows that the consurf sequence retrieved does not
#sufficiently cover the crystal structures queried, OR if the program outputs an error stating that consurfDB retrievals failed by PDB and sequence.

######set inputs to uppercase
def setUpper(object):
    """sets object or group of objects to object.upper()"""
    if not object:
        return object
    else:
        if type(object) == str:
            return object.upper()
        elif type(object) == list:
            return [subObject.upper() for subObject in object]
        elif type(object) == dict:
            return {key.upper(): setUpper(value) for key, value in object.items()}

uniprotIDs      = setUpper(uniprotIDs)
rushes          = setUpper(rushes)
sequence        = setUpper(sequence)
consurfPDB      = setUpper(consurfPDB)
consurfChain    = setUpper(consurfChain)

def main(uniprotID, sequence, email, rushes, forcePDB, consurfPDB, consurfChain, resetSwm, resetConsurf, consurfJob,
         cavitySearch, debug, swmOnly, targetChainName, resetCavities, defaultStructData, getColorCmds):
    def ID_siblings_by_sequence(sequence, email, seqToUniprots = None, projDir = None, cutoff = 80):
        """Run the program using a protein sequence. If the sequence is in recentUniprots, return close uniprotIDs from recentUniprots.
         Otherwise, BLAST the sequence and return closeIDs and tempDir containing blastResults"""

        if not projDir:
            tempDir = make_tempDir()
        else:
            tempDir = projDir

        if seqToUniprots and sequence in seqToUniprots:
            closeUniprotIDs = seqToUniprots[sequence]
            return closeUniprotIDs
        else:
            blastResults    = request_BLAST.search__get_BLAST(tempDir, email, sequence)
            closeUniprotIDs = parse_BLAST.search__find_close_hits(tempDir, blastResults, cutoff)
            return closeUniprotIDs
    def search_by_sequence(sequence, email, seqToUniprots, prevDir):
        """Run the program using a protein sequence. If the sequence is in recentUniprots, return close uniprotIDs from recentUniprots.
         Otherwise, BLAST the sequence and return closeIDs and tempDir containing blastResults"""
        closeUniprotIDs                 = ID_siblings_by_sequence(sequence, email, seqToUniprots)                             # check if sequence is in recent_uniprots, if so return closeUniprotIDs, else BLAST sequence and return closeUniprotIDs
        seqToUniprots[sequence]         = closeUniprotIDs                                                                  # update seqToUniprots with sequence:closeUniprotIDs
        save_to_dir(seqToUniprots, 'recent_uniprots.json', prevDir, force=True)                                     # update seqToUniprots file
        return closeUniprotIDs, seqToUniprots
    def ID_siblings_by_everything(close_uniprotIDs_from_seqSearch, uniprotID, seqToUniprots,
                                  alreadyBlasted, targetSeq, email, prevDir, projDir, cutoff):
        """tries to find closeUniprotIDs in seqToUniprots, if not, runs ID_siblings_by_sequence,
            then, runs fetch_uniprot.search__fetch_entry on closeUniprotIDs and returns closeUniprotEntries"""

        if close_uniprotIDs_from_seqSearch:
            closeUniprotIDs = close_uniprotIDs_from_seqSearch                                                               #if searched by sequence, get closeUniprotIDs from there
        else:
            closeUniprotIDLists = [uniprotIDs for uniprotIDs in seqToUniprots.values() if uniprotID in uniprotIDs]          #else try to get closeUniprotIDs from seqToUniprots
            closeUniprotIDs     = closeUniprotIDLists[0] if closeUniprotIDLists else None

        if not closeUniprotIDs and not alreadyBlasted:                                                                      #if neither worked and you havent already tried blasting once
            closeUniprotIDs = ID_siblings_by_sequence(targetSeq, email, None, projDir=projDir, cutoff=cutoff)                      #run ID_siblings_by_sequence on the sequence

        seqToUniprots[targetSeq]         = closeUniprotIDs                                                                  # update seqToUniprots with sequence:closeUniprotIDs
        save_to_dir(seqToUniprots, 'recent_uniprots.json', prevDir, force=True)                                             # update seqToUniprots file
        return closeUniprotIDs
    try:
        prevDir         = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
        if is_online():
            sibling_IDs_via_seqSearch  = []
            seqToIDs                   = data_from_projDir('recent_uniprots.json', prevDir) or {}
            alreadyBLASTed             = False
            sequence_similarity_cutoff = 80
            #get uniprotID if sequence input
            if sequence and not uniprotID:                                                                                  #if sequence is given and uniprotID is not
                sibling_IDs_via_seqSearch, seqToIDs  = search_by_sequence(                                                  #search for uniprotEntry by sequence
                                                        sequence, email, seqToIDs,
                                                        prevDir, sequence_similarity_cutoff)

                alreadyBLASTed = True                                                                                       #to prevent multiple blasts if no siblings are found
                uniprotID      = sibling_IDs_via_seqSearch[0]                                                               #set uniprotID to closest closeUniprotID
            if uniprotID:

                #get uniprot entry info and project directory
                projDir             = projDir_if_projDir(uniprotID)                                                         #check if a project directory exists for the target UniProt ID
                mainUniprotEntry    = fetch_uniprot     .search__fetch_entry(uniprotID, projDir)                            #get the main_consurfDB target uniprotEntry
                mainUniprotTitle    = mainUniprotEntry  .find('.//{http://uniprot.org/uniprot}name').text                   #Extract title from uniprotEntry
                mainTitle           = f"{uniprotID}_{mainUniprotTitle}"
                targetSeq           = fetch_uniprot.extract_protein_sequence(mainUniprotEntry, projDir)                              # get the target sequence
                if not projDir:                                                                                             #if projDir doesnt already exist
                    projDir         = fetch_uniprot. make_projDir_store_uniprotEntry(uniprotID, mainUniprotEntry)           #create a project directory and store the uniprotEntry
                del mainUniprotEntry
                printWrite('', projDir, firstWrite=True)                                                                                     #create an output file
                printWrite(f'OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO\n'
                      f'project directory:\n{projDir}\n'
                      f'OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO', projDir)

                #get sibling IDs
                allUnipIDs         = ID_siblings_by_everything(sibling_IDs_via_seqSearch,                                  #get related uniprot IDs
                                        uniprotID, seqToIDs, alreadyBLASTed, targetSeq,
                                        email, prevDir, projDir, sequence_similarity_cutoff)

                #print sibling IDs
                if allUnipIDs:
                    IDsOut = '\t'.join(allUnipIDs)
                    printWrite(f'uniprotIDs with >{sequence_similarity_cutoff}% sequence similarity to {uniprotID}:\n{IDsOut}', projDir)
                else:
                    printWrite(f'No other uniprotIDs found with '
                          f'>{sequence_similarity_cutoff}%'
                          f' sequence similarity to {uniprotID}.', projDir)
                allUnipIDs = [uniprotID] + allUnipIDs if uniprotID not in allUnipIDs else allUnipIDs

                #get sibling entries
                def get_sibling_entries(sibling_IDs, projDir):
                    """retrieve uniprot entries from uniprotIDs"""
                    siblingEntries = []
                    for i, sibling_ID in enumerate(sibling_IDs):
                        siblingEntry = fetch_uniprot.search__fetch_entry(sibling_ID, projDir)
                        if siblingEntry not in siblingEntries:
                            siblingEntries.append(siblingEntry)
                            printWrite(f'>{sibling_ID} entry {i+1}/{len(sibling_IDs)}', projDir)
                    return siblingEntries
                siblingEntries = get_sibling_entries(allUnipIDs, projDir)

                #get sibling literature
                fetch_uniprot .lit_if_not_lit(projDir, mainTitle, allUnipIDs, siblingEntries, email, debug)

                #get pdbIDs for all structures from sibling entries
                allPdbIDs = []
                for siblingEntry in siblingEntries:
                    allPdbIDs += fetch_uniprot.extract_pdb_ids(projDir, siblingEntry, targetSeq=targetSeq, rush = rushes)
                del siblingEntries
                allPdbIDs = list(set(allPdbIDs))

                #download all PDB structures
                pdbDir          = fetch_uniprot     .check_subdir(projDir, 'PDBs')
                fetch_uniprot   .search__downloadPDBs(pdbIDs = allPdbIDs, pdbDir=pdbDir, projDir=projDir)

                #get consurf data
                if resetConsurf:
                    delete_from_projDir('consurf_summary.txt', projDir)
                consfData, consfID = fetch_consurf.process_query(consurfJob, allPdbIDs,
                                                                 targetSeq, consurfPDB,
                                                                 consurfChain, uniprotID,
                                                                 email, projDir, debug = debug)

                #process consurf data into dictionaries
                consfDicts = fetch_consurf.get_consfDicts(consfData, targetSeq, projDir) if consfData else None
                if not consfData:
                    issue('failed to get a consurf result. Only steric analysis will be performed.', projDir)
                del consfData

                #get swissmodel structure
                swmGMQE, swmPdbID   = fetch_swissmodel  .main_swissmodel(allUnipIDs, projDir, pdbDir, resetSwm)

                #print PDB and swissmodel results
                if not swmPdbID:
                    issue('failed to get a swissmodel result. no structure will be processed from swissmodel.', projDir)
                    if not allPdbIDs:
                        pass #make a modeling request from swissmodel
                elif not allPdbIDs:
                    noStructMessage = \
                        f"""OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    Not a single valid pdbID found across all siblings.
    however, swissmodel did provide {swmPdbID}, a different protein predicted to have a similar structure.
    All structural and conservation parsing will be done using the model provided,
    with the target sequence fit onto a posssibly foreign crystal structure.
    Ensure that the GMQE and QSQE in the output are above 0.7, as this indicates reliable tertiary and 
    quaternary structure prediction, respectively.
    QSQE may be less critical, as QSQE is often calculated as being below 0.7 even for models developed using
    the correct protein's crystal structure. In other words,these scores rely on the quality of the actual crystal structure,
    and not only on the quality of fitting the target sequence to the retrieved crystal structure.
    OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO"""
                    printWrite(noStructMessage, projDir)
                allPdbIDs. append('SWISSMODEL_'+swmPdbID) if swmPdbID else None

                #process rush request
                swms = [pdbID for pdbID in allPdbIDs if 'SWISSMODEL_' in pdbID]
                if swmOnly and swms and not rushes:
                    allPdbIDs = swms
                if rushes and all([rush in allPdbIDs for rush in rushes]) or forcePDB:
                    allPdbIDs      = rushes                                                                                 #process only 1 sequence when debugging
                elif rushes:
                    badRushes = [rush for rush in rushes if rush not in allPdbIDs]
                    printWrite(f"rush pdbID(s) {', '.join(badRushes)} is/are not structures of {uniprotID}."
                                  f" Did you mean to set rush=False instead?", projDir)
                    raise ValueError('invalid rush selection')
                #finally process structure data

                structureDatas  = parse_PDB.iterate_structures(                                                             #perform structural analyses on all PDBs and swissmodel
                                    allPdbIDs, projDir, pdbDir, consfDicts, targetSeq,
                                    uniprotID, consfID, cavitySearch, defaultStructData,
                                    targetChainName, resetCavities, getColorCmds)
                parse_PDB       . process_structure_datas(structureDatas, defaultStructData, projDir, mainTitle, 'residues')                   #output data aligned by residues
                parse_PDB       . process_structure_datas(structureDatas, defaultStructData, projDir, mainTitle, 'consfScores')                #output data aligned by consfScores
            #if uniprotIDs:
            #    process_results.export_multiple_results(uniprotIDs, defaultStructData)
        else:                                                                                                               #check if online. due to online i/o downstream, probably a cleaner way to do this per-request
            print('terminator requires an internet connection.\n'
                  'Connect to the internet and restart the program.')
    except Exception:
        full_traceback = get_traceback()
        printWrite(f'failed processing {uniprotID} with error:\n{full_traceback}', projDir)

if __name__ == "__main__":
    defaultStructData = {'indexes': [],
                         'structureName': '',
                         'chain ID': '',
                         'N unresolved residues': '',
                         'C unresolved residues': '',
                         'residues': ['failed to find chains'],
                         'consfScores': ['failed to find chains'],
                         'N residue burial': '',
                         'C residue burial': '',
                         'N residue normalized burial': '',
                         'C residue normalized burial': '',
                         'N interface proximity score': '',
                         'C interface proximity score': '',
                         '%mismatch': '',
                         '%gap': '',
                         'N cavities': '',
                         'C cavities': '',
                         'swissmodel GMQE': '',
                         'swissmodel QSQE': ''}

    if uniprotIDs or sequence:
        uniprotIDs = [uniprotIDs] if type(uniprotIDs) == str else uniprotIDs
        for i,uniprotID in enumerate(uniprotIDs):
            print(f'>{uniprotID} uniprot entry {i+1}/{len(uniprotIDs)}')
            def process_rush(uniprotID, rushes):
                if rushes:
                    if isinstance(rushes, str):
                        return [rushes]
                    elif isinstance(rushes, list):
                        return rushes
                    elif isinstance(rushes, dict):
                        pdbs = rushes.get(uniprotID, '')
                        if isinstance(pdbs, str):
                            pdbs = [pdbs] if type(pdbs) == str and pdbs else pdbs
                        return pdbs
                else:
                    return rushes
            theseRushes = process_rush(uniprotID, rushes)
            try:
                main(uniprotID, sequence, email, theseRushes, forcePDB, consurfPDB, consurfChain, resetSwm, resetConsurf, consurfJob,
                    cavitySearch, debug, swmOnly, targetChainName, resetCavities, defaultStructData, getColorCmds)
            except Exception as e:
                full_traceback              = get_traceback()
                print(f'failed processing {uniprotID} with error:\n{full_traceback}')
        process_result.export_multiple_results(uniprotIDs, defaultStructData)
    else:
        print('you must enter a uniprotID or sequence')

