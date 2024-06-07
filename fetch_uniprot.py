from Bio import Align, Entrez, Medline
from QOL_functions import *
from copy import deepcopy
import requests, re

from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import TimeoutException, NoSuchElementException

def search__fetch_entry(uniprotID, projDir):
    """checks projDir for uniprotEntry, if not, calls fetch_entry to get it, saves it to projDir, returns uniprotEntry"""
    def fetch_entry(uniprotID):
        """Fetches a UniProt entry from uniprot.org for a given UniProt ID. Returns the XML entry."""
        url         = f'https://www.uniprot.org/uniprot/{uniprotID}.xml'
        response    = requests.get(url)
        if response.status_code == 200:
            uniprotEntry        = ET.fromstring(response.content)
            return uniprotEntry
        else:
            raise Exception(f"Error fetching UniProt entry for ID {uniprotID}")
    if projDir:
        uniprotEntry = if_not_in_projDir(fetch_entry, f'{uniprotID}_uniprotEntry.xml', projDir, uniprotID = uniprotID)
    else:
        uniprotEntry = fetch_entry(uniprotID)
    return uniprotEntry
def lit_if_not_lit(projDir, uniprotTitle, uniprotIDs, uniprotEntries, email, debug):
    """checks if a literature file exists, if not, extract it from uniprotIDs and NCBI and save it to the project directory"""

    def extract_literature(uniprotID, projDir, debug = False):
        """Extracts literature from uniprot.org for a given UniProt ID. Returns a dictionary of titles and links."""
        printWrite(f'Extracting literature from {uniprotID}', projDir)
        pyName        = os.path.basename(__file__)
        chromeOptions = Options()
        if not debug:
            chromeOptions.add_argument("--headless")

        driver = webdriver.Chrome(options=chromeOptions)

        driver.get(f"https://www.uniprot.org/uniprotkb/{uniprotID}/publications")

        titles__link_directLink = {}

        try:
            WebDriverWait(driver, 10)   .until(EC.presence_of_element_located((By.CSS_SELECTOR, ".card")))
            cards                       = driver.find_elements(By.CSS_SELECTOR, ".card")

            for card in cards:
                try:
                    #find categories of citation to exclude sequence-only citations
                    categories_elements = card.find_elements(By.CSS_SELECTOR, ".decorated-list-item__title.tiny")
                    categories_title    = None
                    for element in categories_elements:
                        if element.text == "Categories":
                            categories_title = element
                            break
                    if categories_title:
                        categories_content = categories_title.find_element(By.XPATH, "./following-sibling::div")

                        if categories_content.text.strip() != "Sequences":

                            # Find the title and link within the card
                            title_element   = card.find_element(By.CSS_SELECTOR, "a[href*='/citations']")
                            links           = card.find_elements(By.CSS_SELECTOR, "a.external-link[rel='noopener'][target='_blank']")
                            if links:
                                link_element = links[0].get_attribute("href")
                                direct_link  = links[-1].get_attribute("href")
                                citation     = ''
                            else:
                                link_element, direct_link = '', ''
                            # If links are not found, extract citation and year
                            citation_element    = card.find_element(By.CSS_SELECTOR, "li > small")
                            citation            = citation_element.text
                            year                = citation.split('(')[-1].split(')')[0]
                            titles__link_directLink[title_element.text] = [link_element, direct_link, citation, year]

                            #add direct link

                            # Store the title text and link href in the dictionary
                except Exception as ex:
                    issue(f'error 1 extracting literature from {uniprotID}:\n{ex}', projDir)
                    pass

        except:
            issue(f'error 2 extracting literature from {uniprotID} (finding cards)', projDir)
            pass
        finally:
            driver.quit()
        return titles__link_directLink
    def get_enzyme_species_names(uniprotEntries):
        """extract species of origin and enzyme names from uniprotEntries"""
        import xml.etree.ElementTree as ET
        speciesNames = set()
        enzymeNames = set()
        for entry in uniprotEntries:
            ns = {'ns0': 'http://uniprot.org/uniprot'}

            scientificName = entry.find(".//ns0:organism/ns0:name[@type='scientific']", namespaces=ns)
            scientificName = scientificName.text if scientificName is not None else None
            speciesNames.add(scientificName)

            commonName = entry.find(".//ns0:organism/ns0:name[@type='common']", namespaces=ns)
            commonName = commonName.text if commonName is not None else None
            speciesNames.add(commonName)

            fullName = entry.findtext(".//ns0:protein/ns0:submittedName/ns0:fullName", namespaces=ns)
            if fullName:
                enzymeNames.add(fullName)
            else:
                recommended_name = entry.findtext(".//ns0:protein/ns0:recommendedName/ns0:fullName", namespaces=ns)
                enzymeNames.add(recommended_name)

            geneName = entry.find(".//ns0:gene/ns0:name", namespaces=ns)
            geneName = geneName.text if geneName is not None else None
            enzymeNames.add(geneName)

            # remove strain details
        speciesNames = {name.split('(')[0] for name in speciesNames if name}
        enzymeNames  = {name.split('(')[0] for name in enzymeNames if name}
        return speciesNames, enzymeNames
    def query_NCBI(projDir, enzyme_names, species_names, email):
        """
        Search for relevant literature for an enzyme based on its name and species of origin.

        Parameters:
        - enzyme (str): The name of the enzyme.
        - species (list): A list of species of origin.

        Returns:
        - list: A list of search results containing literature information.
        """

        def list_to_ORline(aList):
            lineList = [f'({name}[Text Word])' for name in aList]
            ORline = ' OR '.join(lineList)
            return ORline

        Entrez.email = email

        # Construct the search query
        methods       = ['expression', 'purification', 'mechanism', 'activity']
        enzymeORline  = list_to_ORline(enzyme_names)
        speciesORline = list_to_ORline(species_names)
        methodsORline = list_to_ORline(methods)
        query         = f"({enzymeORline}) AND ({speciesORline}) AND ({methodsORline})"

        # Perform literature search using Entrez esearch
        try:
            handle      = Entrez.esearch(db="pubmed", term=query, retmax=100)
            record      = Entrez.read(handle)
            handle      . close()

        except Exception as e:
            e = e.read()
            printWrite(f'failed to retrieve literature from NCBI with error:\n{e}', projDir)

        # Retrieve the list of PubMed IDs from the search results
        pubMedIDs  = record["IdList"]

        # Fetch the details of the retrieved articles using Entrez efetch
        handle      = Entrez    .efetch(db="pubmed", id=",".join(pubMedIDs), rettype="medline", retmode="text")
        records     = Medline   .parse(handle)

        # Process and store lit data
        litInfo = {}
        for record in records:
            if 'TI' in record.keys():
                title  = record.get('TI', '')
                pmid   = record.get('PMID', '')
                pmUrl  = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}"
                try:
                    doi    = record.get("AID", "")
                    if type(doi) == list:
                        doi   = doi[-1].replace(' [doi]', '')  # skip pmi, remove indicator from DOI
                    doiUrl = f'https://doi.org/{doi}'
                except:
                    doiUrl = ''

                year   = record.get('DP', '').split()[0]
                litInfo[title] = [pmUrl,doiUrl, "", year]
        handle.close()
        if not litInfo:
            printWrite('failed to find literature on NCBI', projDir)
            litInfo['no results'] = ['', '', '', '']
        return litInfo
    def process_literature(allLit, projDir, litFileName):
        """processes literature data to a dataframe"""
        literatureDf    = []
        usedTitles      = []
        defaultAssembly = {'Uniprot Entry': '',
                                 'year': '',
                                 'Title': '',
                                 'URL (clickable)': '',
                                 'direct (clickable) (sometimes broken)': '',
                                 'citation (no links)': ''}
        for uniprotEntry, litData in allLit.items():
            assembly                    = deepcopy(defaultAssembly)
            assembly['Uniprot Entry']   = uniprotEntry
            literatureDf                .append(assembly)                         # Append a row for the uniprotEntry

            #sort dictionary by year
            sortedLitList   = sorted(litData.items(), key=lambda x: x[1][3])
            sortedLitData   = dict(sortedLitList)

            for title, data in sortedLitData.items():
                link, directLink, citation, year = data

                linkOut         = f'=HYPERLINK("{link}", "{link}")'                                                            # clickable link
                directLinkOut   = f'=HYPERLINK("{directLink}", "{directLink}")'                                           # clickable direct link
                assembly        = deepcopy(defaultAssembly)

                if title not in usedTitles:
                    usedTitles.append(title)

                    assembly['Title']                                   = title
                    assembly['URL (clickable)']                         = linkOut
                    assembly['direct (clickable) (sometimes broken)']   = directLinkOut
                    assembly['citation (no links)']                     = citation
                    assembly['year']                                    = year

                    literatureDf.append(assembly)

                else:                                                                                                   #to make it year-sortable
                    assembly['Title'] = title
                    literatureDf      .append(assembly)

            literatureDf.append(deepcopy(defaultAssembly))                                                              # Append a blank row for separation
        return literatureDf
    def save_literature(literatureDf, projDir, litFileName):
        """saves literature dataframe to xlsx in project directory"""
        litFilePath = name_from_projDir(litFileName, projDir)
        try:
            save_to_dir(literatureDf, litFileName, projDir, force=True)  # Save the literature data to the project directory
        except:
            input(f'Failed to save {litFileName}. If it is open in another program, close it and press ENTER.')
            save_to_dir(literatureDf, litFileName, projDir, force=True)

    litFileName = f'{uniprotTitle}_literature.xlsx'
    if is_in_Dir(litFileName, projDir):                                                                             #check if literature file exists
        printWrite(f'\nfound {litFileName} in project directory', projDir)
    else:                                                                                                               #else
        printWrite(f'\ngenerating {litFileName}', projDir)
        allLit = {}

        #get literature from each uniprotEntry
        for uniprotID in uniprotIDs:                                                                                    #for uniprotID
            allLit[uniprotID] = extract_literature(uniprotID, projDir, debug)                                                   #extract literature

        #get literature from NCBI
        printWrite('querying NCBI for literature', projDir)
        speciesNames, enzymeNames   = get_enzyme_species_names(uniprotEntries)                                          # get species names of related IDs
        allLit['NCBI']      = query_NCBI(projDir, enzymeNames, speciesNames, email)                                      #query NCBI for literature
        literatureDF                = process_literature(allLit, projDir, litFileName)                          #output the literature data
        save_literature(literatureDF, projDir, litFileName)                                                             #save the literature data to the project directory
def make_projDir_store_uniprotEntry(uniprotID, uniprotEntry):
    """Creates a directory for this project, saves uniprotEntry to it, returns the path to the directory."""
    def return_obj(obj):
        return obj
    uniprotTitle    = uniprotEntry.find('.//{http://uniprot.org/uniprot}name').text     # Extract title from uniprotEntry
    title           = f"{uniprotID} {uniprotTitle}"
    prevDir        = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
    projDir         = os.path.join(prevDir, title)
    if not os.path.exists(projDir):
        os.makedirs(projDir)
    if_not_in_projDir(return_obj, outputFileName = f'{uniprotID}_uniprotEntry.xml', targetDir= projDir, obj = uniprotEntry)
    return projDir
def get_index_consf_maps_percMismatch(oriIdealSeq, oriStructSeq, structIndexes):
    """aligns crystal structure sequence against uniprot sequence and generates
     an index map of ideal indexes to xtal indexes for mapping against consurf
     params:    IdealSeq:       protein sequence from uniprot
                structSeq:      protein sequence directly extracted from xtal structure
                structIndexes:  list of indexes of structure sequence
     returns:
                indexMap:      {struct index: ideal index} with mismatches indexed according to their position
                consfMap:      {struct index: ideal index} with mismatches indexed as None so mutant residues do not receive a consfscore
                """

    pyName = os.path.basename(__file__)

    def replace_atypical_residues(proteinSequence):
        """Replace all content within parentheses with a placeholder character"""
        placeholder = 'X'
        pattern = r"\([^\)]*\)"
        # Replace all occurrences found by the pattern with the placeholder character
        modifiedSequence = re.sub(pattern, placeholder, proteinSequence)
        return modifiedSequence

    structSeq   = replace_atypical_residues(oriStructSeq)
    idealSeq    = replace_atypical_residues(oriIdealSeq)

    def align_sequences(idealSeq, structSeq):
        """aligns sequences using biopython pairwise aligner"""

        aligner     = Align.PairwiseAligner()
        aligner     .match_score            = 1
        aligner     .mismatch_score         = -2
        aligner     .open_gap_score         = -1
        aligner     .extend_gap_score       = 0
        aligner     .target_end_gap_score   = 0
        aligner     .query_end_gap_score    = 0

        aligner.mode    = 'global'  # For global alignment
        alignments      = aligner.align(idealSeq, structSeq)
        alignment       = alignments[0]
        return alignment
    alignment = align_sequences(idealSeq, structSeq)
    def extract_aligned_sequences(alignment):
        """extracts aligned sequences from alignment object"""
        pattern         = r"(target|query)\s+\d+\s+([A-Z\-]+)"
        text            = str(alignment)
        matches         = re.findall(pattern, text)

        idealAlign, structAlign = '', ''
        for match in matches:
            if match[0] == "target":
                idealAlign += match[1]
            elif match[0] == "query":
                structAlign += match[1]
        return idealAlign,  structAlign
    idealAlign, structAlign = extract_aligned_sequences(alignment)
    def _correct_substitutions(struct, ideal):
        """to be able to assign indexkeys, alignment mismatches need to be corrected to allow mismatches without gaps
        a mismatch e.g. GALWA : GTLWA will be aligned as G-ALWA : GT-ALWA. these need to be corrected to a gapless, incorrect overlap"""

        # try incrementing while in twos, since gap-character-gap could result from replacement then insertion,
        # which though rare should not be treated this route

        # consider trimming off c-terminal tags, since they wont have conservation info and likely wont resolve
        while ideal[-1] == '-':  # may help with aligning
            ideal = ideal[:-1]
            struct = struct[:-1]

        i = 0
        while i < len(struct) - 2:  # for every index
            skey, ikey = struct[i:i + 2], ideal[i:i + 2]  # two-index selections of both struct and ideal
            selectRange = 1  # size of selection
            if (skey[1] == '-' and ikey[0] == '-') or (
                    skey[0] == '-' and ikey[1] == '-'):  # if selections contain alternating gaps (e.g. insertions)
                while (skey[-1] == '-' and ikey[-2] == '-') or (
                        skey[-2] == '-' and ikey[-1] == '-'):  # while larger selection still has alternating gaps
                    selectRange += 2  # increase the size of the selection
                    skey, ikey = struct[i:i + 2 + selectRange], ideal[
                                                                i:i + 2 + selectRange]  # generate larger selections
                skey, ikey = struct[i:i + 2 + selectRange - 2], ideal[
                                                                i:i + 2 + selectRange - 2]  # remove selection that failed to show alternating gaps
                snogap, inogap = skey.replace('-', ''), ikey.replace('-',
                                                                     '')  # delete gaps in selection to allow mismatches to overlap
                struct = struct.replace(skey, snogap, 1)  # replace selected regions in original sequences with
                ideal = ideal.replace(ikey, inogap, 1)  # gapless selections
                selectRange -= 2  # decrement selection range to reindex correctly
            i += selectRange  # skip selection to next unexamined index
        return struct, ideal
    def correct_substitutions(struct, ideal):
        """
        Correct alignment mismatches between struct and ideal sequences.
        Mismatches with gaps of any length are corrected to allow mismatches without gaps.
        """
        i = 0
        while i < len(struct):
            if struct[i] == '-' or ideal[i] == '-':
                # Find the end of the gap in both sequences
                gap_end_struct = i
                while gap_end_struct < len(struct) and struct[gap_end_struct] == '-':
                    gap_end_struct += 1

                gap_end_ideal = i
                while gap_end_ideal < len(ideal) and ideal[gap_end_ideal] == '-':
                    gap_end_ideal += 1

                # Check if the gaps in both sequences have the same length
                if gap_end_struct - i == gap_end_ideal - i:
                    # Remove gaps from both sequences
                    struct_nogap = struct[:i] + struct[i:gap_end_struct].replace('-', '') + struct[gap_end_struct:]
                    ideal_nogap = ideal[:i] + ideal[i:gap_end_ideal].replace('-', '') + ideal[gap_end_ideal:]

                    # Update the sequences
                    struct = struct_nogap
                    ideal = ideal_nogap

                    # Move the index to the end of the corrected region
                    i = gap_end_struct
                else:
                    i += 1
            else:
                i += 1

        return struct, ideal
    struct, ideal = correct_substitutions(structAlign, idealAlign)
    def adjust_indexes(structAlign, structIndexes):
        """takes an index map and a structAlign (aligned structSeq), adjust index map with indexes for '-' in structSeq"""
        adjustedIndexes, indexShift = [], 0
        for i, char in enumerate(structAlign):
            if char == '-':
                adjustedIndexes.append(None)
                indexShift += 1
            else:
                adjustedIndexes.append(structIndexes[i - indexShift])

        return adjustedIndexes
    structIndexes = adjust_indexes(struct, structIndexes)
    def consf_index_maps_mismatchNum(struct, ideal, structIndexes):
        """generates a map of indexes for consf scores and an index map for mismatches"""
        consfMap  = {}                                                                                                      #does not index mismatches
        indexMap  = {}                                                                                                      #does index minsmatches
        mismatches = 0
        gaps       = 0
        ideali, structi = 1, 1                                                                                              # begin indexing map from 1
        for i in range(len(struct)):                                                                                        # for every index in mirrored sequences (aligned so same length)
            if struct[i] != '-':                                                                                                # if struct is not a gap
                if ideal[i] != '-':                                                                                                 # if ideal is not a gap
                    if struct[i] == ideal[i]:                                                                                           # if they're the same           #residues match
                        indexMap[structIndexes[i]] = ideali                                                                                # structure index points to ideal index
                        consfMap[structIndexes[i]] = ideali
                        lastMappedIdeali = ideali
                    else:                                                                                                               # if they're different          #mutation in structure
                        mismatches += 1
                        indexMap[structIndexes[i]] = ideali                                                                               # structure index points to None
                        consfMap[structIndexes[i]] = None
                    ideali += 1                                                                                                         # ideal index increases by 1
                else:                                                                                                               # if ideal is a gap                  #insertion in structure
                    gaps += 1
                    indexMap[structIndexes[i]] = None
                    consfMap[structIndexes[i]] = None
            else:                                                                                                               # if struct is a gap
                if ideal[i] != '-':                                                                                                 # if ideal is not a gap             #deletion in structure
                    gaps += 1
                    ideali += 1                                                                                                         # as above, ideal index increases by 1
            structi += 1
        return consfMap, indexMap, lastMappedIdeali, mismatches, gaps

    consfMap, indexMap, lastMappedIdeali, mismatchNum, gapNum = consf_index_maps_mismatchNum(struct, ideal, structIndexes)
    # a key to tell parse_PDB.find_termini_resolutions how many unresolved Cterm residues there are
    idealEnd         = len(idealSeq)
    indexMap[9E100] = idealEnd - lastMappedIdeali  # length of ideal sequence - last ideal index mapped to structure
    percMismatch     = round(100 * (mismatchNum / len(idealSeq)), 1)  # percentage of mismatches in the alignment
    percGaps         = round(100 * gapNum / (len(idealSeq)+len(structSeq)), 1)  # percentage of gaps in the alignment
    return indexMap, consfMap, percMismatch, percGaps
def extract_pdb_ids(projDir, uniprotEntry, targetSeq, rush = False):
    """get pdbIDs from uniprotEntry"""

    uniprotTitle    = uniprotEntry.find('.//{http://uniprot.org/uniprot}name').text     # Extract title from uniprotEntry
    pdbIDs = []

    for db_ref in uniprotEntry.iter('{http://uniprot.org/uniprot}dbReference'):
        if db_ref.get('type') == 'PDB':
            pdbID = db_ref.get('id')
            if pdbID not in pdbIDs:
                if not rush:
                    printWrite(f'\t{pdbID}', projDir)
                pdbIDs.append(pdbID)

    if pdbIDs:
        printWrite(f'found pdbIDs for {uniprotTitle}:', projDir)
    elif not rush and not pdbIDs:
        printWrite(f'no pdbIDs found for {uniprotTitle}', projDir)

    return pdbIDs
def extract_protein_sequence(uniprotEntry, projDir):
    """get protein sequence from uniprotEntry"""
    seqElements    = uniprotEntry.findall('.//{http://uniprot.org/uniprot}sequence')     # Find the sequence element within the XML
    if seqElements:
        targetSeqs            = [element.text for element in seqElements if element.text is not None]
        if targetSeqs:
            sequence  = targetSeqs[0]
            targetSeq = ''.join(sequence.split())        # Remove whitespace
        else:
            raise Exception(f"Failed to find protein sequence in UniProt entry {uniprotEntry}")
    else:
        issue(f'failed to find a sequence in {uniprotEntry}',projDir)

    signalLength = extract_signal_seq_length(uniprotEntry, projDir)
    targetSeq    = targetSeq[signalLength:]  # Remove the signal peptide from the sequence

    return targetSeq


def extract_signal_seq_length(uniprotEntry, projDir):
    """Extract lengths of 'transit peptide' or 'signal peptide' from a UniProt entry."""
    ns = {'ns0': 'http://uniprot.org/uniprot'}  # Define the namespace for the XML
    transitPeptides = uniprotEntry.findall(".//ns0:feature[@type='transit peptide']", namespaces=ns)
    signalPeptides  = uniprotEntry.findall(".//ns0:feature[@type='signal peptide']", namespaces=ns)
    peptideElements = transitPeptides + signalPeptides
    signalLengths = []

    # Iterate over found peptide features
    for feature in peptideElements:
        startIndex = feature.find(".//ns0:begin", namespaces = ns)
        endIndex   = feature.find(".//ns0:end", namespaces = ns)

        if None not in [startIndex, endIndex]:
            begin   = int(startIndex.attrib['position'])
            end     = int(endIndex.attrib['position'])
            length = end - begin + 1
            signalLengths.append(length)
        else:
            issue(f'invalid signal peptide entry in {uniprotEntry}',projDir)

    #issue handle
    if len(signalLengths) > 1:
        issue(f'multiple signal peptide entries in {uniprotEntry}\n'
              f' defaulting to first detected signal peptide.',projDir)
    if signalLengths:
        signalLength = signalLengths[0]
        printWrite(f"protein contains a signal peptide of length: {signalLength}\n"
                   f"this sequence is removed from the protein sequence in Terminator analysis", projDir)
    else:
        signalLength = 0
    return signalLength


def search__downloadPDBs(pdbIDs, pdbDir, projDir):
    """goes through pdbIDs, checks if they're a mutant, saves nonmutants to project directory"""
    def downloadPDB(pdbID, projDir):
        """downloads pdb file from rcsb.org"""
        pdbFileName   = f'{pdbID}.pdb'
        pyName          = os.path.basename(__file__)
        response        = fetchUrl(f'https://files.rcsb.org/download/{pdbID}.pdb', '')
        if response.status_code == 200:
            printWrite(f"Downloaded: {pdbFileName} from rcsb.org", projDir)
            return response.content
        else:
            issue(f"Error ({response.status_code}): Unable to download PDB file for ID {pdbID}", projDir)
            return None
    for pdbID in pdbIDs:
        pdbData     = if_not_in_projDir(downloadPDB, outputFileName = f'{pdbID}.pdb', exactName = True, targetDir=pdbDir, projDir = projDir, pdbID = pdbID)