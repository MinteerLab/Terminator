
from collections import defaultdict
from QOL_functions import *
import copy

### QOL functions
class Hit:
    def __init__(self, attrs = None):
        self.number, self.database, self.ac, self.length, self.start, self.end, self.OS, self.OX, self.GN, self.SV, self.score, self.bits, self.expectation, self.identity, self.positives, self.gaps, self.pattern, self.querySeq, self.matchSeq, self.querySeqStart, self.querySeqEnd, self.matchSeqStart, self.matchSeqEnd, self.scores = [None] * 24
        self.PE = 5 #for hits without a PE score, default to a poor score
        [setattr(self, key, attrs[key]) for key in attrs.keys()]
    def scoreSelf(self, fullQueryLen):
        matchScore  = 1/ self.identity/100     #1/identity (as a decimal) so higher identity contributes lower score

        startGaps       = [{}] *  (self.querySeqStart - 1)                            #empty values for where match wasn't alligned to query
        scores          = [{residue: matchScore} for residue in self.matchSeq]      #score for each residue in the match
        endGaps         = [{}] * (fullQueryLen - len(startGaps + scores))           #empty values for where match wasn't alligned to query
        self.scores     = startGaps + scores + endGaps

    def help(self, attr = None):
        descriptions = {
    "number": "Ordinal number of the hit.",
    "database": "Database where the hit was found.",
    "ac": "Accession number of the hit.",
    "length": "Length of the hit sequence.",
    "start": "Start position of the alignment on the hit sequence.",
    "end": "End position of the alignment.",
    "OS": "Organism species from which the sequence originates.",
    "OX": "Taxonomic identifier for the species.",
    "GN": "Name of the gene associated with the sequence.",
    "PE": "Protein existence score.",
    "SV": "Version number of the sequence in the database.",
    "score": "Alignment score indicating quality.",
    "bits": "Normalized score in a logarithmic scale.",
    "expectation": "Number of hits expected by chance in a database search (lower is better).",
    "identity": "Percentage of identical residues in alignment.",
    "positives": "Percentage of similar residues in alignment.",
    "gaps": "Number of gaps in alignment.",
    "pattern": "Alignment pattern showing matches, mismatches, and gaps.",
    "queryseq": "Aligned sequence of the original requested sequence.",
    "matchseq": "Aligned sequence of THIS hit.",
    "queryseqstart": "Start position of the alignment on the query sequence.",
    "queryseqend": "End position of the alignment on the query sequence.",
    "matchseqstart": "Start position of the alignment on the hit sequence.",
    "matchseqend": "End position of the alignment on the hit sequence.",
    "scores": "A list of 1 element dictionaries, each dictionary is a position in the sequence with the residue as the key and 1/identity as the value. "
    }
        if not attr:
            print(f"request an attribute to define or request 'all' \n available attributes are {[key for key in descriptions.keys()]}")
        elif attr == "all":
            for key in descriptions.keys():
                print(f"{key}\t\t{descriptions[key]}")
        elif attr in descriptions.keys():
            print (f"{attr}\t\t{descriptions[attr]}")
        else:
            print(f"{attr} is not a valid attribute, \n available attributes are {[key for key in descriptions.keys()]}")
def hits_from_blast(blastResultsData):
    """get_consfDicts the BLAST XML file and extracts information from each hit element.
    Each hit's attributes and information from subordinate elements are extracted."""

    def parse_description(description):
        """split description into attributes"""
        parts = description.split()
        parsed = {}
        currentKey = None
        currentValue = []

        for part in parts:
            if "=" in part:
                # Save the previous key-value pair if there was one
                if currentKey is not None:
                    parsed[currentKey] = " ".join(currentValue)
                # Start a new key-value pair
                currentKey, value = part.split("=", 1)
                currentValue = [value]
            else:
                currentValue.append(part)

        # last key-value pair to the dictionary
        if currentKey is not None:
            parsed[currentKey] = " ".join(currentValue)

        return parsed
    namespace = "{http://www.ebi.ac.uk/schema}"
    hits = []
    # Find all 'hit' elements under the 'hits' root
    for hit in blastResultsData.findall(f".//{namespace}hit"):
        hitDict = hit.attrib  # Store hit attributes

        # Assuming there's only one 'alignments' element and one 'alignment' element within each 'hit'
        alignment = hit.find(f"{namespace}alignments/{namespace}alignment")

        # Process each child of 'alignment'
        for child in alignment:
            # Special handling for 'querySeq' and 'matchSeq'
            if child.tag.endswith('querySeq') or child.tag.endswith('matchSeq'):
                baseKey = child.tag.split('}')[1]  # Extract the key name (e.g., 'querySeq', 'matchSeq')
                hitDict[f"{baseKey}Start"] = child.attrib['start']
                hitDict[f"{baseKey}End"] = child.attrib['end']
                hitDict[baseKey] = child.text
            else:
                hitDict[child.tag.split('}')[1]] = child.text  # Store other children with tag name as key

        descAttr      = parse_description(hitDict['description'])
        hitDict      = {**hitDict, **descAttr}
        for key, value in hitDict.items():
            hitDict[key] = str_to_intFloat(value)
        hit = Hit(hitDict)
        hits.append(hit)
        hit.scoreSelf(fullQueryLen= hits[0].length)  # Calculate the score for the hit
    return hits
def search__find_close_hits(projDir, blastResults, minIdentity):
    def find_close_hits(blastResults, minIdentity):
        hits            = hits_from_blast(blastResults)
        closeHits       = [hit for hit in hits if hit.positives > minIdentity and hit.PE < 3]
        closeHitIDs     = [hit.ac for hit in closeHits]
        return closeHitIDs
    closeHitIDs = if_not_in_projDir(find_close_hits, outputFileName='closeHitIDs.json', targetDir=projDir, blastResults=blastResults, minIdentity=minIdentity)
    return closeHitIDs