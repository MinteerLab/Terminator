import copy, os
import pandas as pd
from sortedcontainers import SortedDict
from openpyxl import load_workbook

class AScore():
    def __init__(self, name, value):
        self.name       = name
        self.value      = value

        adjValue, judgement = self.judge()
        self.adjValue   = adjValue
        self.judgement  = judgement

    def judge(self):
        value           = self.value
        name            = self.name
        complexJudgeMags = { 'volume': 200,
                            'consfScore': -1,
                            'activeSiteScore': -5,
                            'conservedSiteScore': -5,
                             'interface proximity score':30,
                             'residue normalized burial':0.4}


        cavJudgeMagKey = name[2:]
        #normalized burial is judged 1 only if its > 0.4
        if value == 'residue normalized burial':
            judgeValue      = complexJudgeMags[cavJudgeMagKey]
            if value > 0.4:
                value *= 5
        #other values are judged at least 1 if nonzero, then additionally according to complexJudgeMags
        elif value != None:
            judgeValue      = value + complexJudgeMags[cavJudgeMagKey] if value else 0
        #none values (no cavities) are judged 0
        else:
            judgeValue      = 0

        judgement       = round(judgeValue,0) // complexJudgeMags[cavJudgeMagKey] if value else 0
        adjValue        = value if judgement else 0
        return adjValue, judgement

    def __repr__(self):
        return self.adjValue

def export_multiple_results(uniprotIDs, defaultStructData):
    """when running multiple UniprotIDs, data are exported to a single file"""
    defaultStructData = copy.deepcopy(defaultStructData)

    del defaultStructData['indexes']
    del defaultStructData['residues']
    del defaultStructData['consfScores']


    #get file path for each uniprotID
    def get_filePaths(uniprotIDs):
        # Get the directory of the current script
        script_dir = os.path.dirname(os.path.abspath(__file__))

        # Get the parent directory of the script directory
        directory = os.path.dirname(script_dir)

        # Create a SortedDict to store the termini scores files
        termini_scores_files = SortedDict({uid: None for uid in uniprotIDs})

        # Iterate over all subdirectories and files in the directory
        for root, dirs, files in os.walk(directory):
            for file in files:
                # Check if the file ends with '_termini_scores.xlsx'
                if file.endswith('_termini_scores.xlsx'):
                    file_path = os.path.join(root, file)

                    # Find the matching UniProt ID for the file
                    for uid in uniprotIDs:
                        if uid in root:
                            termini_scores_files[uid] = termini_scores_files[uid] if termini_scores_files[uid] else []
                            termini_scores_files[uid] .append(file_path)

        duplicateFiles = [key for key,value in termini_scores_files.items() if value and len(value) > 1]
        if duplicateFiles:
            print(f"Duplicate files found for {duplicateFiles}!:")
            for key in duplicateFiles:
                for value in termini_scores_files[key]:
                    print(f"\t{value}")
            input('delete or move the undesired duplicate file out of the Terminator folder, then restart the program.\n')

        incompleteRuns = [key for key in termini_scores_files.keys() if not termini_scores_files[key]]
        if incompleteRuns:
            print(f'the following runs were incomplete (likely awaiting consurfJob) and will be skipped in complete_results file:\n{incompleteRuns}')
        termini_scores_files = {key:value[0] for key, value in termini_scores_files.items() if key not in incompleteRuns}

        return termini_scores_files
    filePaths = get_filePaths(uniprotIDs)

    #extract and process data
    def extract_data_from_excel(file_path, defaultStructData):
        """extracts data from filepath into a {pdbID: {defaultStructData Key: value} of the file}"""

        #extract data from file_path excel
        fileDict = {}
        try:
            data         = pd.read_excel(file_path)
        except FileNotFoundError:
            print(f"File not found: {file_path}")
            return {}

        #index columns according to col1
        rowIndexes = data.set_index(data.columns[0])

        # Find the index of the CONSURF column to exclude it and subsequent columns
        consurfIndex = next(i for i, col in enumerate(data.columns) if col.startswith("CONSURF_"))-1
        validColumns = rowIndexes.iloc[:, :consurfIndex]

        # Create a dictionary with column headers as keys and dictionaries as values
        for columni in validColumns:
            # Populate each nested dictionary using the index keys and column data
            fullColumnDict = validColumns[columni].to_dict()
            columnDict = {key: value for key, value in fullColumnDict.items() if key in defaultStructData}
            fileDict[columni] = columnDict

        return fileDict
    def process_uniprotData(uniprotData):
        """process termini parameters to make a selection"""

        smallerPairs    = [
            ['N interface proximity score', 'C interface proximity score'],
            ['N volume', 'C volume'],
            ['N residue normalized burial','C residue normalized burial']
        ]
        biggerPairs     = [
            ['N consfScore', 'C consfScore'],
            ['N activeSiteScore', 'C activeSiteScore'],
            ['N conservedSiteScore', 'C conservedSiteScore']
        ]

        def bigger_smaller(a,b,relation):
            try:
                if relation == '>':
                    return a > b
                elif relation == '<':
                    return a < b
            except:
                return 'error'
        def compare_pairs(biggerOrSmallerPairs, NCRelation, scores, judgements, entry):
            """compare a parameter from the N and C termini,
            RETURN a judgement for each parameter and a score for the pair"""

            boolScoreDict = {True: -1, False: 1}
            for pair in biggerOrSmallerPairs:

                #extract N and C data for specific parameter
                Ndata, Cdata    = entry[pair[0]], entry[pair[1]]

                #process and judge each parameter
                NScore = AScore(pair[0], Ndata)
                CScore = AScore(pair[1], Cdata)

                #calculate score for pair
                pairJudgements      = [NScore.judgement, CScore.judgement]
                adjNData, adjCData  = [NScore.adjValue, CScore.adjValue]
                if any(pairJudgements):                                                                                 #if either value is significant
                    if adjCData == adjNData:                                                                                #if they're the same
                        score = 0                                                                                               #no preference
                    else:                                                                                                   #if they're different
                        difference      = abs(abs(adjNData) - abs(adjCData))                                            #calculate difference
                        minDifference   = 0.1 * abs(sum([adjNData, adjCData])) / 2                                            #min difference for significance
                        if difference > minDifference:                                                     #if difference > 10% of average
                            score = boolScoreDict[bigger_smaller(adjNData, adjCData, NCRelation)]                                     #prefer more favorable score
                        else:                                                                                                   #if difference < 10% of average
                            score = 0                                                                                               #no preference
                else:                                                                                                   #if neither value is significant
                    score = 0                                                                                               #no preference

                pairKey                            = pair[0][2:]
                scores[f'comparison {pairKey}']    = score
                judgements['N'][f'N J {pairKey}']  = NScore.judgement
                judgements['C'][f'C J {pairKey}']  = CScore.judgement

            return scores, judgements

        for pdbID, pdbData in uniprotData.items():
            if pdbData:

                #initialize scoring dicts
                scores      = {}
                judgements  = {'N':{},'C':{}}

                #extract cavity Data
                NCavities, CCavities = pdbData['N cavities'], pdbData['C cavities']
                NCavities = eval(NCavities) if type(NCavities) == str else None
                CCavities = eval(CCavities) if type(CCavities) == str else None

                #process termini cavities
                def process_cavities(cavities, termKey):
                    """extract cavity data and recalculate parameters using probability"""

                    cavityScores = {'consfScore': 0, 'activeSiteScore': 0, 'conservedSiteScore': 0, 'volume': 0}
                    if cavities:
                        for cavity in cavities:
                            probability = cavity['probability']
                            del cavity['probability']

                            #recalc scores via probability
                            cavityAdjusted              = {key:value*probability for key, value in cavity.items() if key != 'volume'}
                            cavityAdjusted['volume']    = cavity['volume']
                            for key, value in cavityAdjusted.items():
                                cavityScores[key] += value
                        cavityScores = {f'{termKey} {key}':value for key,value  in cavityScores.items()}
                    else:
                        cavityScores = {f'{termKey} {key}':None for key,_       in cavityScores.items()}

                    return cavityScores
                NCavData = process_cavities(NCavities, 'N')
                CCavData = process_cavities(CCavities, 'C')
                pdbData.update(NCavData)
                pdbData.update(CCavData)

                #process scores
                scores, judgements  = compare_pairs(smallerPairs, '<', scores, judgements, pdbData)
                scores, judgements  = compare_pairs(biggerPairs, '>', scores, judgements, pdbData)

                #compare score pairs
                judgements['N']     = {key: -judgements['N'][key] for key, value in judgements['N'].items()}
                numOfBadNScores     = sum(judgements['N'].values())
                numOfBadCScores     = sum(judgements['C'].values())
                pairScoreVals       = list(scores.values())
                pairScoreSum        = sum(pairScoreVals)

                prediction = 0 if pairScoreSum == 0 else -1 if pairScoreSum < 0 else 1

                #update output
                uniprotData[pdbID]['JUDGEMENTS N']    = ''
                uniprotData[pdbID]            .update(judgements['N'])
                uniprotData[pdbID]['JUDGEMENTS C']    = ''
                uniprotData[pdbID]            .update(judgements['C'])
                uniprotData[pdbID]['COMPARISONS']    = ''
                uniprotData[pdbID]            .update(scores)
                uniprotData[pdbID]['final judgement']    = ''


                uniprotData[pdbID]['N risk score']   = abs(numOfBadNScores)
                uniprotData[pdbID]['C risk score']   = numOfBadCScores
                uniprotData[pdbID]['Terminator prefers'] = ['equal', 'C','N'][prediction]
                uniprotData[pdbID]['id again'] = pdbID
        return uniprotData

    processedDatas = {}
    for uniprotID, filePath in filePaths.items():

        #extract data from file
        uniprotData = extract_data_from_excel(filePath, defaultStructData)

        #process termini parameters into a decision
        processedUniprotData = process_uniprotData(uniprotData)
        processedDatas[uniprotID] = processedUniprotData

    #create an excel file with the results
    def save_file(processedDatas, filePaths):
        # Create a list to store the data rows
        allRows = []

        # Iterate over the items in the dictionary
        for uniprotID, uniprotData in processedDatas.items():

            uniprotRows = []
            if uniprotData is None or not uniprotData:
                # If the value is None or an empty dictionary, create a row with only the UniProtID
                row = {'UniprotID': uniprotID}
                allRows.append(row)
                uniprotRows.append(row)
            else:
                for pdbID, pdbData in uniprotData.items():
                    row = {'UniprotID': uniprotID,'pdbID':pdbID, **pdbData}

                    allRows.append(row)
                    uniprotRows.append(row)

            # Create a DataFrame from the uniprot Data
            df = pd.DataFrame(uniprotRows)
            df = df.transpose()

            #write uniprot Data to a new sheet in original file
            oriFilePath = filePaths[uniprotID]
            with pd.ExcelWriter(oriFilePath, engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:
                df.to_excel(writer, sheet_name='analysis')
            print(f"analysis page added to:\n {oriFilePath}")

        #with multiple uniprotIDs, write to compiled_results.xlsx
        if len(processedDatas.keys()) > 1:
            df = pd.DataFrame(allRows)
            output_file = 'compiled_results.xlsx'
            print(f"Excel file '{output_file}' in script directory created successfully.")
            df.to_excel(output_file, index=False)
    save_file(processedDatas, filePaths)
if __name__ == "__main__":
    defaultStructData = {'indexes': [],
                         'structureName': '',
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
    uniprotIDs = ['A0A089Q240', 'P08684', 'P07320', 'P00722', 'Q6QGP7', 'Q8VNN2', 'Q9HV14', 'P20815', 'P40013', 'Q9BHM6', 'Q93VR3', 'P52705', 'P00344', 'P08246', 'P59807', 'A0A0D8BQX7', 'P12493', 'P03973']

    uniprotIDs = [id.upper() for id in uniprotIDs]
    export_multiple_results(uniprotIDs, defaultStructData)
