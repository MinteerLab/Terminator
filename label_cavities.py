import json, pymol, bisect
import numpy as np

json_file = r'C:\Users\16302\University of Utah Dropbox\Minteer Lab\Rokas\main\scripts\231225 terminator\bothterms\D0VWQ0 D0VWQ0_ALCFA\PDBs\SWISSMODEL_2YZ7_cavPAtomLocs.json'

#dirty fix, these are small variations of direct copies from parse_PDB.py
#difficulties with importing from parse_PDB when this script is run in pymol
#difficulties importing libraries because pyMol runs on an older version of Python
class AnAtom:
    def __init__(self, element='', name='', location=None, residue = None, b_factor = None, consfScore = None, colorScore = None, number = None):
        self.element        = element
        self.name           = name
        self.location       = location
        self.residue        = residue
        self.b_factor       = b_factor
        self.consfScore     = consfScore
        self.colorScore     = colorScore
        self.number         = number
        self.burial             = None          #determined in calculate_burials
        self.normBurial         = None          #determined in calculate_burials
        self.isCavity           = None          #determined in crawl_cavities
        self.normBurials        = {}
        self.cavity             = None          #determined in ACavity initialization in crawl_cavities
    def __repr__(self):
        return f"atom {self.name} in {self.residue}"
    def transform(self, matrix= None):
        """ Apply a 3x4 transformation matrix to an atom's location. return new atom. """
        x, y, z = self.location
        new_x = matrix[0][0] * x + matrix[0][1] * y + matrix[0][2] * z + matrix[0][3]
        new_y = matrix[1][0] * x + matrix[1][1] * y + matrix[1][2] * z + matrix[1][3]
        new_z = matrix[2][0] * x + matrix[2][1] * y + matrix[2][2] * z + matrix[2][3]
        transformed_location = np.array([new_x, new_y, new_z])
        transformed_atom = AnAtom(element=self.element, name=self.name,
                                  location=transformed_location, b_factor = self.b_factor,
                                  consfScore = self.consfScore, colorScore = self.colorScore, number = self.number)

        return transformed_atom

def make_positionBins(objects):
    """sort objects into square bins by their position in each dimension, to be called before find_close_objects()"""
    positionBins = []
    for dimension in range(3):
        positionBin = dict()
        for object in objects:
            position = object.location[dimension]
            positionBin[position] = positionBin[position] if position in positionBin else []
            positionBin[position].append(object)
        positionBin = {key: positionBin[key] for key in sorted(positionBin)}
        positionBins.append(positionBin)
    return positionBins
def find_close_objects(reference_object, positionBins, cutoff):
    """ Find atoms within a specified distance of a reference atom. """
    def find_closest_value(targetVal, queryVals):
        """use bisect to find nearest value to targetVal in queryVals"""
        index = bisect.bisect_left(queryVals, targetVal)

        if index == 0:
            return queryVals[0]
        elif index == len(queryVals):
            return queryVals[-1]
        else:
            leftDiff = abs(targetVal - queryVals[index - 1])
            rightDiff = abs(targetVal - queryVals[index])
            return queryVals[index - 1] if leftDiff < rightDiff else queryVals[index]

    atomBins = []
    for dimension, positionBin in enumerate(positionBins):
        positionKeys = list(positionBin.keys())
        minVal       = reference_object.location[dimension] - cutoff
        maxVal       = reference_object.location[dimension] + cutoff
        minReal      = find_closest_value(minVal, positionKeys)
        maxReal      = find_closest_value(maxVal, positionKeys)
        rangeBins    = [positionBin[key] for key in positionKeys if minReal <= key <= maxReal]
        rangeAtoms   = set(atom for atomBin in rangeBins for atom in atomBin)
        atomBins.append(rangeAtoms)
    allDimensionAtoms = set.intersection(*atomBins)

    def get_distance(obj1, obj2):
        """
        Calculate the Euclidean get_distance between two objects.

        :param obj1: First atom with 'location' attribute as a list [x, y, z].
        :param obj2: Second atom with 'location' attribute as a list [x, y, z].
        :return: Distance between obj1 and obj2.
        """

        location1 = np.array(obj1) if isinstance(obj1, list) else obj1.location
        location2 = np.array(obj2) if isinstance(obj2, list) else obj2.location
        return np.linalg.norm(location1 - location2)

    finalAtoms = set()
    for atom in allDimensionAtoms:
        if get_distance(reference_object, atom) <= cutoff:
            finalAtoms.add(atom)
    return finalAtoms

#BODY

def generate_pseudoatoms(json_file):
   # Read the coordinates from the JSON file
    with open(json_file, 'r') as file:
        cavPAtomLocs = json.load(file)


    cavPAtoms = [AnAtom(location=np.array(location), name=f'pAtom {i}') for i, location in enumerate(cavPAtomLocs)]
    cavPAtomLocTree = make_positionBins(cavPAtoms)


    def bin_pAtoms_by_location(cavPAtoms, cavPAtomLocTree, cutoff):
        """collects pAtoms within cutoff of each other into individual cavities"""

        assignedPAtoms = set()
        pAtomCavityBins = []
        for pAtom in cavPAtoms:
            if pAtom not in assignedPAtoms:
                binned = False
                siblingPAtoms = find_close_objects(pAtom, cavPAtomLocTree, cutoff)
                if siblingPAtoms:
                    for pAtom in siblingPAtoms:
                        for bin in pAtomCavityBins:
                            if pAtom in bin:
                                bin.update(siblingPAtoms)
                                binned = True
                                break
                    if not binned:
                        pAtomCavityBins.append(siblingPAtoms)
                    assignedPAtoms.update(siblingPAtoms)

        def join_shared_bins(pAtomCavityBins):
            """Joins bins that share pAtoms"""
            joined_bins = []
            for bin1 in pAtomCavityBins:
                new_bin = bin1.copy()
                for bin2 in joined_bins:
                    if bin1.intersection(bin2):
                        new_bin.update(bin2)
                        joined_bins.remove(bin2)
                joined_bins.append(new_bin)
            return joined_bins

        pAtomCavityBins = join_shared_bins(pAtomCavityBins)

        return pAtomCavityBins
    allCavPAtomBins = bin_pAtoms_by_location(cavPAtoms, cavPAtomLocTree, cutoff=2.5)
    cavPAtomBins    = [bin for bin in allCavPAtomBins if len(bin) > 90]
    colors       = [
    "purple", "teal", "violet", "deepblue", "olive", "hotpink",
    "lightorange", "deepolive", "deepteal", "ruby", "forest",
    "deep purple", "lightred", "darksalmon", "lightpink", "deepgreen",
    "wheat", "grey", "deepgrey", "lightgrey", "darkgrey"
]

    cavPAtomBins = {color: bin for color, bin in zip(colors, cavPAtomBins)}

    for color, cavPAtomBin in cavPAtomBins.items():
        cavName = f'cav_{color}'

        # Iterate over the list of coordinates
        for i, cavPAtom in enumerate(cavPAtomBin):
            pos = list(cavPAtom.location)
            name = f"{cavName}_{i+1}"
            pymol.cmd.pseudoatom(name, pos =pos, vdw = 0.5, chain = cavName, color = color)

# Example usage
generate_pseudoatoms(json_file)
