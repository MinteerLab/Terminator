How to run Terminator:
1.	Install python version 3.9 and Google Chrome on your system
  a.	Terminator was written with python 3.9 on windows 10 and may not function with later or earlier versions of python and will not function on different operating systems. Several of Terminator’s web-scraping functions employ Selenium with chrome dependencies, so google chrome is required for operation.
  b.	Python installation: python.org/downloads/
  c.	Chrome installation: google.com/chrome/
2.	Download the Terminator zip package from our GitHub page:
  a.	github.com/MinteerLab
3.	Extract the downloaded Terminator directory in a convenient location. Terminator outputs are stored as subdirectories in the Terminator folder, titled according to UniProt ID.
4.	Run the installer.py script to install dependency libraries. 
  a.	For totally new python users, scripts can be run by right clicking the script and selecting “edit with Idle”, then with the file open, pressing “F5” or navigating to run> “run module” in the header.
  b.	If the installer script fails for unforeseen reasons, you may need to manually install the external dependencies with pip commands in the command prompt.
  c.	Manual pip installation of dependencies: pip.pypa.io/en/stable/cli/pip_install/
  d.	Dependencies: numpy, pandas, matplotlib, requests, selenium, openpyxl, biopython, scikit-learn
5.	Obtain a UniProt ID or a protein sequence (as long as possible) for your target protein
6.	Open Terminator (“edit with IDLE”), fill in the select input variables, and run the script. Upon completion, ensure that in the output termini_scores file, the ConSurf column has good coverage of the sequence in the UniProt column. If not, input ConsurfPDB and/or consurfChain to specify a structure which has better sequence coverage than the current ConSurf column, and set restConsurf = True. If Consurf fails to retrieve a high-quality sequence, set consurfJob = True, and await an email from Consurf for your job to be completed before running again. 



Input variables:
Terminology:
Users new to python scripting are advised to familiarize themselves with this short section to format Terminator inputs appropriately. 
1.	String: 	a text object, surrounded in single or double quotes, e.g. “cat” or ‘b4’.
2.	List: 	A collection of objects such as strings, surrounded in square brackets, e.g. [‘cat’, ‘dog’]. 
3.	Dictionary:	A dictionary is a collection of key : value pairs (like words : definitions in a real dictionary). Surrounded by curly braces {}, e.g. {‘dog’ : ’a pet’} or {‘cat’: [‘a pet’, ‘a mammal’]}.  
4.	Bool: 	a Boolean value, True or False, first letter capitalized, no surrounding markers.
All objects in python must be properly defined, e.g. rush = “”, an empty string if you do not wish to input a value for rush.

Required variables:
1.	UniProtIDs: A string or list of strings of UniProtIDs. A list results in separate Terminator runs for multiple targets.
2.	Sequence:	A string of a protein sequence if UniProtIDs is not provided.
3.	Email: 	A string. Providing an email address is required for computationally accessing PubMed and BLAST. Although ConSurfDB does not require an email, Consurf Job requests do require one. Terminator does not send the user’s email to any other locations.

Optional variables:
1.	Rushes: 		Specific subset of selected PDBIDs to analyze instead of analyzing every assigned PDBID. This serves only to reduce computation time when cavitySearch = True, see cavitySearch for advice on structure selection.
If a single UniProt ID is provided, rushes can be a single string or list of strings of PDBs. 
If multiple UniProtIDs are provided, rushes must be a dictionary of {string UniProtID : string/list of strings of PDBIDs to focus. Not providing a key for a certain UniProtID simply results in all structures being processed, unless swmOnly = True.
2.	forcePDB:     	a string of a PDB file name. For running a PDB file in the project’s PDB subdirectory even if it isn’t recognized from UniProt entries retrieved for the project, e.g. an AlphaFold structure manually placed there. Only functions with single UniProtID requests. 
3.	consurfPDB:   	a string of a PDB id to specify which PDB should be used to request ConSurf for conservation information. e.g. if the Terminator output reveals many structures lacking coordinate assignment for a significant number of residues (such as structures intended for studying single domains in protein-protein interactions), the user can specify which PDB should be used to retrieve conservation data instead. See consurfJob below. 
4.	consurfChain:  	like consurfPDB, but for specifying the chain of the structure. E.g. “C”
5.	resetConsurf:  	a bool. For deleting old consurfDB data to request new conservation data of a different structure or chain. Because all of Terminator’s database queries are only executed if the expected file retrieved is not present in the project directory, setting a new value for consurfPDB or consurfChain requires setting resetConsurf = True.
6.	consurfJob:    	a bool. Terminator first tries to retrieve conservation information from consurfDB using the target sequence, then the retrieved PDB ids. Because conservation data retrieved from a crystal structure can be confounded by variable residue coordinate assignment (visible in the termini_scores output as gaps in the CONSURF… column), and consurfDB may fail to retrieve a sequence-specific output, consurfJob = True forces Terminator to make a full request for a job from Consurf. This can often take over 20 minutes for ConSurf to complete, compared to instant consurfDB retrievals. If a job is requested, the user will receive an email from Consurf when the job is complete, indicating that they can rerun the Terminator analysis. If consurfDB requests fail by both structure and sequence, a job is requested automatically.
7.	resetSwm:      	a bool. For deleting old SWISSMODEL data to retry if a previous attempt failed, which can occur rarely in spite of a structure being available, for unknown reasons (suspected to occur on the side of the SWISSMODEL server). If Terminator repeatedly fails retrieving SWISSMODELS hours or days apart, the user is advised to try making a manual SWISSMODEL database query on swissmodel.expasy.org, and trying again.
8.	cavitySearch:	a bool. Because cavity calculations are the most time-consuming step in Terminator calculations, the user is advised to first run a target with cavitySearch = False. Then after selecting the best structures to target (by employing the rushes input) based on the number of terminal residue lacking coordinates, mismatch, gaps, and substrate presence (in the structure name), the user should re-run the analysis with cavitySearch = True.
9.	debug:        	a bool. Causes selenium to visibly display browser interactions in certain UniProt and all Consurf interactions. This can help verify if one of these interactions fails because the target database is malfunctioning or has been modified in a way which requires updating the relevant web-scraping functions.
10.	swmOnly: 		a bool. Adds the SWISSMODEL structure to rushes. Although the SWISSMODEL is likely the highest quality structure available for a protein and often corrects issues with the original structure on which the model is based, users should preferentially select structures as discussed in cavitySearch above. If the only structure available is the SWISSMODEL, it is critical for the user to judge the output GMQE and QSQE of the structure to determine its validity. Both metrics should be above 0.7, although occasionally SWISSMODEL may judge that a crystal structure of the actual target protein has an improbable quaternary structure (QSQE < 0.7), as these metrics evaluate crystal structure quality and not merely sequence-structure coherence probability.
11.	targetChainName: a string. like Rush, but for specifying which chain of the protein should be targeted. This is for structures that have a substrate present, but only in one chain. Determining this requires first parsing the structure names in the termini_scores output, followed by a manual analysis of the PDB structure with a protein visualization software such as PyMOL1. Terminator automatically selects the chain with the lowest mismatch and gap percentage vs. the protein sequence retrieved from UniProt.
12.	resetCavities: 	a bool. delete old cavity data for the target to recalculate their positions. Mainly for advanced users to modify cavity searching properties.
13.	getColorCmds: 	a bool. Print a series of commands to execute in PyMOL to color the relevant crystal structure according to ConSurf-derived conservation.
