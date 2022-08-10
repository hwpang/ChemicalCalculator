#!/usr/bin/env python3
#-*- coding: utf-8 -*-

import urllib
import logging
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

# logging.getLogger().setLevel(logging.INFO)
logging.getLogger().setLevel(logging.WARNING)
PEROXIDE_RDKITMOL = Chem.MolFromSmiles('OO')
HYDROGEN_PEROXIDE_RDKITMOL = Chem.rdmolops.AddHs(Chem.MolFromSmiles('OO'))
HYDROGEN_PEROXIDE_MW = Descriptors.MolWt(HYDROGEN_PEROXIDE_RDKITMOL)

class ChemicalSpecies(object):
    """
    A class to hold and process information associated with a chemical species
    """

    def __init__(self,
                 CAS: str,
                 chemical_name: str,
                 smiles: str = None):
        """
        Generate a ChemicalSpecies instance and process necessary chemical information
        - Takes CAS number and chemical name as input
        - Uses CAS number to search for smiles string of the chemical species from cactus.nci.nih.gov
            - Throws warning if smiles string cannot be found
        - Uses smiles string to construct a RDKit molecule
        - Uses RDKit's substructure match function to obtain
            - Number of peroxide group structure (-O-O-) in the chemical species
            - Number of hydrogen peroxide structure (H-O-O-H) in the chemical species
        - Uses RDKit's Descriptors to obtain molecular weight of the chemical species
        
        Args:
            CAS (str): the CAS number for the chemical species
            chemical_name (str): the chemical name for the chemical species
            smiles (str, optional): the smiles for the chemical species
        """
        
        self.CAS = CAS
        self.chemical_name = chemical_name
        
        self.smiles = None
        self.rdkitmol = None
        self.is_peroxide = False
        self.num_peroxide_groups = 0
        self.num_hydrogen_peroxide_groups = 0
        self.molecular_weight = None
        
        if smiles is not None:
            self.smiles = smiles
        else:
            self.smiles = search_smiles(self.CAS)

        if self.smiles is not None:
            self.rdkitmol = Chem.rdmolops.AddHs(Chem.MolFromSmiles(self.smiles))
            self.num_peroxide_groups = len(self.rdkitmol.GetSubstructMatches(PEROXIDE_RDKITMOL, uniquify=True))
            self.num_hydrogen_peroxide_groups = len(self.rdkitmol.GetSubstructMatches(HYDROGEN_PEROXIDE_RDKITMOL, uniquify=True))
            self.molecular_weight = Descriptors.MolWt(self.rdkitmol)
            self.is_peroxide = self.num_peroxide_groups != 0
            
    def __repr__(self):
        """
        Generate string representation that can be used to replicate the object
        """
        return f"ChemicalSpecies(CAS='{self.CAS}', chemical_name='{self.chemical_name}', smiles='{self.smiles}')"

def search_smiles(identifier):
    """
    Search for smiles based on CAS number or chemical name
    """
    smiles = None

    logging.info("smiles not provided... searching on cactus.nci.nih.gov using CAS...")
    url = f"https://cactus.nci.nih.gov/chemical/structure/{identifier}/smiles"
    try:
        out = urllib.request.urlopen(url, timeout=5)
    except:
        logging.error(f"Unable to obtain smiles for {identifier}...")
        logging.error(f"Assuming {identifier} is not a peroxide...")
        logging.error(f"Result of chemical calculator may be wrong for this case...")
    else:
        smiles = out.read().decode('utf-8')
        logging.info(f"Found smiles {smiles} for {identifier}")

    return smiles
                
class ChemicalComposition(object):
    """
    A class to hold a chemical species and its composition
    """
    
    def __init__(self,
                 chemical_species: ChemicalSpecies,
                 min_conc: float,
                 max_conc: float):
        """
        Generate a ChemicalComposition instance and compute necessary concentration information
        
        Args:
            chemical_species (ChemicalSpecies): ChemicalSpecies
            min_conc (float): minimum concentration (mass percent) of the chemical species
            max_conc(float): maximum concentration (mass percent) of the chemical species
        """
        self.chemical_species = chemical_species
        self.min_conc = min_conc
        self.max_conc = max_conc
        
        self.max_hydrogen_peroxide_conc = self.calculate_max_hydrogen_peroxide_conc()
        
    def calculate_max_hydrogen_peroxide_conc(self) -> float:
        """
        Calculate the maximum hydrogen peroxide concentration (mass percent) contributed by the chemical species using
        hydrogen peroxide mass fraction = (number of hydrogen peroxide groups * hydrogen peroxide molecular mass)/chemical species molecular mass
        max hydrogen peroxide concentration = hydrogen peroxide mass fraction * max concentration
        
        Returns:
            float: the maximum hydrogen peroxide concentration
        """
        if self.chemical_species.molecular_weight:
            hydrogen_peroxide_mass_fraction = self.chemical_species.num_hydrogen_peroxide_groups*HYDROGEN_PEROXIDE_MW/self.chemical_species.molecular_weight
            return hydrogen_peroxide_mass_fraction*self.max_conc
        else:
            return 0.0
        
    def __repr__(self):
        """
        Generate string representation that can be used to replicate the object
        """
        return f"ChemicalComposition(chemical_species=ChemicalSpecies('{self.chemical_species.CAS}', '{self.chemical_species.chemical_name}'), min_conc={self.min_conc}, max_conc={self.max_conc}, smiles='{self.chemical_species.smiles}')"

class FormulatedProduct(object):
    """
    A class to hold the compositions of a formulated product
    """

    def __init__(self,
                 product_id : str,
                 chemical_compositions: list[ChemicalComposition]):
        """
        Generate a FormulatedProduct instance to perform chemical calculation
        
        Args:
            product_id (str): product id
            chemical_compositions (list[ChemicalComposition]): list of ChemicalComposition
        """
        self.product_id = product_id
        self.chemical_compositions = chemical_compositions

    def contains_peroxide(self) -> bool:
        """
        Check if the formulated product contains peroxide in its compositions by checking if any of 
        the chemical composition is peroxide and has maximum concentration larger than 0.0
        """
        return any(comp.chemical_species.is_peroxide and comp.max_conc > 0.0 for comp in self.chemical_compositions)
        
    def calculate_Oa(self) -> float:
        """
        Calculate Oa based on formula under CFR 173.128(a)(4)
        """
        return 16*sum(comp.chemical_species.num_peroxide_groups*comp.max_conc/comp.chemical_species.molecular_weight for comp in self.chemical_compositions if comp.chemical_species.is_peroxide)
    
    def is_organic_peroxide(self) -> tuple[bool, str]:
        """
        Determine whether a formulated product meets the definition of organic peroxide based on CFR 173.128(a)(4)
        - Check if the formulated product contains chemical species that has non-zero number of peroxide group structure
            - If the formulated product doesn't contain peroxides, then the formulated product doesn't meet the definition of organic peroxide
        - CFR 173.128(a)(4): If the formulated product contains peroxides, check the conc of hydrogen peroxide and the Oa
            - If the hydrogen peroxide conc and the Oa falls into either CFR 173.128(a)(4)(i) or CFR 173.128(a)(4)(ii), then the formulated product doesn't meet the definition of organic peroxide
            - Otherwise, the formulated product meets the definition of organic peroxide
        
        Returns:
            tuple[bool, str]: a bool to indicate whether the formulated product meets the definition of organic peroxide
                              and the corresponding subsection under CFR 173.128(a)(4)
        """
        if not self.contains_peroxide():
            return False, "Doesn't contain peroxide"
        
        hydrogen_peroxide_conc = sum(comp.max_hydrogen_peroxide_conc for comp in self.chemical_compositions)
        Oa = self.calculate_Oa()
        
        if hydrogen_peroxide_conc <= 1.0 and Oa <= 1.0:
            return False, "Meets CFR 173.128(a)(4)(i)"
        if hydrogen_peroxide_conc > 1.0 and hydrogen_peroxide_conc <= 7.0 and Oa <= 0.5:
            return False, "Meets CFR 173.128(a)(4)(ii)"
        
        return True, "Contains peroxide and doesn't meet CFR 173.128(a)(4)"

def determine_organic_peroxide(file_path: str) -> pd.DataFrame:
    """
    A helper function to use chemical calculator for batch process.
    
    Args:
        file_path (str): file path to a csv file containing a list of formulated products, their product ID, and
                         CAS number, chemical name, min/max conc in mass percent of each chemical ingredients.
                         See data/Product_List.csv for suggested format.
    Returns:
        DataFrame: modified dataframe containing the chemical calculator results
    """
    df = pd.read_csv(file_path)

    all_chemical_compositions = [ChemicalComposition(ChemicalSpecies(df["CAS"][ind], df["Chemical Name"][ind]), df["min (%)"][ind], df["max (%)"][ind]) for ind in df.index]
    df["smiles"] = [comp.chemical_species.smiles for comp in all_chemical_compositions]
    df["Molecular weight"] = [comp.chemical_species.molecular_weight for comp in all_chemical_compositions]
    df["Number of peroxide groups"] = [int(comp.chemical_species.num_peroxide_groups) for comp in all_chemical_compositions]
    df["Number of hydrogen peroxide groups"] = [int(comp.chemical_species.num_hydrogen_peroxide_groups) for comp in all_chemical_compositions]
    df["Max hydrogen peroxide conc (%)"] = [comp.max_hydrogen_peroxide_conc for comp in all_chemical_compositions]
    
    product_ids = list(df["Product ID"].value_counts().keys())
    for product_id in product_ids:
        inds = df.index[df["Product ID"]==product_id]
        product = FormulatedProduct(product_id, [all_chemical_compositions[ind] for ind in inds])
        is_organic_peroxide, reason = product.is_organic_peroxide()
        
        df.loc[inds, "Max product Oa (%)"] = product.calculate_Oa()
        df.loc[inds, "Is product organic peroxide?"] = is_organic_peroxide
        df.loc[inds, "Reason"] = reason

    return df

