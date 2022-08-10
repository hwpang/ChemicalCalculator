# ChemicalCalculator

A chemical calculator designed to determine if a product meets the definition of an Organic Peroxide under CFR 173.128(a)(4)


- Which of the provided formulated products meet the definition Organic Peroxide under CFR 173.128(a)(4)?
    - Product 2 and 3 meet the definition of organic peroxide
- Which of the provided formulated products do not meet the definition?
    - Product 1 and 4 do not meet the definition of organic peroxide
- Which subsection under CFR 173.128(a)(4) was used to determine whether the product meets or does not meet the definition of an Organic Peroxide?
    - Product 2 and 3 are determined as organic peroxide as they contain peroxide but doesn't meet conditions described under CFR 173.128(a)(4). Product 1 and 4 are not determined as organic peroxide as they meet conditions described under CFR 173.128(a)(4). See table below.

        | Product ID      | Is organic peroxide? | Reason                                                |
        | --------------- | -------------------- | ----------------------------------------------------- |
        | Product 1       | No                   | Meets CFR 173.128(a)(4)(i)                            |
        | Product 2       | Yes                  | Contains peroxide and doesn't meet CFR 173.128(a)(4)  |
        | Product 3       | Yes                  | Contains peroxide and doesn't meet CFR 173.128(a)(4)  |
        | Product 4       | No                   | Meets CFR 173.128(a)(4)(i)                            |

- Are any of the products eligible for an exemption?
    - Based on the compositions and ingredients, product 1 seems to be an acne product, while product 2 to 4 seem to be cleaning products. Product 1 may be eligible for an exemption if they fall into a different category of products, e.g. cosmetic products.

## Incremental commits


- Designed a class `ChemicalSpecies` to hold and process information associated with a chemical species
    - Takes CAS number and chemical name as inputs
    - Uses CAS number to search for smiles string of the chemical species from cactus.nci.nih.gov
        - Throws warning if smiles string cannot be found
    - Uses smiles string to construct a RDKit molecule
    - Uses RDKit's substructure match function to obtain
        - Number of peroxide group structure (-O-O-) in the chemical species
        - Number of hydrogen peroxide structure (H-O-O-H) in the chemical species
    - Uses RDKit's Descriptors to obtain molecular weight of the chemical species


- Designed a class `ChemicalComposition` to hold and process information associated with a chemical composition
    - Takes a `ChemicalSpecies`, min conc, and max conc as inputs
    - Uses max conc and information contained in chemical species to calculate max conc of hydrogen peroxide


- Designed a class `FormulatedProduct` to hold and process information associated with a formulated product
    - Takes a list of `ChemicalComposition` and product ID as inputs
    - Check if the formulated product meets the definition of organic peroxide based on CFR 173.128(a)(4)
        - Check if the formulated product contains chemical species that has non-zero number of peroxide group structure
            - If the formulated product doesn't contain peroxides, then the formulated product doesn't meet the definition of organic peroxide
        - CFR 173.128(a)(4): If the formulated product contains peroxides, check the conc of hydrogen peroxide and the Oa
            - If the hydrogen peroxide conc and the Oa falls into either CFR 173.128(a)(4)(i) or CFR 173.128(a)(4)(ii), then the formulated product doesn't meet the definition of organic peroxide
            - Otherwise, the formulated product meets the definition of organic peroxide

- Designed a function to help easy usage of the chemical calculator
    - Takes the path to the product list csv file
    - Returns the modified csv file containing chemical calculator results

- Other miscellaneous tasks
    - Listed out required packages in `environment.yml` file to create conda env
    - Created `setup.py` file to use pip install
    - Created `.gitignore` file to ignore files generated by pip install
    - Created ijupyter notebook for example
    - Created `main.py` for example

## Instructions for installing the package

- Clone the repo by `git clone https://github.com/hwpang/ChemicalCalculator.git`

- There are three options to use this package:
    - Use docker
        - Install Docker (https://www.docker.com/get-started/)
        - Use `docker build -t chem_calc .` in the package folder `ChemicalCalculator` to build image chem_calc
        - Use `docker run -p 8888:8888 chem_calc` to create a container
        - On the terminal, a URL similar to `http://127.0.0.1:8888/lab?token=<token>` can be found
        - Copy the URL and paste in a browser
    - Use pip install
        - Run `pip install jupyter numpy pandas rdkit` to install necessary packages
        - Install the package by `pip install -e .` in the package folder `ChemicalCalculator`
    - Use conda environment
        - Install Anaconda 3 (https://www.anaconda.com/products/individual#Downloads)
        - Use `conda env create -f environment.yml` to create a conda environment `chem_calc`
        - Activate the conda environment by `conda activate chem_calc`
        - Install the package by `pip install -e .` in the package folder `ChemicalCalculator`

## Instructions for running the test suite

- Example code to run the test suite can be found in the `Running_test_suite.ipynb` jupyter notebook in `ipython` folder
    - `cd` to the `ipython` folder
    - `jupyter notebook` to open up the jupyter notebook