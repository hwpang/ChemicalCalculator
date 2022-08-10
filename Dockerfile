FROM jupyter/minimal-notebook
USER root
RUN pip install --no-cache-dir numpy pandas rdkit
ADD . .
RUN pip install --no-cache-dir .
RUN echo "Make sure ChemicalCalculator is installed:"
RUN python -c "from chemical_calculator import *"