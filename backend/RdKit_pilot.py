from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors, Crippen, Lipinski, rdMolDescriptors
import pubchempy as pcp

cpd = input("what molecule do you want to visualize\n"
            "pleace check on the spelling ")
compound = pcp.get_compounds(cpd, 'name')[0]

smiles = compound.connectivity_smiles
mol = Chem.MolFromSmiles(smiles)
def Analyze_chem(smiles):

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ("enter a valid molecule young king")
    props = {
        "MolecularWeight": round(Descriptors.MolWt(mol), 2),
        "LogP": round(Crippen.MolLogP(mol), 2),
        "HBD": Lipinski.NumHDonors(mol),
        "HBA": Lipinski.NumHAcceptors(mol),
        "TPSA": round(rdMolDescriptors.CalcTPSA(mol), 2),
        "RotatableBonds": Lipinski.NumRotatableBonds(mol),
        "NumAtoms": mol.GetNumAtoms(),
    }
    print("molecule properties:", props)

    return props

Draw.MolToImage(mol).show()
Analyze_chem(smiles)