"""
Calculates salt bridges or electrostatic between two chains.

Arguments:
    saltBridges - use it to get stats for pair of 2 charges in distance of less than 4A.
    periphery - use it to look for periphery charges only

See also:
DRIEDING: A Generic Force field for molecular simulations
"""
from __future__ import print_function
import math

from . import DBConfig
from .kdtree import KDTree

VERBOSE = True

def getHbonds(pdb, pdbName):
    conn = DBConfig.get_connection()
    cursor = conn.cursor()
    cursor.execute("""select DonorChain,DonorResId,DonorSymbol,AccChain,AccResId,AccSymbol from
                        Ndrieding
                        inner join donors2
                        on DonorSymbol=donors2.Symbol
                        where PDB='%s'""" % pdbName)

    hbonds = []
    for dChain, dResId, dSymbol, aChain, aResId, aSymbol in cursor.fetchall():
        donorAtom, accAtom = None, None
        for a in pdb.atoms:
            if a.chain == dChain and a.resId == dResId and a.symbol == dSymbol:
                donorAtom = a
            elif a.chain == aChain and a.resId == aResId and a.symbol == aSymbol:
                accAtom = a
        hbonds.append((donorAtom, accAtom))
    return hbonds


def eInteraction(Qi, Qj, R):
    kcal_mol_constant = 322.0637
    return kcal_mol_constant * Qi * Qj / (R ** 2)


def assignCharge(atom, pH=7):
    # how do we assign?

    if atom.residue in ['ASP', 'GLU'] and atom.atomType == 'O' and atom.symbol != 'O':
        return -0.5  # -1
    if atom.symbol == 'OXT':
        return -0.5  # -1

    # arg deloclalized
    # ARG - NH1 NH2 (not NE and N)
    # LYS NZ
    # HIS ND1 NE2
    if atom.residue == 'LYS' and atom.symbol == 'NZ':
        return 1.0
    posRes = ['ARG']
    if 0.1 < pH <= 6:
        posRes.append('HIS')
    if atom.residue in posRes and atom.atomType == 'N' and atom.symbol not in ['N', 'NE']:
        return 0.5  # 1

    return 0


def calcElectrostatic(pdb, interface, exclude_hbonds=False, count_cutoff_distance=5):
    """
    Calculated possible electro interactions, excluding 1-2 and 1-3 interactions
    (already included in angle and bond interactions
     """
    CUTOFF_DISTANCE = 7  # we could have 6?
    if exclude_hbonds:
        hHbonds = getHbonds(pdb, pdb.name)
    else:
        hHbonds = []

    components = []
    oxtAtoms = [(a.chain, a.resId) for a in pdb.atoms if a.symbol == 'OXT']  # C ter
    for a in pdb.atoms:
        if a.symbol == 'O' and (a.chain, a.resId) in oxtAtoms:
            a.symbol = 'OXT'

    for part in pdb.interfaceParts:
        components.append([a for a in interface if a.chain in part and assignCharge(a) != 0])

    # fill kd tree with charged atoms from second component
    comp2 = [atom for atom in pdb.atoms if atom.chain in pdb.interfaceParts[1] and assignCharge(atom) != 0]

    ktree = KDTree.construct_from_data(comp2)
    electroStat = 0.0
    pp, pm, mm = 0, 0, 0
    for atom in components[0]:
        # interactions are not calculated between atoms bonded to each other (1,2 and 1,3 [hbonds])
        Qi = assignCharge(atom)
        if Qi == 0:
            continue

        nearAtoms = list(ktree.findByDistance(query_point=atom.coord, distance=CUTOFF_DISTANCE ** 2))
        contact = [con for con in nearAtoms if not (((con, atom) in hHbonds) or ((atom, con) in hHbonds))]
        for con in contact:
            Qj = assignCharge(con)
            if Qj == 0:
                continue
            R = math.sqrt(atom.distance(con))
            electroStat += eInteraction(Qi, Qj, R)

            if R < count_cutoff_distance:
                if Qi > 0 and Qj > 0:
                    pp += 1
                elif Qi < 0 and Qj < 0:
                    mm += 1
                else:
                    pm += 1

    return electroStat, pp, mm, pm


def is_hydrophobic(atom):
    """
    Checks whether an atom belongs to hydrophobic residue
    :param atom: atom
    :return: True if the atom belongs to hydrophobic residue, otherwise false
    """
    hydrophobic_res = ['VAL', 'ILE', 'LEU', 'PHE', 'TRP', 'CYS', 'ALA', 'PRO']
    return atom.residue in hydrophobic_res


def calcElectroHydrophobic(pdb, interface):
    """
    Calculated possible electro interactions, excluding 1-2 and 1-3 interactions
    (already included in angle and bond interactions
     """
    HYDROPHOBIC_CHARED_CUTOFF_DISTANCE = 4  # we could have 6?

    oxtAtoms = [(a.chain, a.resId) for a in pdb.atoms if a.symbol == 'OXT']  # C ter
    for a in pdb.atoms:
        if a.symbol == 'O' and (a.chain, a.resId) in oxtAtoms:
            a.symbol = 'OXT'

    hydrophobic_positive = set()
    hydrophobic_negative = set()
    hydroElectroOutput = pdb.getFile('.hydrophobicElectro.txt') if VERBOSE else None
    for part in pdb.interfaceParts:
        charged_atoms = [a for a in interface if a.chain in part and assignCharge(a) != 0]
        if not any(charged_atoms):
            continue
        electro_kdtree = KDTree.construct_from_data(charged_atoms)
        other_parts = ''.join([partb for partb in pdb.interfaceParts if partb != part])
        hydrophobic_partners = [a for a in interface if a.chain in other_parts and is_hydrophobic(a)]
        for atom in hydrophobic_partners:
            nearAtoms = list(
                electro_kdtree.findByDistance(query_point=atom.coord, distance=HYDROPHOBIC_CHARED_CUTOFF_DISTANCE ** 2))
            for con in nearAtoms:
                Qi = assignCharge(con)
                R = math.sqrt(atom.distance(con))
                if R < HYDROPHOBIC_CHARED_CUTOFF_DISTANCE:

                    if Qi > 0:
                        hydrophobic_positive.add(
                            (con.chain, con.resId, con.residue, atom.chain, atom.resId, atom.residue))
                    elif Qi < 0:
                        hydrophobic_negative.add((con.chain, con.resId, con.residue, atom.chain, atom.resId,
                                                  atom.residue))
                    if hydroElectroOutput:
                        print(','.join((con.chain, str(con.resId), con.residue, atom.chain, str(atom.resId),
                                        atom.residue, '%.3f' % R)), file=hydroElectroOutput)
    if hydroElectroOutput:
        hydroElectroOutput.close()
    return len(hydrophobic_positive), len(hydrophobic_negative)
