import math

import pandas as pd
import mdtraj as md
import numpy as np

import warnings
import glob

warnings.simplefilter('ignore', DeprecationWarning)

from src.adjust_indices import adjust_index
from src.atom_search import atom_search

########################################################################################################################
########################################################################################################################

class PLConfAnalysis:

    def __init__(self, topology_file, trajectory_file, resi_list, lig_file, criteria=1.3):
        self.topology = md.load_topology(topology_file)
        self.trajectory = trajectory_file
        self.trj_watO = self.topology.select("water and name O")
        self.resi = adjust_index(topology_file, resi_list)
        self.lig = md.load(lig_file)
        self.criteria = criteria
        
        self.lig_don_h = {}      # {H-index: don_N_index}
        self.lig_acc = []        # {index}
        self.prot_don_h = {}     # {H-index: don_N_index}
        self.prot_acc = []       # {index}
        
        self.PL_hbond = []       # Pairs of [P, L] within the criteria of distance and angle.



    ### find H bound to atom (O, N, S) in protein
    def findH(self, trj_frame, pair):

        pair_list = np.array([[pair[0], pair[1]]])
        dist = md.compute_distances(trj_frame, pair_list)
        roundup = round(dist[0][0] * 10, 1)

        # if dist * 10 <= criteria:
        if roundup <= self.criteria:
            return pair

    def ligand_fcn(self):
        
        lig_top = self.lig.topoogy
        lig_at = lig_top.select("element N or element O or element S or element P")
        lig_H = lig_top.select("element H")

        # any possible pair
        lig_at_H = [[x, y] for x in lig_at for y in lig_H]

        for pair in lig_at_H:
            bond_pair = self.findH(self.lig, pair)
            if bond_pair:
                self.lig_don_h[pair[0]] = pair[1]

        for at in lig_at:
            if at not in self.lig_don_h.values():
                self.lig_acc.append(at)

        return


    ### find pairs to be involved in hbonds (idx: residue index)
    def protein_fcn(self, idx):

        sol_acc = []
        sol_don = []

        first_frame = md.load(self.trajectory, frame=0, top=self.topology)

        resn = self.topology.residue(idx).name
        atom_list = atom_search(resn)

        for at in atom_list:
            res_at = self.topology.select("residue " + str(idx) + " and name " + str(at))
            res_H = self.topology.select("residue " + str(idx) + " and element H")
            all_pairs = [[x, y] for x in res_at for y in res_H]

            for pair in all_pairs:
                if self.findH(first_frame, pair):
                    if pair not in sol_don:
                        sol_don.append(pair)
                else:
                    if (pair[0] not in sol_acc) and (self.topology.atom(pair[0]).element.name == "oxygen"):
                        sol_acc.append(pair[0])
                    if (pair[0] not in sol_acc) and (self.topology.atom(pair[0]).element.name == "sulfur"):
                        sol_acc.append(pair[0])

        solute_hbond_pairs = {"sol_acc": sol_acc, "sol_don": sol_don}

        return solute_hbond_pairs


    ### return angle triplet (acc, don_at, don_H)
    def get_angle(self, trj_frame, triplet_list):
        angle_list = []

        for triplet in triplet_list:
            arr_triplet = np.array([triplet])
            angle = round(math.degrees(md.compute_angles(trj_frame, arr_triplet)[0][0]), 3)
            angle_list.append(angle)

        return angle_list


    def prot_lig_interaction(self, distance=3.6, angle=30):

        traj = md.load(self.trajectory, frame=i, top=self.topology)
        num_sol_don = {}

        for idx in self.resi:
            sol_don = self.protein_fcn(idx)["sol_don"]

            for pair in sol_don:
                atom = self.topology.atom(pair[0])

                within_d = md.compute_neighbors(traj, 0.35, query_indices=[pair[0]])

                wat_within_d = []
                for at in within_d[0]:
                    if at in self.trj_watO:
                        wat_within_d.append(at)


                for watO in wat_within_d:
                    ndon = 0
                    triplet_list = [[watO, pair[0], pair[1]]]
                    angle_list = self.get_angle(traj, triplet_list)
                    for angle in angle_list:
                        if angle <= 30:
                            ndon += 1


                    if atom in num_sol_don.keys():
                        num_sol_don[atom] += ndon
                    if atom not in num_sol_don.keys():
                        num_sol_don[atom] = ndon


        return num_sol_don




        return