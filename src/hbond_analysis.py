import mdtraj as md
import numpy as np
import math

from src.adjust_indices import adjust_index
from src.atom_search import atom_search


class HbondAnalysis:

    def __init__(self, topology_file, trajectory_file, resi_list, criteria=1.3):
        self.topology = md.load_topology(topology_file)
        self.trajectory = trajectory_file
        self.trj_watO = self.topology.select("water and name O")
        self.resi = adjust_index(topology_file, resi_list)
        self.criteria = criteria


    ### find H bound to atom (O, N, S) in protein
    def findH(self, trj_frame, pair):

        pair_list = np.array([[pair[0], pair[1]]])
        dist = md.compute_distances(trj_frame, pair_list)
        roundup = round(dist[0][0] * 10, 1)

        # if dist * 10 <= criteria:
        if roundup <= self.criteria:
            return pair


    ### return angle triplet (acc, don_at, don_H)
    def get_angle(self, trj_frame, triplet_list):
        angle_list = []

        for triplet in triplet_list:
            arr_triplet = np.array([triplet])
            angle = round(math.degrees(md.compute_angles(trj_frame, arr_triplet)[0][0]), 3)
            angle_list.append(angle)

        return angle_list


    ### find pairs to be involved in hbonds (idx: residue index)
    def find_pair(self, idx):

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



    ### count number of HB donate (i: frame)
    def num_of_donate(self, i):

        traj = md.load(self.trajectory, frame=i, top=self.topology)

        num_sol_don = {}

        for idx in self.resi:
            sol_don = self.find_pair(idx)["sol_don"]

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



    ### count number of accepting HB (i: frame)
    def num_of_accept(self, i):

        traj = md.load(self.trajectory, frame=i, top=self.topology)

        num_sol_acc = {}

        for idx in self.resi:
            sol_acc = self.find_pair(idx)["sol_acc"]

            for acc in sol_acc:
                atom = self.topology.atom(acc)
                within_d = md.compute_neighbors(traj, 0.35, query_indices=np.asarray([acc]))

                wat_within_d = []
                for at in within_d[0]:
                    if at in self.trj_watO:
                        wat_within_d.append(at)

                for watO in wat_within_d:
                    nacc = 0
                    triplet_list = [[acc, watO, watO + 1], [acc, watO, watO + 2]]
                    angle_list = self.get_angle(traj, triplet_list)
                    for angle in angle_list:
                        if angle <= 30:
                            nacc += 1


                    if atom in num_sol_acc.keys():
                        num_sol_acc[atom] += nacc
                    if atom not in num_sol_acc.keys():
                        num_sol_acc[atom] = nacc

                    

        return num_sol_acc


    def water_neighbors(self, i):

        wat_neighbors = {}

        traj = md.load(self.trajectory, frame=i, top=self.topology)
        for idx in self.resi:
            resn = self.topology.residue(idx).name      # eg) GLU, GLN
            atom_list = atom_search(resn)       # residue atom name list: [N, O, NE2, OE1]

            for at in atom_list:
                nwat = 0
                at_idx = self.topology.select("residue " + str(idx) + " and name " + str(at))[0]
                atom = self.topology.atom(at_idx)       # atom name eg) GLU67-OE1
                within_d = md.compute_neighbors(traj, 0.36, query_indices=np.asarray([at_idx]))

                for wat in within_d[0]:
                    if wat in self.trj_watO:
                        nwat += 1

                if atom in wat_neighbors.keys():
                    wat_neighbors[atom] += nwat
                if atom not in wat_neighbors.keys():
                    wat_neighbors[atom] = nwat
        
        return wat_neighbors

