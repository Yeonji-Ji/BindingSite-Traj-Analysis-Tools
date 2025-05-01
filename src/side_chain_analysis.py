import mdtraj as md
import numpy as np

from src.atom_search import sc_atom_search
from src.adjust_indices import adjust_index


class SideChainAnalysis:

    def __init__(self, topology_file, trajectory_file, rigid_avg_file, flexible_avg_file, resi_list):
        self.topology = md.load_topology(topology_file)
        self.trajectory = trajectory_file
        self.rigid_avg = md.load(rigid_avg_file)
        self.flexible_avg = md.load(flexible_avg_file)

        self.trj_resi = adjust_index(topology_file, resi_list)
        self.avg_resi = adjust_index(rigid_avg_file, resi_list)


    ### return xyz coordinates of at_top of traj/frame
    def list_coords(self, traj, at_top, resn):
        coords_list = {}

        if sc_atom_search(resn) is not None:
            atom_name_list = sorted(sc_atom_search(resn))

            for idx in at_top:
                name = traj.top.atom(idx).name
                if name in atom_name_list:
                    crds = traj.xyz[0, idx, :]*10
                    coords_list[name] = crds

            return coords_list
        else:
            return


    ### return time average coordinates of rigid/flexible structure
    def avg_coords(self, idx):

        ahr_resn = self.rigid_avg.top.residue(idx).name
        ahr_atom_top = self.rigid_avg.top.select("residue " + str(idx + 1))
        ahr_coords = self.list_coords(self.rigid_avg, ahr_atom_top, ahr_resn)

        bbr_resn = self.flexible_avg.top.residue(idx).name
        bbr_atom_top = self.flexible_avg.top.select("residue " + str(idx + 1))
        bbr_coords = self.list_coords(self.flexible_avg, bbr_atom_top, bbr_resn,)

        avg_dic = {"rigid": ahr_coords, "flexible": bbr_coords}

        return avg_dic


    ### return the atom coordinates of the specified residue at the input frame
    def trj_coords(self, i, idx):

        traj = md.load(self.trajectory, frame=i, top=self.topology)
        res = self.topology.select("residue " + str(idx))
        resn = self.topology.residue(idx).name
        trj_coords = self.list_coords(traj, res, resn)

        return trj_coords


    ### calculate distance between atoms and determine the conformation
    def conformation(self, i):
        conf_dic = {}

        for idx in range(len(self.trj_resi)):
            # print(idx)
            trj_idx = self.trj_resi[idx]
            avg_idx = self.avg_resi[idx]

            # print(self.topology.residue(trj_idx))
            # print(self.rigid_avg.topology.residue(avg_idx))


            ahr_coords = self.avg_coords(avg_idx)["rigid"]
            bbr_coords = self.avg_coords(avg_idx)["flexible"]
            trj_coords = self.trj_coords(i, trj_idx)

            # print(ahr_coords)
            # print(trj_coords)

            if trj_coords is not None:
                keys = list(trj_coords.keys())

                ahr_distance, bbr_distance = 0, 0
                for at in keys:
                    dist_a = np.linalg.norm(ahr_coords[at] - trj_coords[at]) * 10
                    dist_b = np.linalg.norm(bbr_coords[at] - trj_coords[at]) * 10
                    ahr_distance += dist_a
                    bbr_distance += dist_b

                if ahr_distance <= bbr_distance:
                    conf = "ahr"
                else:
                    conf = "bbr"


                residue = self.topology.residue(trj_idx)
                conf_dic[residue] = conf

        return conf_dic
