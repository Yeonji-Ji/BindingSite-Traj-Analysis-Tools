import os
import time
import multiprocessing
import subprocess

from scripts.update_output import update_totals
from src.hbond_analysis import HbondAnalysis
from src.side_chain_analysis import SideChainAnalysis

DIC = {"ace":[482, 485, 473, 345, 346, 315, 318],
       "adrb1":[261, 242, 83, 173],
       "akt1":[88, 92, 86],
       "braf":[17, 86, 55],
       "def":[133, 134, 45, 51, 88, 90, 92],
       "fgfr1":[179, 100, 102],
       "fpps":[232, 101, 246, 189],
       "gcr":[41, 119, 47],
       "grik1":[138, 92, 85, 87],
       "hdac2":[136, 298, 144],
       "hmdh":[668, 674, 675, 296, 120, 126, 316, 574],
       "hxk4":[213, 61],
       "kith":[169, 171, 173, 119, 88],
       "mapk2":[48, 162, 96],
       "mcr":[51, 220, 45],
       "mk10":[25, 106, 109]}

rest = {"ahr": "rigid", "bbr": "flexible"}
num_frames = range(80000, 100000)

for target in DIC.keys():
    resi_list = DIC[target]

    ahr_avg_file = "/home/yeonji/Dropbox/myfolder_data/wbp_last/ahr_eq/avg/" + target + "_ahr_eq_l20_avg.pdb"
    bbr_avg_file = "/home/yeonji/Dropbox/myfolder_data/wbp_last/bbr_avg_l20/" + target + "_bbr_l20_average.pdb"

    ahr_path = "/gibbs/yeonji/DUDE/" + target.upper() + "/simulation/ahr_eq/"
    ahr_top_file = ahr_path + target + "_min.prmtop"
    ahr_traj_file = ahr_path + target + "_ahr_eq.nc"

    bbr_path = "/gibbs/vjay/WBP-DUDE/" + target.upper() + "/simulation/"
    bbr_top_file = bbr_path + target + ".prmtop"
    bbr_traj_file = bbr_path + target + "_bbr.nc"


    AHR_Hbond = HbondAnalysis(ahr_top_file, ahr_traj_file, resi_list)
    BBR_Hbond = HbondAnalysis(bbr_top_file, bbr_traj_file, resi_list)
    # AHR_Conf_Analysis = SideChainAnalysis(ahr_top_file, ahr_traj_file, ahr_avg_file, bbr_avg_file, resi_list)
    # BBR_Conf_Analysis = SideChainAnalysis(bbr_top_file, bbr_traj_file, ahr_avg_file, bbr_avg_file, resi_list)


    def ahr_water_process_frame(i):
        don = AHR_Hbond.num_of_donate(i)
        acc = AHR_Hbond.num_of_accept(i)
        neighbors = AHR_Hbond.water_neighbors(i)

        result = {'don': don, 'acc': acc, 'neighbors': neighbors}

        return result

    def bbr_water_process_frame(i):
        don = BBR_Hbond.num_of_donate(i)
        acc = BBR_Hbond.num_of_accept(i)
        neighbors = BBR_Hbond.water_neighbors(i)

        result = {'don': don, 'acc': acc, 'neighbors': neighbors}

        return result

    # def ahr_conf_process_frame(i):
    #     conf = AHR_Conf_Analysis.conformation(i)
    #
    #     result = {'conf': conf}
    #
    #     return result
    #
    # def bbr_conf_process_frame(i):
    #     conf = BBR_Conf_Analysis.conformation(i)
    #
    #     result = {'conf': conf}
    #
    #     return result

    def main():
        start_time = time.time()
        with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
            ahr_output_water = pool.map(ahr_water_process_frame, num_frames)
            # ahr_output_conf = pool.map(ahr_conf_process_frame, num_frames)

            bbr_output_water = pool.map(bbr_water_process_frame, num_frames)
            # bbr_output_conf = pool.map(bbr_conf_process_frame, num_frames)

            df_water_ahr = update_totals(ahr_output_water)
            df_water_bbr = update_totals(bbr_output_water)
            # df_conf_ahr = update_totals(ahr_output_conf)
            # df_conf_bbr = update_totals(bbr_output_conf)

            outdir = "results/edited/" + target.upper() + "/"

            if not os.path.isdir("results/edited/"+target.upper()):
                subprocess.run("mkdir " + outdir, shell=True)

            df_water_ahr.to_csv(outdir + target+"_"+rest["ahr"]+"_PL_Hbonds_water.csv")
            df_water_bbr.to_csv(outdir + target + "_" + rest["bbr"] + "_PL_Hbonds_water.csv")

            # df_conf_ahr.to_csv(outdir + target+"_"+rest["ahr"]+"_sc_confs.csv")
            # df_conf_bbr.to_csv(outdir + target+"_"+rest["bbr"]+"_sc_confs.csv")


        print(f"'{target}' is done")
        end_time = time.time()
        print(end_time - start_time)



    if __name__ == "__main__":
        main()







