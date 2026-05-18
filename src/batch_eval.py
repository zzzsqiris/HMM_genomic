import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("model_file")
parser.add_argument("data_dir")
args = parser.parse_args()

model_file = args.model_file
data_dir = args.data_dir

sn_sum = 0
sp_sum = 0
acc_sum = 0
num_files = 0

for filename in os.listdir(data_dir):
    if filename.endswith(".fa"):
        base_name = filename.replace(".fa", "")

        fa_file = os.path.join(data_dir, filename)
        gff_file = os.path.join(data_dir, base_name + ".gff3")
        pred_file = os.path.join("build", base_name + ".pred.gff3")
        eval_file = os.path.join("build", base_name + ".eval.txt")

        if not os.path.exists(gff_file):
            continue

        pred_cmd = "python src/HMM_pred.py " + model_file + " " + fa_file
        os.system(pred_cmd)

        eval_cmd = "python src/evaluation.py " + pred_file + " " + gff_file
        os.system(eval_cmd)

        with open(eval_file, "r") as f:
            header = f.readline()
            line = f.readline().strip()
            parts = line.split("\t")

        sn = float(parts[0])
        sp = float(parts[1])
        acc = float(parts[2])

        sn_sum += sn
        sp_sum += sp
        acc_sum += acc

        num_files += 1

        # print(base_name, "Sn:", sn, "Sp:", sp, "Acc:", acc)

if num_files == 0:
    print("No files found")
else:
    ave_sn = sn_sum / num_files
    ave_sp = sp_sum / num_files
    ave_acc = acc_sum / num_files

    print("Number of files:", num_files)
    print("Average")
    print("Sn:", round(ave_sn, 4))
    print("Sp:", round(ave_sp, 4))
    print("Acc:", round(ave_acc, 4))
