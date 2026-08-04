import primer3
import pandas as pd

# define the location
up_loc = "results/qpcr/upreg_primers.tsv"
down_loc = "results/qpcr/downreg_primers.tsv"

# read the data
up_reg = pd.read_csv(up_loc, sep="\t")
down_reg = pd.read_csv(down_loc, sep="\t")

up_reg["hairpin_left"] = up_reg["LeftPrimer"].apply(lambda x: primer3.calc_hairpin_tm(x))
up_reg["hairpin_right"] = up_reg["RightPrimer"].apply(lambda x: primer3.calc_hairpin_tm(x))

down_reg["hairpin_left"] = down_reg["LeftPrimer"].apply(lambda x: primer3.calc_hairpin_tm(x))
down_reg["hairpin_right"] = down_reg["RightPrimer"].apply(lambda x: primer3.calc_hairpin_tm(x))

up_reg.to_csv(path_or_buf="results/qpcr/upreg_primers_hairpin.tsv", sep="\t", index=False)
down_reg.to_csv(path_or_buf="results/qpcr/downreg_primers_hairpin.tsv", sep="\t", index=False)