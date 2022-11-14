

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

import seaborn as sns


def mean(s):
    # mean
    mean = s.mean()
    mean = round(mean,2)
    return mean

def max(s):
    # max
    mean = s.max()
    mean = round(mean,2)
    return max

def min(s):
    # min
    mean = s.min()
    mean = round(mean,2)
    return min

def greater05(s):
    # precentage of scores higher or equal than 0.5
    g05 = s.apply(lambda score: 1 if score >= 0.5 else 0).sum()
    g05 = g05/s.shape[0]
    g05 = g05*100
    g05 = round(g05,1)
    return g05

def greater023(s):
    # precentage of scores higher or equal than 0.23
    g023 = s.apply(lambda score: 1 if score >= 0.23 else 0).sum()
    g023 = g023/s.shape[0]
    g023 = g023*100
    g023 = round(g023,1)
    return g023

def write_metrics(filename,s,label=""):
    with open(filename, "a") as f:
        f.write(label+"\n")
        f.write('\tn= {} \
            \n\tMean: {} \
            \n\tPercentage of scores >= 0.5: {} % \
            \n\tPercentage of scores >= 0.23: {} %\
            '.format(s.shape[0],mean(s), greater05(s), greater023(s)))
        f.write("\n\n")

def main():

    abspath = os.getcwd()
    abspath = abspath.split("GPCR_project")[0]+"GPCR_project"

    csv_file = abspath + "/analysis/data_results.csv"
    df = pd.read_csv(csv_file, dtype="str", sep=',', header=0)

    outpath = abspath + "/analysis"

    # prepare dataframe
    df["relative_ASA_below_025"] = df["relative_ASA_below_025"].astype(float)
    df["normalized_alignment_score"] = df["normalized_alignment_score"].astype(float)
    df["average_plDDT_score"] = df["average_plDDT_score"].astype(float)
    df["plDDT_below_025"] = df["plDDT_below_025"].astype(float)
    df["Antigen_Length"] = df["Antigen_Length"].astype(float)




    logfile = outpath + "/info.md"
    with open(logfile, "w") as logf:
        logf.write("### Metrics\n")
    

    # check number available entries
    with open(logfile, "a") as logf:
        logf.write("# Datapoints\n")
        logf.write("Alignment score greater or equal 0.90: {}/{}".format(df[df['normalized_alignment_score'] > 0.90].shape[0],df.shape[0]))
        logf.write("\n")

    # CUTOFF 
    # ignore all values with normalized_alignment_score < 0.90
    df = df[df["normalized_alignment_score"] >= 0.90]


    # disorderd score 
    ax = df.loc[:,["disordered_score"]].plot(kind="density")
    df[df["On_target"] == "TRUE"].loc[:,["disordered_score"]].plot(kind="density", ax = ax,)
    df[df["On_target"] == "FALSE"].loc[:,["disordered_score"]].plot(kind="density", ax = ax,)
    ax.set_xlabel('disordered_score')
    ax.set_xlim([-0.25, 1.25])
    ax.legend(["all", "On_target = TRUE", "On_target = FALSE"])
    #ax.set_ylim([0, 1.6])
    plt.savefig(outpath+"/on_target_disordered.png")

    # disorderd score 
    ax = df.loc[:,["disordered_score"]].plot(kind="density")
    df[df["Off_target"] == "TRUE"].loc[:,["disordered_score"]].plot(kind="density", ax = ax)
    df[df["Off_target"] == "FALSE"].loc[:,["disordered_score"]].plot(kind="density", ax = ax)
    ax.set_xlabel('disordered_score')
    ax.set_xlim([-0.25, 1.25])
    ax.legend(["all", "Off_target = TRUE", "Off_target = FALSE"])
    #ax.set_ylim([0, 1.6])
    plt.savefig(outpath+"/off_target_disordered.png")

    # disorderd score 
    ax = df.loc[:,["disordered_score"]].plot(kind="density")
    df[df["Specific"] == "TRUE"].loc[:,["disordered_score"]].plot(kind="density", ax = ax)
    df[df["Specific"] == "FALSE"].loc[:,["disordered_score"]].plot(kind="density", ax = ax)
    ax.set_xlabel('disordered_score')
    ax.set_xlim([-0.25, 1.25])
    ax.legend(["all", "Specific = TRUE", "Specific = FALSE"])
    #ax.set_ylim([0, 1.6])
    plt.savefig(outpath+"/Specific_disordered.png")



if __name__ == '__main__':
    main()
