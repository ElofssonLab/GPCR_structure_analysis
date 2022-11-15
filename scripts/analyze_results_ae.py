

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
    df["Antigen_Length"] = df["Antigen_Length"].astype(float)
    #df["Antigen_Length"] = df["Antigen_Length"]/df["Antigen_Length"].max()

    df["normalized_alignment_score"] = df["normalized_alignment_score"].astype(float)

    df["average_plDDT_score"] = df["average_plDDT_score"].astype(float)
    df["average_plDDT_score"] = df["average_plDDT_score"]/100 # normalize values: plDDT scores are between 0-100

    df["plDDT_below_50"] = df["plDDT_below_50"].astype(float)


    df["average_relative_ASA"] = df["average_relative_ASA"].astype(float)
    df["relative_ASA_below_025"] = df["relative_ASA_below_025"].astype(float)
    df["relative_H-dssp"] = df["relative_H-dssp"].astype(float)
    df["relative_E-dssp"] = df["relative_E-dssp"].astype(float)
    df["relative_C-dssp"] = df["relative_C-dssp"].astype(float)


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


    value_columns = [ \
       'average_plDDT_score', 'plDDT_below_50',  \
       'average_relative_ASA', 'relative_ASA_below_025', 'relative_H-dssp', \
       'relative_E-dssp', 'relative_C-dssp']

    # disorderd score 
    #plt.style.use('seaborn-deep')
    for value in value_columns:
        #for type in ["Specific", "On_target", "Off_target"]:
        #ax = df.loc[:,[value]].plot(kind="density", color="blue")
        ax=df[(df.On_target=="TRUE") &  (df.Off_target=="FALSE")].loc[:,[value]].plot(kind="density", color="green",legend="on_target")
        df[ df.Off_target=="TRUE"].loc[:,[value]].plot(kind="density", ax = ax, color="red",label="off_target")
        df[(df.On_target=="FALSE") & (df.Off_target=="FALSE")].loc[:,[value]].plot(kind="density", ax = ax, color="blue",label="no_target")
        #df[df[type] == "FALSE"].loc[:,[value]].plot(kind="density", ax = ax, color="red")
        #ax.set_title("TRUE/FALSE: {}/{}".format(df[df[type] == "TRUE"].shape[0], df[df[type] == "FALSE"].shape[0]))
        ax.set_title(value)
        #ax.text(0.8,0.80, "FALSE: {}".format(df[df[type] == "FALSE"].shape[0]))
        #ax.text(0.8,0.60, "TRUE: {}".format(df[df[type] == "TRUE"].shape[0]))
        ax.set_xlabel(value)
        ax.set_xlim([-0.25, 1.25])
        #ax.legend()
        ax.legend(["On target","Off target","No target"])
        plt.savefig(outpath + "/" + "all" +"_"+value+".png")



    """
    # 
    ax = df.loc[:,["relative_ASA_below_025"]].plot(kind="density")
    df[df["On_target"] == "TRUE"].loc[:,["relative_ASA_below_025"]].plot(kind="density", ax = ax)
    df[df["On_target"] == "FALSE"].loc[:,["relative_ASA_below_025"]].plot(kind="density", ax = ax)
    ax.text(0.8,0.80, "FALSE: {}".format(df[df["On_target"] == "FALSE"].shape[0]))
    ax.text(0.8,0.60, "TRUE: {}".format(df[df["On_target"] == "TRUE"].shape[0]))
    #ax.text(0.3,0.75, "on target (FALSE/TRUE/TOTAL): {}/{}/{}".format(df[df["On_target"] == "FALSE"].shape[0], df[df["On_target"] == "TRUE"].shape[0], df.shape[0]))
    ax.set_xlabel('relative_ASA_below_025')
    ax.set_xlim([-0.25, 1.25])
    ax.legend(["all", "On_target = TRUE", "On_target = FALSE"])
    #ax.set_ylim([0, 1.6])
    plt.savefig(outpath+"/on_target_relative_ASA_below_025.png")

    # disorderd score 
    ax = df.loc[:,["relative_ASA_below_025"]].plot(kind="density")
    df[df["Off_target"] == "TRUE"].loc[:,["relative_ASA_below_025"]].plot(kind="density", ax = ax)
    df[df["Off_target"] == "FALSE"].loc[:,["relative_ASA_below_025"]].plot(kind="density", ax = ax)
    ax.set_xlabel('relative_ASA_below_025')
    ax.set_xlim([-0.25, 1.25])
    ax.legend(["all", "Off_target = TRUE", "Off_target = FALSE"])
    #ax.set_ylim([0, 1.6])
    plt.savefig(outpath+"/off_target_relative_ASA_below_025.png")

    # disorderd score 
    ax = df.loc[:,["relative_ASA_below_025"]].plot(kind="density")
    df[df["Specific"] == "TRUE"].loc[:,["relative_ASA_below_025"]].plot(kind="density", ax = ax)
    df[df["Specific"] == "FALSE"].loc[:,["relative_ASA_below_025"]].plot(kind="density", ax = ax)
    ax.set_xlabel('relative_ASA_below_025')
    ax.set_xlim([-0.25, 1.25])
    ax.legend(["all", "Specific = TRUE", "Specific = FALSE"])
    #ax.set_ylim([0, 1.6])
    plt.savefig(outpath+"/Specific_relative_ASA_below_025.png")
    """


if __name__ == '__main__':
    main()
