

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
    df = pd.read_csv(csv_file,  sep=',', header=0)

    outpath = abspath + "/analysis"




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

    # density
    ax = df.loc[:,["Test_ID","disordered_score"]].plot(kind="density")
    #df.loc[:,["pdb","DockQ_AF_0"]].plot(kind="density")
    ax.set_xlabel('DockQscores')
    ax.set_xlim([-0.25, 1.25])
    #ax.set_ylim([0, 1.6])
    plt.savefig(outpath+"/density.png")

    """
    # scatterplot combined
    sns_plot=sns.scatterplot(data=df,x="DockQ_AF_0",y="DockQ_OMF_100X")
    sns_plot=sns.scatterplot(data=df,x="DockQ_AF_0",y="DockQ_OMF_201G")
    sns_plot.text(0.3,0.85,"Omega 201G: "+str(round(df["DockQ_OMF_201G"].mean(),3)))
    sns_plot.text(0.3,0.8,"Omega 100X: "+str(round(df["DockQ_OMF_100X"].mean(),3)))
    sns_plot.text(0.3,0.75,"AF-0: "+str(round(df["DockQ_AF_0"].mean(),3)))

    sns_plot.set(ylim=(0, 1))
    sns_plot.set(xlim=(0, 1))
    plt.legend(loc='lower right')
    plt.plot([0,1],[0,1])
    plt.savefig(outpath+"/scatter.png")


    # histogram
    nberbins = 30
    ax = df.plot.hist(column='DockQ_AF_0', color="b", bins=nberbins)
    df.plot.hist(column='DockQ_OMF_100X',color="g", ax=ax,bins=nberbins,alpha=0.7)
    df.plot.hist(column=['DockQ_OMF_201G'],color="r", ax=ax,bins=nberbins,alpha=0.7)
    ax.set_xlabel('DockQscores')
    ax.set_title('Histogram with %d bins' % nberbins)
    plt.savefig(outpath+"/hist.png")

    # density
    ax = df.loc[:,["pdb","DockQ_AF_0"]].plot(kind="density")
    ax.set_xlabel('DockQscores')
    ax.set_xlim([-0.25, 1.25])
    #ax.set_ylim([0, 1.6])
    plt.savefig(outpath+"/density.png")

    # scatterplot of DockQscores
    ax = df.plot.scatter(x='pdb', y='DockQ_AF_0', color="b",xticks=[],marker="+")
    df.plot.scatter(x='pdb', y='DockQ_AF_0_antibody',color="g", ax=ax,marker="x",alpha=0.7)
    #df.plot.scatter(x=['pdb'], y=['DockQ_OMF_100X'],color="g", ax=ax,marker="x")
    #df.plot.scatter(x=['pdb'], y=['DockQ_OMF_100X'],color="r", ax=ax,marker="-")
    ax.set_ylabel('DockQscores')
    plt.savefig(outpath+"/range.png")


    # # 
    # ax = df.plot.scatter(x='year', y='DockQ_AF_0', color="b",xticks=[],marker="+")
    # #df.plot.scatter(x=['pdb'], y=['DockQ_OMF_100X'],color="g", ax=ax,marker="x")
    # #df.plot.scatter(x=['pdb'], y=['DockQ_OMF_100X'],color="r", ax=ax,marker="-")
    # ax.set_ylabel('DockQscores')
    # plt.savefig(outpath+"/yearly.png")

    # correlation
    with open(logfile, "a") as logf:
        logf.write("# Correlation\n")
        s = df.corr(method='pearson').to_string()
        logf.write(s)

    """


if __name__ == '__main__':
    main()
