# -*- coding: utf-8 -*-

import numpy as np
import os
import matplotlib.pyplot as plt
import csv



def plot_PCA_top2(paired,mountain):

    paired_array = np.load(paired)
    mountain_array = np.load(mountain)

    figure1 = plt.figure()
            
    ind = np.arange(2)    # the x locations for the groups
    width = 0.5       # the width of the bars: can also be len(x) sequence

    p1 = plt.bar(ind[0], np.mean(paired_array), width, yerr=np.std(paired_array), color = 'red')

    p2 = plt.bar(ind[1], np.mean(mountain_array), width, yerr=np.std(mountain_array), color = 'blue')

    plt.xticks(ind, ('Paired Sites', 'Mountain'))

    plt.xlabel('RNA metric')
    plt.ylabel('Percentage of variance')
    plt.title('Percentage of variance represented in the 2D PCA plot')

               
    plt.savefig("Variance_PCA_vs_metric.png")

plot_PCA_top2(paired="paired_top2.npy", mountain="mountain_top2.npy")
