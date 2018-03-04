#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess
import re
import pexpect
import datetime


T_list = [22,27,32,37,42,47,52]
num_sim = 10


ftl_wt = np.zeros((len(T_list),num_sim))
ftl_A56U = np.zeros((len(T_list),num_sim))
ftl_C10U = np.zeros((len(T_list),num_sim))
ftl_C14G = np.zeros((len(T_list),num_sim))
ftl_U22G = np.zeros((len(T_list),num_sim))

ftl_T1 = np.load('22.0.npy')

#print(ftl_T1)

ftl_T2 = np.load('27.0.npy')
ftl_T3 = np.load('32.0.npy')
ftl_T4 = np.load('37.0.npy')
ftl_T5 = np.load('42.0.npy')
ftl_T6 = np.load('47.0.npy')
ftl_T7 = np.load('52.0.npy')


ftl_wt[0] = ftl_T1[0]
ftl_wt[1] = ftl_T2[0]
ftl_wt[2] = ftl_T3[0]
ftl_wt[3] = ftl_T4[0]
ftl_wt[4] = ftl_T5[0]
ftl_wt[5] = ftl_T6[0]
ftl_wt[6] = ftl_T7[0]

ftl_A56U[0] = ftl_T1[1]
ftl_A56U[1] = ftl_T2[1]
ftl_A56U[2] = ftl_T3[1]
ftl_A56U[3] = ftl_T4[1]
ftl_A56U[4] = ftl_T5[1]
ftl_A56U[5] = ftl_T6[1]
ftl_A56U[6] = ftl_T7[1]

ftl_C10U[0] = ftl_T1[2]
ftl_C10U[1] = ftl_T2[2]
ftl_C10U[2] = ftl_T3[2]
ftl_C10U[3] = ftl_T4[2]
ftl_C10U[4] = ftl_T5[2]
ftl_C10U[5] = ftl_T6[2]
ftl_C10U[6] = ftl_T7[2]

ftl_C14G[0] = ftl_T1[3]
ftl_C14G[1] = ftl_T2[3]
ftl_C14G[2] = ftl_T3[3]
ftl_C14G[3] = ftl_T4[3]
ftl_C14G[4] = ftl_T5[3]
ftl_C14G[5] = ftl_T6[3]
ftl_C14G[6] = ftl_T7[3]

ftl_U22G[0] = ftl_T1[4]
ftl_U22G[1] = ftl_T2[4]
ftl_U22G[2] = ftl_T3[4]
ftl_U22G[3] = ftl_T4[4]
ftl_U22G[4] = ftl_T5[4]
ftl_U22G[5] = ftl_T6[4]
ftl_U22G[6] = ftl_T7[4]

figure1 = plt.figure()

names_list = ['wt', 'A56U', 'C10U', 'C14G', 'U22G']
color_list = ['blue','red', 'green', 'black', 'orange']

plt.errorbar(T_list, np.mean(abs(ftl_wt), axis=1), yerr =np.std(abs(ftl_wt), axis=1) , fmt='o', color = "%s" % color_list[0], label = "%s" % names_list[0])

plt.errorbar(T_list, np.mean(abs(ftl_A56U), axis=1), yerr =np.std(abs(ftl_A56U), axis=1) , fmt='o', color = "%s" % color_list[1], label = "%s" % names_list[1])

plt.errorbar(T_list, np.mean(abs(ftl_C10U), axis=1), yerr =np.std(abs(ftl_C10U), axis=1) , fmt='o', color = "%s" % color_list[2], label = "%s" % names_list[2])

plt.errorbar(T_list, np.mean(abs(ftl_C14G), axis=1), yerr =np.std(abs(ftl_C14G), axis=1) , fmt='o', color = "%s" % color_list[3], label = "%s" % names_list[3])

plt.errorbar(T_list, np.mean(abs(ftl_U22G), axis=1), yerr =np.std(abs(ftl_U22G), axis=1) , fmt='o', color = "%s" % color_list[4], label = "%s" % names_list[4])



plt.legend()
plt.xlabel("Temperature")
plt.ylabel("Entropy")
plt.title("RNA ensemble entropy as a function of temperature")

plt.savefig("Entropy_vs_temp.png")
