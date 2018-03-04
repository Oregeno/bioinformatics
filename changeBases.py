#!/usr/bin/python

import sys
import numpy as np
import os
import subprocess
import re
import pexpect

def getName(path):

     name = re.findall("//([^//].*)\..*",path)

     return name[0]

def changeBases(changeFrom,changeTo,fileGiven):

     current = os.path.dirname(os.path.abspath(__file__))

     if fileGiven != 'all':

          with open(fileGiven,'r') as fasta:

               data = fasta.read()
               fasta.close()

          
          partitionedData = re.findall("(>.*\n)(.*)",data)
          basesChanged = partitionedData[0][1].replace(changeFrom,changeTo)

          with open(partitionedData[0][0]+'.fas','w') as fasta:

               fasta.write(partitionedData[0][0] + basesChanged)

               fasta.close()
          
     else:
          for fileName in os.listdir(current):
               if fileName.endswith(".fas"):
               
                    target = current+"//"+fileName
                    with open(target,'r') as fasta:

                         data = fasta.read()
                         fasta.close()


                    partitionedData = re.findall("(>.*\n)(.*)",data)
                    basesChanged = partitionedData[0][1].replace(changeFrom,changeTo)
                    name = getName(target)
                    with open("RNA_"+name+'.fas','w') as fasta:

                         fasta.write(partitionedData[0][0] + basesChanged)

                         fasta.close()


changeBases('A','G','all')
                          

