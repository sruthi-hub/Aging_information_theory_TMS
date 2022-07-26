import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

########################################################################


def Split_By_Ages(data):
    three = data[0:1102]
    eighteen = data[1102:2623]
    twentyfour = data[2623:3855]
    return three, eighteen, twentyfour

def Return_Prob(TF, TG):    
    size = np.size(TF)
    A_count = 0
    B_count = 0
    C_count = 0
    D_count = 0
    sum = 0
    for i in range(0,size):
        if TF[i] == 0 and TG[i] == 0:
            A_count += 1
            sum += 1
        elif TF[i] != 0 and TG[i] == 0:
            B_count += 1
            sum += 1
        elif TF[i] == 0 and TG[i] != 0:
            C_count += 1
            sum += 1
        else:
            D_count += 1
            sum += 1
    
    a = A_count/sum
    b = B_count/sum
    c = C_count/sum
    d = D_count/sum
    return a, b, c, d


########################################################################################
#Empty arrays to append data
#arrays are named to show half <y>_young is on then half <y>_old is on 
x0_zero_half_zero_half = 0
x0_zero_half_one_half = 0
x0_one_half_zero_half = 0
x0_one_half_one_half = 0

x1_zero_half_zero_half = 0
x1_zero_half_one_half = 0
x1_one_half_zero_half = 0
x1_one_half_one_half = 0

########################################################################################
#Imported data
data = pd.read_csv('clean_ordered_facscountmatrix8dec.csv')
pairs = pd.read_csv('top_genes.csv')
num_of_pairs = pairs.shape[0]

#Looping over data and appending to arrays
for i in range(0,num_of_pairs):

    TF = np.array(data[pairs['TF'][i]])
    TG = np.array(data[pairs['TG'][i]])

    #Splitting up age groups
    TF3, TF18, TF24 = Split_By_Ages(TF)
    TG3, TG18, TG24 = Split_By_Ages(TG)

    a3, b3, c3, d3 = Return_Prob(TF3, TG3)
    a24, b24, c24, d24 = Return_Prob(TF24, TG24)

    #Checking Noise for x = 0:
    y_barx03 = c3/(a3+c3)
    y_barx024 = c24/(a24+c24)
    sigmax03 = np.sqrt(y_barx03*(1-y_barx03))
    sigmax024 = np.sqrt(y_barx024*(1-y_barx024))

    if sigmax024 < sigmax03:
        if y_barx03 <= .5:
            if y_barx024 <= .5:
                x0_zero_half_zero_half += 1
            else: #y_barx024 > .5
                x0_zero_half_one_half += 1
        else: #y_barx03 > .5
            if y_barx024 <= .5:
                x0_one_half_zero_half += 1
            else: #y_barx024 > .5
                x0_one_half_one_half += 1
    
    #Checking Noise for x = 1:
    y_barx13 = d3/(b3+d3)
    y_barx124 = d24/(b24+d24)
    sigmax13 = np.sqrt(y_barx13*(1-y_barx13))
    sigmax124 = np.sqrt(y_barx124*(1-y_barx124))

    if sigmax124 < sigmax13:
        if y_barx13 <= .5:
            if y_barx124 <= .5:
                x1_zero_half_zero_half += 1
            else: #y_barx124 > .5
                x1_zero_half_one_half += 1
        else: #y_barx13 > .5
            if y_barx124 <= .5:
                x1_one_half_zero_half += 1
            else: #y_barx124 > .5
                x1_one_half_one_half += 1


#####################################################
#printing results
print("x0 00", x0_zero_half_zero_half)
print("x0 01", x0_zero_half_one_half)
print("x0 10", x0_one_half_zero_half)
print("x0 11", x0_one_half_one_half)

#printing results
print("x1 00", x1_zero_half_zero_half)
print("x1 01", x1_zero_half_one_half)
print("x1 10", x1_one_half_zero_half)
print("x1 11", x1_one_half_one_half)

