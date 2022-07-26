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
    for i in range(0,size):
        if TF[i] == 0 and TG[i] == 0:
            A_count += 1
        elif TF[i] != 0 and TG[i] == 0:
            B_count += 1
        elif TF[i] == 0 and TG[i] != 0:
            C_count += 1
        else:
            D_count += 1
    
    a = A_count/size
    b = B_count/size
    c = C_count/size
    d = D_count/size
    return a, b, c, d


########################################################################################
#Empty arrays to append data

range_youngs = np.array([])
range_olds = np.array([])
range_diffs = np.array([])
range_ratios = np.array([])
std_0_ratios = np.array([])
std_1_ratios = np.array([])
std_0_diffs = np.array([])
std_1_diffs = np.array([])
difference_Px0s = np.array([])
range_decrease = np.array([])
range_increase = np.array([])
std_0_decrease = np.array([])
std_0_increase = np.array([])
std_1_decrease = np.array([])
std_1_increase = np.array([])


#########################################################################################
#Imported data
data = pd.read_csv('clean_ordered_facscountmatrix8dec.csv')
pairs = pd.read_csv('alltftgpairs.csv')
data2 = np.load("4State_All_Results.npz")
num_of_pairs = pairs.shape[0]
pair_index = data2['pair_index']


#Looping over data and appending to arrays
for i in range(0,num_of_pairs):
    if i in pair_index:

        TF = np.array(data[pairs['TF'][i]])
        TG = np.array(data[pairs['TG'][i]])

        #Splitting up age groups
        TF3, TF18, TF24 = Split_By_Ages(TF)
        TG3, TG18, TG24 = Split_By_Ages(TG)

        a3, b3, c3, d3 = Return_Prob(TF3, TG3)
        a24, b24, c24, d24 = Return_Prob(TF24, TG24)

        #if (b3+d3) != 0 and (a3+c3) != 0 and (b24+d24) != 0 and (a24+c24) != 0:
        #if a3 != 0 and b3 != 0 and c3 != 0 and d3 != 0 and a24 != 0 and b24 != 0 and c24 != 0 and d24 != 0:

        range_young = d3/(b3+d3) - c3/(a3+c3)
        range_old = d24/(b24+d24) - c24/(a24+c24)
        range_diffs = np.append(range_diffs, (np.abs(range_old) - np.abs(range_young)))

        if np.abs(range_young) > np.abs(range_old):
            range_decrease = np.append(range_decrease, 1)
        else:
            range_increase = np.append(range_increase, 1)
        
        if range_young != 0:
            range_youngs = np.append(range_youngs, range_young)
            range_olds = np.append(range_olds, range_old)
            range_ratio = np.abs(range_old/range_young)
            #if range_ratio < 10:
            range_ratios = np.append(range_ratios, range_ratio)

        std_0_3 = np.sqrt(a3*c3)/(a3+c3)
        std_0_24 = np.sqrt(a24*c24)/(a24+c24)
        std_0_diffs = np.append(std_0_diffs, (std_0_24 - std_0_3))

        if std_0_3 < std_0_24:
            std_0_increase = np.append(std_0_increase, 1)
        else:
            std_0_decrease = np.append(std_0_decrease, 1)

        if std_0_3 != 0:
            std_0_ratio = std_0_24/std_0_3
            std_0_ratios = np.append(std_0_ratios, std_0_ratio)

        std_1_3 = np.sqrt(b3*d3)/(b3+d3)
        std_1_24 = np.sqrt(b24*d24)/(b24+d24)
        std_1_diffs = np.append(std_1_diffs, (std_1_24 - std_1_3))

        if std_1_3 < std_1_24:
            std_1_increase = np.append(std_1_increase, 1)
        else:
            std_1_decrease = np.append(std_1_decrease, 1)

        if std_1_3 != 0:
            std_1_ratio = std_1_24/std_1_3
            std_1_ratios = np.append(std_1_ratios, std_1_ratio)


        difference_Px0 = a3+c3-a24-c24
        difference_Px0s = np.append(difference_Px0s, np.abs(difference_Px0))

#######################################################################################

#print(np.size(range_ratios))
#print(np.size(std_0_ratios))
#print(np.size(std_1_ratios))
#print(np.size(difference_Px0s))


#Analyzing data
#ratios
ave_range_ratios = np.average(range_ratios)
ave_std_0_ratios = np.average(std_0_ratios)
ave_std_1_ratios = np.average(std_1_ratios)
ave_std_ratios = .5*(ave_std_0_ratios + ave_std_1_ratios)


#differences
ave_range_diffs = np.average(range_diffs)
ave_std_0_diffs = np.average(std_0_diffs)
ave_std_1_diffs = np.average(std_1_diffs)
ave_std_diffs = .5*(ave_std_0_diffs + ave_std_1_diffs)
ave_difference_Px0s = np.average(difference_Px0s)

#print(ave_range_ratios, np.size(range_ratios))
#print(ave_std_0_ratios, np.size(std_0_ratios))
#print(ave_std_1_ratios, np.size(std_1_ratios))
#print(ave_std_ratios)
#print(ave_difference_Px0s, np.size(difference_Px0s))
#print(" ")
#print(range_ratios)
#print(" ")
#print(std_0_ratios)
#print(" ")
#print(std_1_ratios)

#Table 2

print("Average R_diff: ", ave_range_diffs)
print("Average std_0_diff: ", ave_std_0_diffs)
print("Average std_1_diff: ", ave_std_1_diffs)
#print(ave_std_diffs)

I3s = data2['I3s']
I24s = data2['I24s']
pair_index_cross_check = np.array([])
for i in range(0,np.size(range_ratios)):
    if range_ratios[i] > 1:
        if I24s[i] - I3s[i] > 0:
            pair_index_cross_check = np.append(pair_index_cross_check, i)

#print(pair_index_cross_check)
#print(np.size(pair_index_cross_check))

I_ratios_324 = np.array([])
for i in range(0, np.size(I3s)):
    I_ratio_324 = I24s[i]/I3s[i]
    #if I_ratio_324 < 10:
    I_ratios_324 = np.append(I_ratios_324, I_ratio_324)

#print(np.size(I_ratios_324))
#print(I_ratios_324)
ave_I_ratios_324 = np.average(I_ratios_324)
#print(ave_I_ratios_324)

C3s = data2['C3s']
C24s = data2['C24s']

C_ratios_324 = np.array([])
for i in range(0, np.size(C3s)):
    C_ratio_324 = C24s[i]/C3s[i]
    #if I_ratio_324 < 10:
    C_ratios_324 = np.append(C_ratios_324, C_ratio_324)


ave_C_ratios_324 = np.average(C_ratios_324)
#print(ave_C_ratios_324)


#print(np.std(range_ratios, ddof=1))

#print(np.max(range_diffs))
#print(np.min(range_diffs))
#print(np.std(range_diffs, ddof=1))


# Table 1
#Percentages
num_range_decrease = np.sum(range_decrease)
num_range_increase = np.sum(range_increase)
total_in_range =  num_range_decrease + num_range_increase
print("total pairs used in range: ", total_in_range)
print("total with decreasing range: ", num_range_decrease)
print("total with increasing range: ", num_range_increase)

num_std_0_decrease = np.sum(std_0_decrease)
num_std_0_increase = np.sum(std_0_increase)
total_in_std_0 = num_std_0_decrease + num_std_0_increase
print("total pairs used in std 0: ", total_in_std_0)
print("total with decreasing std 0: ", num_std_0_decrease)
print("total with increasing std 0: ", num_std_0_increase)

num_std_1_decrease = np.sum(std_1_decrease)
num_std_1_increase = np.sum(std_1_increase)
total_in_std_1 = num_std_1_decrease + num_std_1_increase
print("total pairs used in std 1: ", total_in_std_1)
print("total with decreasing std 1: ", num_std_1_decrease)
print("total with increasing std 1: ", num_std_1_increase)




