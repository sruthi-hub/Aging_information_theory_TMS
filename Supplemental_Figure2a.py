import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

########################################################################################


def Plot_PDF(TF,TG, bins):
    plt.hist2d(TF,TG, bins, density = True, cmap = 'Blues')
    plt.colorbar()
    plt.title("TF vs TG")
    plt.show()
    return


def Split_By_Ages(data):
    three = data[0:1102]
    eighteen = data[1102:2623]
    twentyfour = data[2623:3855]
    return three, eighteen, twentyfour

################################################################################################

#Imported data
data = pd.read_csv('clean_ordered_facscountmatrix8dec.csv')
pairs = pd.read_csv('alltftgpairs.csv')
num_of_pairs = pairs.shape[0]


Aatf = data['Aatf']
Trp53 = data['Trp53']

TF = np.array(Aatf)
TG = np.array(Trp53)

TF3, TF18, TF24 = Split_By_Ages(TF)
TG3, TG18, TG24 = Split_By_Ages(TG)



origin = np.array([])
x_axis = np.array([])
y_axis = np.array([])
interior_TF = np.array([])
interior_TG = np.array([])
for i in range(0,np.size(TF3)):
    if TF3[i] == 0 and TG3[i] == 0:
        origin = np.append(origin, 0)
    elif TF3[i] == 0 and TG3[i] != 0:
        y_axis = np.append(y_axis, TG[i])
    elif TF3[i] != 0 and TG3[i] == 0:
        x_axis = np.append(x_axis, TF3[i])
    else:
        interior_TF = np.append(interior_TF, TF3[i])
        interior_TG = np.append(interior_TG, TG3[i])

plt.figure()
plt.scatter(interior_TF, interior_TG, c = 'black')
plt.scatter(x_axis, x_axis*0, c = 'black')
plt.scatter(y_axis*0, y_axis, c = 'black')
plt.scatter(origin, origin, c = 'black')
plt.tick_params(which='both',top=True,right=True)
plt.xlabel("Transcription Factor")
plt.ylabel("Target")
plt.savefig("Aatf_vs_Trp53_black")
plt.show()

plt.figure(dpi = 200)
plt.scatter(interior_TF, interior_TG, c = 'orange', label = "State 4", alpha=0.3)
plt.scatter(y_axis*0, y_axis, c = 'green', label = "State 3", alpha=0.3)
plt.scatter(x_axis, x_axis*0, c = 'blue', label = "State 2", alpha=0.3)
plt.scatter(origin, origin, c = 'red', label = "State 1", alpha=0.3)
plt.tick_params(which='both',top=True,right=True)
#plt.xlabel("Transcription Factor")
plt.xticks(fontsize=15)
#plt.ylabel("Target Gene")
plt.yticks(fontsize=15)
plt.legend()
plt.savefig("Aatf_vs_Trp53_colored")
plt.show()




