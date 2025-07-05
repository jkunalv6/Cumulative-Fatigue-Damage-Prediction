import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#from scipy.interpolate import interp1d

S_N_Data = pd.read_excel('C:/Users/jishn/Quick access/F&F code/S-N Data.xlsx')
Load_data = pd.read_excel('C:/Users/jishn/Quick access/F&F code/Load_data.xlsx')

n_cycles = int(input("Number of Load Cycles"))
growth_rate = .1


S_N_Data['Log'] = S_N_Data['N'].apply(lambda x: np.log10(x))


sum = 0

sorted_S_N_data = S_N_Data.sort_values(by='S')
data_log = []
datalog = []
iteration = 0
max_load = Load_data['load '].max()
max_s = S_N_Data['S'].max()
c_s = max_load / max_s
SN_R = sorted_S_N_data[~np.isinf(sorted_S_N_data['N'])]
#interp_func = interp1d(sorted_S_N_data['S'].values, sorted_S_N_data['N'].values, kind='linear', fill_value='extrapolate') # scipy interpolation useful for more complex data curves


while sum !=1 :
  Load_data['Stress1'] = Load_data['load '] / c_s
  Load_data['Effec3'] = np.interp(Load_data['Stress1'], SN_R['S'], SN_R['N']) # numpy linear interpolation 
  #Load_data['Effec'] = interp_func(Load_data['Stress'].values) # scipy interpolation
  Load_data['11/N'] = 1/Load_data['Effec3']
  Load_data['c1/N'] = Load_data['11/N'].mask(Load_data['11/N'] < 1e-6, 0)*Load_data['cycles '] * n_cycles
  sum = Load_data['c1/N'].sum().item()
  iteration += 1
  sum = round(sum , 3)
  data_log.append({'Iteration': iteration, 'Value': sum, 'Area' : round(c_s,3)})
  if sum == 1 :
    break
  elif sum < 1 :
    c_s = c_s*(1-growth_rate)
  elif sum > 1 :
    c_s = c_s*(1+growth_rate)


c_s1 = round(c_s , 3)
data_log = pd.DataFrame(data_log)
print(f"Suggested Cross Section according Palmgren - Miner Theory is ",c_s1)
sum = 0 
c_s = max_load / max_s

while sum !=1 :

    Load_data['Stress'] = Load_data['load '] / c_s
    Load_data['Effec'] = np.interp(Load_data['Stress'], SN_R['S'], SN_R['N'])
    max_cycles = Load_data['Effec'].min()
    if max_cycles <= 730 : # if min value in vector is below 730 then use direct, else use the 14 * x **0.6 
        Load_data['Effec1'] =Load_data['Effec'] - 14 * np.power(Load_data['Effec'], 0.6)
    else :
        Load_data['Effec1'] =Load_data['Effec'] - 14 * np.power(Load_data['Effec'], 0.6)        
    Load_data['1/N'] = 1/Load_data['Effec1']
    
    threshold = 1/(1e6 -(14*(1e6)**.6))
    Load_data['2/N'] = np.where(Load_data['1/N'] < threshold, 0, Load_data['1/N'])
    Load_data['c/N'] = Load_data['2/N']*Load_data['cycles ']
    Load_data['c/N']
    repeated_load = pd.concat([Load_data] * n_cycles, ignore_index=True)

    repeated_load['cumulative_sum'] = repeated_load['c/N'].cumsum()
    stop_row = (repeated_load['cumulative_sum'] >= 1).idxmax()
    # Check if stop_row is 0 and handle it accordingly
    if stop_row == 0:
        leftover = (1 - 0) * repeated_load['Effec1'][stop_row] # Assume previous cumulative sum is 0 if stop_row is 0
    else:
        leftover = (1-repeated_load['cumulative_sum'][stop_row-1]) * repeated_load['Effec1'][stop_row]
    leftover = int(leftover)
    copy_repeated = repeated_load.copy(deep=True)

    copy_repeated = copy_repeated.iloc[stop_row:]
    copy_repeated =copy_repeated.drop("cumulative_sum" , axis = 'columns')
    copy_repeated =copy_repeated.drop("1/N" , axis = 'columns')
    copy_repeated =copy_repeated.drop("2/N" , axis = 'columns')
    copy_repeated =copy_repeated.drop("c/N" , axis = 'columns')
    copy_repeated['Np'] = copy_repeated['Effec'] - copy_repeated['Effec1']
    copy_repeated.loc[stop_row,'cycles '] = copy_repeated['cycles '][stop_row] - leftover
    copy_repeated['c/Np'] = copy_repeated['cycles ']/copy_repeated['Np']
    copy_repeated['c/Np'] = np.where(copy_repeated['Effec1'] > 1e6, 0, copy_repeated['c/Np'])
    copy_repeated['Sum'] = copy_repeated['c/Np'].cumsum()
    sum = copy_repeated['Sum'].max()
    datalog.append({'Value': sum, 'Area' : c_s , 'Fracture Point in concatenated DF' : stop_row})
    sum = round(sum, 3)
    if sum == 1 :
        break 
    elif sum < 1 :
        c_s = c_s*(1-growth_rate)
    elif sum > 1:
        c_s = c_s*(1+growth_rate)

c_s2 = round(c_s, 3)
print(f"Suggested Cross Section according Manson Double Linear Theory is ",c_s2)
datalog = pd.DataFrame(datalog)


fig, ax = plt.subplots(2, 2, figsize=(12, 14))  

ax[0, 0].scatter(data_log['Area'], data_log['Value'], color='b', label="Area v Iterations - Palmgren-Miner")
ax[0, 0].plot(data_log['Area'], [1]*len(data_log['Area']), 'r-', label="Cumulative sum = 1")


ax[0, 1].plot(S_N_Data['Log'], S_N_Data['S'], color='b', label="S-N")
ax[0, 1].set_ylabel('Stress in PSI')
ax[0, 1].set_xlabel('Cycles in 10^x')
ax[0, 1].set_title('S - N curve')

ax[1, 0].scatter(datalog['Area'], datalog['Value'], color='b', label="Area v Iterations - Manson Double Linear Theory")
ax[1, 0].set_xlim(datalog['Area'].min(), datalog['Area'].max())  
ax[1, 0].set_ylim(datalog['Value'].min(), datalog['Value'].max())
ax[1, 0].plot(datalog['Area'], [1]*len(datalog['Area']), 'r-', label="Cumulative sum = 1")
ax[1, 0].set_xlabel('Area')
ax[1, 0].set_ylabel('Value')
ax[1, 0].set_title('Area v Cumulative Sum')
ax[1, 0].legend()

# You still have one empty subplot: ax[1, 1]
# You can either hide it or use it:
ax[1, 1].axis('off')  # hide unused subplot

plt.show()