import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

# Baseline = pd.read_csv("../DATA_envelope_baseline/BSxxxx1--00000.csv")
ENV = pd.read_csv("../data_envelop/DATA_envelope_baseline/ENVxxx1--00000.csv")
BS = pd.read_csv("../data_envelop/DATA_envelope_baseline/BSxxxx1--00000.csv")
peak_df = pd.DataFrame()

# find HIGH LOW Logic
ENV = ENV.drop('time', axis=1)
BS = BS.drop('time', axis=1)
ENV['data2'] = BS
bit = []
for index, row in ENV.iterrows():
    if(row["data"] < row["data2"]):
        bit.append(0.092)
    elif(row["data"] > row["data2"]):
        bit.append(0.097)
# export to csv
bit_df = pd.DataFrame(bit, columns=['data'])
# bit_df.to_csv(r'C:\Users\nut\PycharmProjects\KeyGen\venv\ENV_data.csv', index = False)
DataLine = bit_df
# DataLine = pd.read_csv("ENV_data.csv")


DataLine_array = DataLine['data'].values
class thresholdline:
    #pos val
    th0 = pd.DataFrame(0, index=range(len(DataLine)), columns=range(1))
    th2 = pd.DataFrame(1.93, index=range(len(DataLine)), columns=range(1))
    th1 = (th2 / 2)
    th3 = (3 / 2) * th2
    th4 = th2 * 2
    th6 = th3 * 2
    th8 = th4 * 2
    #negative val
    nth2 = pd.DataFrame(-1.4, index=range(len(DataLine)), columns=range(1))
    nth1 = (nth2 / 2)
    nth3 = (3 / 2) * nth2
    nth4 = nth2 * 2
    nth6 = nth3 * 2
    nth8 = nth4 * 2
    def plot(self):
        plt.plot(self.th0)
        plt.plot(self.th1)
        plt.plot(self.th2)
        plt.plot(self.th3)
        plt.plot(self.th4)
        plt.plot(self.th6)
        plt.plot(self.nth1)
        plt.plot(self.nth2)
        plt.plot(self.nth3)
        plt.plot(self.nth4)
        plt.plot(self.nth6)


threshold1 = thresholdline()
peak = []
width = []
# temp = []
count = 0
temp_count = 0
cal_count = []
for i in range(1, len(DataLine_array)):
    if DataLine_array[i] > DataLine_array[i - 1]:
        peak.append(0.1)
        temp = [-1 * (count / 100)] * (count + 1)
        cal_count = cal_count + [-1 * (count / 100)]
        width = width + temp
        count = 0
    elif DataLine_array[i] < DataLine_array[i - 1]:
        peak.append(0.1)
        temp = [count / 100] * (count+ 1)
        cal_count = cal_count + [(count / 100)]
        width = width + temp
        count = 0
    else:
        peak.append(0.09)
        count = count + 1

raw_data = []
# print(thresholdline.th1.at[1,0])
# print(thresholdline.th2.at[1,0])
# print(thresholdline.th3.at[1,0])
# print(thresholdline.th4.at[1,0])
# print(thresholdline.th6.at[1,0])

for i in range(len(cal_count)):
    # range of positive bit
    if (cal_count[i] <= thresholdline.th1.at[1, 0] + 0.1) & (cal_count[i] >= thresholdline.th1.at[1, 0] - 0.1):
        raw_data = raw_data + [1]
    elif (cal_count[i] <= thresholdline.th2.at[1, 0] + 0.1) & (cal_count[i] >= thresholdline.th2.at[1, 0] - 0.1):
        raw_data = raw_data + [1,1]
    elif (cal_count[i] <= thresholdline.th3.at[1, 0] + 0.1) & (cal_count[i] >= thresholdline.th3.at[1, 0] - 0.1):
        raw_data = raw_data + [1, 1, 1]
    elif (cal_count[i] <= thresholdline.th4.at[1, 0] + 0.1) & (cal_count[i] >= thresholdline.th4.at[1, 0] - 0.1):
        raw_data = raw_data + [1, 1, 1, 1]
    elif (cal_count[i] <= thresholdline.th6.at[1, 0] + 0.1) & (cal_count[i] >= thresholdline.th6.at[1, 0] - 0.1):
        raw_data = raw_data + [1, 1, 1, 1, 1, 1]
    # range of negative bit
    elif (cal_count[i] <= thresholdline.nth1.at[1, 0] + 0.1) & (cal_count[i] >= thresholdline.nth1.at[1, 0] - 0.1):
        raw_data = raw_data + [0]
    elif (cal_count[i] <= thresholdline.nth2.at[1, 0] + 0.1) & (cal_count[i] >= thresholdline.nth2.at[1, 0] - 0.1):
        raw_data = raw_data + [0, 0]
    elif (cal_count[i] <= thresholdline.nth3.at[1, 0] + 0.1) & (cal_count[i] >= thresholdline.nth3.at[1, 0] - 0.1):
        raw_data = raw_data + [0, 0, 0]
    elif (cal_count[i] <= thresholdline.nth4.at[1, 0] + 0.1) & (cal_count[i] >= thresholdline.nth4.at[1, 0] - 0.1):
        raw_data = raw_data + [0, 0, 0, 0]
    elif (cal_count[i] <= thresholdline.nth6.at[1, 0] + 0.1) & (cal_count[i] >= thresholdline.nth6.at[1, 0] - 0.1):
        raw_data = raw_data + [0, 0, 0, 0, 0, 0]
    else:
        raw_data = raw_data + ['start']

# print(width)
print(cal_count)
print(raw_data)
plt.plot(width)
threshold1.plot()
# plt.plot(peak)
# plt.plot(BS)
# plt.plot(ENV)
DataLine['data'] = 100 * DataLine['data']
DataLine['data'] = DataLine['data'] - 8
# plt.plot(DataLine)
plt.show()



# find HIGH LOW Logic
# bit = []
# for index, row in ENV.iterrows():
#     if(row["data"] < row["data2"]):
#         bit.append(0.092)
#     elif(row["data"] > row["data2"]):
#         bit.append(0.097)
# # export to csv
# bit_df = pd.DataFrame(bit, columns=['data'])
# bit_df.to_csv(r'C:\Users\nut\PycharmProjects\KeyGen\venv\ENV_data.csv', index = False)