''' --- Biochemie LDH Enzyme Kinetics Experiment: Data Analysis --- '''

'''Written by Camillo Mizaikoff'''
'''Language used: Python (v3.6)'''

'''Disclaimer: This code was written in the course of one evening and serves its purpose as a means to an end.'''
'''			   It is not suited to handle large amounts of data efficiently and could probably use some editing to make it more efficient!'''
'''			   If you have questions/comments/improvements please send me an email at camillodavid@aol.com '''

'''In its current state the code is meant to loop through 2 folders 'pyruvat' and 'oxamat' which are stored in a common directory defined'''
'''by the variable 'path' below. The range() in the for/if loops defines how many files are in each folder (currently 6 + 6 = 12). '''
'''dE values are calculated by (i) plotting linear regressions through the absorbtion values of the first 50s aswell as (ii) plotting '''
'''Two-Point linear functions through t=0s and t=50s.'''

'''Filetype used is .txt (Exported as 'Data Point Table' from software running on lab computer)'''

import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd 


'''Class for Calculating Linear Regressions'''

class Simple_Linear_Regression:

    def __init__(self, predictor, response):

        self.data_1 = predictor
        self.data_2 = response
        

        self.mean_data_1 = np.average(self.data_1)
        self.mean_data_2 = np.average(self.data_2)


        self.data_1_minus_mean = np.array([])
        self.data_2_minus_mean = np.array([])
        self.x_times_y = np.array([])
        self.data_1_squared = np.array([])
        self.runa()
        self.sum_data_1_squared = np.sum(self.data_1_squared)
        self.sum_x_times_y = np.sum(self.x_times_y)
        self.sum_data_1 = np.sum(self.data_1)
        self.sum_data_2 = np.sum(self.data_2)
        self.a = 0
        self.b = 0
        self.runb()

    def data_minus_mean(self):
        self.data_1_minus_mean = np.subtract(self.data_1, self.mean_data_1)
        self.data_2_minus_mean = np.subtract(self.data_2, self.mean_data_2)
        return self.data_1_minus_mean
        return self.data_2_minus_mean

    def multiply_x_y(self):
        self.x_times_y = self.data_1_minus_mean * self.data_2_minus_mean
        return self.x_times_y

    def square_data_1_minus_mean(self):
        self.data_1_squared = self.data_1_minus_mean ** 2
        return self.data_1_squared

    def calculate_a(self):
        self.a = self.sum_x_times_y / self.sum_data_1_squared
        return self.a


    def calculate_b(self):
        self.b = self.mean_data_2 - self.a * self.mean_data_1
        return self.b

    def ouput_regression(self):
        print ('Linear Regression Function: y = a * x + b')
        print ('with a = ' + str(self.a))
        print ('and b = ' + str(self.b))

    def runa(self):
        
        self.data_minus_mean()
        self.multiply_x_y()
        self.square_data_1_minus_mean()

    def runb(self):
        self.calculate_a()
        self.calculate_b()
        self.ouput_regression()

'''Path to directory where files are stored'''
path = '/Users/camillo/Desktop/BC_Praktikum_Raw_Data/'

'''Variables for loop'''
number = 0

'''Data lists'''
incline_uninhibited_reg = []
v0_uninhibited_reg = []
incline_inhibited_reg = []
v0_inhibited_reg = []
incline_uninhibited_2p = []
v0_uninhibited_2p = []
incline_inhibited_2p = []
v0_inhibited_2p = []

'''Rounded Data Lists'''
incline_uninhibited_reg_round = []
v0_uninhibited_reg_round = []
incline_inhibited_reg_round = []
v0_inhibited_reg_round = []
incline_uninhibited_2p_round = []
v0_uninhibited_2p_round = []
incline_inhibited_2p_round = []
v0_inhibited_2p_round = []


'''Begin Looping Through Directories'''
'''13 = total number of files + 1'''
for i in range(1, 13):

	'''Pyruvate Directory'''
	'''i = number of files in directory'''
	if i < 7:
		number = str(i)
		experiment = 'pyruvat'


		'''Readout TXT Data'''
		df = pd.read_csv(path + experiment + '/' + experiment + '-00' + number + '.txt', sep='\t', names = ['time', 'abs'])



		time = []
		absorbtion = []

		count = 0

		'''Extract Time Values'''
		for i in df['time']:
			if type(i) == str:
				count += 1
			if count in range(4, 54):
				time.append(float(i.replace(',', '.')))



		count = 0

		'''Extract Absorbtion Values'''
		for i in df['abs']:
			if type(i) == str:
				count += 1
			if count in range(4, 54):
				absorbtion.append(float(i.replace(',', '.')))

		'''Loop to collect data from (i) Linear and (ii) Two-Point Regression'''
		for i in range(1, 3):
			if i == 1:

				'''Define Lists to Append Data to'''
				incline_list = incline_uninhibited_reg
				v0_list = v0_uninhibited_reg
				incline_list_round = incline_uninhibited_reg_round
				v0_list_round = v0_uninhibited_reg_round

				'''Determine linear function through linear regression of points abs(t=0) - abs(t=50)'''
				incline = Simple_Linear_Regression(time, absorbtion).a
				y_start = Simple_Linear_Regression(time, absorbtion).b


				'''Absorbtioncoefficient for results in [Micro Molar]'''
				absorbtionskoeff = 0.00622

				'''Calculate Absorbtion after 4 minutes'''
				abs_after_240s = incline*240 + y_start

				'''Calculate Difference in absorbtion from t=0 to t=4min'''
				deltaE = y_start - abs_after_240s
	
				'''Calculate initial change rate of absorbtion'''
				anstieg = deltaE/4

				'''Calculate initial speed'''
				v0 = anstieg/absorbtionskoeff

				'''Append Values to appropriate List'''
				incline_list.append(anstieg)
				v0_list.append(v0)
				incline_list_round.append(round(anstieg, 5))
				v0_list_round.append(round(v0, 2))

			else:

				'''Define Lists to Append Data to'''
				incline_list = incline_uninhibited_2p
				v0_list = v0_uninhibited_2p
				incline_list_round = incline_uninhibited_2p_round
				v0_list_round = v0_uninhibited_2p_round

				'''Determine linear function through abs(t=0) and abs(t=50)'''
				y_start = absorbtion[0]
				incline = (absorbtion[49] - absorbtion[0])/50

				'''Absorbtioncoefficient for results in [Micro Molar]'''
				absorbtionskoeff = 0.00622

				'''Calculate Absorbtion after 4 minutes'''
				abs_after_240s = incline*240 + y_start

				'''Calculate Difference in absorbtion from t=0 to t=4min'''
				deltaE = y_start - abs_after_240s
	
				'''Calculate initial change rate of absorbtion'''
				anstieg = deltaE/4

				'''Calculate initial speed'''
				v0 = anstieg/absorbtionskoeff

				'''Append Values to appropriate List'''
				incline_list.append(anstieg)
				v0_list.append(v0)
				incline_list_round.append(round(anstieg, 5))
				v0_list_round.append(round(v0, 2))

		'''Oxamate Directory'''
		'''i - 6 = total number of files - number of files in pyruvat directory'''
	else:  
		number = str(i - 6)
		experiment = 'oxamat'


		df = pd.read_csv(path + experiment + '/' + experiment + '-00' + number + '.txt', sep='\t', names = ['time', 'abs'])



		time = []
		absorbtion = []

		count = 0

		'''Extract Time Values'''
		for i in df['time']:
			if type(i) == str:
				count += 1
			if count in range(4, 54):
				time.append(float(i.replace(',', '.')))



		count = 0

		'''Extract Absorbtion Values'''
		for i in df['abs']:
			if type(i) == str:
				count += 1
			if count in range(4, 54):
				absorbtion.append(float(i.replace(',', '.')))

		'''Loop to collect data from (i) Linear and (ii) Two-Point Regression'''
		for i in range(1, 3):
			if i == 1:

				'''Define Lists to Append Data to'''
				incline_list = incline_inhibited_reg
				v0_list = v0_inhibited_reg
				incline_list_round = incline_inhibited_reg_round
				v0_list_round = v0_inhibited_reg_round

				'''Determine linear function through linear regression of points abs(t=0) - abs(t=50)'''
				incline = Simple_Linear_Regression(time, absorbtion).a
				y_start = Simple_Linear_Regression(time, absorbtion).b


				'''Absorbtioncoefficient for results in [Micro Molar]'''
				absorbtionskoeff = 0.00622

				'''Calculate Absorbtion after 4 minutes'''
				abs_after_240s = incline*240 + y_start

				'''Calculate Difference in absorbtion from t=0 to t=4min'''
				deltaE = y_start - abs_after_240s
	
				'''Calculate initial change rate of absorbtion'''
				anstieg = deltaE/4

				'''Calculate initial speed'''
				v0 = anstieg/absorbtionskoeff

				'''Append Values to appropriate List'''
				incline_list.append(anstieg)
				v0_list.append(v0)
				incline_list_round.append(round(anstieg, 5))
				v0_list_round.append(round(v0, 2))

			else:

				'''Define Lists to Append Data to'''
				incline_list = incline_inhibited_2p
				v0_list = v0_inhibited_2p
				incline_list_round = incline_inhibited_2p_round
				v0_list_round = v0_inhibited_2p_round

				'''Determine linear function through abs(t=0) and abs(t=50)'''
				y_start = absorbtion[0]
				incline = (absorbtion[49] - absorbtion[0])/50

				'''Absorbtioncoefficient for results in [Micro Molar]'''
				absorbtionskoeff = 0.00622

				'''Calculate Absorbtion after 4 minutes'''
				abs_after_240s = incline*240 + y_start

				'''Calculate Difference in absorbtion from t=0 to t=4min'''
				deltaE = y_start - abs_after_240s
	
				'''Calculate initial change rate of absorbtion'''
				anstieg = deltaE/4

				'''Calculate initial speed'''
				v0 = anstieg/absorbtionskoeff

				'''Append Values to appropriate List'''
				incline_list.append(anstieg)
				v0_list.append(v0)
				incline_list_round.append(round(anstieg, 5))
				v0_list_round.append(round(v0, 2))

'''Print Rounded Values for Linear- and Two-Point Regression'''

print('')
print('')
print('')
print('')
print('')
print('')
print('')
print('')
print('')
print('With Linear Regression:')
print('dE/dt without Inhibitor [1/min]' + str(incline_uninhibited_reg_round))
print('V0 without Inhibitor [µM]' + str(v0_uninhibited_reg_round))
print('')
print('dE/dt with Inhibitor [1/min]' + str(incline_inhibited_reg_round))
print('V0 with Inhibitor [µM]' + str(v0_inhibited_reg_round))
print('')
print('')
print('')
print('With Two-Point-Regression:')
print('dE/dt without Inhibitor [1/min]' + str(incline_uninhibited_2p_round))
print('V0 without Inhibitor [µM]' + str(v0_uninhibited_2p_round))
print('')
print('dE/dt with Inhibitor [1/min]' + str(incline_inhibited_2p_round))
print('V0 with Inhibitor [µM]' + str(v0_inhibited_2p_round))

'''Plot Figures'''

concentration_pyruvat = [66.67, 166.67, 333.33, 666.67, 1333.33, 2333.33]
minor_ticks_x = []
minor_ticks_y = []

for i in range(1, 3001, 100):
	minor_ticks_x.append(i)

for i in range(1, 25):
	minor_ticks_y.append(i)

'''Define Figure Matrix'''
fig = plt.figure(figsize = (25,7))
fig.suptitle('LDH Enzyme Kinetics - Pyruvate Hydrogenation', fontsize = 15)

'''Plot Graph Without Oxamate'''
ax1 = fig.add_subplot(1, 2, 1)
ax1.scatter(concentration_pyruvat, v0_uninhibited_reg_round, s=10, c='red')
ax1.set_ylabel('Reaction Speed V [µM/min]', fontsize = 15)
ax1.set_xlabel('Concentration (Pyruvate) [µM]', fontsize = 15)
ax1.set_title('Without Oxamate Solution', fontsize = 10)
ax1.tick_params('both', which = 'major', labelsize = 10)
ax1.tick_params('both', which = 'minor', labelsize = 10)
ax1.set_yticks(minor_ticks_y, minor = 'true')
ax1.set_xticks(minor_ticks_x, minor = 'true')
ax1.legend(fontsize = 20)
ax1.set_ylim([0,25])
ax1.set_xlim([0,3000])
ax1.grid(which = 'major', color = 'grey', linestyle = '-', linewidth = 0.75)
ax1.grid(which = 'minor', color = 'grey', linestyle = '-', linewidth = 0.5)
ax1.set_facecolor('whitesmoke')

'''Plot Graph With Oxamate'''
ax2 = fig.add_subplot(1, 2, 2)
ax2.scatter(concentration_pyruvat, v0_inhibited_reg_round, s=10, c='red')
ax2.set_ylabel('Reaction Speed V [µM/min]', fontsize = 15)
ax2.set_xlabel('Concentration (Pyruvate) [µM]', fontsize = 15)
ax2.set_title('With Oxamate Solution', fontsize = 10)
ax2.tick_params('both', which = 'major', labelsize = 10)
ax2.tick_params('both', which = 'minor', labelsize = 10)
ax2.set_yticks(minor_ticks_y, minor = 'true')
ax2.set_xticks(minor_ticks_x, minor = 'true')
ax2.legend(fontsize = 20)
ax2.set_ylim([0,25])
ax2.set_xlim([0,3000])
ax2.grid(which = 'major', color = 'grey', linestyle = '-', linewidth = 0.75)
ax2.grid(which = 'minor', color = 'grey', linestyle = '-', linewidth = 0.5)
ax2.set_facecolor('whitesmoke')

'''Display Graphs'''
plt.show()
