
# import matplotlib.pyplot as plt
# import numpy as np
#
# # Fixing random state for reproducibility
# np.random.seed(19680801)
#
# fig, ax = plt.subplots()
#
# # Example data
# people = ('Tom', 'Dick', 'Harry', 'Slim', 'Jim')
# y_pos = np.arange(len(people))
# performance = 3 + 10 * np.random.rand(len(people))
# error = np.random.rand(len(people))
#
# ax.barh(y_pos, performance, xerr=error, align='center')
# ax.set_yticks(y_pos, labels=people)
# ax.invert_yaxis()  # labels read top-to-bottom
# ax.set_xlabel('Performance')
# ax.set_title('How fast do you want to go today?')
#
# plt.show()

import matplotlib.pyplot as plt
import numpy as np
#
# data = {'Barton LLC': 109438.50,
#         'Frami, Hills and Schmidt': 103569.59,
#         'Fritsch, Russel and Anderson': 112214.71,
#         'Jerde-Hilpert': 112591.43,
#         'Keeling LLC': 100934.30,
#         'Koepp Ltd': 103660.54,
#         'Kulas Inc': 137351.96,
#         'Trantow-Barrows': 123381.38,
#         'White-Trantow': 135841.99,
#         'Will LLC': 104437.60}
# data={
#         'R-NO3': 0,
#         'R-N=C=O': 1,
#         'R-N=C=S': 5,
#         'R-O-PO3H2': 21,
#         'R-S-S-R': 31,
#         'R-S(=O)-R': 54,
#         'R-SO3H': 68,
#         'R3P': 76,
#         'C≡C': 93,
#         'R-SH': 188,
#         'R-SO2-R': 190,
#         'R-NO2': 192,
#         'R-COOH': 396,
#         'R-S-R': 411,
#         'R-CHO': 562,
#         'R-COO-R': 709,
#         'R-CONH2': 762,
#         'C6H5OH': 840,
#         'C=C': 1115,
#         'R-O-R': 1467,
#         'X (F, Cl, Br, I)': 1619,
#         'R-OH': 2221,
#         'R-CO-R': 2315,
#         'C6H5': 2702,
#         'R-NH2': 2825,
#         'C-H' : 7495,
#       }
#
# group_data = list(data.values())
# group_names = list(data.keys())
# group_mean = np.mean(group_data)
#
# fig, ax = plt.subplots()
# ax.barh(group_names, group_data)
# labels = ax.get_xticklabels()
# plt.setp(labels, rotation=45, horizontalalignment='right')
# ax.set(xlim=[0, 8000], xlabel='functional grounp number', ylabel='functional group',
#        title='707 molecule attention functional group ')
# plt.show()


#fragemnts
data={
   'C':	11119,
   'O': 4057,
   'c' : 3005,
   'CC':2165,
   'N': 1821,
   'Cl': 1640,
   'CO':1067,
   'C=O':1045,
   'cc': 1018,
   'CCC':748,
   'CN':581,
   'ccc':565,
   'Br':475,
   'CCO':439,
   'CC=O':413,
   'F':400,
    'cccc':387,
    '[O-]':385,
    'CCCC': 310,
   'O=[N+][O-]':302,
    'C=C':281,
    'ccccc':279,
    'c1ccccc1':275,
    'CCN': 273,
    'CCl':268,
    '[N+]=O':209,
     '[N+]':192,
     'O=CO':191,
     'Cc':186,
     'C#N':177,
      }

group_data = list(data.values())
group_names = list(data.keys())
group_mean = np.mean(group_data)

fig, ax = plt.subplots()
ax.barh(group_names, group_data)
labels = ax.get_xticklabels()
plt.setp(labels, rotation=45, horizontalalignment='right')
ax.set(xlim=[0, 12000], xlabel='number', ylabel='fragment',
       title='84992 attention/1328 molecules fragments ')
plt.tight_layout()
plt.show()


# ## functional group and fragments
# data={
# 'C - H' :30880,
# 'R - NH2' :8184,
# 'C6H5' : 6629,
# 'R - OH' : 5945,
# 'R - CO - R' : 5659,
# 'Cl': 5442,
# 'X(F':5442,
# 'Br':5442,
# 'I) ' :5442,
# 'O':4057,
# 'R - O - R':3651,
# 'c':3005,
# 'C = C':2899,
# 'C6H5OH':2007,
# 'R - CONH2':1878,
# 'N':1821,
# 'R - COO - R':1559,
# 'R - CHO':1395,
# 'C = O':1045,
# 'cc':1018,
# 'R - COOH':1015,
# 'R - S - R':719,
# 'ccc':565,
# 'R - NO2':540,
# 'R - SO2 - R':437,
# 'cccc':387,
# '[O -]':385,
# 'R - SH':313,
# 'O = [N +][O -]':302,
# 'ccccc':279,
# }
#
# group_data = list(data.values())
# group_names = list(data.keys())
# group_mean = np.mean(group_data)
#
# fig, ax = plt.subplots()
# ax.barh(group_names, group_data)
# labels = ax.get_xticklabels()
# plt.setp(labels, rotation=45, horizontalalignment='right')
# ax.set(xlim=[0, 31000], xlabel='number', ylabel='functional group and fragments',
#        title='84992 attention/1328 molecules functional group and fragments  ')
# plt.tight_layout()
# plt.show()
#


