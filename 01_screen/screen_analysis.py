#!/usr/bin/env python
# coding: utf-8

# In[4]:


import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import json
from matplotlib.patches import Rectangle
import scipy.stats as stats
from fxns import *

path = "Incucyte_raw_data_v2"

df1 = pd.DataFrame(index=['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P'], columns=['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24'])

df2 = pd.read_csv("44pairs_echo_instructions.csv")

measure_type = ['cell body cluster area.txt', 'cell body clusters.txt', 'neurite branch points per body cluster area.txt', 'neurite branch points per body cluster.txt', 'neurite branch points.txt', 'neurite length per body cluster area.txt', 'neurite length per body cluster.txt', 'neurite length.txt']

plate_list = ['plate 1-1', 'plate 1-2', 'plate 1-3', 'plate 2-1', 'plate 2-2', 'plate 2-3', 'plate 3-1', 'plate 3-2', 'plate 3-3', 'plate 4-1', 'plate 4-2', 'plate 4-3']


# # Pull in raw data files

# In[5]:


# load incucyte txt files into dfs with well names mapped back to echo instructions, stored in dictionary with tuple keys
replicates = {}
for folder in plate_list: #plate names
    for file in measure_type: #measurement types
        df = pd.read_table(path + '/' + folder + '/' + file, header = 1)
        df, cols = mapper(df, folder[6], df2)    
        df.drop('Date Time', axis=1, inplace=True)
        tupl = folder[6], folder[8], file[:-4]
        replicates[tupl] = df
print(f'The number of plates X metrics in dict(replicates) is {len(replicates.keys())}')


# # Clean the data

# In[6]:


# remove unnecessary rows, cast data as float, drops last two time points which were missed by the instrument on certain plates, etc.

for key in list(replicates.keys()):
        df = replicates[key]
        df = df.T
        df.columns = df.loc['Elapsed'].map(lambda x: round(x))
        df.drop(['Elapsed'], axis=0, inplace=True)
        df.index = df.index.str.replace('\.[0-9]+', '', regex=True)
        try:
            df.drop(138, axis=1, inplace=True)
        except KeyError:
            #print(f'{key} doesn\'t have 138hr time point')
            pass
        try:
            df.drop(144, axis=1, inplace=True)
        except KeyError:
            #print(f'{key} doesn\'t have 144hr time point')
            pass

        # Fill hidden spaces
        df = df.replace(r'\s+', 0, regex=True)
        # Remove nan values
        df = df.dropna()

        try:
            df = df.astype('float64')
        except ValueError:
            print(f'unable to cast as flt {key}')

        replicates[key] = df


# # Normalize the data

# In[7]:


#plot_replicate_swarm(replicates, 'miRNA-124', 120, kw='neurite length', title='before_normalization', savefig=False)
#plot_replicate_swarm(replicates, 'miRNA mimic pool (10 nM)', 120, kw='neurite', title='before_normalization', savefig=True)

replicates = cellmetric_normalizePI('miRNA mimic pool (10 nM)', 'siPLK1 (10 nM)', 120, replicates)
replicates = neuritemetric_normalizePA('miRNA mimic pool (10 nM)', 'ATRA (25 uM)', 120, replicates)

#plot_replicate_swarm(replicates, 'ATRA (25 uM)', 120, kw='neurite length', title='after_normalization', savefig=False)
#plot_replicate_swarm(replicates, 'miRNA mimic pool (10 nM)', 120, kw='neurite', title='after_normalization', savefig=True)


# # Agreement of replicates for single time point

# In[2]:


import numpy as np
from scipy import stats

#https://zhiyzuo.github.io/Pearson-Correlation-CI-in-Python/

def pearsonr_ci(x,y,alpha=0.001):
    ''' calculate Pearson correlation along with the confidence interval using scipy and numpy
    Parameters
    ----------
    x, y : iterable object such as a list or np.array
      Input for correlation calculation
    alpha : float
      Significance level. 0.05 by default
    Returns
    -------
    r : float
      Pearson's correlation coefficient
    pval : float
      The corresponding p value
    lo, hi : float
      The lower and upper bound of confidence intervals
    '''

    r, p = stats.pearsonr(x,y)
    r_z = np.arctanh(r)
    se = 1/np.sqrt(x.size-3)
    z = stats.norm.ppf(1-alpha/2)
    lo_z, hi_z = r_z-z*se, r_z+z*se
    lo, hi = np.tanh((lo_z, hi_z))
    return r, p, lo, hi



# Make a 3D scatter plot of the values in replicates[('1','1','neurite length')].loc[:, 120] versus replicates[('1','2','neurite length')].loc[:, 120] versus replicates[('1','3','neurite length')].loc[:, 120]
import scipy
def plot_neurite_length_comparison(replicates, plate, timepoint):
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111, projection='3d')
    # Only plot the values from the dataframes if the index only appears once in the dataframe
    rep1 = replicates[(plate,'1','neurite length')].loc[replicates[(plate,'1','neurite length')].index.duplicated() == False, timepoint]
    rep2 = replicates[(plate,'2','neurite length')].loc[replicates[(plate,'2','neurite length')].index.duplicated() == False, timepoint]
    rep3 = replicates[(plate,'3','neurite length')].loc[replicates[(plate,'3','neurite length')].index.duplicated() == False, timepoint]
    ax.scatter(rep1, rep2, rep3)
    ax.set_xlabel('Replicate 1')
    ax.set_ylabel('Replicate 2')
    ax.set_zlabel('Replicate 3')
    ax.set_title(f'Plate {plate} Neurite Length at Timepoint {timepoint}')
    ax.set_xlim(0, ax.get_xlim()[1])
    ax.set_ylim(0, ax.get_ylim()[1])
    # Draw a line of best fit and report the R value of the line
    x = np.array(rep1)
    y = np.array(rep2)
    z = np.array(rep3)
    A = np.c_[x, y, np.ones(x.shape)]
    C, _, _, _ = scipy.linalg.lstsq(A, z)    # coefficients
    xm = np.linspace(min(x), max(x), 20)
    ym = np.linspace(min(y), max(y), 20)
    zm = C[0]*xm + C[1]*ym + C[2]
    ax.plot(xm, ym, zm, color='red')
        # Print the correlation coefficients and p values, previously with scipy.stats.pearsonr but now with fxn for custom CIs
    data = {'Replicate 1 w/ 2': pearsonr_ci(x, y),
            'Replicate 1 w/ 3': pearsonr_ci(x, z),
            'Replicate 2 w/ 3': pearsonr_ci(y, z)}
    df = pd.DataFrame(data, index=['Correlation Coefficient', 'p Value', 'Low interval', 'High interval'])
    print(df)
    plt.show()

for i in [str(x) for x in range(1,5)]:
    plot_neurite_length_comparison(replicates, i, 120)


# # Agreement of replicates for time slice

# In[9]:


# Make a 3D scatter plot of the values in replicates[('1','1','neurite length')].loc[:, 120] versus replicates[('1','2','neurite length')].loc[:, 120] versus replicates[('1','3','neurite length')].loc[:, 120]
import scipy
def plot_neurite_length_comparison(replicates, plate, timepoint):
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111, projection='3d')
    # Only plot the values from the dataframes if the index only appears once in the dataframe
    rep1 = replicates[(plate,'1','neurite length')].loc[replicates[(plate,'1','neurite length')].index.duplicated() == False, timepoint].mean(axis=1)
    rep2 = replicates[(plate,'2','neurite length')].loc[replicates[(plate,'2','neurite length')].index.duplicated() == False, timepoint].mean(axis=1)
    rep3 = replicates[(plate,'3','neurite length')].loc[replicates[(plate,'3','neurite length')].index.duplicated() == False, timepoint].mean(axis=1)
    ax.scatter(rep1, rep2, rep3)
    ax.set_xlabel('Replicate 1')
    ax.set_ylabel('Replicate 2')
    ax.set_zlabel('Replicate 3')
    ax.set_title(f'Plate {plate} Neurite Length at Timepoint {timepoint}')
    ax.set_xlim(0, ax.get_xlim()[1])
    ax.set_ylim(0, ax.get_ylim()[1])
    # Draw a line of best fit and report the R value of the line
    x = np.array(rep1)
    y = np.array(rep2)
    z = np.array(rep3)
    A = np.c_[x, y, np.ones(x.shape)]
    C, _, _, _ = scipy.linalg.lstsq(A, z)    # coefficients
    xm = np.linspace(min(x), max(x), 20)
    ym = np.linspace(min(y), max(y), 20)
    zm = C[0]*xm + C[1]*ym + C[2]
    ax.plot(xm, ym, zm, color='red')
        # Print the correlation coefficients and p values, previously with scipy.stats.pearsonr but now with fxn for custom CIs
    data = {'Replicate 1 w/ 2': pearsonr_ci(x, y),
            'Replicate 1 w/ 3': pearsonr_ci(x, z),
            'Replicate 2 w/ 3': pearsonr_ci(y, z)}
    df = pd.DataFrame(data, index=['Correlation Coefficient', 'p Value', 'Low interval', 'High interval'])
    print(df)
    plt.show()

for i in [str(x) for x in range(1,5)]:
    plot_neurite_length_comparison(replicates, i, slice(72, 126))
    plot_neurite_length_comparison(replicates, i, slice(96, 126))


# # Collect all samples from all plates into a single complete dictionary of metrics

# In[8]:


# group all replicate wells across any and all plates for each measurement
complete = {}
for m in measure_type:
    m = m[:-4]
    temp = []
    for i in range(1, 5):
        for j in range(1, 4):
            temp.append(replicates[(str(i), str(j), m)])
    complete[m] = pd.concat(temp)

# obtain clean group names
def collect_groups():
    df = complete['cell body cluster area']
    groups = [i for n, i in enumerate(list(df.index[~df.index.str.contains('\.[0-9]')])) if i not in list(df.index[~df.index.str.contains('\.[0-9]')])[:n]]
    return groups

groups = collect_groups()
complete.keys()


# # Synchrony of all controls

# In[8]:


for i in [str(i) for i in [1, 2, 3, 4]]:# plate ID, NOT replicate
    for m in ['neurite length', 'cell body cluster area']: #complete.keys()
        rep1 = replicates[(i, '1', m)]
        rep2 = replicates[(i, '2', m)]
        rep3 = replicates[(i, '3', m)]

        rep1 = rep1[rep1.index.value_counts() == 8]
        rep2 = rep2[rep2.index.value_counts() == 8]
        rep3 = rep3[rep3.index.value_counts() == 8]
        
        plot_synchrony([rep1, rep2, rep3], f'Plate {i}, {m}')
        '''if m in ['neurite length', 'neurite branch points']:
            plt.savefig(f'QC_images/Pearson over time controls on plate {i}, {m}.png')
            plt.savefig(f'QC_images/Pearson over time controls on plate {i}, {m}.svg')'''
        plt.show()
        plt.close()


# # Synchrony of individual controls, note that sorting was required

# In[9]:


def plot_synchrony(reps, title=None, ax=None):
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    synchrony12 = []
    for i in range(len(reps[0].columns)):
        synchrony12.append(np.corrcoef(sorted(reps[0].iloc[:, i].values), sorted(reps[1].iloc[:, i].values))[0, 1])
    #print(stats.ttest_ind(reps[0].iloc[:, i], reps[1].iloc[:, i], equal_var=False)[1])
    ax.plot(reps[0].columns, synchrony12, label='replicate 1 with 2')

    synchrony13 = []
    for i in range(len(reps[0].columns)):
        synchrony13.append(np.corrcoef(sorted(reps[0].iloc[:, i].values), sorted(reps[2].iloc[:, i].values))[0, 1])
    ax.plot(reps[0].columns, synchrony13, label='replicate 1 with 3')

    synchrony23 = []
    for i in range(len(reps[0].columns)):
        synchrony23.append(np.corrcoef(sorted(reps[1].iloc[:, i].values), sorted(reps[2].iloc[:, i].values))[0, 1])
    ax.plot(reps[0].columns, synchrony23, label='replicate 2 with 3')

    ax.set_ylim(ymin=-.4, ymax=1)
    ax.set_xlabel('Time')
    ax.set_ylabel('Pearson correlation')
    ax.legend()
    if title is not None:
        ax.set_title(title)
    return ax

for ctrl in ['ATRA (25 uM)']:
    for i in [str(i) for i in [1, 2, 3, 4]]:# plate ID, NOT replicate
        for m in ['neurite length']: #complete.keys()
            rep1 = replicates[(i, '1', m)]
            rep2 = replicates[(i, '2', m)]
            rep3 = replicates[(i, '3', m)]

            rep1 = rep1[rep1.index == ctrl]
            rep2 = rep2[rep2.index == ctrl]
            rep3 = rep3[rep3.index == ctrl]
            plot_synchrony([rep1, rep2, rep3], f'{ctrl} on plate {i}, {m}')
            plt.show()
            plt.close()


# # Mean of controls over time
# ### list of plates and replicates can be modified to plot specific plates as well

# In[10]:


for m in ['neurite length', 'cell body cluster area']: #complete.keys()
    controls ={}
    dft = replicates[('1','1','neurite length')]
    #print(dft)
    idx_count = dft.index.value_counts()
    ctrl_keys = idx_count[idx_count == 8].index.tolist()
    ctrl_keys.remove('DMSO')
    for key in ctrl_keys:
        for i in [str(i) for i in [1, 2, 3, 4]]:
            for j in [str(j) for j in [1, 2, 3]]:
                rep1 = replicates[(i, j, m)]
                rep1 = rep1[rep1.index == key]
                try:
                    controls[key].append(rep1)
                except KeyError:
                    controls[key] = [rep1]
    plot_controls(controls, measure = m, path='High_level_images')
    print('\n\n\n')


# ### As a control, compare HSA/CI of true miR combinations vs those that are combined from the same family

# In[9]:


from math import log10
def plus_sign(list):
    return [x for x in list if '+' in x]

def find_nl(condition):
    return complete['neurite length'].loc[condition].loc[:, 96:126].mean().mean()

def return_max(one, two):
    one_val = complete['neurite length'].loc[one].loc[:, 96:126].mean().mean()
    two_val = complete['neurite length'].loc[two].loc[:, 96:126].mean().mean()
    if max(one_val, two_val) == one_val:
        return one
    elif max(one_val, two_val) == two_val:
        return two
    else:
        print('ERROR')

def find_nl_pvalue(condition):
    replicate_values_combination = complete['neurite length'].loc[condition].loc[:, 96:126].mean(axis=1)
    one, two = condition.split(' + ')
    max = return_max(one, two)
    single_agent_values = complete['neurite length'].loc[max].loc[:, 96:126].mean(axis=1)
    return stats.ttest_ind(replicate_values_combination, single_agent_values, equal_var=False)[1]


combos = plus_sign(groups)

fams_subset = pd.read_csv('mirna_family_info.csv')

def dictionary_brewer(combo_list):
    family_HSA_data = {}
    for combo in combo_list:
        one, two = combo.split(' + ')
        #calculate delta HSA
        one_nl = find_nl(one)
        two_nl = find_nl(two)
        max_nl = max(one_nl, two_nl)
        combo_nl = find_nl(combo)
        #delta_HSA = combo_nl - max_nl
        delta_HSA = max_nl/combo_nl
        # calculate p-value
        p_val = find_nl_pvalue(combo)
        
        # determine if same families or different
        if one == 'hsa-miR-137':
            one = 'hsa-miR-137-3p'
        if two == 'hsa-miR-137':
            two = 'hsa-miR-137-3p'
        try:
            one_fam = fams_subset.loc[fams_subset['MiRBase ID'] == one, 'miR family'].values[0]
            two_fam = fams_subset.loc[fams_subset['MiRBase ID'] == two, 'miR family'].values[0]
        except IndexError:
            print(one, two)

        if one_fam == two_fam:
            family_HSA_data[combo] = [1, delta_HSA, -log10(p_val)]
        else:
            family_HSA_data[combo] = [0, delta_HSA, -log10(p_val)]
    return family_HSA_data

d = dictionary_brewer(combos)

# first separate the key,value pairs by group
d_0 = {k:v for k,v in d.items() if v[0] == 0}
d_1 = {k:v for k,v in d.items() if v[0] == 1}
topleft = {k:v for k,v in d.items() if v[1] < 1 and v[2]>-log10(0.05)}


plt.figure(figsize=(10,10))
# plot the data
plt.scatter([x[1] for x in list(topleft.values())], [x[2] for x in list(topleft.values())], label='p < 0.05 difference between NL of combination and HSA', color='red', zorder=2)
plt.scatter([x[1] for x in list(d_0.values())], [x[2] for x in list(d_0.values())], color='gray', alpha=.5)
plt.scatter([x[1] for x in list(d_1.values())], [x[2] for x in list(d_1.values())], label='Same families', color='blue')
# add a line at y=0.05 with slope=0
plt.axhline(y=-log10(0.05), color='black', ls='--')
plt.axvline(x=1, color='black', ls='--')
for key, value in topleft.items():
    if value[2]>2.5:
        plt.annotate(key, (value[1], value[2]))
plt.xlabel('Naive CI')
plt.ylabel('-log10 p-value')
plt.legend()
plt.savefig('QC_images/neurite length CI volcano plot.png')
plt.savefig('QC_images/neurite length CI volcano plot.svg')
plt.show()
plt.close()
print(f'number of homofamily combos plotted:\n{len(d_1.values())}')

print(f'number of significant NL synergies:\n{len(topleft.keys())}')

'''q = 'hsa-miR-137'
for h in topleft.keys():
    one, two = h.split(' + ')
    if one == q or two == q:
        print(h)'''


# # Holistic analysis
# ### Find agents that are better than ATRA (lower confluence, longer neurites) which also outperform their HSA

# In[10]:


# modifying for narrower time slice

def plus_sign(list):
    return [x for x in list if '+' in x]

def find_cbca(condition):
    return complete['cell body cluster area'].loc[condition].loc[:, 96:126].mean().mean()
def find_nl(condition):
    return complete['neurite length'].loc[condition].loc[:, 96:126].mean().mean()
def find_nbp(condition):
    return complete['neurite branch points'].loc[condition].loc[:, 96:126].mean().mean()


cytostatic = find_cbca('ATRA (25 uM)')
differentiated = find_nl('ATRA (25 uM)') #('miRNA-124')

combos = plus_sign(groups)

superhits = pd.DataFrame(columns=['combo_cbca', 'combo_nl', 'nl distance from HSA', 'one', 'one_nl', 'two', 'two_nl'])


for combo in combos:
    one, two = combo.split(' + ')
    one_nl = find_nl(one)
    two_nl = find_nl(two)

    max_nl = max(one_nl, two_nl)
    combo_nl = find_nl(combo)
    combo_cbca = find_cbca(combo)
    combo_nbp = find_nbp(combo)

    diff_nl = combo_nl-max_nl

    if combo_cbca < cytostatic and diff_nl>0: # and combo_nl > differentiated 
        superhits.loc[combo] = [combo_cbca, combo_nl, diff_nl, one, one_nl, two, two_nl]

print(f'Cytostatic threshold of {cytostatic:.3f}% confluence, differentiated threshold of {differentiated:.3f} neurite length.\n')
print(f'{superhits.shape[0]} hits identified:\n')

superhits = superhits.sort_values(by='nl distance from HSA', ascending=False)
print(superhits)


# # Trim hits based on p values of CBCA rel ATRA and NL rel HSA

# In[11]:


def return_max(one, two):
    one_val = complete['neurite length'].loc[one].loc[:, 96:126].mean().mean()
    two_val = complete['neurite length'].loc[two].loc[:, 96:126].mean().mean()
    if max(one_val, two_val) == one_val:
        return one
    elif max(one_val, two_val) == two_val:
        return two
    else:
        print('ERROR')


def find_nl_pvalue(condition):
    replicate_values_combination = complete['neurite length'].loc[condition].loc[:, 96:126].mean(axis=1)
    one, two = condition.split(' + ')
    max = return_max(one, two)
    single_agent_values = complete['neurite length'].loc[max].loc[:, 96:126].mean(axis=1)
    return stats.ttest_ind(replicate_values_combination, single_agent_values, equal_var=False)[1]


def find_cbca_pvalue_rel_ATRA(condition):
    replicate_values_combination = complete['cell body cluster area'].loc[condition].loc[:, 96:126].mean(axis=1)
    ATRA = complete['cell body cluster area'].loc['ATRA (25 uM)'].loc[:, 96:126].mean(axis=1)
    return stats.ttest_ind(replicate_values_combination, ATRA, equal_var=False)[1]

def find_nl_pvalue_rel_ATRA(condition):
    replicate_values_combination = complete['neurite length'].loc[condition].loc[:, 96:126].mean(axis=1)
    ATRA = complete['neurite length'].loc['ATRA (25 uM)'].loc[:, 96:126].mean(axis=1)
    return stats.ttest_ind(replicate_values_combination, ATRA, equal_var=False)[1]

for i in superhits.index:
    superhits.loc[i, 'nl_pvalue_over_HSA'] = find_nl_pvalue(i)
    superhits.loc[i, 'cbca_pvalue_rel_ATRA'] = find_cbca_pvalue_rel_ATRA(i)
    superhits.loc[i, 'nl_pvalue_rel_ATRA'] = find_nl_pvalue_rel_ATRA(i)


superhits = superhits[superhits['nl_pvalue_over_HSA']<0.05]
superhits = superhits[superhits['cbca_pvalue_rel_ATRA']<0.05]

#superhits.to_csv('superhits.csv')
print(superhits.shape)


# # Generate heatmaps 
# ### (time intensive, only need to run once)

# In[13]:


# create csvs for HSA heatmaps

for m in ['neurite length', 'cell body cluster area']:
    for t in [slice(96,126)]:
        d = HSA_df(m, t, complete)
        d.to_csv(f'HSA_dfs/{initials(m)}_{t}.csv')
        d = abs_heatmap_df(m, t, 'ATRA (25 uM)', complete)
        d.to_csv(f'ABS_dfs/{initials(m)}_{t}.csv')
'''
for m in complete.keys(): #complete.keys
    for t in range(0, 138, 6):
        d = HSA_df(m, t, complete)
        d.to_csv(f'HSA_dfs/{initials(m)}_{t}.csv')
        d = abs_heatmap_df(m, t, 'ATRA (25 uM)', complete)
        d.to_csv(f'ABS_dfs/{initials(m)}_{t}.csv')'''


# In[33]:





# In[48]:


import matplotlib.pyplot as plt
import matplotlib.patches as patches

m = 'neurite length'
t = slice(96,126)


dft = replicates[('1','1','cell body cluster area')]
idx_count = dft.index.value_counts()
ctrl_keys = idx_count[idx_count == 8].index.tolist()
df = complete[m].loc[ctrl_keys, t]
df = df.mean(axis=1)
df = df.groupby(df.index).mean()

print(m)
d1 = pd.read_csv(f'ABS_dfs/{initials(m)}_{t}.csv', index_col=0)
d1.min().min()

plt.figure(figsize=(10,10))
plt.imshow(df.values.reshape(1,-1), vmin=d1.min().min(), vmax=d1.max().max(), cmap='viridis')
plt.yticks([])
plt.xticks(range(len(df)), df.index, rotation=90)
plt.colorbar()
'''
for i, v in enumerate(df.values):
    plt.text(i, 0, '{:.1f}'.format(v), color='white', fontsize=12)
'''
plt.savefig(f'High_level_images/{m} control strip for heatmap addendum.png')
plt.show()


# # Plot heatmaps

# In[50]:


#switch for annotating miR names vs family names
family_index = False

#d1 is bot left, d2 is top right
def trim_hsa(df):
    df.columns = df.columns.map(lambda x: x[4:])
    df.index = df.index.map(lambda x: x[4:])
    return df

hits_to_advance = ['hsa-miR-124-3p + hsa-miR-363-3p',
'hsa-miR-124-3p + hsa-miR-34b-5p',
'hsa-miR-137 + hsa-miR-450b-3p',
'hsa-miR-137 + hsa-miR-449b-5p',
'hsa-miR-137 + hsa-miR-17-5p',
'hsa-miR-19b-3p + hsa-miR-2110']

for m in ['neurite length', 'cell body cluster area']:
    t = slice(96,126)
    plt.figure(figsize=(15,12))
    if initials(m) == 'nl':
        best = 'highest'
        plt.text(11, 51, f'{m} in absolute units', fontsize=18)
        plt.text(46, 39, f'{m} relative to {best} single agent', fontsize=18, rotation=-90)
    elif initials(m) == 'cbca':
        best = 'lowest'
        plt.text(8, 51, f'{m} in absolute units', fontsize=18)
        plt.text(46, 41, f'{m} relative to {best} single agent', fontsize=18, rotation=-90)
    

    d1 = pd.read_csv(f'ABS_dfs/{initials(m)}_{t}.csv', index_col=0)
    d2 = pd.read_csv(f'HSA_dfs/{initials(m)}_{t}.csv', index_col=0)
    d1 = trim_hsa(d1)
    d2 = trim_hsa(d2)

    plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = False
    plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True

    # convert ticks into families if family_index=True
    ordered_fams = []
    for miR in d2.index:
        if miR == 'miR-137':
            miR = 'miR-137-3p'
        if family_index == True:
            miR = fams_subset.loc[fams_subset['MiRBase ID'] == 'hsa-' + miR, 'miR family'].values[0]
        else:
            pass
        ordered_fams.append(miR)
    
    sns.heatmap(d2, xticklabels=ordered_fams, yticklabels=ordered_fams, square=True, linewidths=0, cbar_kws={"shrink": .4, "use_gridspec":False, "location":"right", 'pad':.01}, mask=np.tril(d2), center=0, cmap="RdGy_r")
    sns.heatmap(d1, xticklabels=ordered_fams, yticklabels=ordered_fams, square=True, linewidths=0, cbar_kws={"shrink": .4, "use_gridspec":False, "location":"bottom", "orientation":"horizontal", "pad":.05}, mask=np.triu(d1),cmap="viridis")
    
    for combo in superhits.index:
        miRs = combo.split(' + ')
        x = d1.columns.get_loc(miRs[0][4:])
        y = d1.index.get_loc(miRs[1][4:])
        if combo in hits_to_advance:
            if x < y:
                plt.gca().add_patch(Rectangle((x, y), 1, 1, fill=False, edgecolor='white', lw=3))
                plt.gca().add_patch(Rectangle((y, x), 1, 1, fill=False, edgecolor='black', lw=3))
            else:
                plt.gca().add_patch(Rectangle((x, y), 1, 1, fill=False, edgecolor='black', lw=3))
                plt.gca().add_patch(Rectangle((y, x), 1, 1, fill=False, edgecolor='white', lw=3))
        else:
            if x < y:
                plt.gca().add_patch(Rectangle((x, y), 1, 1, fill=False, edgecolor='white', lw=1))
                plt.gca().add_patch(Rectangle((y, x), 1, 1, fill=False, edgecolor='black', lw=1))
            else:
                plt.gca().add_patch(Rectangle((x, y), 1, 1, fill=False, edgecolor='black', lw=1))
                plt.gca().add_patch(Rectangle((y, x), 1, 1, fill=False, edgecolor='white', lw=1))
    
    title = f'High_level_images/HSA_ABS_heatmap_d4d5_{initials(m)}'
    if family_index == True:
        title = title + '_by_family'
    plt.savefig(f'{title}.png')
    plt.savefig(f'{title}.svg')
    plt.show()
    plt.close()


# In[28]:


def combo_multiplotter(combo_str, bag_of_wells, restrict_plot=None, path=None):
    if restrict_plot != None:
        plotter = restrict_plot
    else:
        plotter = bag_of_wells.keys()
    for measure in plotter:
        original = bag_of_wells[measure]
        df = original.loc[original.index.str.contains('\+')]
        exp = df.loc[combo_str]
        one, two = combo_str.split(' + ')
        one = original.loc[one]
        two = original.loc[two]
        if 'neurite' in measure:
            ctrl = original.loc['miRNA-124']
            ctrl2 = original.loc['ATRA (25 uM)']
            d = exp, one, two, ctrl, ctrl2
        else:
            ctrl = original.loc['miRNA mimic pool (10 nM)']
            ctrl2 = original.loc['ATRA (25 uM)']
            d = exp, one, two, ctrl, ctrl2

        palette = sns.color_palette("husl", len(d))
        plt.figure()
        for e, j in enumerate(d):
            plt.errorbar(j.columns, j.mean(axis=0), yerr=j.std(axis=0), color=palette[e], label=j.index[0])
        plt.xlabel('Time')
        plt.ylabel(f'{measure}')
        plt.title(combo_str)
        plt.legend()
        if path is not None:
            plt.savefig(f'{path}/{combo_str}_{measure}.png')
            plt.savefig(f'{path}/{combo_str}_{measure}.svg')
        plt.show()
        plt.close()

for i in hits_to_advance:
    combo_multiplotter(i, complete, restrict_plot=['neurite length', 'cell body cluster area'], path='HSAhits_p05_cytostatic_plots')


# In[20]:


from math import log10

pvalues = pd.Series()
for combo in combos:
    pvalues.loc[combo] = find_nl_pvalue(combo)
uniform_array = np.random.uniform(low=0.0, high=1.0, size=946)
pvalues = pvalues.apply(lambda x: -log10(x))

neglog = np.vectorize(log10)
uniform_array = neglog(uniform_array)
uniform_array = np.multiply(uniform_array, -1)

plt.plot(np.sort(uniform_array), np.sort(pvalues), 'ro')
plt.plot([0, 4], [0, 4], 'k-', lw=2)

plt.title('Hypothetical -log(p-values) vs. actual -log(p-values) of neurite length VS HSA')
plt.xlabel('Hypothetical -log(p-values) from uniform distribution')
plt.ylabel('Actual -log(p-values)')
plt.savefig('QC_images/qqplot_NL.png')
plt.show()


# In[21]:


from math import log10

pvalues = pd.Series()
for combo in combos:
    pvalues.loc[combo] = find_cbca_pvalue_rel_ATRA(combo)
uniform_array = np.random.uniform(low=0.0, high=1.0, size=946)
pvalues = pvalues.apply(lambda x: -log10(x))

neglog = np.vectorize(log10)
uniform_array = neglog(uniform_array)
uniform_array = np.multiply(uniform_array, -1)

plt.plot(np.sort(uniform_array), np.sort(pvalues), 'ro')
plt.plot([0, 4], [0, 4], 'k-', lw=2)

plt.title('Hypothetical -log(p-values) vs. actual -log(p-values) of CBCA VS ATRA')
plt.xlabel('Hypothetical -log(p-values) from uniform distribution')
plt.ylabel('Actual -log(p-values)')
plt.savefig('QC_images/qqplot_CBCA.png')
plt.show()


# # For recovering specific incycute images

# In[17]:


hits_to_advance = ['hsa-miR-124-3p + hsa-miR-363-3p',
'hsa-miR-124-3p + hsa-miR-34b-5p',
'hsa-miR-137 + hsa-miR-450b-3p',
'hsa-miR-137 + hsa-miR-449b-5p',
'hsa-miR-137 + hsa-miR-17-5p',
'hsa-miR-19b-3p + hsa-miR-2110']

def sample_locator(combination):
    mir1, mir2 = combination.split(' + ')
    coordinates = []
    mir1_well = df2[df2['Sample Group'] == mir1]['Destination well'].values
    mir1_plate = df2[df2['Sample Group'] == mir1]['Destination Plate Barcode'].values
    for k in range(len(mir1_well)):
        coordinates.append([mir1_well[k], mir1_plate[k]])
    
    mir2_well = df2[df2['Sample Group'] == mir2]['Destination well'].values
    mir2_plate = df2[df2['Sample Group'] == mir2]['Destination Plate Barcode'].values

    for j in range(len(mir2_well)):
        if [mir2_well[j], mir2_plate[j]] in coordinates:
            well = mir2_well[j]
            plate = mir2_plate[j]
            return plate, well


ATRAwells = [key for key, value in nonvarwells.items() if value == 'ATRA (25 uM)']
mir124wells = [key for key, value in nonvarwells.items() if value == 'miRNA-124']
print(nonvarwells.items())
#4HPRwells = [key for key, value in nonvarwells.items() if value == '4HPR (2.5 uM)']
ISOTRETINOINwells = [key for key, value in nonvarwells.items() if value == '13-cis RA (25 uM)']
siPLK1wells = [key for key, value in nonvarwells.items() if value == 'siPLK1 (10 nM)']
mimicpoolwells = [key for key, value in nonvarwells.items() if value == 'miRNA mimic pool (10 nM)']


print(f'ATRA wells are :{ATRAwells}')
print(f'miR-124  wells are :{mir124wells}')
#print(f'4HPR wells are :{4HPRwells}')
print(f'ISOTRETINOIN wells are :{ISOTRETINOINwells}')
print(f'siPLK1 wells are :{siPLK1wells}')
print(f'miRNA mimic pool wells are :{mimicpoolwells}')


for i in hits_to_advance:
    mir1, mir2 = i.split(' + ')
    coordinates = []
    mir1_well = df2[df2['Sample Group'] == mir1]['Destination well'].values
    mir1_plate = df2[df2['Sample Group'] == mir1]['Destination Plate Barcode'].values
    for k in range(len(mir1_well)):
        coordinates.append([mir1_well[k], mir1_plate[k]])
    mir2_well = df2[df2['Sample Group'] == mir2]['Destination well'].values
    mir2_plate = df2[df2['Sample Group'] == mir2]['Destination Plate Barcode'].values

    for j in range(len(mir2_well)):
        if [mir2_well[j], mir2_plate[j]] in coordinates:
            print(f'{i} : well {mir2_well[j]}, plate {mir2_plate[j]}')
            mir1_clip = df2[df2['Sample Group'] == mir1].iloc[:, -3:-1]
            print(f'{mir1} : {mir1_clip[mir1_clip.duplicated()].values}')
            mir2_clip = df2[df2['Sample Group'] == mir2].iloc[:, -3:-1]
            print(f'{mir2} : {mir2_clip[mir2_clip.duplicated()].values}')


# In[21]:


# favorite looking wells are 3-F10, 3-H13, 3-E12, E20
def reverse_lookup(plate, well):
    d = df2[df2['Destination Plate Barcode']==plate]
    return d[d['Destination well']==well]
    
reverse_lookup(3, 'E20')


# 
# # END OF SCRIPT

# In[24]:


'''
### for making CSVs of all NL data points
def plus_sign(list):
    return [x for x in list if '+' in x]

combos = plus_sign(groups)


def return_max(one, two):
    one_val = complete['neurite length'].loc[one].loc[:, 72:126].mean().mean()
    two_val = complete['neurite length'].loc[two].loc[:, 72:126].mean().mean()
    if max(one_val, two_val) == one_val:
        return one
    elif max(one_val, two_val) == two_val:
        return two
    else:
        print('ERROR')


def grab_d3_d5_nl(combination, HSA):
    pdrow = pd.concat([complete['neurite length'].loc[combination].loc[:, 72:126].mean(axis=1).reset_index(), complete['neurite length'].loc[HSA].loc[:, 72:126].mean(axis=1).reset_index()])
    return pdrow.reset_index()


allcombos = []


for combo in combos:
    one, two = combo.split(' + ')
    HSA = return_max(one, two)
    #print(grab_d3_d5_nl(combo, HSA).shape)
    allcombos.append(grab_d3_d5_nl(combo, HSA))
    
    
df = pd.DataFrame(columns=['combo_rep_1', 'combo_rep_2', 'combo_rep_3', 'HSA_rep_1', 'HSA_rep_2', 'HSA_rep_3'])

for c in allcombos:
    
    df.loc[c['index'][0]] = c[0].values
#df.to_csv('nl_hsa_scores.csv')'''


# In[25]:


'''
### for making CSVs of all CBCA data points
def plus_sign(list):
    return [x for x in list if '+' in x]

combos = plus_sign(groups)

def add_columns(df, v):
    for i in range(len(v)):
        df[i] = v[i]
    return df

def grab_d3_d5_cbca(combination):
    pdrow = complete['cell body cluster area'].loc[combination].loc[:, 72:126].mean(axis=1)
    return pdrow.reset_index()

atra_cbca_reps = complete['cell body cluster area'].loc['ATRA (25 uM)'].loc[:, 72:126].mean(axis=1).values


allcombos = []


for combo in combos:
    allcombos.append(grab_d3_d5_cbca(combo))


df = pd.DataFrame(columns=['combo_rep_1', 'combo_rep_2', 'combo_rep_3'])

for c in allcombos:
    df.loc[c['index'][0]] = c[0].values

df = add_columns(df, atra_cbca_reps)
#df.to_csv('cbca_scores.csv')'''


# In[26]:


'''# Pearson correlation of different timepoints (probably not useful)
##### Probably better to do something similiar but for only replicates within the same plate rather than looking at all 96 replicates together
import warnings
warnings.simplefilter(action='ignore', category=UserWarning)

def plot_corr_of_timepoints(dictionary, title=None, ax=None):
    for control in dictionary.keys():
        list_of_plates = dictionary[control]
        all_plates = pd.concat(list_of_plates, axis=0)

        fig, ax = plt.subplots(1, 1, figsize=(6, 4))
        sns.heatmap(all_plates.corr(), vmin=-1, vmax=1)
        if title is not None:
            ax.set_title(f'Correlation of {control}' + f', {title}')
        #plt.savefig('QC_images/' + f'Correlation of {control}' + f', {title}.png')
        #plt.savefig('QC_images/' + f'Correlation of {control}' + f', {title}.svg')
        plt.show()
        plt.close()

for m in ['neurite length', 'cell body cluster area']: # complete.keys()
    controls ={}
    dft = replicates['1','1','neurite length']
    idx_count = dft.index.value_counts()
    ctrl_keys = idx_count[idx_count == 8].index.tolist()
    for key in ctrl_keys:
        for i in [str(i) for i in [1, 2, 3, 4]]:
            for j in [str(j) for j in [1, 2, 3]]:
                rep1 = replicates[(i, j, m)]#.drop('Untreate
                rep1 = rep1[rep1.index == key]
                try:
                    controls[key].append(rep1)
                except KeyError:
                    controls[key] = [rep1]
    print(f'Correlation plots for {m}')
    plot_corr_of_timepoints(controls, title = m)
    print('\n\n\n')'''


# # Freeform search

# In[27]:


'''def time_range_to_df_list(start, end, measurement):
    """
    Takes a start time, end time, and measurement from complete.keys() as keyword arguments to return a list of dataframes read from HSA_dfs/. Note that time points are every 6 hours.
    """
    return [pd.read_csv(f'HSA_dfs/{initials(measurement)}_{t}.csv', index_col=0) for t in range(start, end, 6)]


def highest_value(df):
    # First find the highest value cell in dataframe
    highest_value = df.max().max()
    # Then find the coordinates of the highest value cell
    coordinates = np.where(df == highest_value)
    # Then find the index and column names of the highest value cell
    index = df.index[coordinates[0][0]]
    col = df.columns[coordinates[1][0]]
    return ' + '.join(sorted([index, col]))

def lowest_value(df):
    # First find the highest value cell in dataframe
    lowest_value = df.min().min()
    # Then find the coordinates of the highest value cell
    coordinates = np.where(df == lowest_value)
    # Then find the index and column names of the highest value cell
    index = df.index[coordinates[0][0]]
    col = df.columns[coordinates[1][0]]
    return ' + '.join(sorted([index, col]))

            
window = 4
t = 108
hit_thresh = 1

measures = ['neurite length', 'neurite branch points']
dfs = []
for m in measures:
    dfs = dfs + (time_range_to_df_list(t, t+(window*6), m))

superhits = find_higher_than_one(dfs, higher_than=hit_thresh)
print(f'Hits count:{len(superhits)}\n{superhits}')

print('\nTop performers:')
for j in range(len(measures)):
    for i in range(window):
        unit = (j*window)+ i
        print(f'{highest_value(dfs[unit])} : {measures[j]} @ t={t+(i*6)}')
        l = find_pos([dfs[unit]], threshold=.5)
        print(l)
        highlight_hits(dfs[unit], measures[j], t+(i*6), l)'''


# In[ ]:




