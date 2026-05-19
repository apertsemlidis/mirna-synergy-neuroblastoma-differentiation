import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import json
from matplotlib.patches import Rectangle
import scipy.stats as stats
import warnings
warnings.simplefilter(action='ignore', category=UserWarning)


# for making heatmaps with the miRs clustered by family/sequence/target-space
# to to this I ran an sns.clustermap of mirDB target space of the 44 miRs and grabbed the reordered dendrogram
unordered = ['hsa-miR-103a-3p', 'hsa-miR-106a-5p', 'hsa-miR-106b-5p', 'hsa-miR-107',
       'hsa-miR-10a-5p', 'hsa-miR-10b-5p', 'hsa-miR-124-3p', 'hsa-miR-128-3p',
       'hsa-miR-135b-5p', 'hsa-miR-137', 'hsa-miR-143-5p', 'hsa-miR-17-5p',
       'hsa-miR-18a-5p', 'hsa-miR-18b-5p', 'hsa-miR-193b-3p',
       'hsa-miR-196a-5p', 'hsa-miR-19a-3p', 'hsa-miR-19b-3p',
       'hsa-miR-200a-3p', 'hsa-miR-200b-3p', 'hsa-miR-200c-3p',
       'hsa-miR-204-5p', 'hsa-miR-20a-5p', 'hsa-miR-20b-5p', 'hsa-miR-211-5p',
       'hsa-miR-2110', 'hsa-miR-25-3p', 'hsa-miR-27a-3p', 'hsa-miR-27b-3p',
       'hsa-miR-340-5p', 'hsa-miR-34a-5p', 'hsa-miR-34b-5p', 'hsa-miR-363-3p',
       'hsa-miR-3714', 'hsa-miR-3937', 'hsa-miR-429', 'hsa-miR-449a',
       'hsa-miR-449b-5p', 'hsa-miR-450b-3p', 'hsa-miR-452-5p',
       'hsa-miR-506-3p', 'hsa-miR-873-5p', 'hsa-miR-92a-3p', 'hsa-miR-93-5p']

order = [29,6, 40, 16, 17, 2, 22, 1, 11, 23, 43, 27, 28, 9, 35, 19, 20, 42, 26, 32, 7, 33, 21, 24, 18, 8, 37, 30, 36, 0, 3, 39, 25, 31, 10, 12, 13, 41, 4, 5, 15, 14, 34, 38]
# order the unordered list, according to order
ordered = [unordered[i] for i in order]

#for mapping the controls which were done manually rather than by acoustic liquid handler
with open('nonvarwells.txt') as json_file:
    nonvarwells = json.load(json_file)
with open('plate4nonvarwells.txt') as json_file:
    plate_4_fixed = json.load(json_file)

    
def mapper(df, plate, instructions):
    clip = instructions.loc[instructions['Destination Plate Barcode'] == int(plate)]
    cols = {}
    for column in df.columns:
        try:
            mirs = clip.loc[clip['Destination well'] == column]['Sample Group']
            if mirs.empty==False:
                cols[column] = ' + '.join(sorted(list(set(mirs.values))))
        except ValueError:
            print(f'we had an error with {column}')
    df = df.rename(columns=cols)
    df = df.rename(columns=nonvarwells)
    if int(plate) == 4:
        df = df.rename(columns=plate_4_fixed)
    return df, cols

def plot_replicate_violins(replicates, sample, time, keyword='neurite length'):
    measurements = {}
    plot_labels = {}
    for key in replicates.keys():
        if key[2] not in measurements.keys():
            measurements[key[2]] = []
            plot_labels[key[2]] = []
        measurements[key[2]].append(replicates[key].loc[sample, time])
        plot_labels[key[2]].append(str(key[0]) + '-' + str(key[1]))
    for measurement in [i for i in list(measurements.keys()) if keyword in i]:
        plt.figure()
        plt.title(f'{sample} {measurement} at {time}hrs')
        sns.violinplot(data=measurements[measurement], showmeans=True, showmedians=True)
        plt.xticks(range(0, len(plot_labels[measurement])), plot_labels[measurement], rotation=90)
        plt.show()

def plot_replicate_swarm(replicates, sample, time, kw='neurite length', savefig=False, title=None):
    measurements = {}
    plot_labels = {}
    for key in replicates.keys():
        if key[2] not in measurements.keys():
            measurements[key[2]] = []
            plot_labels[key[2]] = []
        try:
            measurements[key[2]].append(replicates[key].loc[sample, time])
            plot_labels[key[2]].append(str(key[0]) + '-' + str(key[1]))
        except KeyError:
            print(f'Key error at {key}')
        
    for measurement in [i for i in list(measurements.keys()) if kw in i]:
        plt.figure()
        plt.title(f'{sample} {measurement} at {time}hrs')
        sns.boxplot(data=measurements[measurement])
        sns.swarmplot(data=measurements[measurement], color=".2")
        plt.xticks(range(0, len(plot_labels[measurement])), plot_labels[measurement], rotation=90)
        if savefig == True:
            plt.savefig('QC_images/' + title + f' {sample} {measurement} at {time}hrs.svg')
            plt.savefig('QC_images/' + title + f' {sample} {measurement} at {time}hrs.png')
        plt.show()
        plt.close()


            
def cellmetric_normalizePI(neg, pos, time, source):
    for plate in source.keys():
        if plate[2] == 'cell body cluster area':
            p = source[plate].loc[pos].median()[time]
            n = source[plate].loc[neg].median()[time]
            source[plate] = source[plate].apply(lambda x: (p-x)/(p-n))
    return source

def neuritemetric_normalizePA(neg, pos, time, source):
    for plate in source.keys():
        if plate[2] not in ['cell body cluster area', 'cell body clusters']:
            p = source[plate].loc[pos].median()[time]
            n = source[plate].loc[neg].median()[time]
            source[plate] = source[plate].apply(lambda x: (x-n)/(p-n))
    return source

def plot_synchrony(reps, title=None, ax=None):
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    synchrony12 = []
    for i in range(len(reps[0].columns)):
        synchrony12.append(np.corrcoef(reps[0].iloc[:, i], reps[1].iloc[:, i])[0, 1])
    ax.plot(reps[0].columns, synchrony12, label='replicate 1 with 2')

    synchrony13 = []
    for i in range(len(reps[0].columns)):
        synchrony13.append(np.corrcoef(reps[0].iloc[:, i], reps[2].iloc[:, i])[0, 1])
    ax.plot(reps[0].columns, synchrony13, label='replicate 1 with 3')

    synchrony23 = []
    for i in range(len(reps[0].columns)):
        synchrony23.append(np.corrcoef(reps[1].iloc[:, i], reps[2].iloc[:, i])[0, 1])
    ax.plot(reps[0].columns, synchrony23, label='replicate 2 with 3')

    ax.set_ylim(ymin=-.4, ymax=1)
    ax.set_xlabel('Time')
    ax.set_ylabel('Pearson correlation')
    ax.legend()
    if title is not None:
        ax.set_title(title)
    return ax

# avg of controls over time

def plot_controls(dictionary, measure=None, path=None, ax=None):
    for control in dictionary.keys():
        list_of_plates = dictionary[control]
        all_plates = pd.concat(list_of_plates, axis=0)
        dictionary[control] = all_plates
    if ax is None:
        fig, ax = plt.subplots()
    for key, value in dictionary.items():
        mean = value.mean(axis=0)
        std = value.std(axis=0)

        #ax.errorbar(x=mean.index, y=mean, yerr=std, label=key)
        ax.plot(mean, label=key)
        ax.fill_between(mean.index, mean - std, mean + std, alpha=0.1)
    ax.legend()
    title = f'Mean of controls, {measure}'
    ax.set_title(title)
    ax.set_xlabel('Time')
    ax.set_ylabel(title)
    if path != None:
        plt.savefig(f'{path}/' + title + '.png')
        plt.savefig(f'{path}/' + title + '.svg')
    plt.show()
    plt.close()

            
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
        

# make heatmap for each measurement, note that two point normalization has already occured by this point in script and we simply center the measurements at 0

def abs_heatmap_df(measure, time, ctrl, bag_of_wells):
    S = bag_of_wells[measure].loc[:, time]    
   
    combos = S.index[S.index.str.contains('\+')]
    singlets = []
    for i in combos:
        one, two = i.split(' + ')
        if one not in singlets:
            singlets.append(one)
        if two not in singlets:
            singlets.append(two)

    df = pd.DataFrame(index = ordered, columns=ordered)
    # create diagonal of monoagents
    for single in singlets:
        if type(time) == int:
            df.loc[single, single] = S.loc[single].mean()
        else:
            df.loc[single, single] = S.loc[single].mean().mean()
    # create both triangles of combinations
    for comb in combos:
        one, two = comb.split(' + ')
        if type(time) == int:
            df.loc[one, two] = S.loc[comb].mean()
            df.loc[two, one] = S.loc[comb].mean()
        else:
            df.loc[one, two] = S.loc[comb].mean().mean()
            df.loc[two, one] = S.loc[comb].mean().mean()
    df = df.astype(np.float32)
    return df


def abs_heatmapper(df, measure, time, ctrl, path=None):
    plt.figure(figsize=(12,10))
    plt.title(f'{measure} at {time} hours\ncombinations relative to {ctrl}', fontsize=18)
    sns.heatmap(df, cmap='RdBu_r', center=0, square=True, linewidths=.5, cbar_kws={"shrink": .5})
    if path != None:
        plt.savefig(rf'{path}/heatmap_rel_to_{ctrl}_{measure}_at_{time}.png')
        plt.savefig(rf'{path}/heatmap_rel_to_{ctrl}_{measure}_at_{time}.svg')
    plt.show()
    plt.close()

# make heatmap for each measurement,now normalized per highest single agent
def HSA_df(measure, time, bag_of_wells):
    S = bag_of_wells[measure].loc[:, time]
    combos = S.index[S.index.str.contains('\+')]
    singlets = []
    for i in combos:
        one, two = i.split(' + ')
        if one not in singlets:
            singlets.append(one)
        if two not in singlets:
            singlets.append(two)

    df = pd.DataFrame(index = ordered, columns=ordered)
    for single in singlets:
        df.loc[single, single] = 0
    #normalize per highest single agent if measure ! cell body cluster area or cell body clusters
    if measure not in ['cell body cluster area', 'cell body clusters']:
        for comb in combos:
            one, two = comb.split(' + ')
            if type(time) == int:
                df.loc[one, two] = S.loc[comb].mean()-max(S.loc[one].mean(), S.loc[two].mean())
                df.loc[two, one] = S.loc[comb].mean()-max(S.loc[one].mean(), S.loc[two].mean())
            else:
                df.loc[one, two] = S.loc[comb].mean().mean()-max(S.loc[one].mean().mean(), S.loc[two].mean().mean())
                df.loc[two, one] = S.loc[comb].mean().mean()-max(S.loc[one].mean().mean(), S.loc[two].mean().mean())
        df = df.astype(np.float32)
    else:
        for comb in combos:
            one, two = comb.split(' + ')
            if type(time) == int:
                df.loc[one, two] = S.loc[comb].mean()-min(S.loc[one].mean(), S.loc[two].mean())
                df.loc[two, one] = S.loc[comb].mean()-min(S.loc[one].mean(), S.loc[two].mean())
            else:
                df.loc[one, two] = S.loc[comb].mean().mean()-min(S.loc[one].mean().mean(), S.loc[two].mean().mean())
                df.loc[two, one] = S.loc[comb].mean().mean()-min(S.loc[one].mean().mean(), S.loc[two].mean().mean())
        df = df.astype(np.float32)
    return df
'''
def HSA_df(measure, time, bag_of_wells):
    S = bag_of_wells[measure].loc[:, time]
    combos = S.index[S.index.str.contains('\+')]
    singlets = []
    for i in combos:
        one, two = i.split(' + ')
        if one not in singlets:
            singlets.append(one)
        if two not in singlets:
            singlets.append(two)

    df = pd.DataFrame(index = ordered, columns=ordered)
    for single in singlets:
        df.loc[single, single] = 1
    #normalize per highest single agent if measure ! cell body cluster area or cell body clusters
    if measure not in ['cell body cluster area', 'cell body clusters']:
        for comb in combos:
            one, two = comb.split(' + ')
            if type(time) == int:
                df.loc[one, two] = max(S.loc[one].mean(), S.loc[two].mean())/S.loc[comb].mean()
                df.loc[two, one] = max(S.loc[one].mean(), S.loc[two].mean())/S.loc[comb].mean()
            else:
                df.loc[one, two] = max(S.loc[one].mean().mean(), S.loc[two].mean().mean())/S.loc[comb].mean().mean()
                df.loc[two, one] = max(S.loc[one].mean().mean(), S.loc[two].mean().mean())/S.loc[comb].mean().mean()
        df = df.astype(np.float32)
    else:
        for comb in combos:
            one, two = comb.split(' + ')
            if type(time) == int:
                df.loc[one, two] = S.loc[comb].mean()/min(S.loc[one].mean(), S.loc[two].mean())
                df.loc[two, one] = S.loc[comb].mean()/min(S.loc[one].mean(), S.loc[two].mean())
            else:
                df.loc[one, two] = S.loc[comb].mean().mean()/min(S.loc[one].mean().mean(), S.loc[two].mean().mean())
                df.loc[two, one] = S.loc[comb].mean().mean()/min(S.loc[one].mean().mean(), S.loc[two].mean().mean())
        df = df.astype(np.float32)
    return df
'''


def HSA_heatmap_plotter(df, measure, time, path=None):
    plt.figure(figsize=(12,10))
    plt.title(f'{measure} at {time} hours\ncombinations relative to highest single agent of the combo', fontsize=18)
    sns.heatmap(df, cmap='RdBu_r', center=0, square=True, linewidths=.5, cbar_kws={"shrink": .5})  
    if path != None:
        plt.savefig(rf'{path}/HSAnormed{measure}_at_{time}.png')
        plt.savefig(rf'{path}/HSAnormed{measure}_at_{time}.svg')

    plt.show()
    plt.close()

# for cleaning the combinations up, so we don't return both 'A + B' and 'B + A'
def remove_duplicates(list_with_duplicates):
    list_without_duplicates = []
    # strip ' + ' from the strings in the list
    for i in list_with_duplicates:
        i = i.replace(' + ', ' ')
        list_without_duplicates.append(i)
    # split the strings into lists
    for i in range(len(list_without_duplicates)):
        list_without_duplicates[i] = list_without_duplicates[i].split()
    # sort the lists
    for i in range(len(list_without_duplicates)):
        list_without_duplicates[i].sort()
    # join the lists back into strings
    for i in range(len(list_without_duplicates)):
        list_without_duplicates[i] = ' '.join(list_without_duplicates[i])
    # add ' + ' back to the strings
    for i in range(len(list_without_duplicates)):
        list_without_duplicates[i] = list_without_duplicates[i].replace(' ', ' + ')
    # remove duplicates
    list_without_duplicates = list(set(list_without_duplicates))
    return list_without_duplicates

# Identify combinations based on HSA-heatmaps. This function takes a list of heatmaps K and returns a list of tuples corresponding to the index,col coordinate names (a miR combo) which remain positive throughout all members of the heatmap list. Provided a list of HSA normalized heatmaps for end time point measurements, it will return combinations that exceeded their HSA (i.e. HSA set to 0, >0 exceeding HSA)
def find_pos(K, threshold=0):
    best = []
    for row in K[0].index:
        for col in K[0].columns:
            if all(K[i].loc[row, col] > threshold for i in range(len(K))):
                best.append((row + ' + ' + col))
    return remove_duplicates(best)

def find_higher_than_one(K, higher_than=1):
    best = []
    for row in K[0].index:
        for col in K[0].columns:
            if all(K[i].loc[row, col] > higher_than for i in range(len(K))):
                best.append((row + ' + ' + col))
    return remove_duplicates(best)

def find_neg(K, thresh=0):
    worst = []
    for row in K[0].index:
        for col in K[0].columns:
            if all(K[i].loc[row, col] < thresh for i in range(len(K))):
                worst.append((row + ' + ' + col))
    return remove_duplicates(worst)


def initials(string):
    words = string.split()
    result = ""
    for word in words:
        result += word[0]
    return result

def HSA_heatmap_lister(start, end, query):
    k = []
    if isinstance(query, list):
        for m in query:
            for t in range(start, end, 6):
                d=pd.read_csv(f'HSA_dfs/{initials(m)}_{t}.csv', index_col=0)
                k.append(d)
    elif isinstance(query, str):
        for m in [i for i in list(complete.keys()) if query in i]:
            for t in range(start, end, 6):
                d=pd.read_csv(f'HSA_dfs/{initials(m)}_{t}.csv', index_col=0)
                k.append(d)
    else:
        print('Error, query must be string contained in complete.keys() or list like')
    return k


def highlight_hits(d, m, t, h, path=None):
    plt.figure(figsize=(12,10))
    plt.title(f'{m} at {t} hours\ncombinations relative to highest single agent of the combo', fontsize=18)
    sns.heatmap(d, cmap='RdBu_r', center=0, square=True, linewidths=.5, cbar_kws={"shrink": .5}) 
    for combo in h:
        miRs = combo.split(' + ')
        x = d.columns.get_loc(miRs[0])
        y = d.index.get_loc(miRs[1])
        plt.gca().add_patch(Rectangle((x, y), 1, 1, fill=False, edgecolor='black', lw=3))
        plt.gca().add_patch(Rectangle((y, x), 1, 1, fill=False, edgecolor='black', lw=3))

    if path != None:
        plt.savefig(rf'{path}/Superhits highlighted HSA heatmap {m}_at_{t}.png')
        plt.savefig(rf'{path}/Superhits highlighted HSA heatmap {m}_at_{t}.svg')
    plt.show()
    plt.close()
