# -*- coding: utf-8 -*-
"""
Created on Wed Jan 19 15:12:19 2022

@author: xmste
"""

import scipy.io as sio
import seaborn as sns
import pandas as pd

fControl = sio.loadmat('samples.mat')
fDCM = sio.loadmat('samplesDCMFemale.mat')

groups = [fControl, fDCM]
group_names = ['fControl', 'fDCM']

fluxes = []
reactions = []
disease_groups = []

for group, group_name in zip(groups, group_names):
    for idx, sample in enumerate(group['samples'][0:3]):
        for flux in sample:
            fluxes.append(flux)
            reactions.append(idx)
            disease_groups.append(group_name)
            
df = pd.DataFrame({'reaction': reactions, 'fluxes': fluxes, 'group': disease_groups})

g = sns.catplot(x = 'reaction', y = 'fluxes', hue = 'group', data = df)