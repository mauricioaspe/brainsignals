#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 21:42:23 2020

@author: diogenes
"""

#%%
import os
import pickle
import pyjags
import pandas as pd

import statsmodels.api as sm
import statsmodels.formula.api as smf

import matplotlib.pyplot as plt


#%%
path_to_files = '/home/diogenes/Projects/brainsignals/results/info_files/'
path_to_figs = '/home/diogenes/Projects/brainsignals/figs/'

IDs = pickle.load(open('/home/diogenes/Projects/brainsignals/IDs.dict', 'rb'))
WTs = IDs['list_WT']
KOs = IDs['list_KO']


mice_list = sorted(os.listdir(path_to_files))

all_WTs = pd.DataFrame()
all_KOs = pd.DataFrame()
for mouse in mice_list:
    temp = pickle.load(open(path_to_files + mouse, 'rb'))
    temp['path']['ID'] = temp['ID']
    
    if mouse[:4] in WTs:
        temp['path']['condition'] = 'WT'
        all_WTs = all_WTs.append(temp['path'], ignore_index=True)
        color = 'b'
    else:
        temp['path']['condition'] = 'KO'
        all_KOs = all_KOs.append(temp['path'], ignore_index=True)
        color = 'r'
        
        
    plt.plot(temp['path']['entrances'], temp['path']['dwell'], 'o', color=color)
    plt.xlim([0, 600])
    plt.xlabel('time (s)', fontsize=14)
    plt.ylabel('dwell (s)', fontsize=14)
    plt.title(temp['ID'], fontsize=18)
    plt.savefig(path_to_figs + temp['ID'] + '.png', dpi=100)
    plt.close()        
        
    #print(temp['path'])

all_WTs = all_WTs.dropna()   
all_KOs = all_KOs.dropna()   

all_mice = all_WTs.append(all_KOs, ignore_index=True)

#print(all_mice)
#print("Shape {}".format(all_mice.shape))    


#%%

plt.hist(all_WTs['dwell'], bins=20, color='b', alpha=0.7)
plt.hist(all_KOs['dwell'], bins=20, color='r', alpha=0.7)
plt.xlabel('time (s)', fontsize=14)

plt.title('Dwell times', fontsize=16)

plt.savefig(path_to_figs + 'dewells.png', dpi=100)
plt.close()

#%%
plt.plot(all_WTs['entrances'], all_WTs['dwell'], 'ob', alpha=0.7)
plt.plot(all_KOs['entrances'], all_KOs['dwell'], 'or', alpha=0.7)
plt.xlabel('time (s)', fontsize=14)
plt.ylabel('dwell (s)', fontsize=14)

plt.title('Dwell versus entrance times', fontsize=16)

plt.savefig(path_to_figs + 'dwell_vs_time.png', dpi=100)
plt.close()



#%% 
#######################
### Regression analysis

models = {
    'model1': 'dwell ~ condition + entrances'
}

results = dict()
for model in models.keys():
    md = smf.mixedlm(models[model], all_mice, groups=all_mice['ID'])
    mdf = md.fit()
    results[model] = mdf


# In[7]:


print(results['model1'].summary()) # 'model1': 'na2 ~ H + na1 + na1H'







#%%
#####################
### Bayesian analysis


code = '''

model {

    # Priors
    beta0_mu ~ dnorm(0, epsilon)   
    beta1_mu ~ dnorm(0, epsilon) 
    beta2_mu ~ dnorm(0, epsilon)   
    beta3_mu ~ dnorm(0, epsilon)
   
    beta0_sigma ~ dunif(epsilon, 10)
    beta1_sigma ~ dunif(epsilon, 10)
    beta2_sigma ~ dunif(epsilon, 10)
    beta3_sigma ~ dunif(epsilon, 10)
    
    epsilon <- 1e-5

    nu <- nuMinusOne+1
    nuMinusOne ~ dexp(1/29)
    
    
    # Loop for subjects
    for (s in 1:n_subjects) {
    
       beta0[s] ~ dt(beta0_mu, 1/beta0_sigma^2, nu)
       beta1[s] ~ dt(beta1_mu, 1/beta1_sigma^2, nu)
       beta2[s] ~ dt(beta2_mu, 1/beta2_sigma^2, nu) 
       beta3[s] ~ dt(beta3_mu, 1/beta3_sigma^2, nu)
       
    }
    
    
    # Loop for all trials
    for (t in 1:n_trials) {
    
        dwell[t] ~ beta0[S[t]] + beta1[S[t]]*condition[t] + beta2[S[t]]*entrances[t] + beta3[S[t]]*entrances[t]*condition[t]    
    
    }

}

'''

#%%

condition = all_mice['condition'].values
entrances = all_mice['entrances'].values
dwell = all_mice['dwell'].values

n_trials = dwell.shape[0]

subjects   = pd.factorize(all_mice['ID'].index)[0] + 1
n_subjects = pd.unique(subjects).shape[0]


chains = 5
adapt = 5000
iterations = 5000

init_values = None #dict(theta_COM=0.5, theta_HUM=0.5, delta=0) #dict(theta=0.0001, 0.25, 0.75, 1.0])

#print('Initialising model for {}...'.format(group))
model = pyjags.Model(code=code,
                     data=dict(dwell=dwell, entrances=entrances,
                               condition=condition, n_trials=n_trials,
                               n_subjects=n_subjects, S=subjects),
                     init=init_values, chains=chains, adapt=adapt, generate_data=True,
                     progress_bar=False, refresh_seconds=None,
                     threads=1, chains_per_thread=1)


print('Sampling...')
variables = ['beta0', 'beta1', 'beta2', 'beta3']
samples = model.sample(iterations=iterations, vars=variables, thin=1,
                       monitor_type='trace')
