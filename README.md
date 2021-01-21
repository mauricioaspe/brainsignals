# brainsignals
Analysis of LFPs spectral and spiking activity, spiking/LFPs coupling, and other neurophysiological data.

## Coherence analysis
#### Usage
```
python3 4_iCoh.py
```

Filters for different frequency bands in `filter_parameters = {'theta': {'N': 4, 'lowcut': 3.5,  'highcut': 7.5}}` and obtains imaginary part of coherence (32 channels x 32 channels x n events matrix)


## Spectrograms
#### Usage
```
python3 3_spectrograms.py
```
Get a list of IDs and obtains the spectrograms
```
IDs = ['ID1597', 'ID1659'] #, 'ID1678', 'ID1908', 'ID1984', 'ID1985', 'ID2014', 'ID1668', 'ID1665', 'ID2018', 'ID2024', 'ID2013']
spectra.transform(IDs, n_freq=80, substract_baseline=False, epoch_length=4, final_fs=1000.0, save_data=False)
```

## Epoching
##### Usage
```
python3 2_epoching.py
```
Get a list of IDs and make the epochs according to 
```
seconds_pre = 2
seconds_post = 2
```
