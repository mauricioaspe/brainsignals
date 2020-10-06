#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 14 15:39:09 2019

@author: mauro
"""

import os
import numpy as np

# constants for pre-allocating matrices:
MAX_NUMBER_OF_SPIKES = int(1e6)
MAX_NUMBER_OF_EVENTS = 1e6


def loadContinuous(filepath, dtype=float, verbose=True, 
    start_record=None, stop_record=None, ignore_last_record=True):
    """Load continuous data from a single channel in the file `filepath`.
    
    This is intended to be mostly compatible with the previous version.
    The differences are:
    - Ability to specify start and stop records
    - Converts numeric data in the header from string to numeric data types
    - Does not rely on a predefined maximum data size
    - Does not necessarily drop the last record, which is usually incomplete
    - Uses the block length that is specified in the header, instead of
        hardcoding it.
    - Returns timestamps and recordNumbers as int instead of float
    - Tests the record metadata (N and record marker) for internal consistency
    The OpenEphys file format breaks the data stream into "records", 
    typically of length 1024 samples. There is only one timestamp per record.
    Args:
        filepath : string, path to file to load
        dtype : float or np.int16
            If float, then the data will be multiplied by bitVolts to convert
            to microvolts. This increases the memory required by 4 times.
        verbose : whether to print debugging messages
        start_record, stop_record : indices that control how much data
            is read and returned. Pythonic indexing is used,
            so `stop_record` is not inclusive. If `start` is None, reading
            begins at the beginning; if `stop` is None, reading continues
            until the end.
        ignore_last_record : The last record in the file is almost always
            incomplete (padded with zeros). By default it is ignored, for
            compatibility with the old version of this function.
    Returns: dict, with following keys
        data : array of samples of data
        header : the header info, as returned by readHeader
        timestamps : the timestamps of each record of data that was read
        recordingNumber : the recording number of each record of data that
            was read. The length is the same as `timestamps`.
    """
    if dtype not in [float, np.int16]:
        raise ValueError("Invalid data type. Must be float or np.int16")

    if verbose:
        print("Loading continuous data from " + filepath)

    """Here is the OpenEphys file format:
    'each record contains one 64-bit timestamp, one 16-bit sample 
    count (N), 1 uint16 recordingNumber, N 16-bit samples, and 
    one 10-byte record marker (0 1 2 3 4 5 6 7 8 255)'
    Thus each record has size 2*N + 22 bytes.
    """
    # This is what the record marker should look like
    spec_record_marker = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 255])

    # Lists for data that's read
    timestamps = []
    recordingNumbers = []
    samples = []
    samples_read = 0
    records_read = 0
    
    # Open the file
    with open(filepath, 'rb') as f:
        # Read header info, file length, and number of records
        header = readHeader(f)
        record_length_bytes = 2 * header['blockLength'] + 22
        fileLength = os.fstat(f.fileno()).st_size
        n_records = get_number_of_records(filepath)
        
        # Use this to set start and stop records if not specified
        if start_record is None:
            start_record = 0
        if stop_record is None:
            stop_record = n_records
        
        # We'll stop reading after this many records are read
        n_records_to_read = stop_record - start_record
        
        # Seek to the start location, relative to the current position
        # right after the header.
        f.seek(record_length_bytes * start_record, 1)
        
        # Keep reading till the file is finished
        while f.tell() < fileLength and records_read < n_records_to_read:
            # Skip the last record if requested, which usually contains
            # incomplete data
            if ignore_last_record and f.tell() == (
                fileLength - record_length_bytes):
                break
            
            # Read the timestamp for this record
            # litte-endian 64-bit signed integer
            timestamps.append(np.fromfile(f, np.dtype('<i8'), 1))
        
            # Read the number of samples in this record
            # little-endian 16-bit unsigned integer
            N = np.fromfile(f, np.dtype('<u2'), 1).item() 
            if N != header['blockLength']:
                raise IOError('Found corrupted record in block ' + 
                    str(recordNumber))
            
            # Read and store the recording numbers
            # big-endian 16-bit unsigned integer
            recordingNumbers.append(np.fromfile(f, np.dtype('>u2'), 1))
            
            # Read the data
            # big-endian 16-bit signed integer
            data = np.fromfile(f, np.dtype('>i2'), N)
            if len(data) != N:
                raise IOError("could not load the right number of samples")
            
            # Optionally convert dtype
            if dtype == float: 
                data = data * header['bitVolts']
                        
            # Store the data
            samples.append(data)

            # Extract and test the record marker
            record_marker = np.fromfile(f, np.dtype('<u1'), 10)
            if np.any(record_marker != spec_record_marker):
                raise IOError("corrupted record marker at record %d" %
                    records_read)
            
            # Update the count
            samples_read += len(samples)            
            records_read += 1

    # Concatenate results, or empty arrays if no data read (which happens
    # if start_sample is after the end of the data stream)
    res = {'header': header}
    if samples_read > 0:
        res['timestamps'] = np.concatenate(timestamps)
        res['data'] = np.concatenate(samples)
        res['recordingNumber'] = np.concatenate(recordingNumbers)
    else:
        res['timestamps'] = np.array([], dtype=np.int)
        res['data'] = np.array([], dtype=dtype)
        res['recordingNumber'] = np.array([], dtype=np.int)
    return res


def loadSpikes(filepath):
    
    data = { }
    
    print('loading spikes...')
    
    f = open(filepath,'rb')
    header = readHeader(f)
    
    # if float(header[' version']) < 0.4:
    #    raise Exception('Loader is only compatible with .spikes files with version 0.4 or higher')
     
    data['header'] = header 
    numChannels = int(header['num_channels'])
    numSamples = 40 # **NOT CURRENTLY WRITTEN TO HEADER**
    
    spikes = np.zeros((MAX_NUMBER_OF_SPIKES, numSamples, numChannels))
    timestamps = np.zeros(MAX_NUMBER_OF_SPIKES)
    source = np.zeros(MAX_NUMBER_OF_SPIKES)
    gain = np.zeros((MAX_NUMBER_OF_SPIKES, numChannels))
    thresh = np.zeros((MAX_NUMBER_OF_SPIKES, numChannels))
    sortedId = np.zeros((MAX_NUMBER_OF_SPIKES, numChannels))
    recNum = np.zeros(MAX_NUMBER_OF_SPIKES)
    
    currentSpike = 0
    
    # The method tell() returns the current position of the file read/write pointer within the file.

    # The method fstat() returns information about a file associated with the fd.
    # Here is the structure returned by fstat method
    # st_size − total size, in bytes

    # Signature: os.fstat(fd)
    # Docstring:
    # Perform a stat system call on the given file descriptor.

    # Like stat(), but for an open file descriptor.
    # Equivalent to os.stat(fd).
    # Type:      builtin_function_or_method
    
    while f.tell() < os.fstat(f.fileno()).st_size:
        
    # 'description': "'Each record contains 1 uint8 eventType, 1 int64 timestamp,
    # 1 int64 software timestamp, 1 uint16 sourceID, 1 uint16 numChannels (n), 1 uint16 numSamples (m),
    # 1 uint16 sortedID, 1 uint16 electrodeID, 1 uint16 channel, 3 uint8 color codes,
    # 2 float32 component projections, n*m uint16 samples, n float32 channelGains, n uint16 thresholds,
    # and 1 uint16 recordingNumber'",
    #  'date_created': "'24-Apr-2018 165627'",
    #  'electrode': "'SE0'",
    #  'num_channels': '1',
    #  'sampleRate': 30000.0},
        
        eventType = np.fromfile(f, np.dtype('<u1'), 1) # 1 int8 eventType #always equal to 4, discard
        timestamps[currentSpike] = np.fromfile(f, np.dtype('<i8'), 1) # 1 int64 timestamp
        software_timestamp = np.fromfile(f, np.dtype('<i8'), 1) # 1 int64 software timestamp
        source[currentSpike] = np.fromfile(f, np.dtype('<u2'), 1) # 1 uint16 sourceID
        numChannels = np.fromfile(f, np.dtype('<u2'), 1)[0] # 1 uint16 numChannels (n)
        numSamples = np.fromfile(f, np.dtype('<u2'), 1)[0] # 1 uint16 numSamples (m)
        sortedId[currentSpike] = np.fromfile(f, np.dtype('<u2'),1) # 1 uint16 sortedID
        electrodeId = np.fromfile(f, np.dtype('<u2'),1) # 1 uint16 electrodeID
        channel = np.fromfile(f, np.dtype('<u2'),1) # 1 uint16 channel
        color = np.fromfile(f, np.dtype('<u1'), 3) # 3 uint8 color codes
        pcProj = np.fromfile(f, np.float32, 2) # 2 float32 component projections
        sampleFreq = np.fromfile(f, np.dtype('<u2'),1)
        
        waveforms = np.fromfile(f, np.dtype('<u2'), numChannels*numSamples) # n*m uint16 samples
        wv = np.reshape(waveforms, (numChannels, numSamples)) 

        gain[currentSpike,:] = np.fromfile(f, np.float32, numChannels) # n float32 channelGains
        thresh[currentSpike,:] = np.fromfile(f, np.dtype('<u2'), numChannels) # n uint16 thresholds
        
        recNum[currentSpike] = np.fromfile(f, np.dtype('<u2'), 1) # 1 uint16 recordingNumber       
        
        for ch in range(numChannels):
            spikes[currentSpike,:,ch] = (np.float64(wv[ch])-32768)/(gain[currentSpike,ch]/1000)
        
        currentSpike += 1
        
    data['spikes'] = spikes[:currentSpike,:,:]
    data['timestamps'] = timestamps[:currentSpike]
    data['source'] = source[:currentSpike]
    data['gain'] = gain[:currentSpike,:]
    data['thresh'] = thresh[:currentSpike,:]
    data['recordingNumber'] = recNum[:currentSpike]
    data['sortedId'] = sortedId[:currentSpike]
    data['numChannels'] = numChannels
    data['numSamples'] = numSamples
    
    return data


def readHeader(f):
    """Read header information from the first 1024 bytes of an OpenEphys file.
    
    Args:
        f: An open file handle to an OpenEphys file
    
    Returns: dict with the following keys.
        - bitVolts : float, scaling factor, microvolts per bit
        - blockLength : int, e.g. 1024, length of each record (see 
            loadContinuous)
        - bufferSize : int, e.g. 1024
        - channel : the channel, eg "'CH1'"
        - channelType : eg "'Continuous'"
        - date_created : eg "'15-Jun-2016 21212'" (What are these numbers?)
        - description : description of the file format
        - format : "'Open Ephys Data Format'"
        - header_bytes : int, e.g. 1024
        - sampleRate : float, e.g. 30000.
        - version: eg '0.4'
        Note that every value is a string, even numeric data like bitVolts.
        Some strings have extra, redundant single apostrophes.
    """
    header = {}
    
    # Read the data as a string
    # Remove newlines and redundant "header." prefixes
    # The result should be a series of "key = value" strings, separated
    # by semicolons.
    header_string = f.read(1024).decode().replace('\n','').replace('header.','')
    
    # Parse each key = value string separately
    for pair in header_string.split(';'):
        if '=' in pair:
            key, value = pair.split(' = ')
            key = key.strip()
            value = value.strip()
            
            # Convert some values to numeric
            if key in ['bitVolts', 'sampleRate']:
                header[key] = float(value)
            elif key in ['blockLength', 'bufferSize', 'header_bytes']:
                header[key] = int(value)
            else:
                # Keep as string
                header[key] = value

    return header


def get_number_of_records(filepath):
    # Open the file
    with open(filepath, 'rb') as f:
        # Read header info
        header = readHeader(f)
        
        # Get file length
        fileLength = os.fstat(f.fileno()).st_size
        
        # Determine the number of records
        record_length_bytes = 2 * header['blockLength'] + 22
        n_records = int((fileLength - 1024) / record_length_bytes)
        if (n_records * record_length_bytes + 1024) != fileLength:
            raise IOError("file does not divide evenly into full records")
    
    return n_records


