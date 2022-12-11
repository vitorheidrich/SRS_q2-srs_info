import pandas as pd
import numpy as np
import random as rd

def SRS(data, Cmin, set_seed = True, seed = 1):
    
    #set the seed that allows reproducible results
    if set_seed is True:
        np.random.seed(seed)
    
    #check whether samples will have to be discarded due to low library size (sequencing depth)
    if Cmin > min(data.sum(axis = 0)):
        max_lib_size = max(data.sum(axis = 0))
        if Cmin > max_lib_size:
            print("Please select a Cmin <= ", max_lib_size)
            raise ValueError("All samples discarded due to low number of counts.")
        else:
            discarded = data.columns.values[data.sum() < Cmin]
            print("Sample(s) discarded due to low number of counts (number of counts < Cmin): ", end='')
            print(*discarded, sep=", ")
            data = data[data.columns[data.sum() >= Cmin]]
            samples = data.columns
            SRS(data, Cmin, set_seed, seed)
            
    #proceed with samples with library size >= Cmin
    else:
        
        #throw error in case of invalid Cmin
        if Cmin < 0:
            raise ValueError("Cmin < 0. Please select Cmin >= 0.")
        if float(Cmin).is_integer() is False:
            raise ValueError("SRS accepts only integers for Cmin.")
        
        #start normalization
        else:
            #scaling
            data_scaled = data.div(data.sum(axis = 0) / Cmin, axis = 1)
            #picking the integer part of the scaled data
            data_scaled_floor = np.floor(data_scaled)
            #picking the fractional part of the scaled data
            data_scaled_frac = data_scaled - data_scaled_floor
            #n of counts missing for each sample in the integer
            #part of the scaled data in order to reach Cmin
            counts_missing = data_scaled_frac.sum(axis = 0)
            data_normalized = data_scaled_floor
            
            #start distribution of missing counts per sample
            for sample in data.columns:
                counts_missing_s = int(counts_missing[sample])
                data_scaled_frac_s = data_scaled_frac[sample]
                data_scaled_frac_s_rank = data_scaled_frac_s.rank(method = 'min').sort_values(ascending = False)
                data_scaled_floor_s = data_scaled_floor[sample]
                while counts_missing_s > 0:
                    for asv in data_scaled_frac_s_rank.index:
                        if counts_missing_s == 0:
                            break
                        rank_current_asv = data_scaled_frac_s_rank.loc[asv]
                        n_asvs_sharing_frac_rank = data_scaled_frac_s_rank.value_counts()[rank_current_asv]
                        
                        #distribute counts based on the frac rank (descending)
                        if n_asvs_sharing_frac_rank <= counts_missing_s:
                            data_normalized.loc[asv, sample] = data_normalized.loc[asv, sample] + 1
                            counts_missing_s -= 1
                            data_scaled_frac_s_rank = data_scaled_frac_s_rank.drop(asv)
                        
                        #if there are asvs sharing the same frac rank and there are less
                        #missing counts than the number of asvs sharing the same frac rank
                        else:
                            asvs_sharing_frac_rank = data_scaled_frac_s_rank.head(n_asvs_sharing_frac_rank).index
                            data_scaled_floor_s_asfr = data_scaled_floor_s.filter(asvs_sharing_frac_rank)
                            data_scaled_floor_s_asfr_rank = data_scaled_floor_s_asfr.rank(method = 'min').sort_values(ascending = False)
                            for asv in data_scaled_floor_s_asfr_rank.index:
                                if counts_missing_s == 0:
                                    break
                                rank_current_asv = data_scaled_floor_s_asfr_rank.loc[asv]
                                n_asvs_sharing_floor_rank = data_scaled_floor_s_asfr_rank.value_counts()[rank_current_asv]
                                
                                #distribute counts based on the int rank (descending)
                                if n_asvs_sharing_floor_rank <= counts_missing_s:
                                    data_normalized.loc[asv, sample] = data_normalized.loc[asv, sample] + 1
                                    counts_missing_s -= 1
                                    data_scaled_floor_s_asfr_rank = data_scaled_floor_s_asfr_rank.drop(asv)
                                
                                #if these asvs also share the same int rank, then
                                #distribute missing counts randomly (without replacement)
                                else:
                                    while counts_missing_s > 0:
                                        random_asv = rd.choice(data_scaled_floor_s_asfr_rank.index)
                                        data_normalized.loc[random_asv, sample] = data_normalized.loc[random_asv, sample] + 1
                                        counts_missing_s -= 1
                                        data_scaled_floor_s_asfr_rank = data_scaled_floor_s_asfr_rank.drop(random_asv)
    
    #return df with SRS-normalized data (counts)
    return(data_normalized.astype(int))
