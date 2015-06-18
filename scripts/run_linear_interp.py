import time
t0 = time.time()
import sys; sys.path.append('../slip_rate_tools/')

#import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
import slip_rate_tools as srt

'''
Test data.  This will be removed when the adding data interface is built.
'''

offset_df = pd.read_csv('../test_data/offsets.csv')
offset_df['offset_m'] = offset_df.offset_in * 200.

t1 = offset_df[offset_df.unit == 'T1']
qa = offset_df[offset_df.unit == 'Qa']
qao = offset_df[offset_df.unit == 'Qao']

#qa['offset_m'] += 200.

t1_age = {'mean': 24., 'sd':8.}
qa_age = {'mean': 50., 'sd':20.}
qao_age = {'mean':100., 'sd':32.}

#qao_age['mean'] += 200

T1 = srt.OffsetMarker(age_mean=t1_age['mean'], age_sd=t1_age['sd'],
                      offset_vals=t1.offset_m, offset_probs=t1.rel_prob)

Qa = srt.OffsetMarker(age_mean=qa_age['mean'], age_sd=qa_age['sd'],
                      offset_vals=qa.offset_m, offset_probs=qa.rel_prob)

Qao = srt.OffsetMarker(age_mean=qao_age['mean'], age_sd=qao_age['sd'],
                      offset_vals=qao.offset_m, offset_probs=qao.rel_prob)


'''
Real function
'''

def run_linear_interp(offset_list, n_iters, zero_offset_age=0.,
                      check_increasing=False, check_rate_change=False):
    '''
    Main linear interpolation function. Runs both 

    Arguments:
    offset_list = list of OffsetMarkers, in increasing age order.
    n_iters = integer specifying number of Monte Carlo iterations.
    zero_offset_age = float specifying the age at which no more offset occurs;
                      i.e. the last age of faulting.
    check_increasing = Boolean value to ensure increasing offsets with
                       increasing ages.
    check_rate_change = Boolean value to perform a piecewise linear
                        interpolation (2 pieces) 
    '''

    srt.check_unit_consistency(offset_list)

    n_pts = len(offset_list) + 1 # accounting for 0 offset

    print('sampling offset markers')
    age_arr, off_arr = srt.make_age_offset_arrays(offset_list, n_iters,
                                             check_increasing=check_increasing)

    print('doing fits')
    results_df = srt.do_linear_fits(age_arr, off_arr, 
                                    check_rate_change=check_rate_change,
                                    trim_results=True)

    results_df['log_like_1'] = srt.log_likelihood(results_df.sumsq1, n_pts)

    if check_rate_change==True:
        num_1_odds = srt.rate_change_test(results_df, 4, print_res=True)


    elif check_rate_change==False:
        print('\nbest fit slip rate results:')
        print(results_df.m.describe())

    return results_df, age_arr


if __name__ == '__main__':
    res_df, age_arr= run_linear_interp([T1, Qa, Qao], 500, 
                                       check_increasing=True, 
                                       check_rate_change=True)

    rate_hist_df = srt.make_rate_hist_array(res_df, age_arr, n_segments=2)

    cum_hist_df = srt.make_cum_hist_array(rate_hist_df)

    t1 = time.time()
    print(t1-t0)

    #plt.figure()
    #for i in res_df.index




