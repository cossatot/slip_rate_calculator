# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 14:00:54 2012

@author: Richard
"""

import time
from collections import OrderedDict
import numpy as np
import numpy.ma as ma 
#import Splines
#from Splines import spline1d
from scipy.interpolate import interp1d
from scipy.optimize import minimize, curve_fit, leastsq
from scipy.stats import gaussian_kde  # do this for some lists
import pandas as pd
#import matplotlib.pyplot as plt

try:  #consider moving this to separate file
    from matplotlib.backends import qt_compat
    from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.figure import Figure
    from matplotlib import collections as mc
    from PyQt4 import QtGui, QtCore

except:
    pass

# TODO: Use Bayes to refine offset estimates given slip rate constraints


def tspline_interpolate():
    pass


def fit_history_spline(age_array, offset_array):

    return interp1d(age_array, offset_array)


def sample_slip_history(age_array, offset_array, time_array, 
                        extend_time=False):

    history_spline = fit_history_spline(age_array, offset_array)

    if extend_time == False:
        time_array = time_array[time_array <= np.max(age_array)]
    elif extend_time == True:
        raise Exception('extrapolating not supported yet')

    return history_spline(time_array)


def inverse_transform_sample(vals, probs, n_samps, n_interp=1000, seed=False,
                            seed_val=69):
    pdf_range, pdf_vals = make_pdf(vals, probs, n_interp)
    cdf_range, cdf_vals = make_cdf(pdf_range, pdf_vals)

    cdf_interp = interp1d(cdf_vals, cdf_range)

    if seed == True:
        np.random.seed(seed_val)
    
    samps = np.random.rand(n_samps)

    return cdf_interp(samps)


def make_pdf(vals, probs, n_interp=1000):
    val_min = np.min(vals)
    val_max = np.max(vals)
    
    pdf_range = np.linspace(val_min, val_max, n_interp)

    pmf = interp1d(vals, probs)
    pmf_samples = pmf(pdf_range)
    pdf_vals = pmf_samples / np.sum(pmf_samples)

    return pdf_range, pdf_vals


def make_cdf(pdf_range, pdf_vals):
    return (pdf_range, np.cumsum(pdf_vals))


class OffsetMarker:
    """Represents an offset geologic marker.
    
    Attributes:
        offsets: list of possible offset distances for the given marker.
            If offset_type = normal, offsets = [mean, sd]
        offset_probs: list of probabilities of corresponding offset distances
        offset_dist_type: offset prob distribution (normal, uniform, arbitrary)
        ages: list of possible ages for the given marker
        age_probs: list of probabilities of corresponding ages
        age_dist_type: age prob. distribution (normal, uniform, arbitrary)
        source: Source for information (e.g., what article, field campaign)
    
    """
    # TODO: Need to make a random.choice setting for large arrays of vals

    def __init__(self, offsets=np.array([]), offset_probs=None,
                 offset_vals=None, offset_mean=None, offset_median=None,
                 offset_sd=None, offset_mad=None,
                 offset_min=None, offset_max=None,
                 offset_seed=None,
                 offset_dist_type='unspecified', offset_units='unspecified',
                 ages=np.array([]),
                 age_probs=None, age_vals=None, 
                 age_mean=None, age_median=None, age_sd=None, age_mad=None,
                 age_min=None, age_max=None,
                 age_seed=None,
                 age_dist_type='unspecified', age_units='unspecified',
                 source='None'):

            self.offsets = offsets
            self.offset_probs = offset_probs          
            self.offset_vals = offset_vals
            self.offset_mean = offset_mean
            self.offset_median = offset_median
            self.offset_sd = offset_sd
            self.offset_mad = offset_mad
            self.offset_min = offset_min
            self.offset_max = offset_max
            self.offset_units = offset_units

            if offset_dist_type != 'unspecified':
                self.offset_dist_type = offset_dist_type
            elif offset_dist_type == 'unspecified':
                if offset_mean is not None and offset_sd is not None:
                    self.offset_dist_type = 'normal'
                elif (offset_min is not None and offset_max is not None
                      and offset_sd == None):
                    self.offset_dist_type = 'uniform'
                elif offset_probs is not None and offset_vals is not None:
                    self.offset_dist_type = 'arbitrary'

            self.ages = ages
            self.age_probs = age_probs          
            self.age_vals = age_vals
            self.age_mean = age_mean
            self.age_median = age_median
            self.age_sd = age_sd
            self.age_mad = age_mad
            self.age_min = age_min
            self.age_max = age_max
            self.age_units = age_units

            if age_dist_type != 'unspecified':
                self.age_dist_type = age_dist_type
            elif age_dist_type == 'unspecified':
                if age_mean is not None and age_sd is not None:
                    self.age_dist_type = 'normal'
                elif (age_min is not None and age_max is not None
                      and age_sd == None):
                    self.age_dist_type = 'uniform'
                elif age_probs is not None and age_vals is not None:
                    self.age_dist_type == 'arbitrary'
                     
            self.source = source
            
    
    def sample_offset_from_normal(self, n):
        """Generates n-length sample from normal distribution of offsets"""

        return sample_from_bounded_normal(self.offset_mean, self.offset_sd, n,
                                          self.offset_min, self.offset_max)
    
    def sample_offset_from_uniform(self, n):
        """Generates n-length sample from uniform distribution of ages"""
        
        return np.random.uniform(self.offset_min, self.offset_max, n)
    
    def sample_offset_from_arbitrary(self, n):
        """not supported yet"""
        offset_sample = inverse_transform_sample(self.offset_vals,
                                                 self.offset_probs, n)
        return offset_sample
    
    def sample_offset(self, n):
        """Generates n-length array of samples from distribution"""
        if self.offset_dist_type == 'normal':
            offset_sample = self.sample_offset_from_normal(n)
        
        elif self.offset_dist_type == 'uniform':
            offset_sample = self.sample_offset_from_uniform(n)
        
        elif self.offset_dist_type == 'arbitrary':
            offset_sample = self.sample_offset_from_arbitrary(n)
        
        else:
            print('What is the offset distribution type?')
        
        return offset_sample    
    
    def sample_age_from_normal(self, n):
        """Generates n-length sample from normal distribution of ages"""
        if self.age_min:
            age_min = self.age_min
        else:
            age_min = 0.

        age_sample = sample_from_bounded_normal(self.age_mean, self.age_sd, n,
                                           age_min, self.age_max)

        return age_sample
    
    def sample_age_from_uniform(self, n):
        """Generates n-length sample from uniform distribution of ages"""
        age_sample = (np.random.rand(n) * 2 -1) * self.age_err + self.age_mean
        #TODO: change mean, err to min, max
        
        return age_sample
        
    def sample_age_from_arbitrary(self, n):
        """not supported yet"""
        age_sample = inverse_transform_sample(self.age_vals, self.age_probs, n)

        return age_sample
    
    def sample_age(self, n):
        """Generates n-length array of samples from distribution"""
        if self.age_dist_type == 'normal':
            age_sample = self.sample_age_from_normal(n)
        
        elif self.age_dist_type == 'uniform':
            age_sample = self.sample_age_from_uniform(n)
        
        elif self.age_dist_type == 'arbitrary':
            pass
        
        else:
            print('What is the age distribution type?')
        
        return age_sample

    def sample(self, n):
        age_sample = self.sample_age(n)
        offset_sample = self.sample_offset(n)
        
        asl = len(age_sample)
        osl = len(offset_sample)
        
        if asl > osl:
            age_sample = age_sample[0:osl]
        elif osl > asl:
            offset_sample = offset_sample[0:asl]
        
        return age_sample, offset_sample


def offset_list_from_gui(tabledata, table_header):
    offsets_d = offset_markers_from_gui(tabledata, table_header)

    return list(offsets_d.values())

    
def offset_markers_from_gui(tabledata, table_header):
    offsets_d = OrderedDict()

    for row in tabledata:
        off_mark_d = offset_marker_dict_from_row(row, table_header)
        offsets_d[off_mark_d['Name']] = offset_marker_from_dict(off_mark_d)

    return offsets_d


def offset_marker_dict_from_row(row, table_header):
    # header_table should be passed from gui
    off_mark_d = OrderedDict()

    for i, key in enumerate(table_header):
        off_mark_d[key] = row[i]

    return off_mark_d


def offset_marker_from_dict(off_row_d):
    or_d = off_row_d

    args = {'offset_units': or_d['Offset_Units'],
            'age_units': or_d['Age_Units']}

    # get offset arguments
    if or_d['Offset_Type'] == 'mean':
        if not np.isscalar(or_d['Offset']):
            raise Exception('Mean Offset has to be a scalar!')
        args['offset_mean'] = or_d['Offset']

    elif or_d['Offset_Type'] == 'median':
        if not np.isscalar(or_d['Offset']):
            raise Exception('Median Offset has to be a scalar!')
        args['offset_median'] = or_d['Offset']

    elif or_d['Offset_Type'] == 'list':
        if len(or_d['Offset']) < 2:
            raise Exception('List Offsets have to be longer than 1!')
        args['offset_vals'] = or_d['Offset']
    
    else:
        raise Exception('Offset_Type must be mean, median or list!')

    # get offset err arguments
    # TODO: More consistency checking between arg types
    if or_d['Offset_Err_Type'] == 'sd':
        if not np.isscalar(or_d['Offset_Err']):
            raise Exception('sd Offset_Err must be a scalar!')
        args['offset_sd'] = or_d['Offset_Err']

    elif or_d['Offset_Err_Type'] == 'mad':
        if not np.isscalar(or_d['Offset_Err']):
            raise Exception('mad Offset_Err must be a scalar!')
        args['offset_mad'] = or_d['Offset_Err']

    elif or_d['Offset_Err_Type'] == 'minmax':
        if not np.isscalar(or_d['Offset_Err']):
            raise Exception('minmax Offset_Err must be a scalar!')
        if not np.isscalar(or_d['Offset']):
            raise Exception('Mean Offset has to be a scalar!')
        args['offset_min'] = or_d['Offset'] - or_d['Offset_Err']
        args['offset_max'] = or_d['Offset'] + or_d['Offset_Err']
        args['offset_sd'] = None # just to make sure the class inits right
   
    elif or_d['Offset_Err_Type'] == 'probs':
        if len(or_d['Offset_Err']) < 2:
            raise Exception('probs Offset_Err have to be longer than 1!')
        args['offset_probs'] = or_d['Offset_Err']
        # check to make sure offset vals are set too?

    elif or_d['Offset_Err_Type'] == 'kde':
        if len(or_d['Offset_Err']) < 2:
            raise Exception('kde Offset_Err have to be longer than 1!')
        args['offset_probs'] = kde(or_d['Offset'])

    else:
        raise Exception('Offset_Err_Type must be sd, mad, minmax, probs, '
                        +'or kde!')

    # get age arguments
    if or_d['Age_Type'] == 'mean':
        if not np.isscalar(or_d['Age']):
            raise Exception('Mean Age has to be a scalar!')
        args['age_mean'] = or_d['Age']

    elif or_d['Age_Type'] == 'median':
        if not np.isscalar(or_d['Age']):
            raise Exception('Median Age has to be a scalar!')
        args['age_median'] = or_d['Age']

    elif or_d['Age_Type'] == 'list':
        if len(or_d['Age']) < 2:
            raise Exception('List Ages have to be longer than 1!')
        args['age_vals'] = or_d['Age']
    
    else:
        raise Exception('Age_Type must be mean, median or list!')

    # get age err arguments
    # TODO: More consistency checking between arg types
    if or_d['Age_Err_Type'] == 'sd':
        if not np.isscalar(or_d['Age_Err']):
            raise Exception('sd Age_Err must be a scalar!')
        args['age_sd'] = or_d['Age_Err']

    elif or_d['Age_Err_Type'] == 'mad':
        if not np.isscalar(or_d['Age_Err']):
            raise Exception('mad Age_Err must be a scalar!')
        args['age_mad'] = or_d['Age_Err']

    elif or_d['Age_Err_Type'] == 'minmax':
        if not np.isscalar(or_d['Age_Err']):
            raise Exception('minmax Age_Err must be a scalar!')
        if not np.isscalar(or_d['Age']):
            raise Exception('Mean Age has to be a scalar!')
        args['age_min'] = or_d['Age'] - or_d['Age_Err']
        args['age_max'] = or_d['Age'] + or_d['Age_Err']
        args['age_sd'] = None # just to make sure the class inits right
   
    elif or_d['Age_Err_Type'] == 'probs':
        if len(or_d['Age_Err']):
            raise Exception('probs Age_Err have to be longer than 1!')
        args['age_probs'] = or_d['Age_Err']
        # check to make sure age vals are set too?

    elif or_d['Age_Err_Type'] == 'kde':
        if len(or_d['Age_Err']) < 2:
            raise Exception('kde Age_Err have to be longer than 1!')
        args['age_probs'] = kde(or_d['Age'])
        
    else:
        raise Exception('Age_Err_Type must be sd, mad, minmax, probs, '
                        +'or kde!')

    return OffsetMarker(**args)


def kde(vals):
    # not sure how to do this yet
    # need to match input length?  or just resample? pass resampling to class?
    # will need to re-do vals too!
    
    raise Exception('Not Implemented Yet')


######
### stats functions
#####

def sample_from_bounded_normal(mean, sd, n, sample_min=None, sample_max=None):

    sample = np.random.normal(mean, sd, n)
    sample = trim_distribution(sample, sample_min=sample_min, 
                                       sample_max=sample_max)

    while len(sample) < n:
        next_sample = np.random.normal(mean, sd, n)
        next_sample = trim_distribution(next_sample, sample_min, sample_max)
        sample = np.hstack([sample, next_sample])

    return sample[:n]


def trim_distribution(sample, sample_min=None, sample_max=None):

    if sample_min is not None and sample_max is not None:
        if sample_min >= sample_max:
            raise Exception('min must be less than max!')

    if sample_min is not None:
        sample = sample[sample >= sample_min]

    if sample_max is not None:
        sample = sample[sample <= sample_max]

    return sample


def check_monot_increasing(in_array):
    """Checks to see if array is monotonically increasing, returns bool value
    """
    dx = np.diff(in_array)

    return np.all(dx >= 0)


def check_unit_consistency(offset_list):
    off_unit_list = [om.offset_units for om in offset_list]
    age_unit_list = [om.age_units for om in offset_list]

    for off_u in off_unit_list:
        if off_u != off_unit_list[0]:
            raise Exception('OffsetMarker units not consistent.')
    
    for age_u in age_unit_list:
        if age_u != age_unit_list[0]:
            raise Exception('OffsetMarker units not consistent.')

    return


def get_log_pts(p_min, p_max, n_pts=50, base=np.e):
    """Generates n_pts length array of logarithmically spaced points"""
    if p_min == 0:
        pts_array = np.hstack([0, np.logspace(np.log(1e-5), np.log(p_max),
                                                num=n_pts-1, base=base)])
    else:
        pts_array = np.logspace(p_min, p_max, num=n_pts, base=base)

    return pts_array


def make_age_offset_arrays(offset_list, n, force_increasing=False, 
                           zero_offset_age=0., seed=False, seed_value=None,
                           sample_chunks=1):

    # TODO: implement sample chunking (using n samples per marker per fit)
    
    if seed == True:
        np.random.seed(seed_value)
    
    age_array = np.zeros((n, len(offset_list)+1 * sample_chunks))
    off_array = np.zeros((n, len(offset_list)+1 * sample_chunks))
    
    age_array[:,0] = zero_offset_age
    
    for i, off_mark in enumerate(offset_list):
        col = i+1
        age_array[:,col], off_array[:,col] = off_mark.sample(n)
        
    if force_increasing == True:
        
        def make_inc_bool(age_array, off_array, n):
        
            inc_bool = np.ones((age_array.shape[0]), dtype=int)
            for row in range(n):
                age_inc = check_monot_increasing(age_array[row,:])
                off_inc = check_monot_increasing(off_array[row,:])
                
                if not (age_inc and off_inc):
                    inc_bool[row] = 0
                    
            inc_bool = np.array(inc_bool, dtype=bool)
                
            return inc_bool
    
        inc_bool = make_inc_bool(age_array, off_array, n)
                    
        age_array = age_array[inc_bool, :]
        off_array = off_array[inc_bool, :]
        
        while age_array.shape[0] < n:
            
            next_age_array, next_off_array = make_age_offset_arrays(
                                                offset_list, n,
                                                force_increasing=False,
                                                zero_offset_age=zero_offset_age)
            
            next_inc_bool = make_inc_bool(next_age_array, next_off_array, n)
            
            next_age_array = next_age_array[next_inc_bool, :]
            next_off_array = next_off_array[next_inc_bool, :]
           
            off_array = np.vstack([off_array, next_off_array])
            age_array = np.vstack([age_array, next_age_array])
            
    return age_array[:n,:], off_array[:n,:]


def piece_lin_objective(params, x_data, y_data): 
    '''docs

    Modified from a function by Andreas Hillboll on the StatsModels
    mailing list.
    '''
    y1 = 0.
    y2, y3, x2 = params
    x1, x3 = x_data[0], x_data[-1] 
    Xbefore = y1 + (x_data - x1) * (y2 - y1) / (x2 - x1) 
    Xafter = y2 + (x_data - x2) * (y3 - y2) / (x3 - x2) 
    Xbreak = np.where(x_data <= x2, Xbefore, Xafter) 
    return (ma.masked_invalid(Xbreak - y_data)**2).sum()


def piece_lin_opt(x_data, y_data):
    
    init_guesses = (np.mean(y_data), np.mean(y_data), np.mean(x_data))
    bounds = ((0, np.max(y_data)), (0., np.max(y_data)), (0., np.max(y_data)))
    
    
    res = minimize(piece_lin_objective, init_guesses, (x_data, y_data),
                   method="TNC", bounds=bounds)
    
    sum_sq_err = piece_lin_objective(res.x, x_data, y_data)
    
    y2, y3, x2 = res.x
    
    slope1 = y2 / x2
    slope2 = ((y3 - y2) / (np.max(x_data) - x2))
    breakpoint = x2
    
    return slope1, slope2, breakpoint, sum_sq_err


def piecewise_linear(x, breakpt, m1, m2):
    return np.piecewise(x, [x < breakpt], [lambda x: m1 * x,
                             lambda x: m2 * x + (m1 * breakpt) - m2 * breakpt])


def piecewise_linear_objective(params, x_data, y_data):

    return ( (y_data - piecewise_linear(x_data, *params))**2).sum()   


def penalized_piecewise_linear_objective(params, x_data, y_data, weight=0.1):
    
    breakpt, m1, m2 = params
    
    resids = np.array( (y_data - piecewise_linear(x_data, *params)) )
    
    rate_change_penalization = np.sum(np.abs(resids)) * np.abs(m1 - m2) * weight
    
    new_resids = np.append(resids, rate_change_penalization)

    return new_resids


def piecewise_linear_opt(x_data, y_data):

    breakpt_guess = np.median(x_data)
    m1_guess = x_data.max() / y_data.max()
    m2_guess = x_data.max() / y_data.max()
    init_vals = [breakpt_guess, m1_guess, m2_guess]

    try:
        params, cov_matrix = curve_fit(piecewise_linear, x_data, y_data, 
                                       init_vals)
    except RuntimeError:
        results = minimize(piecewise_linear_objective, init_vals,
                           (x_data, y_data), method='SLSQP')
        #print('slsqp')
        params = results.x

    # params = 

    breakpt, m1, m2 = params

    errs = y_data - piecewise_linear(x_data, breakpt, m1, m2)

    sum_sq_err = np.sum(errs**2)

    return m1, m2, breakpt, sum_sq_err


def penalized_piecewise_linear_opt(x_data, y_data, weight=0.3):
    breakpt_guess = np.median(x_data)
    m1_guess = x_data.max() / y_data.max()
    m2_guess = x_data.max() / y_data.max()
    init_vals = (breakpt_guess, m1_guess, m2_guess)
    
    params, success = leastsq(penalized_piecewise_linear_objective, init_vals, 
                              args=(x_data, y_data, weight))

    breakpt, m1, m2 = params
    
    errs = y_data - piecewise_linear(x_data, breakpt, m1, m2)

    sum_sq_err = np.sum(errs**2)

    return m1, m2, breakpt, sum_sq_err


def lin_fit(x_data, y_data):
    x = x_data[:,np.newaxis]
    m, _, _, _ = np.linalg.lstsq(x, y_data)
    m = m[0]
    
    sum_sq_err = ((y_data - (m * x_data))**2).sum()
    
    return m, sum_sq_err


def do_linear_fits(age_arr, off_arr, fit_type=None, trim_results=True,
                   n_linear_pieces=None):

    n_iters = age_arr.shape[0]
    
    if fit_type == 'piecewise':
        if n_linear_pieces == 2:
            results_columns = ['m1', 'm2', 'breakpt', 'sumsq2', 'm', 'sumsq1']
        else:
            raise Exception('Only 2 piece piecewise-linear fits supported')
    elif fit_type == 'linear':
        results_columns = ['m', 'sumsq1']

    results_arr = np.zeros( (n_iters, len(results_columns) ) )

    if fit_type == 'linear':
        for i in range(n_iters):
            xd = age_arr[i,:]
            yd = off_arr[i,:]

            results_arr[i,:] = lin_fit(xd, yd)

    elif fit_type == 'piecewise':
        for i in range(n_iters):
            xd = age_arr[i,:]
            yd = off_arr[i,:]

            results_arr[i, 4:6] = lin_fit(xd, yd)

            #results_arr[i, 0:4] = piece_lin_opt(xd, yd)
            #results_arr[i, 0:4] = piecewise_linear_opt(xd, yd)
            results_arr[i, 0:4] = penalized_piecewise_linear_opt(xd, yd)

    results_df = pd.DataFrame(results_arr, columns=results_columns)
     
    if fit_type == 'piecewise':
       if trim_results==True:
           # option will be set in the GUI
           
           results_df = results_df[(results_df.breakpt > age_arr[:,0])
                                   &(results_df.breakpt < age_arr[:,-1])]
           
           m1_75 = results_df.m1.describe()['75%']
           m2_75 = results_df.m2.describe()['75%']

           m1_max = 5 * m1_75
           m2_max = 5 * m2_75

           results_df = results_df[(results_df.m1 < m1_max)] 
           results_df = results_df[(results_df.m2 < m2_max)] 


    return results_df


def log_likelihood(sum_sq, n):
    
    return -n / 2 * np.log(sum_sq)


def BIC(log_likelihood, n, p):
    
    return log_likelihood - ( 0.5 * p * np.log(n / 2 * np.pi))


def find_nearest_index(array, value):
    idx = (np.abs(array-value)).argmin()
    return idx


def rate_change_test(results_df, n_offsets, print_res=False):
    results_df['log_like_2'] = log_likelihood(results_df.sumsq2, n_offsets)
    n_iters_out = results_df.shape[0]

    p1 = 1 # number of parameters for single linear fit
    p2 = 3 # number of parameters for 2 part piecewise fit
    results_df['bic_1'] = BIC(results_df.log_like_1, n_offsets, p1)
    results_df['bic_2'] = BIC(results_df.log_like_2, n_offsets, p2)

    num_1_count = results_df[results_df.bic_1 > results_df.bic_2].shape[0]
    num_2_count = n_iters_out - num_1_count
    num_1_odds = num_1_count / n_iters_out
    num_2_odds = num_2_count / n_iters_out
    
    if print_res==True:
        if num_1_odds > num_2_odds:
            print('1 line fits best.  {}/{} ({}% chance)'.format(num_1_count,
                                                               n_iters_out,
                                                               num_1_odds*100))
            print('\nbest fit slip rate results:')
            print(results_df.m.describe())


        else:
            print('2 lines fit best.  {}/{} ({}% chance)'.format(num_2_count,
                                                               n_iters_out,
                                                               num_2_odds*100))
            print('\nbest fit slip rate results:')
            print('rate 1 (younger):')
            print(results_df.m1.describe())
            print('rate change timing:')
            print(results_df.breakpt.describe())
            print('rate 2 (older):')
            print(results_df.m2.describe())
            print('rate_change:')
            print((results_df.m2 - results_df.m1).describe())
   
    return num_1_odds


def linear_rate_interp(rate, run_time_max, sim_time_max, zero_offset_age=0.,
                       num_pts=1000):
    ''' Makes a history array of slip rates.  In this case, the slip
    rate is a constant from zero_offset_age to run_time_max, and is
    zero outside of those boundaries.  Returns a Pandas Series.
    
    Arguments:
    rate (float): slip rate.
    run_time_max (float): Maximum age of slip rate for this MC iteration,
                          i.e. age of oldest offset feature. Times older
                          than this will have zero slip rate.
    sim_time_max (float): Maximum age of oldest feature in the whole MC
                          simulation. This determines the length of the 
                          array.
    zero_offset_age (float): Youngest age of faulting.  Times younger than
                             this time will have zero rate.
    num_pts (int): Number of points in the array.
    '''

    times = np.linspace(zero_offset_age, sim_time_max, num_pts)
    slip_rate_history = pd.Series(index=times, data=np.zeros(num_pts))

    slip_rate_history.ix[zero_offset_age : run_time_max] = rate

    return slip_rate_history


def piecewise_rate_interp(rate1, rate2, breakpt, run_time_max, sim_time_max,
                          zero_offset_age=0., num_pts=1000):

    times = np.linspace(zero_offset_age, sim_time_max, num_pts)
    slip_rate_history = np.zeros(num_pts)
    
    
    zero_offset_idx = find_nearest_index(times, zero_offset_age)
    run_time_max_idx = find_nearest_index(times, run_time_max)
    breakpt_idx = find_nearest_index(times, breakpt)
    
    slip_rate_history[zero_offset_idx : breakpt_idx] = rate1
    slip_rate_history[breakpt_idx : run_time_max_idx] = rate2


    return slip_rate_history


def make_rate_hist_array(results_df, age_arr, n_segments=1, num_pts=1000,
                         zero_offset_age=0., return_array=False,
                         sim_time_max='mc_age_max'):

    if sim_time_max == 'mc_age_max':
        sim_time_max = np.max(age_arr)

    times = np.linspace(zero_offset_age, sim_time_max, num_pts)

    rate_hist_df = pd.DataFrame(columns=times, index=results_df.index)
    rate_hist_ar = np.zeros((len(results_df.index), num_pts))

    if n_segments == 1:
        for i in rate_hist_df.index:
            rate = results_df.ix[i, 'm']
            run_time_max = age_arr[i, -1]
            rate_hist_df.ix[i, :] = linear_rate_interp(rate, run_time_max,
                                                       sim_time_max,
                                                       zero_offset_age,
                                                       num_pts)
    elif n_segments == 2:
        for i, row in enumerate(results_df.index):
            rate1 = results_df.ix[row, 'm1']
            rate2 = results_df.ix[row, 'm2']
            breakpt = results_df.ix[row, 'breakpt'] 
            run_time_max = age_arr[row, -1]
            #rate_hist_df.ix[i, :] = piecewise_rate_interp(rate1, rate2,
            #rate_hist_df.loc[i, :] = piecewise_rate_interp(rate1, rate2,
            rate_hist_ar[i, :] = piecewise_rate_interp(rate1, rate2,
                                                          breakpt, 
                                                          run_time_max,
                                                          sim_time_max,
                                                          zero_offset_age,
                                                          num_pts)
    else:
        raise Exception('Only 1 or 2 rates supported now.')
    
    #return rate_hist_df if return_array == True else rate_hist_df.values
    return rate_hist_ar

def make_cum_hist_array(rate_hist_array):

    return np.cumsum(rate_hist_array, axis=0)


def run_interp_from_gui(offset_list, run_config_dict):
    t0 = time.time()

    rc = run_config_dict

    check_unit_consistency(offset_list)

    n_offsets = len(offset_list) + 1

    print('sampling offset markers')
    age_arr, off_arr = make_age_offset_arrays(offset_list, rc['n_iters'],
                                       force_increasing=rc['force_increasing'],
                                       seed=rc['random_seed'],
                                       seed_value=rc['random_seed_value'])

    print('doing fits')
    if rc['fit_type'] in ['linear', 'piecewise']:
        results_df = do_linear_fits(age_arr, off_arr, fit_type=rc['fit_type'],
                                    n_linear_pieces=rc['n_linear_pieces'])
    else:
        raise Exception('fit type not implemented yet')

    results_df['log_like_1'] = log_likelihood(results_df.sumsq1, n_offsets)

    if rc['fit_type']  == 'linear':
        print(results_df.m.describe())

    elif rc['fit_type'] == 'piecewise':
        one_rate_odds = rate_change_test(results_df, n_offsets, print_res=True)
    
    print("\ndone in {:.2f} seconds".format(time.time() - t0))

    return results_df, age_arr, off_arr


def trim_age_offset_arrays(res_df, age_arr, off_arr=None):
    """
    Trims age and offset arrays based on retained values from the results_df.
    """
    good_inds = res_df.index.values

    age_arr_trim = age_arr[good_inds, :]

    if off_arr is not None:
        off_arr_trim = off_arr[good_inds, :]

        return age_arr_trim, off_arr_trim

    else:
        return age_arr_trim


def cumulative_offsets(prev_age, prev_rate, new_age, new_rate):
    return prev_age * prev_rate + (new_age - prev_age) * new_rate


def get_line_pts_from_results(res_df, age_array, n_pieces=1):
    
    if age_array.shape[0] != res_df.shape[0]:
        age_array = trim_age_offset_arrays(res_df, age_array)

    n_pts = n_pieces + 1
    n_rows = res_df.shape[0]

    x_array = np.zeros((res_df.shape[0], n_pts))
    y_array = x_array.copy()
    
    if n_pieces == 1:
        x_array[:,0] = age_array[:,0]
        x_array[:,1] = age_array.max(axis=1)

        y_array[:,0] = 0.
        y_array[:,1] = age_array.max(axis=1) * res_df.m

    elif n_pieces == 2:
        x_array[:,0] = age_array[:,0]
        x_array[:,1] = res_df.breakpt.values
        x_array[:,2] = age_array.max(axis=1)

        y_array[:,0] = 0.
        y_array[:,1] = res_df.breakpt.values * res_df.m1
        y_array[:,2] = cumulative_offsets(x_array[:,1], res_df.m1, 
                                          x_array[:,2], res_df.m2)
    else:
        raise Exception('only 1 or 2 piece lines currently implemented.')

    line_pts = [tuple(zip(x_array[i,:], y_array[i,:])) for i in range(n_rows)]

    return line_pts


def plot_slip_histories_from_gui(res_df, age_arr, run_config_dict, 
                                 offset_arr=None, offset_list=None,
                                 show_data=False, show_samples=False):

    if run_config_dict['fit_type'] == 'linear':
        n_pieces = 1

    elif run_config_dict['fit_type'] == 'piecewise':
        n_offsets = age_arr.shape[1]
        one_rate_odds = rate_change_test(res_df, n_offsets, print_res=False)

        if one_rate_odds >= 0.5:
            n_pieces = 1
        elif one_rate_odds < 0.5:
            # TODO: Use best value if n_linear_pieces > 2, when implemented
            n_pieces = run_config_dict['n_linear_pieces']
    
    elif run_config_dict['fit_type'] == 'cubic':
        raise Exception('cubic spline plotting not implemented yet')

    else:
        raise Exception('fit type needs to be linear, piecewise, or cubic')

    line_pts = get_line_pts_from_results(res_df, age_arr, n_pieces)

    line_coll = mc.LineCollection(line_pts, linewidths=0.1, colors='k')

    canvas = MplCanvas()
    ax = canvas.axes

    if show_data == True:
        ax.errorbar()

    if show_samples == True:

        n_iters = run_config_dict['n_iters']

        if n_iters < 10:
            sym = 'o'
        elif 10 <= n_iters < 100:
            sym = '.'
        elif 100 <= n_iters:
            sym = ','

        ax.plot(age_arr.ravel(), offset_arr.ravel(), sym)
    
    ax.add_collection(line_coll)
    ax.autoscale()

    canvas.show()


def plot_histograms_from_gui(run_config_dict, results_df):
    rc = run_config_dict


    canvas = MplCanvas()


    if rc['fit_type'] == 'linear':
        canvas.axes.hist(results_df.m, bins=50)

    elif rc['fit_type'] == 'piecewise':
        pass

    canvas.show()
    #return


class MplCanvas(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        # We want the axes cleared every time plot() is called
        self.axes.hold(False)

        #self.compute_initial_figure()

        #
        FigureCanvas.__init__(self, fig)
        self.setParent(parent)

        #FigureCanvas.setSizePolicy(self,
        #                           QtGui.QSizePolicy.Expanding,
        #                           QtGui.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
