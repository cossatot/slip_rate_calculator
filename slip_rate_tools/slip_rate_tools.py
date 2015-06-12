# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 14:00:54 2012

@author: Richard
"""

import numpy as np
#import Splines
from Splines import spline1d
from scipy.interpolate import interp1d
from scipy.optimize import minimize


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
        offset_dist_type: offset prob. distribution (normal, uniform, arbitrary)
        ages: list of possible ages for the given marker
        age_probs: list of probabilities of corresponding ages
        age_dist_type: age prob. distribution (normal, uniform, arbitrary)
        source: Source for information (e.g., what article, field campaign)
    
    """
    def __init__(self, offsets=np.array([]), offset_probs=None,
                 offset_vals=None, offset_mean=None, offset_median=None,
                 offset_sd=None, offset_mad=None,
                 offset_min=None, offset_max=None,
                 offset_seed=None,
                 offset_dist_type='unspecified', ages=np.array([]),
                 age_probs=None, age_vals=None, 
                 age_mean=None, age_median=None, age_sd=None, age_mad=None,
                 age_min=None, age_max=None,
                 age_seed=None,
                 age_dist_type='unspecified', 
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

        age_sample = sample_bounded_normal(self.age_mean, self.age_sd, n,
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


def sample_bounded_normal(mean, sd, n, sample_min=None, sample_max=None):

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


def check_increasing(in_array):
    """Checks to see if array is monotonically increasing, returns bool value
    """

    dx = np.diff(in_array)

    return np.all(dx >= 0)


def get_log_pts(p_min, p_max, n_pts=50, base=np.e):
    """Generates n_pts length array of logarithmically spaced points"""
    if p_min == 0:
        pts_array = np.hstack([0, np.logspace(np.log(1e-5), np.log(p_max),
                                                num=n_pts-1, base=base)])
    else:
        pts_array = np.logspace(p_min, p_max, num=n_pts, base=base)

    return pts_array


def tspline_interpolate():
    pass


def make_age_offset_arrays(offset_list, n, check_increasing=False, 
                           zero_offset_age=0.):
    
    #TODO: use random seeding
    
    
    age_array = np.zeros((n, len(offset_list)+1))
    off_array = np.zeros((n, len(offset_list)+1))
    
    age_array[:,0] = zero_offset_age
    
    for i, off_mark in enumerate(offset_list):
        col = i+1
        age_array[:,col], off_array[:,col] = off_mark.sample(n)
        
    if check_increasing == True:
        
        def make_inc_bool(age_array, off_array, n):
        
            inc_bool = np.ones((age_array.shape[0]), dtype=int)
            for row in range(n):
                age_inc = check_increasing(age_array[row,:])
                off_inc = check_increasing(off_array[row,:])
                
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
                                                check_increasing=False,
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
                   method="L-BFGS-B", bounds=bounds)
    
    sum_sq_err = piece_lin_objective(res.x, x_data, y_data)
    
    y2, y3, x2 = res.x
    
    slope1 = y2 / x2
    slope2 = ((y3 - y2) / (np.max(x_data) - x2))
    breakpoint = x2
    
    return slope1, slope2, breakpoint, sum_sq_err


def lin_fit(x_data, y_data):
    x = x_data[:,np.newaxis]
    m, _, _, _ = np.linalg.lstsq(x, y_data)
    m = m[0]
    
    sum_sq_err = ((y_data - (m * x_data))**2).sum()
    
    return m, sum_sq_erri


def log_like(sum_sq, n):
    
    return -n / 2 * np.log(sum_sq)


def BIC(log_like, n, p):
    
    return log_like - ( 0.5 * p * np.log(n / 2 * np.pi))
