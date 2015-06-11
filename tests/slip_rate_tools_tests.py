# -*- coding: utf-8 -*-
"""
Created on Sat Sep 01 17:34:29 2012

@author: Richard
"""

#import sys
#sys.path.append('C:\\python_modules')
#sys.path.append('C:\\itchy\code_repos\\src-working')
import numpy as nmp
import slip_rate_tools as srt

#om = srt.OffsetMarker(offsets = [])


def truncate_samples_test():
	ar1 = nmp.array([[1, 3, 4, 5], [2, 3, 4, 5]])
	ar2 = nmp.array([[7, 8, 9], [11, 12, 13]])
	sample_list = [ar1, ar2]
	ar1_short, ar2_short = srt.truncate_samples(sample_list)
	print 'ar1 shape =', ar1.shape, 'and ar1_short shape =', ar1_short.shape

def check_increasing_test():
    inc_array = nmp.array([0, 2, 4, 5, 6])
    dec_array = nmp.array([3, 3, 5, 1, 4])
    inc_check = srt.check_increasing(inc_array)
    dec_check = srt.check_increasing(dec_array)
    print 'good array is', inc_check
    print 'bad array is', dec_check

truncate_samples_test()
check_increasing_test()
