import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from slip_rate_tools import inverse_transform_sample, make_pdf, make_cdf

off = pd.read_csv('/home/itchy/research/kf_cosmo/data/offsets.csv')

off['offset_m'] = off.offset_in * 200
t1 = off[off.unit == 'T1']

pdf_range, pdf_vals = make_pdf(t1.offset_m, t1.rel_prob)
cdf_range, cdf_vals = make_cdf(pdf_range, pdf_vals)


t1_samps = inverse_transform_sample(t1.offset_m, t1.rel_prob, 5000)


plt.figure(1)
plt.subplot(311)
plt.plot(pdf_range, pdf_vals)

plt.subplot(312)
plt.plot(cdf_range, cdf_vals)

plt.subplot(313)
plt.hist(t1_samps, bins=50)

plt.show()

