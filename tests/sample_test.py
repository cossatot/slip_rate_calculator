import slip_rate_tools as srt
import matplotlib.pyplot as plt

samp_min = 0
samp_max = 4

dist = srt.sample_bounded_normal(3, 2, 1000, sample_min=samp_min, 
                                 sample_max=samp_max)

print(len(dist))

plt.hist(dist)

plt.show()
