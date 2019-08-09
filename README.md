# Sinc resampler

A public domain sinc resampler implementation.

Should build on almost anything, tested on Solus, Ubuntu 18.04 and MSYS2.

Includes a Python script for visualizing the result. Useful for seeing (and hearing) the effect of different parameters and tradeoffs (table size, order). Requires Python 3, scipy, numpy and simpleaudio.

Just run `make && bin/test_resamp -c 0.8 -o 32 && python3 show.py` to get a result that might look something like the one below.

![alt text](https://github.com/intealls/sinc_resampler/blob/master/python.png "Figure 1")