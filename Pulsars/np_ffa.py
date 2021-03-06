# -*- coding: utf-8 -*-
from math import log
import numpy as NP

diag = False

def shift_and_add(a,shift):
  """
  Combine two rows with and without a left shift

  This is the inner loop of the fast folding algorithm. It operates
  on a matrix of 2 rows by N columns.  The new upper row is the vector
  sum of the upper and lower row.  The new lower row is the result of
  shifting the lower row by 'shift' and vector adding it to the upper
  row to create the new lower row.
  """
  n_rows,period = a.shape
  if n_rows != 2:
    print("shift_and_add: shift_and_add requires 2 rows")
    print("shift_and_add: received:\n",a)
    return [[],[]]
  # If you try to shift more than the array length you get no shift
  shift = shift % period
  new_upper = a[0] + NP.append(a[1,shift:  ],a[1,:shift])
  new_lower = a[0] + NP.append(a[1,shift+1:],a[1,:shift+1])
  new = NP.append(new_upper,new_lower).reshape((n_rows,period))
  return new
  
def process_column(a,step):
  """
  Process one column in the Staelin (1969) diagram.

  Step 1 processes the left column, step 2 the one right of that, etc.
  """
  n_rows,period = a.shape
  radix = 2**step
  upper_rows = list(range(0,n_rows,radix))
  for upper_row in upper_rows:
    pair_uppers = list(range(upper_row,upper_row+radix/2))
    for pair_upper in pair_uppers:
      pair_lower = pair_upper + radix/2
      offset = pair_upper % radix
      new_pair = shift_and_add(a[pair_upper:pair_lower+1:radix/2],offset)
      if upper_row == 0 and pair_upper == 0:
        new = new_pair
      else:
        new = NP.append(new,new_pair,axis=0)
  return new

def reshape_data(samples,period):
  """
  Reshape 1-D data into a 2-D array with row length of period.

  The array is conceptually similar to a column in Fig. 1 of
  Staelin (1969).  This format makes it easy to implement the
  algorithm.
  
  Note
  ====
  If there are enough samples to justify it, the array is padded with
  zeros to make period*2**N data.  Otherwise, the samples are truncated
  to make perios*2**N data, to make the computation faster.  What is
  'enough' needs to be properly analyzed.  For now, the rule is that
  the number of zeros should not exceed half of the number of samples.
  """
  n_samp = len(samples)
  if diag:
    print("reshape_data:",n_samp,"samples")
  n_rows = len(samples)/period
  if diag:
    print("reshape_data:",n_rows,"rows")
  power_of_2 = int(log(n_rows)/log(2))+1
  if diag:
    print("reshape_data: Exponent =",power_of_2)
  n_rows = 2**power_of_2
  if diag:
    print("reshape_data: New number of rows =",n_rows)
  num_pad = n_rows*period-len(samples)
  if diag:
    print("reshape_data:",num_pad,"samples of 0 for padding")
  if num_pad < n_samp/2:
    print("reshape_data:",num_pad,"zeros to be added")
    padded = NP.append(samples,NP.zeros(num_pad))
    print("reshape_data: Padded data shape",padded.shape)
    a = padded.reshape((n_rows,period))
  else:
    power_of_2 -= 1
    print("reshape_data: New exponent =",power_of_2)
    n_rows = 2**power_of_2
    print("reshape_data: New number of rows =",n_rows)
    a = samples[:n_rows*period].reshape((n_rows,period))
  if diag:
    print("reshape_data: Returned shape:",a.shape)
  return a, power_of_2

def np_ffa(samples,period):
  """
  Numpy FFA - Fast Folding Algorithm
  
  Search 'samples' for pulses with period between 'period' and 'period'+1
  """
  a, power_of_2 = reshape_data(samples,period)
  num_steps = power_of_2 + 1
  for step in range(1,num_steps):
    radix = 2**step
    a = process_column(a,step)
  result = compute_periods(a)
  return result

def compute_periods(a):
  """
  Computes the period corresponding to each row of the final array

  This could be off by one time point.
  """
  n_rows,n_columns = a.shape
  periods = []
  for row in range(n_rows):
    period = float(n_columns) \
              + float(n_columns+1)*float(row)/(n_rows*n_columns)
    periods.append(period)
  return periods,a

def best_period(periods,pulses):
  """
  Find the best period for these data and the summed pulse.
  """
  imax = pulses.argmax()/pulses.shape[1]
  return imax,pulses[imax]
  
