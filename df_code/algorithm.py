#!/usr/bin/python3 
# algorithm namespace 

import os 
import sys 

#_______________________________________________________________________________
def binarySearch(x,key): 
   # search array x for value key and return bounding indices 
   # input 
   # - key = value to search for in array x 
   # - x   = array to be searched 
   # output 
   # - lo  = index such that x[lo] < key 
   # - hi  = index such that x[hi] > key
   N = len(x) 
   comparisonCount = 1;    # count the number of comparisons (optional)
   lo = 0;
   hi = N-1;
   
   # To start, find the subscript of the middle position.
   position = int( (lo + hi)/2 )
   
   while (x[position]!=key) and (lo<=hi): 
      comparisonCount = comparisonCount + 1;
      if x[position]>key:
         # decrease position by one.
         hi = position - 1;
      else:
         # increase position by one.
         lo = position + 1;
      position = int( (lo + hi)/2 );
   
   dump = lo;
   if lo<=hi:
      # Here we have an exact match to the key
      lo = position;
      hi = position+1;
   else:
      # Here the lower bound surpassed the upper
      lo = hi;
      hi = dump;
   # to safeguard against values that are outside the boundaries of the grid 
   if hi>=N: 
      hi = N-1;
      lo = N-2;
   if hi==0:
      lo = 0;
      hi = 1;

   return lo,hi 
