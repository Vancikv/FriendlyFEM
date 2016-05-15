'''
Created on May 12, 2016

@author: Werner
'''

base_optional = ['nodes', 'stiffdef']
base_compulsory = ['E', 'nu', 'density']
beam_lin_optional = base_optional + []
beam_lin_compulsory = base_compulsory + ['A', 'Iy', 'Ik', 'k']
