'''
Created on May 15, 2016

@author: Werner

Contains argument lists for all subclasses of CommonInit.

optional arguments - opt_classname
compulsory arguments - cmp_classname

Every new subclass of Common init should add its argument lists here.
'''

# Elements

## Level 1:
opt_Element = ['nodes', 'stiffdef']
cmp_Element = ['E', 'nu', 'density']
##

### Level 2
opt_ElemBeamGrid = opt_Element + []
cmp_ElemBeamGrid = cmp_Element + ['A', 'Ik', 'Iy', 'k']
###

# /Elements

# Nodes

## Level 1:
opt_Node = []
cmp_Node = ['x','y','nnodedofs','F_ext','supports']
##

# /Nodes