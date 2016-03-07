# Sort_Forests

# Description
Sorts [Consistent Tree](https://bitbucket.org/pbehroozi/consistent-trees) outputs into contiguous order.
*Written in C, probably only works on Linux*. 

# Motivation
Consistent tree output is spread out over various ``tree_*_*_*.dat`` files with a ``locations.dat`` file
specifying the ``tree_root_id``, ``fileid`` and ``offset`` to read in the data. ``Sort_forests`` resorts
the entire data-set such that all forests are grouped contiguously in the same file -- speeding up IO.

The code is written to be very fast and uses lower level system calls. 

# Installation

## Pre-requisites
1. A C compiler (gcc/icc/clang). gcc is default.

## Compiling
```
$ git clone https://github.com/manodeep/sort_forests/
$ make
```

## Running

```
$ ./sort_forests <input directory> <output directory>
```

Input directory should contain ``forests.list``, ``locations.dat`` and ``tree_*_*_*.dat`` files.
Output directory will contain the new ``forests.list``, ``locations.dat`` and ``tree_*_*_*.dat`` files.
In addition, the output directory also contains ``forests_and_locations_new.dat`` file that contains
``forest id``, ``tree_root_id``, ``fileid``, ``offset``, ``filename`` and a new column _``bytes``_ that
gives the number of (ASCII) bytes in each tree. This ``bytes`` column will let you pre-allocated a
buffer for a tree and read in the entire tree in one go (rather than line by line). 

# Author

``Sort_Forests`` is written/maintained by Manodeep Sinha. Please contact the [author](mailto:manodeep@gmail.com) in
case of any issues. 

# LICENSE

``Sort_Forests`` is released under the MIT license. Basically, do what you want
with the code including using it in commercial application.

# Project URL

* version control (https://github.com/manodeep/sort_forests/)
                               






