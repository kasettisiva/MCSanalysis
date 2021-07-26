#!/bin/bash
 
cd $MRB_BUILDDIR
nice ninja install $@
cd -

