#!/bin/bash
cd /home
wget https://s3.amazonaws.com/preview-ica-expression-data/ica_bone_marrow_h5.h5
python example_ica.py
