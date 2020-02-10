#!/bin/bash
convert Ypred01_Global.png Ypred01err_Global.png +append Ypredcombined_Global.png
convert Ypred01_Malvinas.png Ypred01err_Malvinas.png +append Ypredcombined_Malvinas.png
convert Ypred0_Malvinas.png Ypred1_Malvinas.png +append Ypred_sep_scales_Malvinas.png
