#-----------------------------------------------------------------------------------------------------
# This function uses code from netmeta to ensure that treatment codes are organized in ascending order
# As written, it takes results from the pairwise function.Ultimately will need to think of how this works for winbugs
#-----------------------------------------------------------------------------------------------------

wo = as.character(pa_reac_contrast$treat1) > as.character(pa_reac_contrast$treat2)

pa_reac_contrast$TE[wo] = -pa_reac_contrast$TE[wo]
ttreat1 = as.character(pa_reac_contrast$treat1)
pa_reac_contrast$treat1[wo] = pa_reac_contrast$treat2[wo]
pa_reac_contrast$treat2[wo] = ttreat1[wo]


