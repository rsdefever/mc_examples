import numpy as np
from scipy.stats import linregress
import unyt as u
import pandas as pd

t_start = 100  # 10 ps
t_end = 20000  # 20 ns

sizes = [1.0, 1.5, 2.0]
npairs = [0, 4, 8]

data = []
for size in sizes:
    for npair in npairs:
        new_row = {}
        filen = f"../{size}nm_{npair}pairs/msd_water.xvg"
        # Load in the data and throw out the part we don't use to fit
        msd = np.genfromtxt(filen, skip_header=21)
        locs = np.where(np.logical_and(msd[:,0] > t_start, msd[:,0] < t_end))[0]
        msd = msd[locs]
        result = linregress(msd[:,0], msd[:,1])
        print(f"Fit MSD: R^2 = {result.rvalue}")
        diffusivity = (result.slope / 4.0) * u.nm**2 / u.ps
        
        new_row["size_nm"] = size
        new_row["npairs"] = npair
        new_row["diffusivity_water_cm^2/s"] = diffusivity.to_value(u.cm**2 / u.s)

        if npair > 0:
            filen = f"../{size}nm_{npair}pairs/msd_Cl.xvg"
            msd = np.genfromtxt(filen, skip_header=21)
            locs = np.where(np.logical_and(msd[:,0] > t_start, msd[:,0] < t_end))[0]
            msd = msd[locs]
            result = linregress(msd[:,0], msd[:,1])
            print(f"Fit MSD: R^2 = {result.rvalue}")
            diffusivity = (result.slope / 4.0) * u.nm**2 / u.ps
 
            new_row["diffusivity_Cl_cm^2/s"] = diffusivity.to_value(u.cm**2 / u.s)

            filen = f"../{size}nm_{npair}pairs/msd_Na.xvg"
            msd = np.genfromtxt(filen, skip_header=21)
            locs = np.where(np.logical_and(msd[:,0] > t_start, msd[:,0] < t_end))[0]
            msd = msd[locs]
            result = linregress(msd[:,0], msd[:,1])
            print(f"Fit MSD: R^2 = {result.rvalue}")
            diffusivity = (result.slope / 4.0) * u.nm**2 / u.ps
 
            new_row["diffusivity_Na_cm^2/s"] = diffusivity.to_value(u.cm**2 / u.s)
        else:
            new_row["diffusivity_Na_cm^2/s"] = np.nan
            new_row["diffusivity_Cl_cm^2/s"] = np.nan

        data.append(new_row)

df = pd.DataFrame(data)
df.to_csv("results.csv")

