import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn
import pandas

from matplotlib import rcParams

rcParams['font.sans-serif'] = "Arial"
rcParams['font.family'] = "sans-serif"
seaborn.set_palette("husl", 3)


def main():

    results = pandas.read_csv("results.csv", index_col=0)
    
    fig, ax = plt.subplots()
    
    ax.scatter(
        results["size_nm"],
        results["diffusivity_water_cm^2/s"], #*1e5,
        c=results["npairs"],
        marker="o",
        alpha=0.6,
        s=120,
        norm=matplotlib.colors.Normalize(vmin=0, vmax=10),
        label="Water",
    )
    ax.scatter(
        results["size_nm"],
        results["diffusivity_Na_cm^2/s"], #*1e5,
        c=results["npairs"],
        marker="s",
        alpha=0.6,
        norm=matplotlib.colors.Normalize(vmin=0, vmax=10),
        s=120,
        label=r"Na$^+$",
    )
    ax.scatter(
        results["size_nm"],
        results["diffusivity_Cl_cm^2/s"], #*1e5,
        c=results["npairs"],
        marker="^",
        alpha=0.6,
        norm=matplotlib.colors.Normalize(vmin=0, vmax=10),
        s=120,
        label=r"Cl$^-$",
    )
    
    
    ax.set_yscale("log")
    ax.set_xlabel("Pore width (nm)", fontsize=24, labelpad=10)
    ax.set_ylabel(r"D (cm$^2$/s)", fontsize=24, labelpad=10)
    
    ax.xaxis.set_major_locator(ticker.AutoLocator())
    #ax.yaxis.set_major_locator(ticker.AutoLocator())
    #ax.yaxis.set_major_locator([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0])
    #ax.set_yticks([0.1, 0.5, 1.0, 5.0])
    #ax.set_yticks([0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9, 2.0, 3.0, 4.0], minor=True)
    #ax.get_yaxis().get_major_formatter().labelOnlyBase = False
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    #ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    
    ax.tick_params(axis="both", direction="in", which="both", length=4, labelsize=18, pad=12, top=True, right=True)
    ax.tick_params(axis="both", which="major", length=8)
    ax.legend(fontsize=18)
    ax.set_ylim(10**(-6), 10**(-4))
    leg = ax.get_legend()
    leg.legendHandles[0].set_color("black")
    leg.legendHandles[1].set_color("black")
    leg.legendHandles[2].set_color("black")
    
    fig.tight_layout()
    fig.savefig("diffusivity.pdf")


if __name__ == "__main__":
    main()

