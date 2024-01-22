# cheops_visibility

`cheops_visibility.py` is a short Python script to compute whether a target is visible or not with the space telescope CHEOPS (Characterising Exoplanet Satellite; [Benz et al. 2021](https://ui.adsabs.harvard.edu/abs/2021ExA....51..109B)).

## Dependencies
`numpy`, `matplotlib`, `astropy`  
optional: `tqdm`

## How to run the code
In a terminal window, run one of the following lines:
- `python cheops_visibility.py TARGET_NAME_1 [TARGET_NAME_2 ...]` to print out the visibility range of your targets
- `python cheops_visibility.py -p TARGET_NAME_1 [TARGET_NAME_2 ...]` to display a plot with the visibility range of your targets

*Note that target names with spaces should be written within quotes.*

In a Python shell or notebook:
```
from matplotlib import pyplot as plt
from astropy import coordinates as coord, units as u
from cheops_visibility import compute_visibility_ranges

# Generate a list of targets
## the list can contain target names (str)
target_list = ["KELT-11", "WASP-189"]
## or target coordinates (astropy.coordinates.SkyCoord)
target_list += [coord.SkyCoord(ra=0*u.deg, dec=0*u.deg), ]
## or tuples in the format (target_name, target_coord)
target_list += [("My star", coord.SkyCoord(ra=266.4168371*u.deg, dec=29.0078106*u.deg)), ]

print_results = True    # to print out the visibility ranges of your targets
plot_results = True     # to plot the visibility ranges of your targets
visibility_ranges = compute_visibility_ranges(target_list, print_results=print_results,
                                              plot_results=plot_results)
if plot_results:
    plt.show()
```
