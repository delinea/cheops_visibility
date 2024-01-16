# cheops_visibility

`cheops_visibility.py` is a short Python script to compute wether a target is visible or not with the space telescope CHEOPS (Characterising Exoplanet Satellite; [Benz et al. 2021](https://ui.adsabs.harvard.edu/abs/2021ExA....51..109B)).

## Dependencies
`numpy`, `matplotlib`, `astropy`  
optional: `tqdm`

## How to run the code
In a terminal window:
- `python cheops_visibility.py TARGET_NAME` to print out the visibility range of your target `TARGET_NAME`
- `python cheops_visibility.py TARGET_NAME -p` to display a plot with the visibility range of your target `TARGET_NAME`

For several targets:
1. Open a Python environment
2. Generate a list of targets `target_list` containing target names (`str`), target coordinates (`astropy.coordinates.SkyCoord`) or tuples `(target_name, target_coord)`
3. `from cheops_visibility import get_visibility_ranges, plot_visibility`
4. `get_visibility_ranges(target_list)` to print out the visibility ranges of your targets
5. `plot_visibility(target_list)` to plot the visibility ranges of your targets
