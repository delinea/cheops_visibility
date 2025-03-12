import numpy as np
from matplotlib import pyplot as plt, ticker
from astropy import time, units as u, coordinates as coord
try:
    import tqdm
except ImportError:
    tqdm = None


REF_DATE = time.Time("2001-01-01")


def get_Sun_separation(target_coord, date):
    sun = coord.get_sun(date)
    sep = sun.separation(target_coord)
    return sep


def get_cheops_visibility(target_list, sun_exclusion_angle_in_deg=120,
                          time_res_in_days=1/24):
    year = REF_DATE + np.arange(0, 365, time_res_in_days)*u.day
    visibility = dict(time=year)
    if isinstance(target_list, str):
        target_list = [target_list, ]
    elif isinstance(target_list, coord.SkyCoord):
        if len(target_list.shape) == 0:
            target_list = [target_list, ]
        target_list = [target_list, ]
    if tqdm is not None:
        target_list = tqdm.tqdm(target_list, leave=False,
                                desc="Computing target visibility",)
    for i, target in enumerate(target_list):
        if isinstance(target, coord.SkyCoord):
            target_name = "Target #{}".format(i+1)
        elif isinstance(target, tuple):
            try:
                target_name, target_coord = target
                target_name = str(target_name)
                assert isinstance(target_coord, coord.SkyCoord)
            except Exception:
                err_txt = ("targets provided as a tuple must have the "
                           "following format: "
                           "(target_name, astropy.coordinates.SkyCoord)")
                raise ValueError(err_txt)
        else:
            target_name = target
            try:
                target_coord = coord.get_icrs_coordinates(target)
            except coord.name_resolve.NameResolveError:
                err_txt = ("'{}' cannot be resolved by Sesame: you can try "
                           "using astropy 'SkyCoord' coordinates instead"
                           .format(target))
                raise coord.name_resolve.NameResolveError(err_txt)
        sep = get_Sun_separation(target_coord, year)
        target_vis = sep > sun_exclusion_angle_in_deg*u.deg
        visibility[target_name] = target_vis
    return visibility


def compute_visibility_ranges(target_list, sun_exclusion_angle_in_deg=120,
                              time_res_in_days=1/24, print_results=True,
                              plot_results=False):
    kwargs = dict(sun_exclusion_angle_in_deg=sun_exclusion_angle_in_deg,
                  time_res_in_days=time_res_in_days)
    visibility = get_cheops_visibility(target_list, **kwargs)
    time = visibility.pop("time")
    vis_ranges = dict()
    for name, vis in visibility.items():
        dvis = np.diff(vis.astype(int))
        t_start = time[np.concatenate((vis[:1], dvis == 1))]
        t_end = time[np.concatenate((dvis == -1, vis[-1:]))]
        assert len(t_start) == len(t_end)
        vis_ranges[name] = np.transpose(np.stack((t_start, t_end)))

    if print_results:
        reduced_vis_ranges = dict()
        name_len = 10
        for name, vr in vis_ranges.items():
            rvr_i = []
            for vr_i in vr:
                t_start = vr_i[0].datetime
                t_end = vr_i[1].datetime
                rvr_j = "{}  <->  {}".format(t_start.strftime("%b %d"),
                                             t_end.strftime("%b %d"))
                rvr_i.append(rvr_j)
            reduced_vis_ranges[name] = rvr_i
            if len(name) > name_len:
                name_len = len(name)
        print("-"*80)
        for name, vr in reduced_vis_ranges.items():
            if len(vr) > 0:
                vr = "    |    ".join(vr)
            else:
                vr = "{:^19}".format("NOT VISIBLE")
            print("  >  {{:{}}}:    {{}}".format(name_len).format(name, vr))
        print("-"*80)

    if plot_results:
        n_targets = len(vis_ranges)
        fig, ax = plt.subplots(figsize=(8, 2+0.5*n_targets),
                               gridspec_kw=dict(left=0.25))
        ax.set_title("CHEOPS visibility", fontsize=16)
        for i, (name, vr) in enumerate(vis_ranges.items()):
            kwargs = dict(ls="-", lw=3, marker="|", ms=15, mew=3,
                          c="C{}".format(i))
            ax.plot(np.nan, np.nan, label=name, **kwargs)
            for vr_i in vr:
                t_i = [vr_j.jd for vr_j in vr_i]
                ax.plot(t_i-REF_DATE.jd, np.full(2, n_targets-1-i), **kwargs)
            ax.set_ylim([-0.5, n_targets-0.5])
        # ax.legend()
        ax.set_xlim([0, 365])
        ax.set_yticks(range(n_targets), list(vis_ranges.keys())[::-1])
        days_per_month = {"Jan": 31, "Feb": 28, "Mar": 31, "Apr": 30,
                          "May": 31, "Jun": 30, "Jul": 31, "Aug": 31,
                          "Sep": 30, "Oct": 31, "Nov": 30, "Dec": 31}
        months = list(days_per_month.keys())
        days = np.array(list(days_per_month.values()))
        cum_days = np.cumsum(days)
        major_ticks = ticker.FixedLocator(np.concatenate(((0, ), cum_days)))
        minor_ticks = ticker.FixedLocator(cum_days - 0.5*days)
        ax.xaxis.set_major_locator(major_ticks)
        ax.xaxis.set_minor_locator(minor_ticks)
        ax.xaxis.set_major_formatter(ticker.NullFormatter())
        ax.xaxis.set_minor_formatter(ticker.FixedFormatter(months))
        ax.tick_params(axis="x", which="minor", length=0)

    return vis_ranges


if __name__ == "__main__":
    # Importing modules
    import argparse

    # Input parameters
    argparser_kwargs = dict(prog="CHEOPSVisibility",
                            description="Is my target visible with CHEOPS?")
    parser = argparse.ArgumentParser(**argparser_kwargs)
    parser.add_argument("target_name", type=str, nargs="+",
                        help="Name(s) of the target(s)")
    parser.add_argument("-p", "--plot", action="store_true",
                        help="Plot the visibility")
    parser.add_argument("-sea", "--sun_exclusion_angle", default=117,
                        type=float, help="Sun Exclusion Angle [deg]")
    parser.add_argument("-dt", "--time_resolution", default=1/24, type=float,
                        help="Temporal resolution [days]")
    args = parser.parse_args()
    target_name = args.target_name
    plot = args.plot
    sea = args.sun_exclusion_angle
    dt = args.time_resolution

    compute_visibility_ranges(target_name, sun_exclusion_angle_in_deg=sea,
                              time_res_in_days=dt, print_results=True,
                              plot_results=plot)
    if plot:
        plt.show()
