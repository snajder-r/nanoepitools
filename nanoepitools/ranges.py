import pandas as pd


def intersect_locations_with_ranges(locations, ranges, nonevalue=-1):
    """
    Takes a list of sites, and a list of genomic regions, and maps
    the sites to the regions in linear time (O(n+m), where n is the number
    of locations and m the number of regions). Requires both locations and
    ranges to be sorted.

    A few important prerequisites on the input:
      * Both dataframes need to have the columns in the order
        specified below, since the function uses itertuples to iterate
        quickly. Column names are irrelevant.
      * Both dataframes need to be sorted by the chromosome column, such
        that we can compare the two. So, if they are both using strings
        for the chromosome, they should both be sorted lexographically.
        If you are unsure about this, set is_sorted to false, so this
        function will perform the correct sorting for you.
      * Both dataframes are allowed to have further columns, as long
        as the first few columns are as specified.

    :param locations: pandas dataframe of the format: (location, ...)
    :param ranges: pandas dataframe of the format:
                   (start_location, end_location,...)
    :param nonevalue: what value to use for a "None" index. Since not all
        datatypes are nullable, you should choose whatever works best
        for your datatype. Default is -1 since this works with int64.
        CAVEAT: If you provide "None" but the index datatype is not
        nullable, pandas will just convert it (to float, for example)
        and you will end up with indices that are not comparable.
    :return: A panda series where the index is the same index as
             in locations, and the value is either the index of the matching
             region, or None if the location is in no region
    """

    region_membership = pd.Series(
        data=nonevalue, index=locations.index, dtype=ranges.index.dtype
    )

    en_ranges = enumerate(ranges.itertuples())
    en_loc = enumerate(locations.itertuples())

    _, region = next(en_ranges)
    _, loc = next(en_loc)

    try:
        """
        When accessing the tuples, remember:
            loc[0]: cpg index
            loc[1]: location
            region[0]: region index
            region[1]: start site
            region[2]: end site
        """
        while True:
            """ If location is behind region, spool location """
            while loc[1] < region[1]:
                _, loc = next(en_loc)

            """ If location is past region, spool region """
            while loc[1] > region[2]:
                _, region = next(en_ranges)

            """ Check all the constraints """
            if region[1] <= loc[1] <= region[2]:
                region_membership[loc[0]] = region[0]
            else:
                """ This happens if we spooled the region past the location """
                pass

            """ Get next cpg location """
            _, loc = next(en_loc)

    except StopIteration:
        pass

    return region_membership
