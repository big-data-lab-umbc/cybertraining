
#Example for downloading ERA-5 dataset from ECMWF website
#!/usr/bin/env python

import cdsapi

c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-single-levels-monthly-means',
    {
        'format': 'netcdf', #specify data format
        'product_type': 'monthly_averaged_reanalysis', #specify temporal resolution
         'variable': [
            'total_cloud_cover', 'total_column_cloud_ice_water', 'total_column_cloud_liquid_water', #specify variable names
        ],
        'year': [
            '1979', '1980', '1981',
            '1982', '1983', '1984',
            '1985', '1986', '1987',
            '1988', '1989', '1990',
            '1991', '1992', '1993',
            '1994', '1995', '1996',
            '1997', '1998', '1999',
            '2000', '2001', '2002',
            '2003', '2004', '2005',
            '2006', '2007', '2008',
            '2009', '2010', '2011',
            '2012', '2013', '2014',
            '2015', '2016', '2017',
            '2018', '2019'
        ],
        'month': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12'
        ],
        'time': '00:00',
        'grid': [1.0, 1.0], #specify spatial resolution
    },
    'ERA5_monthly_clouds_1979_2019_1deg.nc') #save to file
