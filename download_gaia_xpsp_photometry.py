__author__ = "Felipe Almeida-Fernandes"
__email__ = "felipefer42@gmail.com"

import pandas as pd

from astropy import units as u
from astropy.coordinates import SkyCoord

from astroquery.vizier import Vizier

from gaiaxpy import generate as gaiaxpy_generate
from gaiaxpy import PhotometricSystem, load_additional_systems


def download_gaia3_XPSP_ids(ra0, dec0, width = 1, height = 1):
    """
    Downloads a list of Gaia DR3 IDs with spectra from the vizier database
    in the region covered by a splus image. Changes column names for the
    format used by the pipeline and saves results to a fits file

    Parameters
    ----------
    ra0: float
        center ra coordinate (deg)
        
    dec0: float
        center dec coordinate (deg)
        
    width: float
        width of the square region to query (deg)
        
    height: float
        height of the square region to query (deg)
    
    Returns
    -------
    Returns catalog with gaia IDs for objects with XPcont > 0.5 (selects flag 1)
    """

    # Vizier Table
    cat_id = "I/355/gaiadr3"

    columns = ['RA_ICRS', 'DE_ICRS', 'Source']

    column_filters = {'Plx': '>0', 'Gmag': '<30',
                      'BPmag': '<30', 'RPmag': '<30',
                      'XPcont': ">0.5"}

    # Get data from Vizier
    v = Vizier()

    # Change row limit to get whole catalog
    v.ROW_LIMIT = -1

    # List of columns to download
    v.columns = columns
    
    # Apply restrictions
    v.column_filters = column_filters
    
    # Convert center coordinate to astropy coord
    c = SkyCoord(ra=ra0, dec=dec0, unit=(u.deg, u.deg))
    
    # Send query
    query = v.query_region(c, width = width*u.deg, height = height*u.deg,
                           catalog = cat_id, cache = False)
    
    # Retrieve catalog
    catalog = query[cat_id]
    
    # Convert to dataframe
    catalog = catalog.to_pandas()
    
    return catalog
        

def generate_gaiaxpsp_splus_photometry(gaia3_xpsp_id_catalog,
                                       phot_systems_list):

    """
    Uses gaiaxpy to download gaia synthetic spectra for the sources
    listed in the gaia_source_catalog and generate synthetic photometry
    for the filters in filter_dict.keys()

    gaia3_xpsp_id_catalog: str
        gaia3 XPSP catalog containing the column 'Source'
        
    phot_systems_list: list of PhotometricSystem
        lists the photometric systems to generate synthetic magnitudes

    Returns
    -------
    Returns gaia3_xpsp_id_catalog catalog with added synthetic magnitudes 
    """

    # Rename Source column to match gaiaxpy output
    gaia3_xpsp_id_catalog = gaia3_xpsp_id_catalog.rename({"Source":"source_id"},
                                                         axis = "columns")

    # Get list of sources
    source_list = list(gaia3_xpsp_id_catalog["source_id"])

    # Run gaiaxpy
    print(f"Running gaiaxpy for {len(source_list)} sources.")
    synthetic_photometry = gaiaxpy_generate(source_list,
                                           photometric_system=phot_systems_list)

    # Get RA and DEC from source catalog
    output_pd = gaia3_xpsp_id_catalog.merge(synthetic_photometry, 
                                            how='inner', 
                                            on="source_id")

    # return gaia synthetic photometry
    return output_pd



if __name__ == '__main__':
    
    # Load additional photometric systems to gaiaxpy
    # see https://gaia-dpci.github.io/GaiaXPy-website/tutorials/Additional%20systems%20tutorial.html

    # Include additional photometric systems
    
    # File CTIO_S_PLUS.gaiaXPy_dr3_v1.xml from http://svo2.cab.inta-csic.es/theory/fps/
    # must be in the directory add_filters_xml_path
    
    add_filters_xml_path = "/home/felipe/add_filter_systems"
    PhotometricSystem = load_additional_systems(add_filters_xml_path)
    
    # Running example for S-PLUS
    phot_systems_list = [PhotometricSystem.USER_CTIO_S_PLUS]
    
    # Get source ID catalog
    source_id_catalog = download_gaia3_XPSP_ids(ra0 = 0.0, dec0 = 0.0,
                                                width = 1.0, height = 1.0)
    
    # Generate photometry
    phot_cat = generate_gaiaxpsp_splus_photometry(source_id_catalog,
                                                  phot_systems_list)
                                                  
    print(phot_cat)
                                                  
    # Saving result to a file
    phot_cat.to_csv("example_gaia_photometry.csv", index = False)
