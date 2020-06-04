# Standard library imports.
import os
import urllib
#import xml.etree.ElementTree as ET
import lxml.etree as ET

# Third party imports.
import requests

# https://opendap.larc.nasa.gov/opendap/hyrax/CALIPSO/LID_L1-Standard-V4-10/catalog.xml

# https://opendap.larc.nasa.gov/opendap/hyrax/CALIPSO/LID_L1-Standard-V4-10/contents.html

# https://www-calipso.larc.nasa.gov/data/BROWSE/production/V4-10/
# 2020-01-31/2020-01-31_00-51-51_V4.10_1_6.png

# https://www-calipso.larc.nasa.gov/products/lidar/browse_images/exp_showdate.php?browse_date=2020-01-14

# https://www-calipso.larc.nasa.gov/products/lidar/browse_images/
# show_detail.php?s=production&v=V4-10&browse_date=2020-01-31&
# orbit_time=00-51-51&page=1&granule_name=CAL_LID_L1-Standard-V4-10.2020-01-31T00-51-51ZN.hdf


class Sitemap:
    
    catalog_url = "catalog.xml"
    
    def __init__(self, base_url):
        
        self.base_url  = base_url
        self.site_tree = {base_url : {}}
        
        self.crawl(self.base_url)
        
    def crawl(self, url):
        
        print(url)
        url = url + "/" if url[-1] != "/" else url
        url = urllib.parse.urljoin(url, self.catalog_url)
        
        print(url)
        
        r    = requests.get(url)
        root = ET.fromstring(r.content)


def remove_namespace(tag):
    
    tag = tag.split("}")[1] if "}" in tag else tag
    
    return tag


def remove_namespaces(tags):
    
    tags = [remove_namespace(tag) for tag in tags]
    
    return tags


def remove_namespaces_from_attrib(attrib):
    
    keys   = attrib.keys()
    values = attrib.values()
    
    keys = remove_namespaces(keys)
    
    attrib = dict(zip(keys, values))
    
    return attrib
        
        
def xml_to_dict(element):
    
    result = {}

    for child in element.iterchildren():
                
        # Removes namespace.
        tag = child.tag.split("}")[1] if "}" in child.tag else child.tag
        
        if tag == "dataset":
            
            name = child.attrib["name"]
            
            result[name] = xml_to_dict(child)
            
        elif tag == "catalogRef":
            
            attrib = remove_namespaces_from_attrib(child.attrib)
            
            name = attrib["name"]
            url  = "https://opendap.larc.nasa.gov/" + attrib["ID"] + "catalog.xml"
            
            print(url)
            
            result[name] = xml_to_dict(ET.fromstring(requests.get(url).content))
            
    return result
    
    
if __name__ == "__main__":
    
    #sitemap = Sitemap("https://opendap.larc.nasa.gov/opendap/hyrax/CALIPSO")
    
    url  = "https://opendap.larc.nasa.gov/opendap/CALIPSO/LID_L1-Standard-V4-10/catalog.xml"
    root = ET.fromstring(requests.get(url).content)
    
    x = xml_to_dict(root)

