1. Goto https://esgf-node.llnl.gov/search/cmip6/
2. Specify following search criteria in the left tab filters:

  - Source ID: GFDL-CM4
  - ExperimentID: amip
  - Variant Label: r1i1p1f1
  - Frequency: mon
  - Variable: cllcalipso, clmcalipso, clhcalipso, cltcalipso, rsut, rsdt, rsutcs
3. Press Search
4. Download 'wget script' against each of the datasets returned. 
5. Upload the script to taki workspace where data needs to be downloaded.
6. A downloaded version of the script has been placed for you. 
6. Run the script using the following command:
    wget-XXXXXX -s
7. You can refer to the FAQs for more assistance: https://esgf.github.io/esgf-user-support/faq.html#error-http-request-to-openid-relying-party-service-failed
