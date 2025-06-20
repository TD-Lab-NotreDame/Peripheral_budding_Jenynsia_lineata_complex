import ipyrad.analysis as ipa
import pandas as pd

# init a conversion tool
converter = ipa.vcf_to_hdf5(
    name="filtered_final_94.recode",
    data="inputs/filtered_final_94.recode.vcf.gz",
    ld_block_size=500,
)
# run the converter
converter.run()
