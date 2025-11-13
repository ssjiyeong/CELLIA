# run_cellia_web.py

import os
import scanpy as sc
from cellia import *
from cellia_web import *


def main():
    print("Running CELLIA...")

    # Load input data
    adata = sc.read_h5ad("dataset/CRC.h5ad")

    # Run processing pipeline
    result = cellia_run(
        adata=adata,
        tissue_db="crc|colon",
        tissue_type="human colorectal cancer (CRC)",
        api_key=os.getenv("API_KEY"),
        n_top_markers=15
    )

    adata.write_h5ad("dataset/CRC_cellia_annotated.h5ad")

    print(result)

    print("CELLIA annotation completed successfully.")
    print("Save to dataset/CRC_cellia_annotated.h5ad")

    rationale_path = "../cellia_output/gpt_explanations_db.json"


    launch_cap_style_app(
        adata=adata,
        port=8051,
        debug=True,
        rationale_json_path=rationale_path,
    )



if __name__ == "__main__":
    main()