
def GM_adjust_map_core_raw(map_):
    if map_.RunPars.charges_model == "HF-CIS":
        map_.rawcore["frequency_data_file_linear"] = [
            "diff_charges_HF-CIS.txt"]
    # else: B3LYP is present in core by default.
