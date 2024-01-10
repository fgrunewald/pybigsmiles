def _parse_non_st_obj_gsmile(gsmile_str):
    gmsile_str = gmsile_str.replace("|", "").strip()
    if "(" in gsmile_str and ")" in gsmile_str:
        key, value = gsmile_str.split("(")
        value = value[:-1]
    else:
        key = "molecular_weigth"
        value = float(gsmile_str)
    return key, value
