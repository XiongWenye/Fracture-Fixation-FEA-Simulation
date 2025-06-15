def get_material_type(x, y, bone_width, fixator_thickness, fracture_params):
    """
    Determine the material type.
    """

    fixator_length, bone_length, fracture_width = fracture_params
    if y < fixator_thickness or y > bone_width + fixator_thickness:
        return "fixator"
    fracture_start = (
        (fixator_length - bone_length) / 2 + bone_length / 2 - fracture_width / 2
    )
    fracture_end = fracture_start + fracture_width
    if fracture_start <= x <= fracture_end:
        return "callus"
    else:
        return "bone"
