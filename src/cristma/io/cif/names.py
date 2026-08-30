"""Ordered CIF 1.1 aliases used by structural mapping."""

CELL_A = ("_cell_length_a",)
CELL_B = ("_cell_length_b",)
CELL_C = ("_cell_length_c",)
CELL_ALPHA = ("_cell_angle_alpha",)
CELL_BETA = ("_cell_angle_beta",)
CELL_GAMMA = ("_cell_angle_gamma",)

HM_SYMBOL = ("_space_group_name_h-m_alt", "_symmetry_space_group_name_h-m")
HALL_SYMBOL = ("_space_group_name_hall", "_symmetry_space_group_name_hall")
IT_NUMBER = ("_space_group_it_number", "_symmetry_int_tables_number")
SETTING = ("_space_group_it_coordinate_system_code", "_symmetry_cell_setting")
ORIGIN_CHOICE = ("_space_group_origin_choice",)
SYMMETRY_OPERATION = (
    "_space_group_symop_operation_xyz",
    "_symmetry_equiv_pos_as_xyz",
)

ATOM_LABEL = ("_atom_site_label",)
ATOM_TYPE = ("_atom_site_type_symbol",)
FRACT_X = ("_atom_site_fract_x",)
FRACT_Y = ("_atom_site_fract_y",)
FRACT_Z = ("_atom_site_fract_z",)
OCCUPANCY = ("_atom_site_occupancy",)
OXIDATION = ("_atom_site_oxidation_number",)
WYCKOFF = ("_atom_site_wyckoff_symbol",)
MULTIPLICITY = ("_atom_site_symmetry_multiplicity",)
DISORDER_ASSEMBLY = ("_atom_site_disorder_assembly",)
DISORDER_GROUP = ("_atom_site_disorder_group",)
U_ISO = ("_atom_site_u_iso_or_equiv",)
B_ISO = ("_atom_site_b_iso_or_equiv",)

ANISO_LABEL = ("_atom_site_aniso_label",)
ANISO_U11 = ("_atom_site_aniso_u_11",)
ANISO_U22 = ("_atom_site_aniso_u_22",)
ANISO_U33 = ("_atom_site_aniso_u_33",)
ANISO_U12 = ("_atom_site_aniso_u_12",)
ANISO_U13 = ("_atom_site_aniso_u_13",)
ANISO_U23 = ("_atom_site_aniso_u_23",)

FORMULA = ("_chemical_formula_sum", "_chemical_formula_structural")
METADATA = {
    "mineral_name": ("_chemical_name_mineral",),
    "common_name": ("_chemical_name_common",),
    "systematic_name": ("_chemical_name_systematic",),
    "publication_title": ("_publ_section_title",),
    "journal": ("_journal_name_full",),
    "year": ("_journal_year",),
    "volume": ("_journal_volume",),
    "page_first": ("_journal_page_first",),
    "page_last": ("_journal_page_last",),
    "doi": ("_journal_paper_doi",),
}
