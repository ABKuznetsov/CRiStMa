from cristma.core.values import MissingKind, parse_measured_value


def test_parses_standard_uncertainty_at_last_digits():
    value = parse_measured_value("7.6959(2)", unit="angstrom")

    assert value.value == 7.6959
    assert value.uncertainty == 0.0002
    assert value.raw == "7.6959(2)"
    assert value.unit == "angstrom"


def test_preserves_distinct_cif_missing_states():
    assert parse_measured_value("?").missing is MissingKind.UNKNOWN
    assert parse_measured_value(".").missing is MissingKind.INAPPLICABLE
    assert parse_measured_value(None).missing is MissingKind.ABSENT
