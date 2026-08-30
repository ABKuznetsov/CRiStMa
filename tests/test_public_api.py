from cristma import __version__


def test_public_package_has_semantic_version():
    assert __version__ == "0.1.0"
