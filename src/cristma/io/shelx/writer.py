"""Writers for native SHELX documents."""

from .document import ShelxDocument


def write_shelx_document(
    document: ShelxDocument,
    *,
    mode: str = "preserve",
) -> str:
    """Render a SHELX document without altering untouched source text."""

    if mode != "preserve":
        raise ValueError("ShelxDocument supports only preserve-mode writing")
    return document.render_preserved()


__all__ = ["write_shelx_document"]
