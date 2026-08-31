"""Count-driven frame indexing for XYZ and extXYZ sources."""

from __future__ import annotations

from cristma.diagnostics import Diagnostic, Severity

from .document import XyzDocument, XyzFrameSpan


def _diagnostic(severity: Severity, code: str, message: str) -> Diagnostic:
    return Diagnostic(severity=severity, code=code, message=message)


def index_xyz(
    source: str,
    source_name: str | None = None,
) -> tuple[XyzDocument, tuple[Diagnostic, ...]]:
    """Index complete frames without parsing their numerical atom columns."""

    lines = source.splitlines(keepends=True)
    offsets = [0]
    for line in lines:
        offsets.append(offsets[-1] + len(line))

    diagnostics: list[Diagnostic] = []
    frames: list[XyzFrameSpan] = []
    line_index = 0
    while line_index < len(lines):
        if not lines[line_index].strip():
            if all(not line.strip() for line in lines[line_index:]):
                break
            diagnostics.append(
                _diagnostic(
                    Severity.WARNING,
                    "xyz.frame.blank_between_frames",
                    "Blank line between XYZ frames was ignored",
                )
            )
            line_index += 1
            continue

        count_token = lines[line_index].strip()
        try:
            if len(count_token.split()) != 1:
                raise ValueError
            atom_count = int(count_token)
            if atom_count < 0:
                raise ValueError
        except ValueError:
            diagnostics.append(
                _diagnostic(
                    Severity.ERROR,
                    "xyz.frame.count_invalid",
                    f"Invalid XYZ atom count: {count_token!r}",
                )
            )
            break

        comment_index = line_index + 1
        end_line = comment_index + 1 + atom_count
        if comment_index >= len(lines) or end_line > len(lines):
            diagnostics.append(
                _diagnostic(
                    Severity.WARNING,
                    "xyz.frame.incomplete",
                    f"XYZ frame {len(frames)} is incomplete and was ignored",
                )
            )
            break

        frames.append(
            XyzFrameSpan(
                index=len(frames),
                atom_count=atom_count,
                start_offset=offsets[line_index],
                end_offset=offsets[end_line],
                comment_start_offset=offsets[comment_index],
                comment_end_offset=offsets[comment_index + 1],
                atom_rows_start_offset=offsets[comment_index + 1],
            )
        )
        line_index = end_line

    return XyzDocument(source, source_name, tuple(frames)), tuple(diagnostics)


__all__ = ["index_xyz"]
